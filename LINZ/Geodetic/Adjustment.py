#!/usr/bin/python
# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import math
import inspect
import os.path
from collections import namedtuple
import re
import numpy as np
from numpy import linalg

from .Network import Network
from .Station import Station
from .Observation import Observation, ObservationValue

'''
Module to adjust network coordinates.  
'''

# Error classes used by the adjustment

class ConvergenceError( RuntimeError ): pass
class SingularityError( RuntimeError ): pass

# Options controlling an adjustment

class Options( object  ):
    '''
    Class to define the options used in an adjustment.

    '''
    @staticmethod
    def boolOption(value):
        return 'yes'.startswith(value.lower()) or 'true'.startswith(value.lower())

    def __init__(self,config_file=None):
        self.__dict__['_values']={}

    def __getattr__( self, name ):
        if not name.startswith('_') and name in self._values:
            return self.__dict__['_values'][name]
        raise KeyError(name)

    def __setattr__( self, name, value ):
        if name.startswith("_"):
            raise KeyError(name)
        if name in self._values:
            self._values[name]=value

    def set( self, name, value ):
        self.__setattr__(name,value)

    def setupOptions( self, **options ):
        '''
        Adds option and default values to the options array
        '''
        self._values.update(options)

    def update( self, options, addoptions=False ):
        '''
        Update options from options list, but only if the option is already defined.
        '''
        if isinstance (options,Options):
            options=options.values
        for k,v in options.iteritems():
            if k in self._values or addoptions:
                self._values[k]=v
            else:
                raise KeyError(k)

class ObsEqn( object ):

    def __init__( self, nrow, nparam, diagonalCovar=True, haveSchreiber=False ):
        self.obseq=np.zeros((nrow,nparam))
        self.obsres=np.zeros((nrow,1))
        self.obscovar=np.zeros((nrow,1)) if diagonalCovar else np.zeros((nrow,nrow))
        self.diagonal=diagonalCovar
        self.schreiber=np.ones((nrow,1)) if haveSchreiber else None

    def setCovar( self, obscovar ):
        self.obscovar=np.array(obscovar)
        self.diagonal=obscovar.shape[1] == 1

    def setSchreiber( self, schreiber ):
        self.schreiber=None if schreiber is None else np.array(schreiber)

class Plugin( object ):

    def __init__( self, adjustment ):
        self.adjustment=adjustment
        self.options=adjustment.options
        self.initPlugin()

    def initPlugin( self ):
        pass

    def pluginOptions( self ):
        return {}

class Adjustment( object ):

    '''
    Class to adjust a network.  Can take parameters:

       stations:  Iterator over Station objects
       observations: Iterator over observations
       options: Options used in adjustment
       config_file: Name of a config file to read options from
       config: Additional configuration options as a dictionary (see Options.loadConfig)
       output_file: The name of the listing file

    The following entry points can be overridden in subclasses:

    High level functions:
        setup*
        calculateSolution*!
        report*
        writeOutputFiles*

    Functions to include additional parameters.  setupParameters should call
    addParameter and addParameterUpdate. observationEquation should call Adjustment.observationEquation
    and then update the resultant equations (may not be necessary if parameters are
    applied indirectly via the coordParamMapping).

        setupParameters*
        setupStationCoordMapping*!
        calcStationOffsets*
        observationEquation

    The adjustment can also be modified by adding a plugin as a class which implements 
    any the functions marked with *.  Note that for the calculateSolution function only the last defined
    plugin function will be used.  

    Additional plugin hooks are
        setConfigOption!
        preSetup!
        postSetup
        preSetupParameters!
        postSetupParameters
        sumNormalEquations
        updateObservationEquation

    The order in which plugins are loaded is defines the order in which the plugin hooks will be run.
    The Adjustment class itself is the first plugin loaded.  Those functions marked with ! are run in
    reverse order of loading.
    
    The parameters in the plugin function will be the same as
    for the corresponding Adjustment function.  

        class MyPlugin( Adjustment.Plugin ):

            # Optional plugin configuration 

            def initPlugin( self ):
                ....

            def pluginOptions( self ):
                return dict(
                    item=value,
                    ...
                    )

            # Plugin hook functions ...
            # Reading configuration items

            def setConfigOption( self, item, value ):
                if item=='opt1':
                    self.options.item=value
                elif:
                    ...
                else:
                    return False
                return True

            # Other functions implemented by plugin

            def setupParameters( self )
                ....

    '''

    def __init__( self, 
                 stations=Network(), 
                 observations=[], 
                 options=None,
                 config_file=None,
                 output_file=None,
                 plugins=[],
                 verbose=False):

        self.stations=stations
        self.observations=observations
        self.originalXyz={}

        self.options=Options()
        self.setupOptions()
        if verbose:
            self.options.verbose=True
        outputset=False
        if output_file is not None:
            self.setOutputFile( output_file )
            outputset=True
        self.parameters=[]
        self.plugins=[self]
        self.pluginFuncs={}

        # Each plugin is added before reading configuration
        if plugins is not None:
            for plugin in plugins:
                self.addPluginClass( plugin )
        self.solved=0

        # Setup configuration
        if config_file is not None:
            self.loadConfigFile( config_file, sys.stderr.write )
        if options is not None:
            self.options.update(options)
        if verbose:
            self.options.verbose=True
        if output_file is None:
            self.setOutputFile()

    def setOutputFile( self, output_file=None ):
        if output_file is None:
            output_file = self.options.listingFile
        if isinstance(output_file,basestring):
            output_file=open(output_file,"w")
        self.output=output_file

    def addPluginClass( self, pluginClass ):
        if not inspect.isclass(pluginClass) or not issubclass( pluginClass, Plugin ):
            raise RuntimeError('Invalid Adjustment plugin class '+str(pluginClass))
        for plugin in self.plugins:
            if type(plugin) == pluginClass:
                return
        plugin=pluginClass(self)
        plugin.initPlugin()
        pluginopts=plugin.pluginOptions()
        self.options.update(pluginopts,addoptions=True)
        self.addPlugin( plugin )

    def addPlugin( self, plugin ):
        self.plugins.append(plugin)
        self.pluginFuncs={}

    def addPluginClassByName( self, pluginClassName ):
        # Two possible names for plugin module, one as entered and one
        # converted to mixed case.  To support consistency in configuration
        # file parameters using underscore separator.

        origname=pluginClassName
        altname=re.sub(r'(?:_|^)(\w)',lambda m: m.group(1).upper(),pluginClassName)

        # Try first for module in the path of the script
        path0=os.path.dirname(os.path.abspath(sys.argv[0]))

        # Then in the Adjustment directory.  Search relative to the
        # base directory
        path1=os.path.dirname(os.path.abspath(__file__))
        packages=Adjustment.__module__.split('.')[:-1]
        for package in packages:
            path1=os.path.dirname(path1)
        packages='.'.join(packages)
        module=None

        modnames=set()
        for name in [origname,altname]:
            oldpath=list(sys.path)
            try:
                sys.path[:]=[path0]
                try:
                    module=__import__(name,fromlist=['*'])
                except ImportError:
                    modnames.add(name)
                    pass
                if module is not None:
                    break
                sys.path[:]=[path1]
                modname=name
                if packages:
                    modname=packages+'.'+modname
                try:
                    module=__import__(modname,fromlist=['*'])
                except ImportError:
                    modnames.add(modname)
                    pass
                if module is not None:
                    break
            finally:
                sys.path[:]=oldpath
        if module is None:
            raise RuntimeError('Cannot load Adjustment plugin module '+name+'('+','.join(modnames)+')')
        added=False
        for item in dir(module):
            if item.startswith('_'):
                continue
            value=getattr(module,item)
            if not inspect.isclass(value):
                continue
            if not issubclass(value,Plugin):
                continue
            if value == Plugin:
                continue
            self.addPluginClass(value)
            added=True
        if not added:
            raise RuntimeError('No Adjustment plugin found in module '+pluginClassName)

    def getPluginFunction( self, funcname, **options ):
        # Default options
        if funcname not in self.pluginFuncs:
            prepost=options.pop('runPrePostFunctions',False)
            reverse=options.get('reverse',False)
            prefunc=None
            postfunc=None
            if prepost:
                prefuncname=re.sub('^(.)',lambda m: 'pre'+m.group(1).upper(),funcname)
                options['reverse']=True
                prefunc=self.getPluginFunction(prefuncname,**options)
                postfuncname=re.sub('^(.)',lambda m: 'post'+m.group(1).upper(),funcname)
                options['reverse']=False
                postfunc=self.getPluginFunction(postfuncname,**options)

            firstOnly=options.get('firstOnly',False)
            firstTrue=options.get('firstTrue',False)
            funcs=[
                getattr(p,funcname,None) for p in self.plugins
                if getattr(p,funcname,None) is not None and callable(getattr(p,funcname))
                ]
            if reverse:
                funcs.reverse()
            if firstOnly:
                funcdef=funcs[0]
            elif firstTrue:
                def funcdef( *params ):
                    for f in funcs:
                        result=f(*params)
                        if result is not None and result is not False:
                            return result
                    return
            else:
                funcdef=lambda *params: [f(*params) for f in funcs]
            if prepost:
                impfunc=funcdef
                funcdef=lambda *params: [prefunc(*params),impfunc(*params),postfunc(*params)]
            self.pluginFuncs[funcname]=funcdef

        return self.pluginFuncs[funcname]

    def runPluginFunction( self, funcname, *params, **options ):
        funcdef=self.getPluginFunction(funcname,**options)
        return funcdef(*params)

    def setupOptions(self):
        '''
        Setup options.  Installs adjustment options into an Options class
        '''
        self.options.setupOptions(
            # File names
            listingFile=None,
            stationFiles=[],
            dataFiles=[],
            outputStationFile=None,
            outputStationFileGeodetic=None,
            outputStationCovariances=False,
            outputStationErrorEllipses=False,
            outputStationOffsets=False,
            outputResidualFile=None,
            outputResidualFileOptions={},

            # Station options
            fixedStations=[],
            acceptStations=[],
            floatStations=[],
            reweightObsType={},
            ignoreMissingStations=False,

            # Adjustment options
            minIterations=0,
            maxIterations=10,
            convergenceTolerance=0.0001,
            adjustENU=False,
            refractionCoefficient=0.075,

            # Output options
            verbose=False,

            # Specific outputs - only apply if verbose is True
            debugStationOffsets=False,
            debugObservationEquations=False,
            debugFloatStations=False,
            )

    def splitConfigValue( self, value ):
        '''
        Split a configuration value based on white space.  Values
        may be quoted with ".." to include whitespace 
        '''
        fieldre=r'\s+(\"[^\"]*\"|[^\s\"]\S*)'
        value=' '+value
        if not re.match('^(?:'+fieldre+')+$',value):
            return None
        parts=[]
        for m in re.finditer(fieldre,value):
            item=m.group(1)
            if item.startswith('"'):
                item=item[1:-1]
            parts.append(item)
        return parts

    def _splitStationList( self, stnlist ):
        for s in stnlist.split():
            if s.startswith('@'):
                n=Network()
                n.readCsv(s[1:])
                for s in n.stations():
                    yield s.code()
            else:
                yield s

    def setConfigOption( self, item, value ):
        item=item.lower()
        value=value.strip()
        # Plugins
        if item == 'use_plugin':
            self.addPluginClassByName( value )
        # File names
        elif item == 'listing_file':
            self.options.listingFile=value
        elif item == 'coordinate_file':
            self.options.stationFiles.append(value)
        elif item == 'data_file':
            parts=self.splitConfigValue(value)
            if parts is None:
                raise RuntimeError("Invalid data_file definition "+value)
            filename=parts[0]
            attributes={}
            for p in parts[1:]:
                m = re.match(r'^(\w+)\=(.+)$',p)
                if m:
                    attributes[m.group(1)]=m.group(2)
                else:
                    raise RuntimeError("Invalid data_file attribute "+p)
            self.options.dataFiles.append(
                {'filename':filename,'attributes':attributes})
        elif item == 'output_coordinate_file':
            parts=self.splitConfigValue(value)
            self.options.outputStationFile=parts[0]
            for option in parts[1:]:
                option=option.lower()
                if option == 'geodetic':
                    self.options.outputStationFileGeodetic=True
                elif option == 'xyz':
                    self.options.outputStationFileGeodetic=False
                elif option == 'covariances':
                    self.options.outputStationCovariances=True
                elif option == 'ellipses' or option == 'error_ellipses':
                    self.options.outputStationErrorEllipses=True
                elif option == 'offsets':
                    self.options.outputStationOffsets=True
                else:
                    raise RuntimeError('Invalid output_coordinate_file option '+option)
        elif item == 'residual_csv_file':
            parts=self.splitConfigValue(value)
            self.options.outputResidualFile=parts[0]
            options={}
            for p in parts[1:]:
                m = re.match(r'^(\w+)(?:\=(.+))?$',p)
                if m:
                    options[m.group(1)]=m.group(2)
                else:
                    raise RuntimeError("Invalid data_file attribute "+p)
            self.options.outputResidualFileOptions=options
        # Station selection
        elif item == 'fix':
            self.options.fixedStations.extend(((True,v) for v in self._splitStationList(value)))
        elif item == 'free':
            self.options.fixedStations.extend(((False,v) for v in self._splitStationList(value)))
        elif item == 'accept':
            self.options.acceptStations.extend(((True,v) for v in self._splitStationList(value)))
        elif item == 'reject':
            self.options.acceptStations.extend(((False,v) for v in self._splitStationList(value)))
        elif item == 'float':
            try:
                hfs,vfs,slist=value.split(None,2)
                hfloat=float(hfs)
                vfloat=float(vfs)
                self.options.floatStations.extend((((hfloat,vfloat),v) for v in self._splitStationList(value)))
            except:
                raise RuntimeError("Invalid float option: "+value)
                
        elif item == 'ignore_missing_stations':
            self.options.ignoreMissingStations=Options.boolOption(value)

        # Observation options
        elif item == 'reweight_observation_type':
            if value != '':
                match=re.match(r'([a-z]{2})\s+(\d+(?:\.\d+)?)$',value,re.I)
                if not match:
                    raise RuntimeError('Invalid reweight_observation_type option: '+value)
                typecode=match.group(1).upper()
                factor=float(match.group(2))
                self.options.reweightObsType[typecode]=factor

        # Adjustment options
        elif item == 'convergence_tolerance':
            self.options.convergenceTolerance=float(value)
        elif item == 'max_iterations':
            self.options.maxIterations=int(value)
        elif item == 'min_iterations':
            self.options.minIterations=int(value)
        elif item == 'adjust_enu':
            self.options.adjustENU=Options.boolOption(value)
        elif item == 'refraction_coefficient':
            self.options.refractionCoefficient=float(value)

        # Output options
        elif item == 'verbose':
            self.options.verbose=Options.boolOption(value)

        # Debug options
        elif item == 'debug':
            for debugopt in value.lower().split():
                if debugopt == 'observation_equations':
                    self.options.debugObservationEquations=True
                elif debugopt == 'station_offsets':
                    self.options.debugStationOffsets=True
                elif debugopt == 'float_stations':
                    self.options.debugFloatStations=True
                else:
                    raise RuntimeError('Invalid debug option {0}'.format(debugopt))
        elif item == 'debug_observation_equations':
            self.options.debugObservationEquations=Options.boolOption(value)
        elif item == 'debug_station_offsets':
            self.options.debugStationOffsets=Options.boolOption(value)
        elif item == 'debug_float_stations':
            self.options.debugFloatStations=Options.boolOption(value)
        else:
            raise RuntimeError('Unrecognized configuration item: '+item+': '+value)

    def setConfig( self, item, value ):
        '''
        Set a configuration option based on the text item/value read from a 
        configuration file
        '''
        self.runPluginFunction('setConfigOption',item,value,firstTrue=True,reverse=True)

    def loadConfigFile( self, config_file, write=None ):
        cfg={}
        nerrors=0
        with open(config_file) as cfgf:
            for l in cfgf:
                try:
                    l=l.strip()
                    if len(l)==0 or l.startswith('#'):
                        continue
                    try:
                        parts=l.split(None,1)
                        item=parts[0].lower()
                        value=parts[1] if len(parts)==2 else 'y'
                    except:
                        raise RuntimeError("Invalid configuration line: ",l)
                    self.setConfig( item, value )
                except Exception as ex:
                    if write is not None:
                        write(ex.message)
                        write("\n")
                        nerrors+=1
                    else:
                        raise
        if nerrors > 0:
            raise RuntimeError("Stopped with {0} errors in configuration file".format(nerrors))

    def write( self, message, debug=False ):
        if self.output is not None:
            self.output.write( message )
        if self.options.verbose and not debug: 
            sys.stdout.write(message)

    def writeDebugOutput( self, message ):
        self.write( message, True )

    def loadDataFiles( self ):

        # Load the coordinate files
        if len(self.options.stationFiles) > 0:
            self.write("\nInput coordinate file:\n")
            for crdfile in self.options.stationFiles:
                self.write("  {0}\n".format(crdfile))
                self.stations.readCsv(crdfile)

        # Save initial station coordinates
        originalXyz={}
        for stn in self.stations.stations():
            try:
                xyz=stn.xyz()
                originalXyz[stn.code()]=np.array(xyz)
            except:
                pass
        self.originalXyz=originalXyz

        # Load the data files
        reweight=self.options.reweightObsType
        if len(reweight) > 0:
            self.write("\nReweighting input observations:\n")
            for otype in sorted(reweight.keys()):
                self.write("  Observation type {0} errors scaled by {1:6.3f}\n"
                           .format(otype,reweight[otype]))

        first=True
        for obsfile in self.options.dataFiles:
            filename=obsfile['filename']
            attributes=obsfile['attributes']
            if first:
                self.write("\nObservation files:\n")
                first=False
            self.write("  {0}\n".format(filename))
            filetype=''
            filetype=''
            match=re.match(r'^(\w+)\:(.*)$',filename)
            if match:
                filetype=match.group(1).lower()
                filename=match.group(2)
            if filetype == '':
                match=re.search(r'\.(\w+)$',filename)
                if match:
                    filetype=match.group(1).lower()
            if filetype == '' and re.search(r'\.snx(\.gz)?$',filename,re.I):
                filetype = 'snx'
            if filetype == 'msr':
                from . import MsrFile
                reader=MsrFile.read
            elif filetype == 'snx':
                from . import SinexObsFile
                reader=SinexObsFile.read
            elif filetype == 'csvhor':
                from . import CsvObsFile
                reader=CsvObsFile.readConvertToHorDist
            elif filetype == 'csv':
                from . import CsvObsFile
                reader=CsvObsFile.read
            else:
                raise RuntimeError('Invalid observation file type '+filetype)

            for obs in reader(filename,**attributes):
                typecode=obs.obstype.code
                if typecode in reweight:
                    factor=reweight[typecode]
                    for ov in obs.obsvalues:
                        if ov.stderr is not None:
                            ov.stderr *= factor
                    if obs.covariance is not None:
                        obs.covariance *= factor*factor
                self.observations.append(obs)

        if self.options.acceptStations:
            self.filterObsByStation(self.options.acceptStations)

    def usedStations( self, includeFixed=True ):
        '''
        Get dictionary with keys matching the stations used in the list
        of observations.  If includeFixed is False then only adjusted
        stations are included.
        '''
        used={}
        for o in self.observations:
            for obsval in o.obsvalues:
                used[obsval.inststn]=1
                used[obsval.trgtstn]=1
        used={code:1 for code in used if self.stations.get(code) is not None}
        if not includeFixed:
            fixedstn=self.useStationFunc(self.options.fixedStations,False)
            used={code:1 for code in used if not fixedstn(code)}
        return used

    def missingStationList( self, clean=False ):
        missing=set()
        good_station=lambda x: x is None or self.stations.get(x) is not None
        for obs in self.observations:
            for o in obs.obsvalues:
                if not good_station(o.inststn):
                    missing.add((False,o.inststn))
                if not good_station(o.trgtstn):
                    missing.add((False,o.trgtstn))
        if clean:
            self.filterObsByStation(missing)
        missing=list(missing)
        missing.sort( key=lambda x: x[1] )
        return missing

    def useStationFunc( self, uselist, default=True ):
        '''
        Create a function to test a station code against a set of
        criteria.  The criteria are defined in a uselist which is
        a list of pairs (use,code).  
        
          use is boolean, signifying the station is to be used or not.  
        
          code is one of:
            *       meaning all codes
            re:exp  meaning codes matching the regular expression re
            code    a specific station code.  Note - to enter a specific code
                    starting re:, prefix the code with ':'

        The default status applies if there are no criteria in the list.
        Otherwise the default status is the opposite of the first criteria
        in the list (ie if it is True, then the default is False, and vice versa)
        or None if it is not True or False

        All station code matching is case insensitive
        '''
        if not uselist:
            return lambda code: default

        def usefunc( use, code ):
            code=code.lower()
            if code == '*':
                f=lambda c,s: use
            elif code.startswith('re:'):
                codere=re.compile('^'+code[3:]+'$',re.I)
                f=lambda c,s: use if codere.match(c) else s
            else:
                if code.startswith(':'):
                    code=code[1:]
                f=lambda c,s: use if c==code else s
            return f

        # If first action is to set status to True, then initial status must
        # be false
        status0=uselist[0][0]
        startStatus=False if status0 is True else True if status0 is False else None
        funcs=tuple(usefunc(use,code) for use,code in uselist)

        def usestation( code ):
            if code is None:
                return True
            code=code.lower()
            status=startStatus
            for f in funcs: status=f(code,status)
            return status
        return usestation

    def filterObsByStation( self, uselist ):
        '''
        Filter observations by station.  Takes a list of stations
        keep or stations to remove, or both (though that doesn't make 
        much sense!)
        '''
        usestn=self.useStationFunc(uselist)
        useobs=lambda obs,ov: usestn(ov.inststn) and usestn(ov.trgtstn)
        self.filterObservations(useobs)

    def filterObservations( self, select ):
        '''
        Filter observations using a function select.  Select takes
        parameters (obs, obsvalue) and returns True to retain the 
        observation.
        '''
        goodobs=[]
        for obs in self.observations:
            goodvalues=[]
            goodrows=[]
            nvalue=obs.obstype.nvalue
            for i,o in enumerate(obs.obsvalues):
                if select(obs,o):
                    goodvalues.append(o)
                    goodrows.extend(range(i*nvalue,(i+1)*nvalue))
            if len(goodvalues)==len(obs.obsvalues):
                goodobs.append(obs)
                continue
            if len(goodvalues)==0:
                continue
            covariance=obs.covariance
            if covariance is not None:
                covariance=covariance[goodrows,:][:,goodrows]
            o=Observation(obs.obstype.code,obsdate=obs.obsdate,obsvalue=goodvalues,covariance=covariance)
            goodobs.append(o)
        self.write("Filtering {0} observations to {1}\n".format(len(self.observations),len(goodobs)))
        self.observations=goodobs

    def countObservations( self ):
        counts={}
        for obs in self.observations:
            obstype=obs.obstype
            if obstype not in counts:
                counts[obstype]=[0,0]
            counts[obstype][0] += 1
            counts[obstype][1] += len(obs.obsvalues)
        return counts

    def initParameters( self ):
        self.coordParamMapping={}
        self.parameters=[]
        self.updateFuncs=[]
        self.nparam=0

    def addParameter( self, paramname ):
        '''
        Adds a parameter to the adjustment and returns its parameter number
        '''
        self.parameters.append(paramname)
        self.nparam += 1
        return self.nparam-1

    def addParameterUpdate( self, updateFunc ):
        '''
        Function used to update parameters.  Adjustment residual calculation
        assumes that parameters are updated at each iteration, so that at final
        iteration residuals are not changed.

        The updateFunc should take as parameters the list of calculated
        parameter values, and return a boolean converged status
        '''
        self.updateFuncs.append(updateFunc)

    def setupStationCoordMapping( self ):
        '''
        Defines a mapping from station coordinates to parameters.

        Sets up a dictionary coordParamMapping from station id to a tuple of three values
            paramno: parameter number array - parameters used to define stations
            dxyzdp: (n,3) array defining differential of x,y,z against param
            updatexyz: if true then the station coordinate is updated based on
                the parameters in updateStationCoordParameters

        Default implementation simply adjusts x,y,z or e,n,u coord for each station.

        This routine can be overridden to parameterize station coordinates 
        in terms of other parameters.

        Alternatively this can be done using plugin functions.  In either case the 
        functions should define mappings in the coordParamMapping dictionary.  This
        routine will only set up mappings for stations not included in this 
        dictionary.

        '''

        # Get a list of used stations - remove fixed stations from it
        # Result is a list of stations to adjust

        usedStations=self.usedStations(includeFixed=False)

        mapping=self.coordParamMapping

        adjustenu=self.options.adjustENU
        axes=('east','north','up') if adjustenu else ('X','Y','Z')
        for code in sorted(usedStations.keys()):
            if code in mapping:
                continue
            if adjustenu:
               station=self.stations.get(code)
               prmxyz=station.enu()
            else:
               prmxyz=np.identity(3)
            prmnos=[self.addParameter(code+' '+axis) for axis in axes]
            mapping[code]=(prmnos,prmxyz,True)

    # Update the station coordinates and return a parameter indicating the 
    # maximum coordinate change

    def updateStationCoordParameters( self, paramValues ):
        '''
        Applies the calculated updates to the coordinate parameters and
        calculates the maximum offset at the station.
        '''
        if not self.solved:
            raise RuntimeError("Cannot update parameters if equations not solved")
        maxoffset=0.0
        maxcode=None
        printoffsets=self.options.debugStationOffsets
        for code,mapping in self.coordParamMapping.iteritems():
            stn=self.stations.get(code)
            xyz=stn.xyz()
            dxyz = paramValues[mapping[0]].dot(mapping[1])
            if printoffsets:
                denu = stn.enu().dot(dxyz.T).T
                self.writeDebugOutput("  {0} ENU change {1:.4f} {2:.4f} {3:.4f}\n".format(
                    code,denu[0,0],denu[0,1],denu[0,2]))
            if len(mapping) > 2 and mapping[2]:
                xyz += dxyz.reshape((3,))
                stn.setXYZ(xyz)
            offset=linalg.norm(dxyz)
            if offset > maxoffset:
                maxoffset=offset
                maxcode=code
        return maxoffset,maxcode

    def updateParameters( self, paramValues ):
        '''
        Function to update adjustment parameters based on the least squares solution
        vector paramValues.  Returns the maximum coordinate offset and the code of the
        corresponding station.
        '''
        converged=True
        for update in self.updateFuncs:
            if not update(paramValues):
                converged=False
        coordUpdate,code=self.updateStationCoordParameters( paramValues )
        if coordUpdate > self.options.convergenceTolerance:
            converged=False
        return converged,coordUpdate,code


    def setupParameters( self ):
        '''
        Function to set up adjustment parameters.  
        Parameters are added with addParameter(paramname), which returns
        the parameter row number in the adjustment.
        Note parameters also configured by the coordParamMapping dictionary 
        setupStationCoordMapping()
        '''
        pass

    def setupNormalEquations( self ):
        self.initParameters()
        self.runPluginFunction('setupParameters',runPrePostFunctions=True)
        self.runPluginFunction('setupStationCoordMapping',reverse=True)
        nprm=self.nparam
        self.solved=0
        self.N=np.zeros((nprm,nprm))
        self.b=np.zeros((nprm,1))
        self.ssr=0.0
        self.nobs=0
        return nprm

    def getStation( self, code ):
        return self.stations.get(code)

    def calcStationOffsets( self, obs ):
        '''
        Determines the offsets of stations that apply for a particular observation.
        Returns a station offset for the instrument and target station of the observation.

        The offset can account for factors such as instrument heights and deformation.
        Each offset is returned as an ENU offset local to the station.

        Returns a list of entries, one for each observation value.  Each entry in the list
        has five values:
            offsettype     eg Station.OFFSET_H, Station.OFFSET_ENU, ...
            inst_offset    offset for instrument station
            trgt_offset    offset for target station
            inst_params    None or mapping of offset parameter(s) to offset value
            trgt_params    None or mapping of offset parameter(s) to offset value

        Each set of parameters is defined by
            paramno: parameter number array - parameters used to define stations
            denudp: (n,3) array defining differential of E,N,U against n param.
                    For Station.OFFSET_H this is (n)
                    

        Subclasses can override for additional offsets.  Plugin offsets are added
        to form a total offset
        '''

        offsets=[]
        for o in obs.obsvalues:
            offi=None
            offt=None
            if o.insthgt is not None and o.insthgt != 0.0:
                offi=self.getStation(o.inststn).genu()[2]*o.insthgt
            if o.trgthgt is not None and o.trgthgt != 0.0:
                offt=self.getStation(o.trgtstn).genu()[2]*o.trgthgt
            offsets.append((Station.OFFSET_XYZ,offi,offt,None,None))
        return offsets

    def convertOffsetToXYZ( self, obsvalue, offset ):
        offsettype=offset[0]
        if offsettype==Station.OFFSET_XYZ:
            return offset
        result=[Station.OFFSET_XYZ,offset[1],offset[2],offset[3],offset[4]]
        for i,stn in ((1,obsvalue.inststn),(2,obsvalue.trgtstn)):
            offset=result[i]
            if offset is None:
                continue
            doffset=result[i+2]
            if offsettype==Station.OFFSET_H:
                offset=np.array((0.0,0.0,offset))
                if doffset is not None:
                    ddp=doffset[1]
                    zero=np.zeros(ddp.shape)
                    doffset=(doffset[0],np.array([zero,zero,ddp]).T)
            stn=self.getStation(stn)
            enu=stn.enu() if offsettype==Station.OFFSET_ENU else stn.genu()
            offset=offset.dot(enu)
            if doffset is not None:
                doffset=(doffset[0],doffset[1].dot(enu))
            result[i]=offset
            result[i+2]=doffset
        return result

    def compileStationOffsets( self, obs ):
        offsets=None
        for poffsets in self.runPluginFunction('calcStationOffsets',obs):
            if poffsets is None:
                continue
            compiled=[]
            for i,obsvalue in enumerate(obs.obsvalues):
                poffset=self.convertOffsetToXYZ(obsvalue,poffsets[i])
                if offsets is not None:
                    offset=offsets[i]
                    tpo,offi,offt,doffi,dofft=offset
                    ptpo,poffi,pofft,dpoffi,dpofft=poffset
                    if poffi is not None:
                        if offi is not None:
                            offi=np.add(offi,poffi)
                        else:
                            offi=poffi
                    if pofft is not None:
                        if offt is not None:
                            offt=np.add(offt,pofft)
                        else:
                            offt=pofft
                    if dpoffi is not None:
                        if doffi is not None:
                            doffi[0].extend(dpoffi[0])
                            doffi=(doffi[0],np.vstack((doffi[1],dpoffi[1])))
                        else:
                            doffi=dpoffi
                    if dpofft is not None:
                        if dofft is not None:
                            dofft[0].extend(dpofft[0])
                            dofft=(dofft[0],np.vstack((dofft[1],dpofft[1])))
                        else:
                            dofft=dpofft
                    poffset=(Station.OFFSET_XYZ,offi,offt,doffi,dofft)
                compiled.append(poffset)
            offsets=compiled
        return offsets

    def observationEquation( self, obs ):
        '''
        Forms the observation equations for an observation o, which may
        be an array of observation values
        '''
        refcoef=self.options.refractionCoefficient
        obstype=obs.obstype
        nval=obstype.nvalue
        nrow=len(obs.obsvalues)*nval
        obseqn=ObsEqn(nrow,self.nparam)
        diagonal=True
        if obs.covariance is not None:
            obseqn.setCovar(obs.covariance)
            diagonal=False
        obsres=obseqn.obsres
        obseq=obseqn.obseq
        obscovar=obseqn.obscovar

        offsets=self.compileStationOffsets( obs )
        offsettype,offi,offt,doffi,dofft=Station.OFFSET_XYZ,None,None,None,None

        for i,o in enumerate(obs.obsvalues):
            stf=self.stations.get(o.inststn)
            if stf is None:
                raise RuntimeError("Station "+o.inststn+" is not defined")
            stt=None
            if o.trgtstn is not None:
                stt=self.stations.get(o.trgtstn)
                if stt is None:
                    raise RuntimeError("Station "+o.trgtstn+" is not defined")
            # Get XYZ offsets and differential with respect to parameters
            if offsets is not None:
                offsettype,offi,offt,doffi,dofft=offsets[i]

            # Calculate the observed value and its dependence on XYZ coordinates
            calcval,ddxyz0,ddxyz1=obstype.calcobs(
                stf,trgtstn=stt,instofst=offi,trgtofst=offt,refcoef=refcoef,
                ddxyz=True,offsettype=offsettype)

            # Set up the observation residual and covariance
            i0=i*nval
            i1=i0+nval
            if nval == 1:
                obsres[i,0]=o.value-calcval
                if diagonal:
                    obscovar[i,0]=o.stderr*o.stderr
            else:
                obsres[i0:i1,0]=np.array(o.value)-calcval

            # Set up the coordinate parameters
            mapping=self.coordParamMapping.get(o.inststn,None)
            if mapping is not None:
                obseq[i0:i1,mapping[0]]=mapping[1].dot(ddxyz0)
            mapping=self.coordParamMapping.get(o.trgtstn,None)
            if mapping is not None:
                obseq[i0:i1,mapping[0]]+=mapping[1].dot(ddxyz1)

            # Set up the offset parameters
            for stn,doff,dxyz in ((stf,doffi,ddxyz0),(stt,dofft,ddxyz1)):
                if doff is None:
                    continue
                prms,dprm=doff
                obseq[i0:i1,prms]=dprm.dot(dxyz)

        # Need special handling for horizontal angles

        if obstype.code == 'HA':
            assert diagonal
            # Round to multiple of 180.
            obsres[:] = np.remainder(obsres-obsres[0]+180.0,360.0)-180.0
            # Remove weighted mean residual
            wgt=1.0/obscovar
            wgtsum=wgt.sum()
            obsres -= (wgt.T.dot(obsres))/wgtsum
            schreiber=np.ones((nrow,1))
            obseqn.setSchreiber(schreiber)

        # Handle 360 degree wrapping for azimuths

        elif obstype.code == 'AZ':
            obsres[:] = np.remainder(obsres+180.0,360.0)-180.0

        self.runPluginFunction('updateObservationEquation',obs,obseqn)
        return obseqn

    # Note: May be scope for simplification/optimisation here using more
    # sophisticated routines such as numpy.linalg.lstsqu.  Also look at
    # sparse matrix options for observation equations (and possibly normal
    # equations)

    def sumObservation( self, obseqn ):
        assert self.solved==0,"Cannot sum observations after they have been solved"
        if obseqn.diagonal:
            weight=np.diag(1.0/obseqn.obscovar.flatten())
        else:
            weight=linalg.inv(obseqn.obscovar)
        A=obseqn.obseq
        y=obseqn.obsres
        self.N += A.T.dot(weight.dot(A))
        self.b += A.T.dot(weight.dot(y))
        self.ssr += y.T.dot(weight.dot(y))[0,0]
        self.nobs += y.size
        A2=obseqn.schreiber
        if A2 is not None:
            ws=A2.T.dot(weight.dot(A2))
            if ws.size==1:
                ws=1.0/ws
            else:
                ws=np.linalg.inv(ws)
            n12=A.T.dot(weight.dot(A2))
            b2=A2.T.dot(weight.dot(y))
            self.N -= n12.dot(ws.dot(n12.T))
            self.b -= n12.dot(ws.dot(b2))
            self.ssr -= b2.T.dot(ws.dot(b2))[0,0]
            self.nobs -= A2.shape[1]

    def covariance( self ):
        '''
        Return the parameter covariance matrix
        '''
        assert self.solved > 0,"calcResiduals requires the equations to be solved"
        if self.solved==1:
            self.N=linalg.inv(self.N)
            self.solved=2
        return self.N

    def stationENUCovariance( self, code ):
        '''
        Return the ENU coordinate covariance matrix
        '''
        assert self.solved > 0,"calcResiduals requires the equations to be solved"
        N=self.covariance()
        mapping=self.coordParamMapping.get(code)
        stn=self.stations.get(code)
        if stn is None or mapping is None or len(mapping[0]) == 0:
            return np.zeros((3,3))
        obseq=np.zeros((self.nparam,3))
        obseq[mapping[0]]=mapping[1]
        cvr=obseq.T.dot(N.dot(obseq))
        enu=stn.enu()
        cvrenu=enu.dot((enu.dot(cvr)).T)
        return cvrenu

    def calcResidual( self, obs ):
        '''
        Calculate the observation residuals. Returns 

        resvalue, rescovar
        '''
        assert self.solved > 0,"calcResiduals requires the equations to be solved"
        N=self.covariance()
        obseqn=self.observationEquation(obs)
        A=obseqn.obseq
        calccovar=A.dot(N.dot(A.T))
        # Apply affect of Schreiber equations to covariance matrix...?
        A2=obseqn.schreiber
        if A2 is not None:
            if obseqn.diagonal:
                wda2=(1.0/obseqn.obscovar)*A2
            else:
                wda2=linalg.solve(obseqn.obscvr,A2)

            ws=A2.T.dot(wda2)
            if ws.size==1:
                ws=1.0/ws
            else:
                ws=np.linalg.inv(ws)
            n12=A.T.dot(wda2)
            u12=n12.dot(ws)
            tmp1=u12.T.dot(N)
            tmp2=A.dot(tmp1.T).dot(A2.T)
            calccovar -= tmp2
            calccovar -= tmp2.T
            calccovar += A2.dot(ws.dot(A2.T))
            calccovar += A2.dot(tmp1.dot(u12)).dot(A2.T)
            
        rescovar=np.diag(obseqn.obscovar.flatten()) if obseqn.diagonal else obseqn.obscovar
        rescovar-=calccovar
        resvalue=obseqn.obsres

        return resvalue, rescovar


    def writeObservationEquations(self,compactRows=False):
        '''
        Debugging option to print observation equations
        '''
        self.setupNormalEquations()
        self.writeDebugOutput("\nParameters:\n")
        for p in self.parameters:
            self.writeDebugOutput("    {0}\n".format(p))
        self.writeDebugOutput("\nObservation Equations:\n")
        for o in self.observations:
            ov=o.obsvalues[0]
            self.writeDebugOutput("\nObservation: {0} from {1} to {2}{3}\n".format(
                o.obstype.code,ov.inststn,ov.trgtstn,' ...' if len(o.obsvalues) > 1 else ''))
            oe=self.observationEquation(o)
            factor=1.0
            offset=0
            # Convert to compatible format with SNAP for checking 
            # (radians, 0 offset for first HA residuals)
            if o.obstype.code in ('AZ','ZD','HA'):
                factor=math.radians(1)
            if o.obstype.code == 'HA':
                offset=oe.obsres[0]
            self.writeDebugOutput("  ObsRes: {0}\n".format((oe.obsres-offset)*factor))
            self.writeDebugOutput("  ObsEq:  {0}\n".format(oe.obseq*factor))
            self.writeDebugOutput("  ObsCvr: {0}\n".format(oe.obscovar*factor*factor))
            if oe.schreiber != None:
                self.writeDebugOutput("  Schrb: {0}\n".format(oe.schreiber))

    def observationEquations( self ):
        for obs in self.observations:
            yield self.observationEquation(obs)

    def sumNormalEquations( self ):
        for oe in self.observationEquations():
            self.sumObservation(oe)
        if len(self.options.floatStations) > 0:
            self.sumFloatStations()

    def sumFloatStations(self):
        # Sum float stations if requested...
        # Create a GX observation from the station coordinates. 
        # Note that if there is no original coordinate then the 
        # float coordinate changes at each iteration, so the
        # solution will not converge quickly!
        func=self.useStationFunc(self.options.floatStations,None)
        for code in self.coordParamMapping:
            params=self.coordParamMapping[code][0]
            if len(params) == 0:
                continue
            float=func(code)
            if float is None:
                continue
            if self.options.debugFloatStations:
                self.writeDebugOutput(
                    "  Floating {0}: {1:.4f} {2:.4f}\n".format(code,float[0],float[1]))
            enuerr=np.diag((float[0],float[0],float[1]))
            enuerr=enuerr*enuerr
            stn=self.stations.get(code)
            code=stn.code()
            xyz=self.originalXyz.get(code)
            if xyz is None:
                if self.options.debugFloatStations:
                    self.writeDebugOutput("  Using transient coord for float of {0}\n".format(code))
                xyz=stn.xyz()
            xyz=np.array(xyz)
            enu=stn.enu()
            cvr=enu.T.dot(enuerr.dot(enu))
            obs=Observation('GX',covariance=cvr)
            obs.addObservation(ObservationValue(code,value=xyz))
            oe=self.observationEquation(obs)
            self.sumObservation(oe)

    def solveEquations( self ):
        try:
            self.x=linalg.solve(self.N, self.b).flatten()
            self.solved=1
        except linalg.LinAlgError:
            self.solved=-1
            raise SingularityError()

    def runOneIteration( self ):
        self.setupNormalEquations()
        self.runPluginFunction('sumNormalEquations')
        self.solveEquations()
        self.dof=self.nobs-self.nparam
        self.seu=math.sqrt(self.ssr/self.dof) if self.dof > 0 else 1.0
        return self.updateParameters( self.x )

    def writeObservationSummary( self ):
        self.write("\nObservation summary:\n")
        counts=self.countObservations()
        for obstype in sorted(counts.keys(),key=lambda x:x.code):
            count=counts[obstype]
            suffix=''
            if obstype.code == 'HA':
                suffix=' in {0} rounds'.format(count[0])
            self.write("  {0} {1} observations{2}\n".format(count[1],obstype.name,suffix))
        self.write("  Using refraction coefficient {0:.3f}\n".format(self.options.refractionCoefficient)) 

    def residuals( self ):
        '''
        Iterator function to return the observations and residuals.  At each iteration 
        returns an observation, residual vector, and residual covariance vector.
        '''
        for obs in self.observations:
            res,rescovar=self.calcResidual(obs)
            yield obs,res,rescovar
 
    def writeSummaryStats( self ):
        self.write("\nSum of squared residuals:     {0:.6f}\n".format(self.ssr))
        self.write("Number of parameters:         {0}\n".format(self.nparam))
        self.write("Number of observations:       {0}\n".format(self.nobs))
        self.write("Degrees of freedom:           {0}\n".format(self.dof))
        self.write("Standard error of unit weight {0:.4f}\n".format(self.seu))

    def calcResidualSummary( self, keyfunc=None ):
        '''
        Returns a summary of observation residuals classified by key
        which defaults to the observation type code.  Otherwise should be
        a function taking an observation as a parameter and returning key

        The summary is returned as a dictionary based on key, for which 
        each item is a tuple (ssr,nvalue).  SSR is the sum of squared 
        standardised residuals, nvalue is the number of values summed.
        '''
        if keyfunc is None:
            keyfunc=lambda obs: obs.obstype.code

        obsTypeSummary={}
        for obs,res,rescvr in self.residuals():
            key=keyfunc(obs)
            if key not in obsTypeSummary:
                obsTypeSummary[key]=[0.0,0]
            summary=obsTypeSummary[key]

            # What follows is not strictly rigorous ... but I think is indicative of 
            # magnitude of residuals

            if obs.covariance is not None:
                tol=1.0e-12 # Random!
                cvrinv=linalg.pinv(rescvr,tol)
                rank=linalg.matrix_rank(rescvr,tol)
                summary[0] += res.T.dot(cvrinv.dot(res))[0,0]
                summary[1] += rank

            elif obs.obstype.nvalue == 1:
                for obsv,resv,rescvr in zip(obs.obsvalues,res.flatten(),rescvr.diagonal().flatten()):
                    stderr=obsv.stderr
                    if rescvr >= stderr*stderr*1.0e-10:
                        summary[0] += resv*resv/rescvr
                        summary[1] += 1
            else:
                # Currently not handling this (vector data, no covariance)..
                pass 
        return obsTypeSummary

    def writeResidualSummary( self, keyfunc=None, title=None ):
        summary=self.calcResidualSummary(keyfunc)
        if title is None:
            title='Summary of residuals'
            self.write("\n{0}\n\n  {1:<10s} {2:4s} {3:8s}\n"
                       .format(title,'Type','NRes','  RMS'))
        for key in sorted(summary.keys()):
            ssr,nres=summary[key]
            rms=math.sqrt(ssr/nres) if nres > 0 else 1.0
            self.write("  {0:<10} {1:4d} {2:8.4f}\n".format(key,nres,rms))

    def writeResidualCsv( self, csvfile, options={} ):
        '''
        Dumps residuals of scalar observations to a CSV file
        '''
        import csv

        attributes=[]
        if 'attributes' in options:
            attdef=options['attributes']
            attributes=[x for x in re.split(r'\W+',attdef) if len(x) > 0]

        wkt='wkt' in options

        with open(csvfile,"wb") as csvf:
            csvw=csv.writer(csvf)
            cols=['fromstn','tostn','fromhgt','tohgt','obstype','obsset','value','error','calcval','residual','reserr','stdres']
            cols.extend(attributes)
            if wkt:
                cols.append('wkt')
            csvw.writerow(cols)
            lastset=0
            for obs,res,rescvr in self.residuals():
                if obs.obstype.nvalue > 1:
                    continue
                set=None
                if len(obs.obsvalues) > 0:
                    lastset += 1
                    set=lastset
                vformat="{0:.6f}"
                eformat="{0:.6f}"
                seformat="{0:.4f}"
                for obsv,resv,rescvr in zip(obs.obsvalues,res.flatten(),rescvr.diagonal().flatten()):
                    stderr=obsv.stderr
                    if rescvr >= stderr*stderr*1.0e-10:
                        reserr=math.sqrt(rescvr)
                    else:
                        reserr=None

                    data=[
                        obsv.inststn,
                        obsv.trgtstn,
                        obsv.insthgt,
                        obsv.trgthgt,
                        obs.obstype.code,
                        set,
                        vformat.format(obsv.value),
                        eformat.format(stderr),
                        vformat.format(obsv.value-resv),
                        vformat.format(resv),
                        eformat.format(reserr) if reserr is not None else None,
                        seformat.format(abs(resv/reserr)) if reserr is not None else None
                        ]
                    data.extend([obsv.attributes.get(a) for a in attributes])
                    if wkt:
                        obswkt=''
                        try:
                            llh0=self.stations.get(obsv.inststn).llh()
                            if obsv.trgtstn is not None:
                                llh1=self.stations.get(obsv.trgtstn).llh()
                                obswkt='LINESTRING Z({0:.9f} {1:.9f} {2:.4f},{3:.9f} {4:.9f} {5:.4f})'.format(
                                    llh0[0],llh0[1],llh0[2],llh1[0],llh1[1],llh1[2])
                            else:
                                obswkt='POINT Z({0:.9f} {1:.9f} {2:.4f})'.format(
                                    llh0[0],llh0[1],llh0[2])
                        except:
                            obswkt=None
                        data.append(obswkt)
                    csvw.writerow(data)
        
    def _getExtraStationDataFunc( self ):
        writeCovar=self.options.outputStationCovariances
        writeEllipse=self.options.outputStationErrorEllipses
        writeOffsets=self.options.outputStationOffsets
        if not writeCovar and not writeEllipse and not writeOffsets:
            return None
        columns=[]
        if writeCovar:
            columns.extend(('stderr_e','stderr_n','stderr_u','corr_en','corr_eu','corr_nu'))
        if writeEllipse:
            columns.extend(('maxerr_h','minerr_h','azmaxerr'))
            if not writeCovar:
                columns.extend(('stderr_u',))
        if writeOffsets:
            columns.extend(('offset_e','offset_n','offset_u'))

        def func(code):
            if code is None: 
                return columns
            data={}
            if writeCovar or writeEllipse:
                cvrenu=self.stationENUCovariance( code )
                stderr=np.sqrt(cvrenu.diagonal())
                div=stderr*stderr.reshape((3,1))
                data['stderr_e']=stderr[0]
                data['stderr_n']=stderr[1]
                data['stderr_u']=stderr[2]
                data['corr_en']=cvrenu[0,1]/div[0,1] if div[0,1] > 0.0 else 0.0
                data['corr_eu']=cvrenu[0,2]/div[0,2] if div[0,2] > 0.0 else 0.0
                data['corr_nu']=cvrenu[1,2]/div[1,2] if div[1,2] > 0.0 else 0.0
                if writeEllipse:
                    v1=(cvrenu[1,1]+cvrenu[0,0])/2.0
                    v2=(cvrenu[1,1]-cvrenu[0,0])/2.0
                    v3=cvrenu[0,1]
                    v4=math.sqrt(v2*v2+v3*v3)
                    azmax=math.degrees(math.atan2(v3,v2)/2.0) if v4 > 0.0 else 0.0
                    if azmax < 0.0:
                        azmax += 180.0
                    emax=math.sqrt(max(0.0,v1+v4))
                    emin=math.sqrt(max(0.0,v1-v4))
                    data['maxerr_h']=emax
                    data['minerr_h']=emin
                    data['azmaxerr']=azmax
            if writeOffsets:
                if code in self.originalXyz:
                    stn=self.stations.get(code)
                    offset=stn.xyz()-self.originalXyz[code]
                    offset=stn.enu().dot(offset)
                    data['offset_e']=offset[0]
                    data['offset_n']=offset[1]
                    data['offset_u']=offset[2]
            for c in data:
                data[c]="{0:.5f}".format(data[c])

            return [data.get(c) for c in columns]

        return func

    def _getCoordinateCovariance( self, code ):
        if code is None:
            return 
        cvrenu=self.stationENUCovariance( code )
        stderr=np.sqrt(cvrenu.diagonal())
        div=stderr*stderr.reshape((3,1))
        return ["{0:.4f}".format(x) for x in 
                (stderr[0],stderr[1],stderr[2],
                cvrenu[0,1]/div[0,1] if div[0,1] > 0.0 else 0.0,
                cvrenu[0,2]/div[0,2] if div[0,2] > 0.0 else 0.0,
                cvrenu[1,2]/div[1,2] if div[1,2] > 0.0 else 0.0)]

    def writeOutputFiles( self ):
        if self.options.outputResidualFile is not None:
            self.writeResidualCsv(self.options.outputResidualFile,self.options.outputResidualFileOptions)

        if self.options.outputStationFile is not None:
            self.stations.writeCsv(self.options.outputStationFile,
                                   geodetic=self.options.outputStationFileGeodetic,
                                   extradata=self._getExtraStationDataFunc())


    def setup( self ):
        pass

    def ignoreMissingStations( self ):
        # Deal with missing stations
        options=self.options
        if options.ignoreMissingStations:
            missing=self.missingStationList()
            if len(missing) > 0:
                self.write("Ignoring missing stations:\n")
                for use,stn in missing:
                    self.write("  {0}\n".format(stn))
                self.filterObsByStation(missing)
        if options.acceptStations:
            self.filterObsByStation(options.acceptStations)

    def calculateSolution( self ):
        options=self.options

        self.setupNormalEquations()
        self.write("\nCalculating {0} parameters\n".format(self.nparam))

        if self.options.debugObservationEquations:
           self.writeObservationEquations()

        converged=True
        for i in range(options.maxIterations):
            converged,coordUpdate,code=self.runOneIteration()
            self.write("Iteration {0}: max coord change {1:.4f}m at {2}\n".format(i+1,coordUpdate,code))
            if converged and i >= options.minIterations-1 :
                break

        if not converged:
            raise ConvergenceError('Adjustment failed to converge')

    def report(self):
        self.writeSummaryStats()
        self.writeResidualSummary()

    def runSetup( self ):
        self.loadDataFiles()
        if self.options.verbose:
            self.writeObservationSummary()
        self.runPluginFunction('setup',runPrePostFunctions=True)
        self.ignoreMissingStations()
        self.setupNormalEquations()

    def runCalculateSolution(self):
        self.runPluginFunction('calculateSolution',firstOnly=True,reverse=True)

    def runOutputs(self):
        self.runPluginFunction('report')
        self.runPluginFunction('writeOutputFiles')

    def run( self ):
        self.runSetup()
        try:
            self.runCalculateSolution()
        finally:
            self.runOutputs()

    @staticmethod
    def main(plugins=None):
        '''
        Main function to run the adjustment

        Optionally can take a list of plugin classes 
        '''
        import argparse
        import os.path
        parser=argparse.ArgumentParser(description="Adjust network")
        parser.add_argument('config_file',help="Adjustment definition file")
        parser.add_argument('output_file',nargs='?',help="Output file - default is standard output")
        parser.add_argument('-r','--residual_file',help="Residual CSV file")
        parser.add_argument('-s','--output-coordinate-file',help="Output station CSV file")
        parser.add_argument('-v','--verbose',action='store_true',help="Verbose output")
        parser.add_argument('-c','--create',action='store_true',help="Create an example adjustment configuration file")
        args=parser.parse_args()

        # Crude configuration file reading, parsed to a dictionary
        # Multiple options compiled to '\n' separated string
        # Blank options have value 'y'

        plugins=plugins or []
        if args.create:
            if os.path.exists(args.config_file):
                print("Config file "+args.config_file+" already exists - not creating!")
                sys.exit()
            import inspect
            classes=[Adjustment]
            classes.extend(plugins)
            with open(args.config_file,'w') as cf:
                spacer=''
                for cls in classes:
                    source=inspect.getfile(cls)
                    source=os.path.splitext(source)[0]+'.adj'
                    if not os.path.exists(source):
                        if cls==Adjustment:
                            print("Cannot find example configuration file :-(")
                            sys.exit()
                        continue
                    with open(source) as sf:
                        cf.write(spacer)
                        cf.write(sf.read())
                    spacer="\n"
            print("Example configuration file "+args.config_file+" created.")
            sys.exit()

        # Set up and run the adjustment

        adj=Adjustment(plugins=plugins,config_file=args.config_file,output_file=args.output_file,verbose=args.verbose)
        adj.run()

if __name__=="__main__":
    Adjustment.main()
