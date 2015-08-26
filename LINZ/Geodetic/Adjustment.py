#!/usr/bin/python
# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import math
from collections import namedtuple
import re
import numpy as np
from numpy import linalg

from .Network import Network
from .Observation import Observation

'''
Module to adjust network coordinates.  
'''

# Error classes used by the adjustment

class ConvergenceError( RuntimeError ): pass
class SingularityError( RuntimeError ): pass

# Options controlling an adjustment

class Options( object  ):
    '''
    Class to define the options used in an adjustment

    Subclasses can override the __init__ and loadConfig function to read 
    additional information from a configuration.  Note that loadConfig may
    be called more than once with different configuration settings, so items
    should be initiallized in __init__ and only updated in loadConfig.
    '''

    def __init__(self,config_file=None):
        # File names
        self.listingFile=None
        self.stationFile=None
        self.dataFiles=[]
        self.outputStationFile=None
        self.outputResidualFile=None
        # Station options
        self.fixedStations=[]
        self.acceptStations=[]
        self.rejectStations=[]
        self.reweightObsType={}
        self.ignoreMissingStations=False
        self.calcMissingCoords=False
        self.localGeoidModel=None
        # Adjustment options
        self.maxIterations=10
        self.convergenceTolerance=0.0001
        self.adjustENU=False
        self.refractionCoefficient=0.075
        # Output options
        self.verbose=False
        self.debugStationOffsets=False

        # Specific outputs - only apply if verbose is True
        self.debugStationOffsets=False
        self.debugObservationEquations=False
        self.debugCalcMissingCoords=False
        if config_file is not None:
            self.loadConfigFile( config_file )

    def loadConfigFile( self, config_file ):
        cfg={}
        with open(config_file) as cfgf:
            for l in cfgf:
                l=l.strip()
                if len(l)==0 or l.startswith('#'):
                    continue
                try:
                    parts=l.split(None,1)
                    item=parts[0].lower()
                    value=parts[1] if len(parts)==2 else 'y'
                except:
                    print("Invalid configuration line: ",l)
                    continue
                item=item.lower()
                cfg[item]=value if item not in cfg else cfg[item]+'\n'+value
        self.loadConfig( cfg )

    def getbool(self,config,param,current): 
        if param not in config:
            return current
        return 'yes'.startswith((config.get(param) or 'yes').lower())

    def loadConfig( self, cfg ):
        # File names
        self.listingFile=cfg.get('listing_file',None) or self.listingFile
        self.stationFile=cfg.get('coordinate_file',None) or self.stationFile
        self.dataFiles.extend((df for df in cfg.get('data_file','').split("\n") if len(df) > 0))
        self.outputStationFile=cfg.get('output_coordinate_file',None) or self.outputStationFile
        self.outputResidualFile=cfg.get('residual_csv_file',None) or self.outputResidualFile
        # Station selection

        self.addFixedStations(cfg.get('fix','').split())
        self.addAcceptStations(cfg.get('accept','').split())
        self.addRejectStations(cfg.get('reject','').split())
        self.ignoreMissingStations=self.getbool(cfg,'ignore_missing_stations',self.ignoreMissingStations)
        self.calcMissingCoords=self.getbool(cfg,'calculate_missing_stations',self.calcMissingCoords)
        geoidModel=cfg.get('local_geoid','')
        if geoidModel != '':
            try:
                gparts=geoidModel.split()
                code=gparts.pop(0)
                geoidHeight,xi,eta=(float(v) for v in gparts[0:3])
                range=float(gparts[4]) if len(gparts)==5 else None
            except:
                raise RuntimeError("Invalid local_geoid: "+geoidModel)
            self.localGeoidModel={'code':code,'height':geoidHeight,'xi':xi,'eta':eta, 'range':range}

        # Observation options
        for reweight in cfg.get('reweight_observation_type','').split('\n'):
            reweight=reweight.strip()
            if reweight == '':
                continue
            match=re.match(r'([a-z]{2})\s+(\d+(?:\.\d+)?)$',reweight,re.I)
            if not match:
                raise RuntimeError('Invalid reweight_observation_type option: '+reweight)
            typecode=match.group(1).upper()
            factor=float(match.group(2))
            self.reweightObsType[typecode]=factor

        # Adjustment options
        self.maxIterations=10
        self.convergenceTolerance=float(cfg.get('convergence_tolerance',self.convergenceTolerance))
        self.maxIterations=int(cfg.get('max_iterations',self.maxIterations))
        self.adjustENU=self.getbool(cfg,'adjust_enu',self.adjustENU)
        self.refractionCoefficient=float(cfg.get('refraction_coefficient',self.refractionCoefficient))
        # Output options
        self.verbose=self.getbool(cfg,'verbose',self.verbose)
        # Debug options
        self.debugObservationEquations=self.getbool(cfg,'debug_observation_equations',self.debugObservationEquations)
        self.debugStationOffsets=self.getbool(cfg,'debug_station_offsets',self.debugStationOffsets)
        self.debugCalcMissingCoords=self.getbool(cfg,'debug_calculate_missing_stations',self.debugCalcMissingCoords)


    def addFixedStations( self, codelist ):
        '''
        Add a list of additional stations to the fixed stations list
        '''
        self.fixedStations.extend(codelist)

    def addAcceptStations( self, codelist ):
        '''
        Add a list of additional stations to use in the adjustment 
        (other stations will be rejected)
        '''
        self.acceptStations.extend(codelist)

    def addRejectStations( self, codelist ):
        '''
        Add a list of additional stations to be rejected from the adjustment
        '''
        self.rejectStations.extend(codelist)

class ObsEq( object ):

    def __init__( self, obsres, obseq, obscovar, schreiber=None ):
        self.obsres=obsres
        self.obseq=obseq
        self.obscovar=obscovar
        self.diagonal=obscovar.shape[1] == 1
        self.schreiber=schreiber

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
        setup
        run
        report
        writeOutputFiles

    Functions to include additional parameters.  setupParameters should call
    addParameter and addParameterUpdate. observationEquation should call Adjustment.observationEquation
    and then update the resultant equations (may not be necessary if parameters are
    applied indirectly via the coordParamMapping.

        setupParameters
        observationEquation
        updateParameters

    '''

    def __init__( self, 
                 stations=Network(), 
                 observations=[], 
                 options=Options(), 
                 config_file=None,
                 config=None,
                 output_file=None ):

        self.stations=stations
        self.observations=observations

        self.options=options
        if config_file is not None:
            options.loadConfigFile( config_file )
        if config is not None:
            options.loadConfig( config )
        
        if output_file is None:
            output_file = options.listingFile
        if isinstance(output_file,basestring):
            output_file=open(output_file,"w")

        self.output=output_file
        self.parameters=[]
        self.solved=0

    def write( self, message, debug=False ):
        if self.output is not None:
            self.output.write( message )
        if self.options.verbose and not debug: 
            sys.stdout.write(message)

    def writeDebugOutput( self, message ):
        self.write( message, True )

    def loadDataFiles( self ):

        # Load the coordinate file
        if self.options.stationFile is not None:
            self.write("\nInput coordinate file:\n")
            self.write("  {0}\n".format(self.options.stationFile))
            self.stations.loadCsv(self.options.stationFile)

        # Load the data files
        reweight=self.options.reweightObsType
        if len(reweight) > 0:
            self.write("\nReweighting input observations:\n")
            for otype in sorted(reweight.keys()):
                self.write("  Observation type {0} errors scaled by {1:6.3f}\n"
                           .format(otype,reweight[otype]))

        self.write("\nObservation files:\n")
        for obsfile in self.options.dataFiles:
            self.write("  {0}\n".format(obsfile))
            if obsfile.lower().endswith('.msr'):
                from . import MsrFile
                reader=MsrFile.read
            elif re.search(r'\.snx(\.gz)?$',obsfile,re.I):
                from . import SinexObsFile
                reader=SinexObsFile.read
            else:
                from . import CsvObsFile
                reader=CsvObsFile.read
            for obs in reader(obsfile):
                typecode=obs.obstype.code
                if typecode in reweight:
                    factor=reweight[typecode]
                    for ov in obs.obsvalues:
                        if ov.stderr is not None:
                            ov.stderr *= factor
                    if obs.covariance is not None:
                        obs.covariance *= factor*factor
                self.observations.append(obs)

    def usedStations( self ):
        '''
        Get dictionary with keys matching the stations used in the list
        of observations.
        '''
        used={}
        for o in self.observations:
            for obsval in o.obsvalues:
                used[obsval.inststn]=1
                used[obsval.trgtstn]=1
        used={code:1 for code in used.keys() if self.stations.get(code) is not None}
        return used

    def missingStationList( self, clean=False ):
        missing=set()
        good_station=lambda x: x is None or self.stations.get(x) is not None
        for obs in self.observations:
            for o in obs.obsvalues:
                if not good_station(o.inststn):
                    missing.add(o.inststn)
                if not good_station(o.trgtstn):
                    missing.add(o.trgtstn)
        if clean:
            self.filterObsByStation(remove=missing)
        return missing

    def filterObsByStation( self, keep=None, remove=None ):
        '''
        Filter observations by station, can either specify stations to
        keep or stations to remove, or both (though that doesn't make 
        much sense!)
        '''
        kfunc=None
        rfunc=None
        if keep is not None and len(keep) > 0:
            kfunc=lambda obs,ov: ov.inststn in keep and ov.trgtstn in keep
        if remove is not None and len(remove) > 0:
            rfunc=lambda obs,ov: not (ov.inststn in remove or ov.trgtstn in remove)
        func=(rfunc if kfunc is None else
              kfunc if rfunc is None else
              lambda obs,ov: rfunc(obs,ov) and kfunc(obs,ov))
        if func is not None:
            self.filterObservations(func)

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
                covariance=covariance[goodrows,goodrows]
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
        '''
        self.updateFuncs.append(updateFunc)

    def setupStationCoordMapping( self ):
        '''
        Defines a mapping from station coordinates to parameters.

        Sets up a dictionary from station id to a tuple of three values
            paramno: parameter number array - parameters used to define stations
            dxyzdp: (n,3) array defining differential of x,y,z against param
            updatexyz: if true then the station coordinate is updated based on
                the parameters in updateStationCoordParameters

        Default implementation simply adjusts x,y,z coord for each station.

        The routines setStationCoordMapping and updateStationCoordParameters
        can be overridden to parameterize station coordinates in terms of other
        parameters.

        '''

        # Get a list of used stations - remove fixed stations from it
        # Result is a list of stations to adjust

        usedStations=self.usedStations()
        for code in self.options.fixedStations:
            usedStations.pop(code,None)

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
            if mapping[2]:
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
        for update in self.updateFuncs:
            update(paramValues)
        return self.updateStationCoordParameters( paramValues )


    def setupParameters( self ):
        '''
        Function to set up adjustment parameters.  
        Parameters are added with addParameter(paramname), which returns
        the parameter row number in the adjustment.
        Also sets up the coordParamMapping dictionary as described in
        setupStationCoordMapping()

        '''
        return self.setupStationCoordMapping()

    def setupNormalEquations( self ):
        self.initParameters()
        self.setupParameters()
        nprm=self.nparam
        self.solved=0
        self.N=np.zeros((nprm,nprm))
        self.b=np.zeros((nprm,1))
        self.ssr=0.0
        self.nobs=0
        return nprm

    def observationEquation( self, obs ):
        '''
        Forms the observation equations for an observation o, which may
        be an array of observations
        '''
        refcoef=self.options.refractionCoefficient
        obstype=obs.obstype
        nval=obstype.nvalue
        nrow=len(obs.obsvalues)*nval
        obseq=np.zeros((nrow,self.nparam))
        obsres=np.zeros((nrow,1))
        obscovar=obs.covariance
        diagonal=obscovar is None
        schreiber=None
        if diagonal:
            obscovar=np.zeros((nrow,1))
        for i,o in enumerate(obs.obsvalues):
            stf=self.stations.get(o.inststn)
            if stf is None:
                raise RuntimeError("Station "+o.inststn+" is not defined")
            stt=None
            if o.trgtstn is not None:
                stt=self.stations.get(o.trgtstn)
                if stt is None:
                    raise RuntimeError("Station "+o.trgtstn+" is not defined")
            calcval,ddxyz0,ddxyz1=obstype.calcobs(stf,trgtstn=stt,insthgt=o.insthgt,trgthgt=o.trgthgt,refcoef=refcoef,ddxyz=True)
            i0=i*nval
            i1=i0+nval
            if nval == 1:
                obsres[i,0]=o.value-calcval
                if diagonal:
                    obscovar[i,0]=o.stderr*o.stderr
            else:
                obsres[i0:i1,0]=np.array(o.value)-calcval
            mapping=self.coordParamMapping.get(o.inststn,None)
            if mapping is not None:
                obseq[i0:i1,mapping[0]]=mapping[1].dot(ddxyz0)
            mapping=self.coordParamMapping.get(o.trgtstn,None)
            if mapping is not None:
                obseq[i0:i1,mapping[0]]+=mapping[1].dot(ddxyz1)

        # Need special handling for horizontal angles

        if obstype.code == 'HA':
            assert diagonal
            # Round to multiple of 180.
            obsres = np.remainder(obsres-obsres[0]+180.0,360.0)-180.0
            # Remove weighted mean residual
            wgt=1.0/obscovar
            wgtsum=wgt.sum()
            obsres -= (wgt.T.dot(obsres))/wgtsum
            schreiber=np.ones((nrow,1))

        # Handle 360 degree wrapping for azimuths

        elif obstype.code == 'AZ':
            obsres = np.remainder(obsres+180.0,360.0)-180.0

        return ObsEq(obsres,obseq,obscovar,schreiber)

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

    def solveEquations( self ):
        try:
            self.x=linalg.solve(self.N, self.b).flatten()
            self.solved=1
        except linalg.LinAlgError:
            self.solved=-1
            raise SingularityError()

    def runOneIteration( self ):
        self.setupNormalEquations()
        for oe in self.observationEquations():
            self.sumObservation(oe)
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
                summary[0] += res.T.dot(cvrinv.dot(res))
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

    def writeResidualCsv( self, csvfile ):
        '''
        Dumps residuals of scalar observations to a CSV file
        '''
        import csv
        with open(csvfile,"wb") as csvf:
            csvw=csv.writer(csvf)
            csvw.writerow(('fromstn','tostn','fromhgt','tohgt','obstype','obsset','value','error','calcval','residual','reserr','stdres'))
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
                    csvw.writerow((
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
                    ))
        
    def writeOutputFiles( self ):
        if self.options.outputResidualFile is not None:
            self.writeResidualCsv(self.options.outputResidualFile)

        if self.options.outputStationFile is not None:
            self.stations.writeCsv(self.options.outputStationFile)
        
    def setupLocalGeoid( self, geoidModel ):
        refStationCode=geoidModel['code']
        refStation=self.stations.get( geoidModel['code'] )
        if refStation is None:
            raise RuntimeError("Local geoid model reference station {0} is not defined"
                               .format(refStationCode))
        gheight=geoidModel['height']
        xi=geoidModel['xi']/3600.0
        eta=geoidModel['eta']/3600.0
        xirad=math.radians(xi)
        etarad=math.radians(eta)
        range=geoidModel.get('range')
        if range is not None:
            range = range*range

        xyz0=refStation.xyz()
        enu_axes=refStation.enu()

        for s in self.network.stations:
            denu=enu.dot(s.xyz()-xyz0)
            if range is not None and (denu[0]*denu[0]+denu[1]*denu[1]) > range:
                continue
            s.setXiEta([xi,eta])

    def setup( self ):
        options=self.options
        self.loadDataFiles()

        # Calculate missing stations
        if options.calcMissingCoords:
            from . import StationLocator
            debug=options.debugCalcMissingCoords
            write=None
            if debug:
                write=self.write
                write("\nCalculating missing station coordinates\n")
            nupdated=StationLocator.locateStations(self.stations,self.observations,write)
            if nupdated > 0:
                self.write("\nApproximate coordinates calculated for {0} stations\n"
                           .format(nupdated))

        # Apply a geoid model
        if options.localGeoidModel is not None:
            gm=options.localGeoidModel
            code=gm['code']
            geoidHeight=gm['height']
            xieta=[gm['xi']/3600.0,gm['eta']/3600.0]
            range=gm['range']
            self.stations.setLocalGeoid(code,geoidHeight,xieta,range)

        # Deal with missing stations
        if options.ignoreMissingStations:
            missing=self.missingStationList()
            if len(missing) > 0:
                self.write("Ignoring missing stations:\n")
                for stn in missing:
                    self.write("  {0}\n".format(stn))
                self.filterObsByStation(remove=missing)
        self.filterObsByStation(keep=options.acceptStations,remove=options.rejectStations)

        # Set up the normal equations
        self.setupNormalEquations()

    def run( self ):
        options=self.options

        self.setupParameters()
        self.write("\nCalculating {0} parameters\n".format(self.nparam))

        if self.options.debugObservationEquations:
           self.writeObservationEquations()

        converged=False
        for i in range(options.maxIterations):
            coordUpdate,code=self.runOneIteration()
            self.write("Iteration {0}: max coord change {1:.4f}m at {2}\n".format(i,coordUpdate,code))
            if coordUpdate < options.convergenceTolerance:
                converged=True
                break

        if not converged:
            raise ConvergenceError('Adjustment failed to converge')

    def report(self):
        self.writeSummaryStats()
        self.writeResidualSummary()

def main(adjustment_class=Adjustment):
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

    if args.create:
        if os.path.exists(args.config_file):
            print("Config file "+args.config_file+" already exists - not creating!")
            sys.exit()
        import inspect
        source=inspect.getfile(adjustment_class)
        source=os.path.splitext(source)[0]+'.adj'
        if not os.path.exists(source):
            print("Cannot find example configuration file :-(")
            sys.exit()
        with open(source) as sf, open(args.config_file,'w') as cf:
            cf.write(sf.read())
            print("Example configuration file "+args.config_file+" created.")
            sys.exit()

    cfg={}
    if args.verbose:
        cfg['verbose']='yes'

    # Set up and run the adjustment

    adj=adjustment_class(config_file=args.config_file,output_file=args.output_file,config=cfg)
    adj.setup()
    adj.writeObservationSummary()
    adj.run()
    adj.report()
    adj.writeOutputFiles()

if __name__=="__main__":
    main(Adjustment)
