# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from collections import namedtuple
import re
import math
import numpy as np

from .Adjustment import Plugin, Options, ObsEqn
from .Station import Station

class SetupHeightPlugin( Plugin ):
    '''
    Plugin to set/calculate instrument setup height offset based on 
    observation attribute
    '''
    SetupHeightConstraint=namedtuple('SetupValue','code codere value error')
        
    def pluginOptions(self):
        self.setupHeights={}
        self.setupFuncs=None
        return dict(
            calculateSetupHeights=False,
            calculateSetupAttributes=['inst_setup','trgt_setup'],
            validSetupRegex=[],
            invalidSetupRegex=[],
            fixSetupValue={},
            floatSetupValue={},
            )

    def _setupFunc( self, attr ):
        if attr == 'inststn':
            return lambda v: v.inststn
        elif attr == 'trgtstn':
            return lambda v: v.trgtstn
        else:
            return lambda v: v.attributes.get(attr)

    def setConfigOption( self, item, value ):
        if item == 'calculate_setup_heights':
            self.options.calculateSetupHeights=Options.boolOption(value)
        elif item == 'inst_trgt_setup_attributes':
            parts=value.split()
            if len(parts) != 2:
                raise RuntimeError("inst_trgt_setup_attributes must define two attributes")
            self.options.calculateSetupAttributes=parts
            self.setupFuncs=tuple(self._setupFunc(p) for p in parts)
        elif item == 'valid_setup_regex':
            try:
                self.options.validSetupRegex.append(re.compile(value))
            except:
                raise RuntimeError("Invalid valid_setup_regex "+value)
        elif item == 'invalid_setup_regex':
            try:
                self.options.invalidSetupRegex.append(re.compile(value))
            except:
                raise RuntimeError("Invalid invalid_setup_regex "+value)
        elif item == 'fix_setup_height':
            m=re.match(r'^(\S+)\s+([-+]?\d+(?:\.\d+)?)$',value)
            if m:
                code=m.group(1)
                codere=None
                if code.startswith('re:'):
                    try:
                        codere=re.compile(code[3:],re.I)
                    except:
                        raise RuntimeError("Invalid regular expression in fix_setup_height")
                value=float(m.group(2))
                constraint=self.SetupHeightConstraint(code,codere,value,None)
                self.options.fixSetupValue[code]=constraint
            else:
                raise RuntimeError("Invalid fix_setup_height value "+value)
        elif item == 'float_setup_height':
            m=re.match(r'^(\S+)\s+([-+]?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)$',value)
            if m:
                code=m.group(1)
                codere=None
                if code.startswith('re:'):
                    try:
                        codere=re.compile(code[3:],re.I)
                    except:
                        raise RuntimeError("Invalid regular expression in float_setup_height")
                value=float(m.group(2))
                error=float(m.group(3))
                constraint=self.SetupHeightConstraint(code,codere,value,error)
                self.options.floatSetupValue[code]=constraint
            else:
                raise RuntimeError("Invalid fix_setup_height value "+value)
        else:
            return False
        return True

    def setupParameters( self ):
        if not self.options.calculateSetupHeights:
            return
        if self.setupFuncs is None:
            return
        setupFuncs = self.setupFuncs
        oldsetups=self.setupHeights
        self.setupHeights={}
        setups={}
        for o in self.adjustment.observations:
            for v in o.obsvalues:
                for f in setupFuncs:
                    setup=f(v)
                    if setup is not None:
                        if setup not in setups:
                            setups[setup] = 1
                        else:
                            setups[setup] += 1
        fixSetups=self.options.fixSetupValue
        floatSetups=self.options.floatSetupValue
        haveParams=False
        
        for setup in sorted(setups):
            valid=True
            if setup not in fixSetups and setup not in floatSetups:
                if len(self.options.validSetupRegex) > 0:
                    valid=False
                    for vre in self.options.validSetupRegex:
                        if vre.match(setup):
                            valid=True
                            break
                if valid and len(self.options.invalidSetupRegex) > 0:
                    for vre in self.options.invalidSetupRegex:
                        if vre.match(setup):
                            valid=False
                            break

                if not valid:
                    continue
            paramno=-1
            isfixed=False
            error=None
            float=None
            value=0.0
            if setup in fixSetups:
                value=fixSetups[setup].value
                isfixed=True
            else:
                for constraint in fixSetups.values():
                    if constraint.codere is None:
                        continue
                    if constraint.codere.match(setup):
                        value=constraint.value
                        isfixed=True
                        break
            if not isfixed:
                paramno=self.adjustment.addParameter('Setup '+setup+' height offset')
                haveParams=True
                if setup in floatSetups:
                    value=floatSetups[setup].value
                    float=value
                    error=floatSetups[setup].error
                else:
                    for constraint in floatSetups.values():
                        if constraint.codere is not None and constraint.codere.match(setup):
                            value=constraint.value
                            float=value
                            error=constraint.error
                            break
            if setup in oldsetups:
                value=oldsetups[setup]['value']
            self.setupHeights[setup]= {
                 'count': setups[setup],
                 'paramno': paramno,
                 'value': value,
                 'float': float,
                 'error': error,
                }

        if haveParams:
            self.adjustment.addParameterUpdate( self.updateSetupParams )

    def sumNormalEquations( self ):
        for height in self.setupHeights.values():
            if height['paramno'] >= 0 and height['error'] is not None:
                obseqn=ObsEqn(1,self.adjustment.nparam)
                obseqn.obseq[0,height['paramno']]=1.0
                obseqn.obsres[0]=height['float']-height['value']
                obseqn.obscovar[0]=height['error']**2
                self.adjustment.sumObservation(obseqn)

    def updateSetupParams( self, paramValues ):
        maxoffset=0.0
        for sev in self.setupHeights.values():
            paramno=sev['paramno']
            if paramno >= 0:
                offset=paramValues[paramno]
                sev['value'] += offset
                offset=abs(offset)
                if offset > maxoffset:
                    maxoffset=offset
        return maxoffset < self.options.convergenceTolerance

    def calcStationOffsets( self, obs ):
        offsets=[]
        havesetups=False
        setupFuncs=self.setupFuncs
        for v in obs.obsvalues:
            sevalue=[Station.OFFSET_H,0.0,0.0,None,None]
            offsets.append(sevalue)
            for i,f in enumerate(setupFuncs):
                setup=f(v)
                if setup in self.setupHeights:
                    havesetups=True
                    sedata=self.setupHeights[setup]
                    sevalue[i+1]=sedata['value']
                    paramno=sedata['paramno']
                    if paramno >= 0:
                        sevalue[i+3]=([paramno],np.array([1.0]))
        return offsets if havesetups else None

    def report( self ):
        if len(self.setupHeights) == 0:
            return
        write=self.adjustment.write
        write("\nCalculated setup heights:\n")
        write("\n{0:<10s} {1:>8s} {2:>8s}\n".format
              ("Setup","Height","Error"))

        covar=self.adjustment.covariance()
        for s in sorted(self.setupHeights):
            sev=self.setupHeights[s]
            paramno=sev['paramno']
            value=sev['value']
            if paramno < 0:
                write("\n{0:<10s} {1:8.4f} {2:>8s}\n".format
                      (s,value,"-    "))
            else:
                float=' (floated {0:.4f})'.format(sev['error']) if sev['error'] is not None else ''
                try:
                    stderr=math.sqrt(covar[paramno,paramno])
                    write("\n{0:<10s} {1:8.4f} {2:8.4f}{3}\n".format
                        (s,value,stderr,float))
                except ValueError:
                    stderr=covar[paramno,paramno]
                    write("\n{0:<10s} {1:8.4f} ({2:7.4f}{3})\n".format
                        (s,value,stderr,float))
