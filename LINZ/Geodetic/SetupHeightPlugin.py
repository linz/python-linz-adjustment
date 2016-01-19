# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import re
import math
import numpy as np

from .Adjustment import Plugin, Options
from .Station import Station

class SetupHeightPlugin( Plugin ):
    '''
    Plugin to set/calculate instrument setup height offset based on 
    observation attribute
    '''
        
    def pluginOptions(self):
        self.setupHeights={}
        return dict(
            calculateSetupHeights=False,
            calculateSetupAttributes=['inst_setup','trgt_setup'],
            validSetupRegex=[],
            invalidSetupRegex=[],
            fixSetupValue={},
            )

    def setConfigOption( self, item, value ):
        if item == 'calculate_setup_heights':
            self.options.calculateSetupHeights=Options.boolOption(value)
        elif item == 'inst_trgt_setup_attributes':
            parts=value.split()
            if len(parts) != 2:
                raise RuntimeError("inst_trgt_setup_attributes must define two attributes")
            self.options.calculateSetupAttributes=parts
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
            m=re.match(r'^(\w+)\s+([-+]?\d+(?:\.\d+)?)$',value)
            if m:
                self.options.fixSetupValue[m.group(1)]=float(m.group(2))
            else:
                raise RuntimeError("Invalid fix_setup_height value "+value)
        else:
            return False
        return True

    def setupParameters( self ):
        if not self.options.calculateSetupHeights:
            return
        setupAttributes=self.options.calculateSetupAttributes
        oldsetups=self.setupHeights
        self.setupHeights={}
        setups={}
        for o in self.adjustment.observations:
            for v in o.obsvalues:
                for attr in setupAttributes:
                    setup=v.attributes.get(attr,None)
                    if setup is not None:
                        if setup not in setups:
                            setups[setup] = 1
                        else:
                            setups[setup] += 1
        fixed=self.options.fixSetupValue
        haveParams=False
        for setup in sorted(setups):
            valid=True
            if setup not in fixed:
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
            value=oldsetups[setup]['value'] if setup in oldsetups else 0.0
            if setup in fixed:
                value=fixed[setup]
            else:
                paramno=self.adjustment.addParameter('Setup '+setup+' height offset')
                haveParams=True
            self.setupHeights[setup]= {
                 'count': setups[setup],
                 'value': value,
                 'paramno': paramno
                }

        if haveParams:
            self.adjustment.addParameterUpdate( self.updateSetupParams )

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
        setupAttributes=self.options.calculateSetupAttributes
        for v in obs.obsvalues:
            sevalue=[Station.OFFSET_H,0.0,0.0,None,None]
            offsets.append(sevalue)
            for i,attr in enumerate(setupAttributes):
                setup=v.attributes.get(attr,None)
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
                try:
                    stderr=math.sqrt(covar[paramno,paramno])
                    write("\n{0:<10s} {1:8.4f} {2:8.4f}\n".format
                        (s,value,stderr))
                except ValueError:
                    stderr=covar[paramno,paramno]
                    write("\n{0:<10s} {1:8.4f} ({2:7.4f})\n".format
                        (s,value,stderr))
