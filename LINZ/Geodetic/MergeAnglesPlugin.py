# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from collections import namedtuple
import numpy as np
import numpy.linalg as la

from .Adjustment import Plugin
from .Observation import Observation, ObservationValue
from . import StationLocator

_haset=namedtuple('set','obs targets')

def _meanangle( angles ):
    '''
    Derive a mean angle accounting for possible 360 degree offsets
    '''
    if len(angles) == 1:
        return angles[0]
    a0=angles[0]-180
    return np.mean(np.mod(np.array(angles)-a0,360.0))+a0

def _reorient( trialset, mergeset ):
    '''
    Approximately reorients trial set to be aligned with mergeset
    '''
    common={}
    for obs in trialset:
        for o in obs.obsvalues:
            if o.trgtstn not in common:
                common[o.trgtstn]=([],[])
            common[o.trgtstn][0].append(o.value)
    for obs in mergeset:
        for o in obs.obsvalues:
            if o.trgtstn in common:
                common[o.trgtstn][1].append(o.value)
    offset=_meanangle([
        _meanangle(c[0])-_meanangle(c[1])
        for c in common.values()
        if len(c[1]) > 0
        ])
    for obs in trialset:
        for o in obs.obsvalues:
            o.value=np.mod(o.value-offset,360.0)

def _mergedObs( obslist ):
    '''
    Merges an approximately oriented set of angles

    Constructs and solves a least squares problem for the direction to
    each target and the orientation of each set, 
    '''

    nset=len(obslist)
    if nset < 2:
        return obslist[0]
    targets=set()
    nobs=0
    for obs in obslist:
        for o in obs.obsvalues:
            targets.add(o.trgtstn)
            nobs += 1
    ntarget=len(targets)
    targetid={t:i for i,t in enumerate(sorted(targets))}

    # Form normal equations.  Weighted by stderr.
    # Each obs consists of direction + set orientation
    A=np.zeros((nobs+1,ntarget+nset))
    b=np.zeros((nobs+1))
    nrow=0
    sumwgt=0.0
    for iobs,obs in enumerate(obslist):
        iset=iobs+ntarget
        for o in obs.obsvalues:
            wgt=1.0/o.stderr
            A[nrow,iset]=wgt
            A[nrow,targetid[o.trgtstn]]=wgt
            b[nrow]=o.value*wgt
            sumwgt += wgt
            nrow += 1

    # Add constraint on mean set orientation
    A[nrow,ntarget:]=sumwgt*nrow

    # Solve normal equations
    N=A.T.dot(A)
    x=A.T.dot(b)
    N=la.inv(N)
    x=N.dot(x)

    x=np.mod(x-x[0],360.0)

    # Form the observations
    obs=Observation('HA')
    inststn=obslist[0].obsvalues[0].inststn

    for t in sorted(targets):
        i=targetid[t]
        value=x[i]
        stderr=np.sqrt(N[i,i])
        obs.addObservation(ObservationValue(inststn,t,value,stderr))

    return obs

def _mergedSets( haobs ):
    # First merge intersecting sets
    hasets=[_haset([o],set(v.trgtstn for v in o.obsvalues)) for o in haobs]
    while hasets:
        trialset=hasets.pop()
        mergeset=next((ha for ha in hasets if not ha.targets.isdisjoint(trialset.targets)),None)
        if mergeset is not None:
            _reorient(trialset.obs,mergeset.obs)
            mergeset.obs.extend(trialset.obs)
            mergeset.targets.update(trialset.targets)
        else:
            yield _mergedObs(trialset.obs)

def _mergeAngleObservations( observations ):
    newobs=[]
    haobs={}
    for obs in observations:
        if obs.obstype.code != 'HA':
            newobs.append(obs)
        else:
            inststn=obs.obsvalues[0].inststn
            if inststn not in haobs:
                haobs[inststn]=[]
            haobs[inststn].append(obs)

    for inststn in sorted(haobs.keys()):
        newobs.extend(_mergedSets(haobs[inststn]))
    return newobs

class MergeAnglesPlugin( Plugin ):
    '''
    Adjustment plugin class preprocesses observations to combine sets of horizontal angles
    with common targets.  Not really useful except when filtering data to not include 
    common targets - still allows data to be used.  Not strictly equivalent statistically
    to proper calculation of targets, but close!
    '''

    def preSetup( self ):
        self.adjustment.observations=_mergeAngleObservations(self.adjustment.observations)


