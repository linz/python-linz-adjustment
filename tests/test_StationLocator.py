# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os.path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'LINZ'))
import unittest
from Geodetic import Station, Network, Observation, StationLocator
from Geodetic.Observation import Observation as Obs, ObservationValue as ObsVal

def HA( stfrom, targets, offset=0.0, error=0.001 ):
    ha=Obs('HA')
    refha=None
    for stto in targets:
        angle=stfrom.azimuthTo(stto)
        if refha is None:
            refha=angle-offset
        angle=angle-refha
        if angle < 0:
            angle += 360.0
        ha.addObservation(ObsVal(stfrom.code(),stto.code(),angle,error))
    return ha

def AZ(stfrom,stto,error=0.001):
    angle=stfrom.azimuthTo(stto)
    return Obs('AZ',obsvalue=ObsVal(stfrom.code(),stto.code(),angle,error))

def SD(stfrom,stto,error=0.001):
    angle=stfrom.distanceTo(stto)
    return Obs('SD',obsvalue=ObsVal(stfrom.code(),stto.code(),angle,error))

def ZD(stfrom,stto,error=0.001):
    angle=stfrom.zenithDistanceTo(stto)
    return Obs('ZD',obsvalue=ObsVal(stfrom.code(),stto.code(),angle,error))

def LV(stfrom,stto,error=0.001):
    angle=stfrom.heightDifferenceTo(stto)
    return Obs('LV',obsvalue=ObsVal(stfrom.code(),stto.code(),angle,error))

class StationLocatorTestCase( unittest.TestCase ):


    def test_000_basic_offset( self ):
        '''
        Simple offset bearing, distance, zenith distance
        '''
        net=Network.Network()
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        net.addStation(st1)
        obs=[]
        obs.append(AZ(st1,st2))
        obs.append(SD(st1,st2))
        obs.append(ZD(st1,st2))
        nfixed=StationLocator.locateStations(net,obs,sys.stdout.write)
        st2c=net.get('ST2')
        self.assertTrue(nfixed==1,'Failed to add one station calculated ({0} fixed)'.format(nfixed))
        self.assertTrue(st2c is not None,'ST2 not added to network')
        offset=st2c.distanceTo(st2)
        self.assertTrue(offset < 1.0,'ST2 out of position  by {0:.1f} metres'.format(offset))

    def test_001_basic_offset_lv( self ):
        '''
        Simple offset bearing, distance, height difference 
        '''
        net=Network.Network()
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        net.addStation(st1)
        obs=[]
        obs.append(AZ(st1,st2))
        obs.append(SD(st1,st2))
        obs.append(LV(st1,st2))
        nfixed=StationLocator.locateStations(net,obs,sys.stdout.write)
        st2c=net.get('ST2')
        self.assertTrue(nfixed==1,'Failed to add one station calculated ({0} fixed)'.format(nfixed))
        self.assertTrue(st2c is not None,'ST2 not added to network')
        offset=st2c.distanceTo(st2)
        self.assertTrue(offset < 1.0,'ST2 out of position  by {0:.1f} metres'.format(offset))

    def test_003_basic_offset_ha( self ):
        '''
        Simple offset angle, distance, zenith distance
        '''
        net=Network.Network()
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st3=Station.Station('ST3',llh=(171.02,-45.01,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        net.addStation(st1)
        net.addStation(st3)
        obs=[]
        obs.append(HA(st1,[st3,st2]))
        obs.append(SD(st1,st2))
        obs.append(ZD(st1,st2))
        nfixed=StationLocator.locateStations(net,obs,sys.stdout.write)
        st2c=net.get('ST2')
        self.assertTrue(nfixed==1,'Failed to add one station calculated ({0} fixed)'.format(nfixed))
        self.assertTrue(st2c is not None,'ST2 not added to network')
        offset=st2c.distanceTo(st2)
        self.assertTrue(offset < 1.0,'ST2 out of position  by {0:.1f} metres'.format(offset))

if __name__=="__main__":
    unittest.main()
