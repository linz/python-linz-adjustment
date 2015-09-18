# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import numpy as np
import os.path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'LINZ'))
import unittest
from Geodetic import Station, Network, StationLocator

from createobs import HA, AZ, SD, ZD, LV, GX, Traverse

class StationLocatorTestCase( unittest.TestCase ):


    def checkLocator( self, net, obs, expected, tolerance=1.0 ):
        nfixed=StationLocator.locateStations(net,obs,sys.stdout.write)
        nexpected=len(expected)
        self.assertTrue(nfixed==nexpected,'Failed to add {0} calculated stations ({1} fixed)'.format(nexpected,nfixed))
        for s in expected:
            sf=net.get(s.code())
            self.assertTrue(sf is not None,'{0} not added to network'.format(s.code()))
            if sf is not None:
                offset=s.distanceTo(sf)
                print("{0} offset {1:.2f}".format(s.code(),offset))
                if offset > tolerance:
                    print("Offset exceeded!")
                    print(sf.code(),"located at",sf.llh())
                    print(s.code(),"should be",s.llh())
                    self.fail('{0} out of position  by {1:.1f} metres'.format(s.code(),offset))


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
        self.checkLocator(net,obs,[st2])

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
        self.checkLocator(net,obs,[st2])

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
        self.checkLocator(net,obs,[st2])

    def test_004_basic_offset_with_GX( self ):
        '''
        Simple using GX fixed station plus bearing, distance, zenith distance
        '''
        net=Network.Network()
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        obs=[]
        obs.append(GX(st1))
        obs.append(AZ(st1,st2))
        obs.append(SD(st1,st2))
        obs.append(ZD(st1,st2))
        self.checkLocator(net,obs,[st1,st2])

    def test_005_traverse( self ):
        '''
        Extending network along a traversej
        '''
        net=Network.Network()
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        st3=Station.Station('ST3',llh=(171.348,-44.825,40.0))
        st4=Station.Station('ST4',llh=(171.3589,-44.8262,310.0))
        stations=[st1,st2,st3,st4]
        net.addStation(st1)
        obs=[]
        obs.append(AZ(st1,st2))
        obs.extend(Traverse(stations))
        self.checkLocator(net,obs,stations[1:])

    def test_006_unoriented_traverse( self ):
        '''
        Angle connected at each end
        '''
        net=Network.Network()
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.008,-44.995,20.0))
        st3=Station.Station('ST3',llh=(171.004,-44.995,40.0))
        net.addStation(st1)
        net.addStation(st3)
        obs=Traverse([st1,st2,st3])
        self.checkLocator(net,obs,[st2])

    def test_007_unoriented_traverse2( self ):
        '''
        Traverse connected at each end
        '''
        net=Network.Network()
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.004,-44.995,40.0))
        st3=Station.Station('ST3',llh=(171.0,-44.995,20.0))
        st4=Station.Station('ST4',llh=(171.0083,-44.9965,50.0))
        net.addStation(st1)
        net.addStation(st2)
        obs=Traverse([st1,st4,st3,st2])
        self.checkLocator(net,obs,[st3,st4])

if __name__=="__main__":
    unittest.main()
