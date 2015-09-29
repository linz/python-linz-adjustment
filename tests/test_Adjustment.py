
# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import numpy as np
import os.path
import StringIO
from LINZ import fileunittest

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'LINZ'))
from Geodetic import Station, Network, Adjustment

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from createobs import HA, AZ, SD, ZD, LV, GB, GX, Traverse

class AdjustmentTestCase( fileunittest.TestCase ):

    def offsetStation( self, st, dxyz ):
        st.setXYZ( st.xyz() + dxyz )

    def runAdjustment( self, test, adj, checkListing=True, checkGeoid=False ):
        outputfile=StringIO.StringIO()
        try:
            oldio=sys.stdout
            sys.stdout=outputfile
            adj.run()
        finally:
            sys.stdout=oldio
            output=outputfile.getvalue()
            outputfile.close()
        if checkListing:
            self.check(test+': Output',output)
        self.check(test+': Obs',adj.observations)
        for s in adj.stations.stations():
            self.check(test+': Station {0} coordinates'.format(s.code()),s.llh())
            if checkGeoid:
                grav=list(s.xieta())
                grav.append(s.geoidHeight())
                self.check(test+': Station {0} geoid'.format(s.code()),grav)


    def test_001_basic_offset( self ):
        '''
        Simple offset bearing, distance, zenith distance
        '''
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        obs=[]
        obs.append(AZ(st1,st2))
        obs.append(SD(st1,st2))
        obs.append(ZD(st1,st2))
        obs.append(LV(st1,st2))
        self.offsetStation(st2,[0.5,1.0,-0.3])
        net=Network.Network()
        net.addStation(st1)
        net.addStation(st2)
        adj=Adjustment.Adjustment(stations=net,observations=obs,verbose=True)
        adj.setConfig('fix','ST1')
        adj.setConfig('refraction_coefficient','0.0')
        self.runAdjustment('Test 1',adj)

    def test_002_gx( self ):
        '''
        Simple XYZ observation
        '''
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        obs=[]
        obs.append(GX(st1))
        self.offsetStation(st1,[0.5,1.0,-0.3])
        net=Network.Network()
        net.addStation(st1)
        adj=Adjustment.Adjustment(stations=net,observations=obs,verbose=True)
        self.runAdjustment('Test 2',adj)

    def test_003_gb( self ):
        '''
        Simple DXYZ observation
        '''
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        obs=[]
        obs.append(GB(st1,st2))
        self.offsetStation(st2,[0.5,1.0,-0.3])
        net=Network.Network()
        net.addStation(st1)
        net.addStation(st2)
        adj=Adjustment.Adjustment(stations=net,observations=obs,verbose=True)
        adj.setConfig('fix','ST1')
        adj.setConfig('refraction_coefficient','0.0')
        self.runAdjustment('Test 3',adj)

    def test_004_ha( self ):
        '''
        Simple offset bearing, distance, zenith distance
        '''
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.035,-44.982,20.0))
        st3=Station.Station('ST3',llh=(171.012,-44.972,20.0))
        obs=[]
        obs.append(HA(st2,[st1,st3]))
        obs.append(SD(st1,st2))
        obs.append(ZD(st1,st2))
        self.offsetStation(st2,[0.5,1.0,-0.3])
        net=Network.Network()
        net.addStation(st1)
        net.addStation(st2)
        net.addStation(st3)
        adj=Adjustment.Adjustment(stations=net,observations=obs,verbose=True)
        adj.setConfig('fix','ST1')
        adj.setConfig('fix','ST3')
        adj.setConfig('refraction_coefficient','0.0')
        self.runAdjustment('Test 4',adj)

    def test_005_basic_offset( self ):
        '''
        With station heights
        '''
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        obs=[]
        obs.append(AZ(st1,st2,insthgt=2.5,trgthgt=4.1))
        obs.append(SD(st1,st2,insthgt=2.5,trgthgt=-4.1))
        obs.append(ZD(st1,st2,insthgt=-2.5,trgthgt=4.1))
        obs.append(LV(st1,st2,insthgt=-2.5,trgthgt=-4.1))
        self.offsetStation(st2,[0.5,1.0,-0.3])
        net=Network.Network()
        net.addStation(st1)
        net.addStation(st2)
        adj=Adjustment.Adjustment(stations=net,observations=obs,verbose=True)
        adj.setConfig('fix','ST1')
        adj.setConfig('refraction_coefficient','0.0')
        self.runAdjustment('Test 5',adj)

    def test_100_locator_plugin( self ):
        '''
        Station locator plugin
        '''
        from Geodetic.StationLocatorPlugin import StationLocatorPlugin
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        obs=[]
        obs.append(AZ(st1,st2))
        obs.append(SD(st1,st2))
        obs.append(ZD(st1,st2))
        self.offsetStation(st2,[0.5,1.0,-0.3])
        net=Network.Network()
        net.addStation(st1)
        adj=Adjustment.Adjustment(stations=net,observations=obs,
                                  verbose=True,plugins=[StationLocatorPlugin])
        adj.setConfig('fix','ST1')
        adj.setConfig('calculate_missing_stations','yes')
        adj.setConfig('debug_calculate_missing_stations','yes')
        adj.setConfig('refraction_coefficient','0.0')
        self.runAdjustment('Test 100',adj)

    def test_101_locator_plugin_config_load( self ):
        '''
        Station locator plugin via configuration
        '''
        from Geodetic.StationLocatorPlugin import StationLocatorPlugin
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        obs=[]
        obs.append(AZ(st1,st2))
        obs.append(SD(st1,st2))
        obs.append(ZD(st1,st2))
        self.offsetStation(st2,[0.5,1.0,-0.3])
        net=Network.Network()
        net.addStation(st1)
        adj=Adjustment.Adjustment(stations=net,observations=obs,
                                  verbose=True)
        adj.setConfig('use_plugin','station_locator_plugin')
        adj.setConfig('fix','ST1')
        adj.setConfig('calculate_missing_stations','yes')
        adj.setConfig('debug_calculate_missing_stations','yes')
        adj.setConfig('refraction_coefficient','0.0')
        self.runAdjustment('Test 101',adj)

    def test_150_geoid_plugin_config_load( self ):
        '''
        Local geoid plugin via configuration
        '''
        from Geodetic.StationLocatorPlugin import StationLocatorPlugin
        st1=Station.Station('ST1',llh=(171.0,-45.0,100.0))
        st2=Station.Station('ST2',llh=(171.05,-45.0,100.0))
        st3=Station.Station('ST3',llh=(171.0,-44.93,100.0))
        st4=Station.Station('ST4',llh=(171.0,-44.80,100.0))
        obs=Traverse([st1,st2,st3,st4])
        net=Network.Network()
        net.addStation(st1,st2,st3,st4)
        adj=Adjustment.Adjustment(stations=net,observations=obs,
                                  verbose=True)
        adj.setConfig('use_plugin','local_geoid_model_plugin')
        adj.setConfig('fix','ST1 ST2 ST3 ST4')
        adj.setConfig('local_geoid','ST1 25.6 30.0 -45.0 15000')
        adj.setConfig('refraction_coefficient','0.0')
        self.runAdjustment('Test 150',adj,checkListing=True,checkGeoid=True)


if __name__=="__main__":
    fileunittest.main()
