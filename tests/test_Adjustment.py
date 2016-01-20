
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

# Install path relative to LINZ to avoid potentially picking up installed modules...
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'LINZ'))
from Geodetic import Station, Network, Adjustment

basedir=os.path.dirname(os.path.abspath(__file__))
sys.path.append(basedir)
from createobs import HA, AZ, SD, ZD, LV, GB, GX, Traverse

class AdjustmentTestCase( fileunittest.TestCase ):

    def offsetStation( self, st, dxyz ):
        st.setXYZ( st.xyz() + dxyz )

    def runAdjustment( self, test, adj, checkListing=True, checkGeoid=False, outputfiles={} ):
        global basedir
        basepath=basedir
        outputfile=StringIO.StringIO()
        try:
            adj.setOutputFile(outputfile)
            adj.run()
        finally:
            output=outputfile.getvalue()
            outputfile.close()
            output=output.replace(basedir,'')
        if checkListing:
            self.check(test+': Output',output)
        # Clean directory references
        for o in adj.observations:
            for obs in o.obsvalues:
                if 'source' in obs.attributes:
                    obs.attributes['source']=obs.attributes['source'].replace(basedir,'')
        self.check(test+': Obs',adj.observations)
        for s in adj.stations.stations():
            self.check(test+': Station {0} coordinates'.format(s.code()),s.llh())
            if checkGeoid:
                grav=list(s.xieta())
                grav.append(s.geoidHeight())
                self.check(test+': Station {0} geoid'.format(s.code()),grav)
        # Output files
        for f in sorted(outputfiles):
            filename=outputfiles[f]
            try:
                with open(filename) as opf:
                    filedata=opf.read()
            except:
                filedata='missing'
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except:
                pass
            self.check(test+': Output file '+f,filedata)

    def outputFilePath( self, filename ):
        global basedir
        outdir=os.path.join(basedir,'outputfiles')
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        return os.path.join(outdir,filename)

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

    def test_050_output_coords( self ):
        '''
        With station heights
        '''
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.001,-44.995,20.0))
        obs=[]
        obs.append(AZ(st1,st2,0.0015))
        obs.append(SD(st1,st2,0.005))
        obs.append(LV(st1,st2,0.025))
        net=Network.Network()
        net.addStation(st1)
        net.addStation(st2)
        adj=Adjustment.Adjustment(stations=net,observations=obs,verbose=True)
        adj.setConfig('fix','ST1')
        adj.setConfig('refraction_coefficient','0.0')
        cfname=self.outputFilePath('test50_coords.csv')
        adj.setConfig('output_coordinate_file',cfname)
        self.runAdjustment('Test 50',adj,outputfiles={'coords':cfname})

    def test_051_output_coord_covar( self ):
        '''
        With station heights
        '''
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.001,-44.995,20.0))
        obs=[]
        obs.append(AZ(st1,st2,0.0015))
        obs.append(SD(st1,st2,0.005))
        obs.append(LV(st1,st2,0.025))
        net=Network.Network()
        net.addStation(st1)
        net.addStation(st2)
        adj=Adjustment.Adjustment(stations=net,observations=obs,verbose=True)
        adj.setConfig('fix','ST1')
        adj.setConfig('refraction_coefficient','0.0')
        cfname=self.outputFilePath('test51_coords.csv')
        adj.setConfig('output_coordinate_file','covariances '+cfname)
        self.runAdjustment('Test 51',adj,outputfiles={'coords':cfname})

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

    def test_200_csv_attribute( self ):
        '''
        Read CSV file attributes
        '''
        global basedir
        df=os.path.join(basedir,'data','testadj1.csv')
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=(171.35,-44.82,20.0))
        net=Network.Network()
        net.addStation(st1,st2)
        adj=Adjustment.Adjustment(stations=net, verbose=True)
        adj.setConfig('data_file','"'+df+'" attributes=eqpt,setup')
        adj.setConfig('fix','ST1 ST2')
        adj.loadDataFiles()
        i=0
        for o in adj.observations:
            for v in o.obsvalues:
                i+=1
                self.check('Test 200: Observation '+str(i)+' attributes:',[v.attributes.get('eqpt'),v.attributes.get('setup')])


    def test_250_setup_height_calcs( self ):
        '''
        Height setup calculator
        '''
        global basedir
        df=os.path.join(basedir,'data','testadj1.csv')
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=[ 171.0024646771,  -45.000160237 ,   15.9343333328])
        net=Network.Network()
        net.addStation(st1,st2)
        adj=Adjustment.Adjustment(stations=net, verbose=True)
        adj.setConfig('data_file','"'+df+'" attributes=eqpt,setup')
        adj.setConfig('fix','ST1')
        adj.setConfig('use_plugin','setup_height_plugin')
        adj.setConfig('calculate_setup_heights','true')
        adj.setConfig('inst_trgt_setup_attributes','none setup')
        adj.setConfig('valid_setup_regex','[A-Z]')
        adj.setConfig('fix_setup_height','A 0.25')
        adj.setConfig('debug_observation_equations','true')
        self.runAdjustment('Test 250',adj,checkListing=True)

    def test_251_setup_height_re( self ):
        '''
        Height setup calculator
        '''
        global basedir
        df=os.path.join(basedir,'data','testadj1.csv')
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=[ 171.0024646771,  -45.000160237 ,   15.9343333328])
        net=Network.Network()
        net.addStation(st1,st2)
        adj=Adjustment.Adjustment(stations=net, verbose=True)
        adj.setConfig('data_file','"'+df+'" attributes=eqpt,setup')
        adj.setConfig('fix','ST1')
        adj.setConfig('use_plugin','setup_height_plugin')
        adj.setConfig('calculate_setup_heights','true')
        adj.setConfig('inst_trgt_setup_attributes','none setup')
        adj.setConfig('valid_setup_regex','[A-Z]')
        adj.setConfig('fix_setup_height','[AB] 0.25')
        adj.setConfig('debug_observation_equations','true')
        self.runAdjustment('Test 251',adj,checkListing=True)

    def test_252_setup_height_float( self ):
        '''
        Height setup calculator
        '''
        global basedir
        df=os.path.join(basedir,'data','testadj1.csv')
        st1=Station.Station('ST1',llh=(171.0,-45.0,10.0))
        st2=Station.Station('ST2',llh=[ 171.0024646771,  -45.000160237 ,   15.9343333328])
        net=Network.Network()
        net.addStation(st1,st2)
        adj=Adjustment.Adjustment(stations=net, verbose=True)
        adj.setConfig('data_file','"'+df+'" attributes=eqpt,setup')
        adj.setConfig('fix','ST1')
        adj.setConfig('use_plugin','setup_height_plugin')
        adj.setConfig('calculate_setup_heights','true')
        adj.setConfig('inst_trgt_setup_attributes','none setup')
        adj.setConfig('valid_setup_regex','[A-Z]')
        adj.setConfig('fix_setup_height','A 0.25')
        adj.setConfig('fix_setup_height','B 0.75 0.0025')
        adj.setConfig('debug_observation_equations','true')
        self.runAdjustment('Test 252',adj,checkListing=True)

if __name__=="__main__":
    fileunittest.main()
