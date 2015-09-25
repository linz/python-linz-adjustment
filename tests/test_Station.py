# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os.path
import numpy as np
from numpy import array
from LINZ import fileunittest

sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'LINZ'))
from Geodetic import Station

class StationTestCase( fileunittest.TestCase ):

    def checkStation( self, test, station, checkGrav=True, checkEnu=False ):
        self.check(test+' code',station.code())
        self.check(test+' name',station.name())
        self.check(test+' llh',station.llh())
        self.check(test+' xyz',station.xyz())
        if checkGrav:
            self.check(test+' xieta',station.xieta())
            self.check(test+' geoidHeight',station.geoidHeight())
        if checkEnu:
            self.check(test+' enu',station.enu())
            self.check(test+' genu',station.genu())

    def test_001_attributes( self ):
        st1=Station.Station('A23C',llh=(172.345,-41.286,23.5))
        self.checkStation('Test1a',st1,checkGrav=True,checkEnu=True)
        st2=Station.Station('A23-C',xyz=(-4754070.1304,625885.2064,-4191614.0057),name='A23 mark C')
        self.checkStation('Test1b',st2,checkGrav=True,checkEnu=True)
        st3=Station.Station('A23C',llh=(172.345,-41.286,23.5),xieta=(0.1,0.4),geoidhgt=25.0)
        self.checkStation('Test1c',st3,checkGrav=True,checkEnu=True)
        enu=st3.enu()
        genu=st3.genu()
        self.check('Test1_genu',np.degrees(genu[2].dot(enu.T)),'Gravitational vertical not correct',0.00001)

    def test_002_offsets( self ):
        st1=Station.Station('A23C',llh=(172.345,-41.286,23.5),xieta=(10.0,-4.0),geoidhgt=25.0)
        xyz=st1.xyz()
        enu=st1.enu()
        genu=st1.genu()

        xyz2=st1.xyz(offset=15.2)
        self.check('Test2 offset none',(xyz2-xyz).dot(genu.T))
        xyz2=st1.xyz(offset=23.2,offsettype=st1.OFFSET_H)
        self.check('Test2 offset h',(xyz2-xyz).dot(genu.T))
        xyz2=st1.xyz(offset=[-13.2,19.7,23.2],offsettype=st1.OFFSET_GENU)
        self.check('Test2 offset genu',(xyz2-xyz).dot(genu.T))
        xyz2=st1.xyz(offset=[-15.2,29.7,18.2],offsettype=st1.OFFSET_ENU)
        self.check('Test2 offset enu',(xyz2-xyz).dot(enu.T))
        xyz2=st1.xyz(offset=[-18.2,-19.7,8.2],offsettype=st1.OFFSET_XYZ)
        self.check('Test2 offset xyz',(xyz2-xyz))

    def test_003_calculations( self ):
        st1=Station.Station('A23C',llh=(172.345,-41.286,23.5),xieta=(10.0,-4.0),geoidhgt=25.0)
        st2=Station.Station('A23D',llh=(172.3421,-41.287,118.27),xieta=(3.0,1.5),geoidhgt=35.0)
        for f in (st1.distanceTo,st1.azimuthTo,st1.geodeticAzimuthTo,st1.heightDifferenceTo ):
            name=f.__name__
            self.check('Test3 '+name,f(st2))
            fdata=f(st2,ddxyz=True)
            self.check('Test3 '+name+' ddxyz',fdata)
            # Check differentials
            fdiff=[]
            for i in range(3):
                offset=[0.0,0.0,0.0]
                offset[i]=1.0
                fcalc=f(st2,instofst=offset,offsettype=st1.OFFSET_XYZ)
                fdiff.append(fcalc-fdata[0])
            maxerror=np.max(np.abs(np.subtract(fdata[1],fdiff)))
            # Tolerance because currently heightDifferenceTo differentials do not 
            # factor in slope of geoid (which is significant in this extreme test case)
            tolerance = 5.0e-3 if name != 'heightDifferenceTo' else 0.15
            if maxerror > tolerance:
                self.reportFail("Test3 {0} inst station differentials in error by {1}\n  {2} {3}"
                                .format(name,maxerror,fdata[1],fdiff))
            fdiff=[]
            for i in range(3):
                offset=[0.0,0.0,0.0]
                offset[i]=1.0
                fcalc=f(st2,trgtofst=offset,offsettype=st1.OFFSET_XYZ)
                fdiff.append(fcalc-fdata[0])
            maxerror=np.max(np.abs(np.subtract(fdata[2],fdiff)))
            # Tolerance because currently heightDifferenceTo differentials do not 
            # factor in slope of geoid (which is significant in this extreme test case)
            tolerance = 5.0e-3 if name != 'heightDifferenceTo' else 0.15
            if maxerror > tolerance:
                self.reportFail("Test3 {0} trgt station differentials in error by {1}\n  {2} {3}"
                                .format(name,maxerror,fdata[2],fdiff))

if __name__=="__main__":
    fileunittest.main()
