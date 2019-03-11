# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import numpy as np
import os.path
from LINZ.Geodetic.Test import fileunittest

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'LINZ'))
from Geodetic import MergeAnglesPlugin
from Geodetic.Observation import Observation, ObservationValue

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from createobs import HA, AZ, SD, ZD, LV, GX, Traverse

class StationLocatorTestCase( fileunittest.TestCase ):

    def checkMergeAngles( self, test, obslist ):
        for i,o in enumerate(obslist):
            self.check("{0} input obs {1}:\n".format(test,i),o)
        obslist=MergeAnglesPlugin._mergeAngleObservations(obslist)
        for i,o in enumerate(obslist):
            self.check("{0} output obs {1}:\n".format(test,i),o)

    def angleSet( self, fromstn, *angles ):
        o=Observation('HA')
        for tostn,value in zip(angles[0::2],angles[1::2]):
            o.addObservation(ObservationValue(fromstn,tostn,value,0.001))
        return o

    def distanceSet( self, fromstn, *angles ):
        o=Observation('SD')
        for tostn,value in zip(angles[0::2],angles[1::2]):
            o.addObservation(ObservationValue(fromstn,tostn,value,0.001))
        return o

    def test_000_simplesets( self ):
        '''
        Simple offset bearing, distance, zenith distance
        '''
        obslist=[]
        obslist.append(self.angleSet('S0','S1',42.0,'S2',152.0))
        obslist.append(self.angleSet('S1','S0',32.0,'S2',254.0))
        obslist.append(self.distanceSet('S1','S0',1032.0,'S2',254.0))
        obslist.append(self.angleSet('S0','S2',42.0,'S3',252.0))
        obslist.append(self.angleSet('S0','S3',42.0,'S4',252.0,'S4a',255.0))
        obslist.append(self.angleSet('S0','S5',42.0,'S6',252.0))
        obslist.append(self.angleSet('S0','S7',42.0,'S8',252.0))
        obslist.append(self.angleSet('S0','S8',42.0,'S9',252.0))
        self.checkMergeAngles('Test 000',obslist)

if __name__=="__main__":
    fileunittest.main()
