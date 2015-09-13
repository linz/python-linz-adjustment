# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import numpy as np
import StringIO
import os.path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'LINZ'))
import unittest
from Geodetic import CsvObsFile

test1='''
fromstn,tostn,az_value,az_error,sd_value,sd_error
ST1,ST2,95.2365,0.0015,281.322,0.008
ST1,ST2,5.2365,0.0015,,
ST2,ST3,,,181.522,0.008
'''

test1_result='''
SD,ST1,ST2,0.0000,0.0000,281.3220,0.0080

AZ,ST1,ST2,0.0000,0.0000,95.2365,0.0015

AZ,ST1,ST2,0.0000,0.0000,5.2365,0.0015

SD,ST2,ST3,0.0000,0.0000,181.5220,0.0080
'''

test2='''
fromstn,tostn,obsset,ha_value,ha_error,sd_value,sd_error
ST1,ST2,1,95.2365,0.0015,281.322,0.008
ST1,ST3,1,5.2365,0.0015,,
ST1,ST4,1,5.2365,0.0015,,
ST2,ST1,1,95.2365,0.0015,281.322,0.008
ST2,ST3,1,5.2365,0.0015,,
ST2,ST1,2,95.2365,0.0015,281.322,0.008
ST2,ST3,2,5.2365,0.0015,,
'''

test2_result='''
SD,ST1,ST2,0.0000,0.0000,281.3220,0.0080

HA,ST1,ST2,0.0000,0.0000,95.2365,0.0015
HA,ST1,ST3,0.0000,0.0000,5.2365,0.0015
HA,ST1,ST4,0.0000,0.0000,5.2365,0.0015

SD,ST2,ST1,0.0000,0.0000,281.3220,0.0080

HA,ST2,ST1,0.0000,0.0000,95.2365,0.0015
HA,ST2,ST3,0.0000,0.0000,5.2365,0.0015

SD,ST2,ST1,0.0000,0.0000,281.3220,0.0080

HA,ST2,ST1,0.0000,0.0000,95.2365,0.0015
HA,ST2,ST3,0.0000,0.0000,5.2365,0.0015
'''

test3='''
fromstn,tostn,obstype,obsset,value,error
ST1,ST2,HA,1,95.2365,0.0015
ST1,ST2,SD,1,195.2365,0.0025
ST1,ST3,HA,1,5.2365,0.0015
ST1,ST3,SD,1,83.1365,0.0035
ST1,ST4,HA,1,5.2365,0.0015
ST2,ST1,HA,1,95.2365,0.0015
ST2,ST3,HA,1,5.2365,0.0015
ST2,ST1,HA,2,95.2365,0.0015
ST2,ST3,HA,2,5.2365,0.0015
'''
test3_result='''
SD,ST1,ST2,0.0000,0.0000,195.2365,0.0025

SD,ST1,ST3,0.0000,0.0000,83.1365,0.0035

HA,ST1,ST2,0.0000,0.0000,95.2365,0.0015
HA,ST1,ST3,0.0000,0.0000,5.2365,0.0015
HA,ST1,ST4,0.0000,0.0000,5.2365,0.0015

HA,ST2,ST1,0.0000,0.0000,95.2365,0.0015
HA,ST2,ST3,0.0000,0.0000,5.2365,0.0015

HA,ST2,ST1,0.0000,0.0000,95.2365,0.0015
HA,ST2,ST3,0.0000,0.0000,5.2365,0.0015
'''

class CsvObsFileTestCase( unittest.TestCase ):

    def check_csv_output( self, input, expected ):
        source=StringIO.StringIO(input.lstrip())
        result=''
        for obs in CsvObsFile.read(source):
            result=result+"\n"+str(obs)
        result=result.strip()
        expected=expected.strip()
        if result != expected:
            print("=================================")
            print(result)
            print("=================================")
        self.assertEqual(result,expected)

    def test_001_simple_obs_csv( self ):
        '''
        Simple obs csv
        '''
        self.check_csv_output(test1,test1_result)

    def test_002_ha_obs_csv( self ):
        '''
        Obs csv with HA
        '''
        self.check_csv_output(test2,test2_result)

    def test_003_obstype_csv( self ):
        '''
        Obs csv with obstype column
        '''
        self.check_csv_output(test3,test3_result)

if __name__=="__main__":
    unittest.main()
