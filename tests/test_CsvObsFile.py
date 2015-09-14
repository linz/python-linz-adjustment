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

test1az='''
fromstn,tostn,az_value,az_error
ST1,ST2,95.2365,0.0015
ST2,ST3,87.2365,0.0020
'''

test1az_result='''
AZ,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1

AZ,ST2,ST3,0.0000,0.0000,87.2365,0.0020,source=input:2
'''

test1sd='''
fromstn,tostn,sd_value,sd_error
ST1,ST2,95.2365,0.0015
ST2,ST3,87.2365,0.0020
'''
test1sd_result='''
SD,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1

SD,ST2,ST3,0.0000,0.0000,87.2365,0.0020,source=input:2
'''

test1zd='''
fromstn,tostn,zd_value,zd_error
ST1,ST2,95.2365,0.0015
ST2,ST3,87.2365,0.0020
'''

test1zd_result='''
ZD,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1

ZD,ST2,ST3,0.0000,0.0000,87.2365,0.0020,source=input:2
'''

test1lv='''
fromstn,tostn,lv_value,lv_error
ST1,ST2,95.2365,0.0015
ST2,ST3,87.2365,0.0020
'''

test1lv_result='''
LV,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1

LV,ST2,ST3,0.0000,0.0000,87.2365,0.0020,source=input:2
'''

test1ha='''
fromstn,tostn,ha_value,ha_error,obsset
ST1,ST2,0.0,0.0015
ST1,ST3,95.2365,0.0015
ST2,ST3,87.2365,0.0020
ST2,ST1,305.28,0.003
ST2,ST4,23.92,0.002
'''
test1ha_result='''
HA,ST1,ST2,0.0000,0.0000,0.0000,0.0015,source=input:1
HA,ST1,ST3,0.0000,0.0000,95.2365,0.0015,source=input:2

HA,ST2,ST3,0.0000,0.0000,87.2365,0.0020,source=input:3
HA,ST2,ST1,0.0000,0.0000,305.2800,0.0030,source=input:4
HA,ST2,ST4,0.0000,0.0000,23.9200,0.0020,source=input:5
'''

test1sdh='''
fromstn,tostn,fromhgt,tohgt,sd_value,sd_error
ST1,ST2,0.12,1.3,95.2365,0.0015
ST2,ST3,0.15,1.52,87.2365,0.0020
'''
test1sdh_result='''
SD,ST1,ST2,0.1200,1.3000,95.2365,0.0015,source=input:1

SD,ST2,ST3,0.1500,1.5200,87.2365,0.0020,source=input:2
'''

test1='''
fromstn,tostn,az_value,az_error,sd_value,sd_error
ST1,ST2,95.2365,0.0015,281.322,0.008
ST1,ST2,5.2365,0.0015,,
ST2,ST3,,,181.522,0.008
'''

test1_result='''
SD,ST1,ST2,0.0000,0.0000,281.3220,0.0080,source=input:1

AZ,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1

AZ,ST1,ST2,0.0000,0.0000,5.2365,0.0015,source=input:2

SD,ST2,ST3,0.0000,0.0000,181.5220,0.0080,source=input:3
'''

test1_result='''
SD,ST1,ST2,0.0000,0.0000,281.3220,0.0080,source=input:1

AZ,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1

AZ,ST1,ST2,0.0000,0.0000,5.2365,0.0015,source=input:2

SD,ST2,ST3,0.0000,0.0000,181.5220,0.0080,source=input:3
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
SD,ST1,ST2,0.0000,0.0000,281.3220,0.0080,source=input:1

HA,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1
HA,ST1,ST3,0.0000,0.0000,5.2365,0.0015,source=input:2
HA,ST1,ST4,0.0000,0.0000,5.2365,0.0015,source=input:3

SD,ST2,ST1,0.0000,0.0000,281.3220,0.0080,source=input:4

HA,ST2,ST1,0.0000,0.0000,95.2365,0.0015,source=input:4
HA,ST2,ST3,0.0000,0.0000,5.2365,0.0015,source=input:5

SD,ST2,ST1,0.0000,0.0000,281.3220,0.0080,source=input:6

HA,ST2,ST1,0.0000,0.0000,95.2365,0.0015,source=input:6
HA,ST2,ST3,0.0000,0.0000,5.2365,0.0015,source=input:7
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
SD,ST1,ST2,0.0000,0.0000,195.2365,0.0025,source=input:2

SD,ST1,ST3,0.0000,0.0000,83.1365,0.0035,source=input:4

HA,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1
HA,ST1,ST3,0.0000,0.0000,5.2365,0.0015,source=input:3
HA,ST1,ST4,0.0000,0.0000,5.2365,0.0015,source=input:5

HA,ST2,ST1,0.0000,0.0000,95.2365,0.0015,source=input:6
HA,ST2,ST3,0.0000,0.0000,5.2365,0.0015,source=input:7

HA,ST2,ST1,0.0000,0.0000,95.2365,0.0015,source=input:8
HA,ST2,ST3,0.0000,0.0000,5.2365,0.0015,source=input:9
'''

test4='''
fromstn,tostn,daytime,obstype,test_attr,obsset,value,error
ST1,ST2,night,AZ,test,1,95.2365,0.0015
ST1,ST2,day,SD,,1,195.2365,0.0025
'''

test4_result='''
AZ,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1

SD,ST1,ST2,0.0000,0.0000,195.2365,0.0025,source=input:2
'''

test4a_result='''
AZ,ST1,ST2,0.0000,0.0000,95.2365,0.0015,daytime=night source=input:1 test_attr=test

SD,ST1,ST2,0.0000,0.0000,195.2365,0.0025,daytime=day source=input:2 test_attr=None
'''

test5='''
inst,target,insthgt,tgthgt,sd_value,error
ST1,ST2,0.12,1.3,95.2365,0.0015
ST2,ST3,0.15,1.52,87.2365,0.0020
'''

test5_result='''
SD,ST1,ST2,0.0000,0.0000,95.2365,0.0015,source=input:1

SD,ST2,ST3,0.0000,0.0000,87.2365,0.0020,source=input:2
'''



class CsvObsFileTestCase( unittest.TestCase ):

    def check_csv_output( self, input, expected, varname='result', attributes=None, colnames=None ):
        source=StringIO.StringIO(input.lstrip())
        result=''
        for obs in CsvObsFile.read(source,attributes=attributes,colnames=colnames):
            result=result+"\n"+str(obs)
        result=result.strip()
        expected=expected.strip()
        if result != expected:
            print(">"*20)
            print(varname+"='''")
            print(result)
            print("'''")
            print("<"*20)
        self.assertEqual(result,expected)

    def test_001_simple_obs_csv( self ):
        '''
        Simple obs csv
        '''
        self.check_csv_output(test1az,test1az_result,'test1az_result')
        self.check_csv_output(test1sd,test1sd_result,'test1sd_result')
        self.check_csv_output(test1zd,test1zd_result,'test1zd_result')
        self.check_csv_output(test1lv,test1lv_result,'test1lv_result')
        self.check_csv_output(test1ha,test1ha_result,'test1ha_result')
        self.check_csv_output(test1sdh,test1sdh_result,'test1sdh_result')
        self.check_csv_output(test1,test1_result,'test1_result')

    def test_002_ha_obs_csv( self ):
        '''
        Obs csv with HA
        '''
        self.check_csv_output(test2,test2_result,'test2_result')

    def test_003_obstype_csv( self ):
        '''
        Obs csv with obstype column
        '''
        self.check_csv_output(test3,test3_result,'test3_result')

    def test_004_obstype_csv( self ):
        '''
        Obs with attributes
        '''
        self.check_csv_output(test4,test4_result,'test4_result')
        self.check_csv_output(test4,test4a_result,'test4a_result',attributes='daytime test_attr?')

    def test_005_colnames( self ):
        '''
        Data with column renaming
        '''
        self.check_csv_output(test5,test5_result,'test5_result',colnames={
            'inst':'fromstn','target':'tostn',
            'value':'sd_value','error':'sd_error'})

if __name__=="__main__":
    unittest.main()
