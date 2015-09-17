# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import StringIO
import os.path
from LINZ import fileunittest

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'LINZ'))
from Geodetic import CsvObsFile

test1az='''
fromstn,tostn,az_value,az_error
ST1,ST2,95.2365,0.0015
ST2,ST3,87.2365,0.0020
'''
test1sd='''
fromstn,tostn,sd_value,sd_error
ST1,ST2,95.2365,0.0015
ST2,ST3,87.2365,0.0020
'''
test1zd='''
fromstn,tostn,zd_value,zd_error
ST1,ST2,95.2365,0.0015
ST2,ST3,87.2365,0.0020
'''
test1lv='''
fromstn,tostn,lv_value,lv_error
ST1,ST2,95.2365,0.0015
ST2,ST3,87.2365,0.0020
'''
test1ha='''
fromstn,tostn,ha_value,ha_error,obsset
ST1,ST2,0.0,0.0015
ST1,ST3,95.2365,0.0015
ST2,ST3,87.2365,0.0020
ST2,ST1,305.28,0.003
ST2,ST4,23.92,0.002
'''
test1sdh='''
fromstn,tostn,fromhgt,tohgt,sd_value,sd_error
ST1,ST2,0.12,1.3,95.2365,0.0015
ST2,ST3,0.15,1.52,87.2365,0.0020
'''
test1='''
fromstn,tostn,az_value,az_error,sd_value,sd_error
ST1,ST2,95.2365,0.0015,281.322,0.008
ST1,ST2,5.2365,0.0015,,
ST2,ST3,,,181.522,0.008
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

test4='''
fromstn,tostn,daytime,obstype,test_attr,obsset,value,error
ST1,ST2,night,AZ,test,1,95.2365,0.0015
ST1,ST2,day,SD,,1,195.2365,0.0025
'''

test5='''
inst,target,insthgt,tgthgt,sd_value,error
ST1,ST2,0.12,1.3,95.2365,0.0015
ST2,ST3,0.15,1.52,87.2365,0.0020
'''

class CsvObsFileTestCase( fileunittest.TestCase ):

    def check_csv_output( self, input, varname, attributes=None, colnames=None ):
        global printAllResults
        source=StringIO.StringIO(input.lstrip())
        result=''
        for obs in CsvObsFile.read(source,attributes=attributes,colnames=colnames):
            result=result+"\n"+str(obs)
        self.check(varname,result)

    def test_001_simple_obs_csv( self ):
        '''
        Simple obs csv
        '''
        self.check_csv_output(test1az,'test1az_result')
        self.check_csv_output(test1sd,'test1sd_result')
        self.check_csv_output(test1zd,'test1zd_result')
        self.check_csv_output(test1lv,'test1lv_result')
        self.check_csv_output(test1ha,'test1ha_result')
        self.check_csv_output(test1sdh,'test1sdh_result')
        self.check_csv_output(test1,'test1_result')

    def test_002_ha_obs_csv( self ):
        '''
        Obs csv with HA
        '''
        self.check_csv_output(test2,'test2_result')

    def test_003_obstype_csv( self ):
        '''
        Obs csv with obstype column
        '''
        self.check_csv_output(test3,'test3_result')

    def test_004_obstype_csv( self ):
        '''
        Obs with attributes
        '''
        self.check_csv_output(test4,'test4_result')
        self.check_csv_output(test4,'test4a_result',attributes='daytime test_attr?')

    def test_005_colnames( self ):
        '''
        Data with column renaming
        '''
        self.check_csv_output(test5,'test5_result',colnames={
            'inst':'fromstn','target':'tostn',
            'value':'sd_value','error':'sd_error'})

if __name__=="__main__":
    fileunittest.main()
