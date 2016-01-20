
# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import csv
import numpy as np

from .Station import Station
from .CsvRecord import Reader

class Network:
    '''
    Class representing a network of stations
    '''
    
    _defaultColumnsXYZ='code name X Y Z xi eta geoidhgt'.split()
    _defaultColumnsLLH='code name longitude latitude ellheight xi eta geoidhgt'.split()

    def __init__ ( self ):
        self._index={}
        self._geodeticCsv=None

    def addStation( self, *stations ):
        for stn in stations:
            self._index[stn.code()]=stn

    def loadCsv( self, csvfile, colnames=None ):
        '''
        Load stations from a csv file.  
        Column ids (and default column names) are:
            code
            name
            x y z
            lon lat hgt
            longitude latitude (ellheight or height)
            xi eta
            geoid
        Column names are treated as case insensitive.  Either x y z or 
        lon lat hgt should be provided. code is required.  The column names
        can be updated by a colnames parameter, a dictionary from column id 
        to actual column name.

        colnames is an optional dictionary mapping from column ids to actual
        column names.
        '''

        reader=Reader('station','''
                      code name? 
                      x:float? y:float? z:float?
                      lon:float? lat:float? hgt:float?
                      longitude:float? latitude:float? height:float? ellheight:float?
                      xi:float? eta:float?
                      geoidhgt:float?
                      ''')

        for stn in reader.readCsv(csvfile,colnames):
            code=stn.code
            if not code:
                continue
            name=stn.name or code
            xyz=None
            llh=None
            xieta=None
            if stn.x is not None:
                xyz=(stn.x,stn.y,stn.z)
                
            elif stn.lon is not None:
                llh=(stn.lon,stn.lat,stn.hgt)
                self._geodeticCsv=True
            elif stn.longitude is not None:
                hgt=stn.ellheight if stn.ellheight is not None else stn.height
                llh=(stn.longitude,stn.latitude,hgt)
                self._geodeticCsv=True
            if stn.xi is not None:
                xieta=(stn.xi/3600,stn.eta/3600)
            geoidhgt=stn.geoidhgt
            stn=Station(code,name=name,xyz=xyz,llh=llh,xieta=xieta,geoidhgt=geoidhgt)
            self.addStation(stn)

    def _datafunc( self, extradata ):
        extradata = extradata or {}
        coldata={}
        for v in extradata.values():
            coldata.update(v)
        columns=sorted(coldata.keys())
        if not columns:
            return lambda x: []
        def func(code):
            if code is none:
                return columns
            data=extradata.get(code,{})
            return [data.get(c) for c in columns]

    def writeCsv( self, csvfile, geodetic=None, colnames=None, extradata=None ):
        '''
        Write stations to a csv file.  Column names default to 
        Network._defaultColumnsXYZ, or Network._defaultColumnsLLH for 
        if geodetic is True (which means output lat/lon/ellheight parameters)
        If geodetic is not defined then will use XYZ or LLH depending on the type
        of the last file loaded.

        colnames is an optional dictionary mapping from column ids to actual
        column names.

        extradata can be: 
            a dictionary of dictionaries, keyed on station code with values
            a key/pair dictionary
            a function that returns a list of column names if called with None
            as a parameter, or a list of values if called with code as a parameter.
        '''
        
        if not callable(extradata):
            extradata=self._datafunc(extradata)

        geodetic=self._geodeticCsv if geodetic is None else geodetic

        cols=Network._defaultColumnsLLH if geodetic else Network._defaultColumnsXYZ
        colnames=colnames or {}
        cols=[colnames.get(c,c) for c in cols]
        cols.extend(extradata(None))

        with open(csvfile,'wb') as csvfh:
            csvw=csv.writer(csvfh)
            csvw.writerow(cols)
            for stn in self.stations():
                row=[stn.code(),stn.name()]
                if geodetic:
                    row.extend(("{0:.{1}f}".format(c,n) for c,n in zip(stn.llh(),(9,9,4))))
                else:
                    row.extend(("{0:.4f}".format(c) for c in stn.xyz()))
                xieta=stn.xieta()
                row.append(xieta[0]*3600)
                row.append(xieta[1]*3600)
                row.append(stn.geoidHeight())
                row.extend(extradata(stn.code()))
                csvw.writerow(row)

    def get( self, code ):
        return self._index.get(code,None)

    def stations( self ):
        for code in sorted(self._index.keys()):
            yield self._index[code]
