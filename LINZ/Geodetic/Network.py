
import csv
import numpy as np

from Station import Station
from CsvRecord import Reader

class Network:
    '''
    Class representing a network of stations
    '''
    
    _defaultColumnsXYZ='code name X Y Z xi eta geoidhgt'.split()
    _defaultColumnsLLH='code name longitude latitude ellheight xi eta geoidhgt'.split()

    def __init__ ( self ):
        self._index={}
        self._geodeticCsv=None

    def addStation( self, stn ):
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

    def writeCsv( self, csvfile, geodetic=None, colnames=None ):
        '''
        Write stations to a csv file.  Column names default to 
        Network._defaultColumnsXYZ, or Network._defaultColumnsLLH for 
        if geodetic is True (which means output lat/lon/ellheight parameters)
        If geodetic is not defined then will use XYZ or LLH depending on the type
        of the last file loaded.

        colnames is an optional dictionary mapping from column ids to actual
        column names.
        '''

        geodetic=self._geodeticCsv if geodetic is None else geodetic

        cols=Network._defaultColumnsLLH if geodetic else Network._defaultColumnsXYZ
        if colnames:
            cols=list((colnames.get(c,c) for c in cols))

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
                csvw.writerow(row)

    def get( self, code ):
        return self._index.get(code,None)

    def stations( self ):
        for code in sorted(self._index.keys()):
            yield self._index[code]
        
    def setLocalGeoid( self, code, geoidHeight, xieta, range=None ):
        '''
        Applies a local geoid model to the stations.  The model is defined
        by a geoid height and deflection at one of the marks.  Optionally 
        define a range beyond which the model is not defined.
        '''
        refStation=self.get(code)
        if refStation is None:
            raise RuntimeError("Local geoid model reference station {0} is not defined"
                               .format(code))

        xyz0=refStation.xyz()
        enu=refStation.enu()

        range2=None if range is None else range*range
        # Note slope is dh/dn, dh/de as xieta are ordered lat,lon
        slope=-np.radians(xieta)

        for s in self.stations():
            offset=s.xyz()-xyz0
            denu=enu.dot(offset)
            if range is not None and (denu[0]*denu[0]+denu[1]*denu[1]) > range:
                continue
            ghgt=geoidHeight+slope.dot((denu[1],denu[0]))
            s.setXiEta(xieta)
            s.setGeoidHeight(ghgt)
