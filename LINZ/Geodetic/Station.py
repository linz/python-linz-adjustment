
# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import math
from LINZ.Geodetic import Ellipsoid

'''
Station: module defining a station used in an adjustment
'''

class Station( object ):

    _ellipsoid = Ellipsoid.GRS80

    OFFSET_H=1
    OFFSET_ENU=2
    OFFSET_GENU=3
    OFFSET_XYZ=4

    def __init__( self, code, name=None, llh=None, xyz=None, xieta=None, geoidhgt=None ):
        '''
        Initiallize a station.  Can specify either xyz, or llh coordinates.
        Can also specify deflection of the vertical xieta, and geoid height 
        geoidhgt.  Point is identified by a code, but can also specify an 
        name.
        '''
        self._code=code
        self._name=name or code
        self._xyz=None
        self.setXiEta( xieta )
        self.setGeoidHeight( geoidhgt )
        self._xyz=None
        if llh is not None:
            self.setLatLonHgt( llh )
        elif xyz is not None:
            self.setXYZ( xyz )

    def code (self ):
        '''
        Return the code identifying the point
        '''
        return self._code

    def name (self ):
        '''
        Return the name identifying the point
        '''
        return self._name

    def setXiEta( self, xieta ):
        '''
        Set the north (xi) and east (eta)  deflection of the vertical
        (in degrees)
        '''
        self._xieta=[0.0,0.0] if xieta is None else xieta
        if self._xyz is not None:
            self.setXYZ( self._xyz )

    def setGeoidHeight( self, geoidhgt ):
        '''
        Set the geoid height of the point
        '''
        self._geoidhgt=geoidhgt or 0.0
        if self._xyz is not None:
            self.setXYZ( self._xyz )

    def setLatLonHgt( self, llh ):
        '''
        Set the coordinates of the point in terms of longitude,
        latitude, and height
        '''
        xyz=Station._ellipsoid.xyz( llh )
        self.setXYZ( xyz )

    def setXYZ( self, xyz ):
        '''
        Set the coordinates of the point in terms of geocentric XYZ 
        coordinate
        '''
        self._xyz=np.array(xyz)
        lon,lat,hgt=Station._ellipsoid.geodetic(xyz)
        self._llh=[lon,lat,hgt]
        self._enu=Station._ellipsoid.enu_axes(lon,lat)
        self._genu=Station._ellipsoid.enu_axes(
            lon+self._xieta[1]/math.cos(math.radians(lat)),
            lat+self._xieta[0])

    def llh( self ):
        ''' 
        Return the coordinates of the point as lon, lat, height
        (lon/lat in degrees, height in metres)
        '''
        return Station._ellipsoid.geodetic(self.xyz())

    def xyz( self, offset=None, offsettype=None ):
        '''
        Return the coordinates of the point as geocentric X, Y, Z

        Can include an offset instofst, which is either None, 
        a floating point offset height, or a numpy array of E,N,U
        offsets.
        '''
        if self._xyz is None:
            raise RuntimeError('Coordinates not defined for station '+self._code)
        if offset is None:
            return self._xyz
        if offsettype is None:
            offsettype=Station.OFFSET_H
        if offsettype == Station.OFFSET_H:
            offset=[0,0,offset]
        if offsettype == Station.OFFSET_H or offsettype == Station.OFFSET_GENU:
            offset=np.array(offset).dot(self._genu)
        elif offsettype == Station.OFFSET_ENU:
            offset=np.array(offset).dot(self._enu)
        elif offsettype != Station.OFFSET_XYZ:
            raise RuntimeError('Invalid offset type in Station.xyz')
        return self._xyz+offset

    def xieta( self ):
        '''
        Return the deflection of the vertical (north, east) in degrees
        '''
        return self._xieta

    def geoidHeight( self ):
        '''
        Return the geoid height
        '''
        return self._geoidhgt

    def enu( self ):
        '''
        Returns the east/north/up vectors in the geodetic coordinate system
        '''
        # Check coordinates are defined
        if self._enu is None:
            raise RuntimeError('Coordinates not defined for station '+self._code)
        return self._enu

    def genu( self ):
        '''
        Returns the east/north/up vectors in the gravitational coordinate system
        '''
        # Check coordinates are defined
        if self._genu is None:
            raise RuntimeError('Coordinates not defined for station '+self._code)
        return self._genu

    # Calculate functions have a common set of parameters
    #  self (instrument station)
    #  trgtstn (target station)
    #  instofst (instrument station height)
    #  trgtofst (target station height)
    #  refcoef (refraction coefficient)
    #  ddxyz (calculate differentials wrt x,y, z)
    #
    # May be good to work out a tidier way of handling refraction coefficient
    # and other unused parameters.


    def vectorTo( self, trgtstn=None, instofst=None, trgtofst=None, refcoef=None, ddxyz=False, offsettype=None ):
        diff=trgtstn.xyz(trgtofst, offsettype=offsettype) - self.xyz(instofst,offsettype=offsettype)
        if not ddxyz:
            return diff
        return diff, -np.identity(3), np.identity(3)

    def distanceTo( self, trgtstn=None, instofst=None, trgtofst=None, refcoef=None, ddxyz=False, offsettype=None ):
        '''
        Calculate distance to trgtstn point and optionally its 
        differentials wrt X,Y,Z coordinates of base point and target
        point
        '''
        diff=self.vectorTo(trgtstn, instofst, trgtofst, offsettype=offsettype)
        dist=np.sqrt(np.vdot(diff,diff))
        if not ddxyz:
            return dist
        diff /= dist
        return dist, -diff, diff

    def horizontalDistanceTo( self, trgtstn=None, instofst=None, trgtofst=None, refcoef=None, ddxyz=False, offsettype=None ):
        '''
        Calculate distance to trgtstn point and optionally its 
        differentials wrt X,Y,Z coordinates of base point and target
        point
        '''
        diff=self.vectorTo(trgtstn, instofst, trgtofst, offsettype=offsettype)
        # Want the mean vertical vector...
        enu=self.enu()[2]+trgtstn.enu()[2]
        enu /= np.sqrt(np.vdot(diff,diff))
        diff = diff - enu * np.vdot(diff,enu)
        dist=np.sqrt(np.vdot(diff,diff))
        if not ddxyz:
            return dist
        diff /= dist
        return dist, -diff, diff

    def azimuthTo( self, trgtstn=None, instofst=None, trgtofst=None, refcoef=None, ddxyz=False, geodetic=False, offsettype=None ):
        '''
        Calculate bearing (degrees) to trgtstn point and its 
        differentials wrt X,Y,Z coordinates of base point and target
        point
        '''
        diff=self.vectorTo(trgtstn, instofst, trgtofst, offsettype=offsettype)
        enu=self._enu if geodetic else self._genu
        denu=enu.dot(diff)
        brng=math.degrees(math.atan2(denu[0],denu[1]))
        if not ddxyz:
            return brng
        # Calculate the differential in ENU system
        # (wrt self coordinates)
        diff=np.array([denu[1],-denu[0],0])
        diff /= np.vdot(diff,diff)
        # Convert to XYZ
        diff = enu.T.dot(diff)
        # Convert to differential of degrees
        diff *= np.degrees(1)
        return brng, -diff, diff

    def geodeticAzimuthTo( self, trgtstn=None, instofst=None, trgtofst=None, refcoef=None, ddxyz=False, offsettype=None ):
        return self.azimuthTo( trgtstn=trgtstn, instofst=instofst, trgtofst=trgtofst, refcoef=refcoef, ddxyz=ddxyz, geodetic=True, offsettype=offsettype )

    def zenithDistanceTo( self, trgtstn=None, instofst=None, trgtofst=None, refcoef=0.0, ddxyz=False, offsettype=None ):
        '''
        Calculate zenith distance (degrees) to trgtstn point and 
        its differentials wrt X,Y,Z coordinates of base point and target
        point
        '''
        diff=self.vectorTo(trgtstn, instofst, trgtofst, offsettype=offsettype)
        denu=self._genu.dot(diff)
        hordist=math.hypot(denu[0],denu[1])
        brng=math.atan2(hordist,denu[2])
        if refcoef != 0.0:
            subtend=np.sqrt((1.0-self._enu[2].dot(trgtstn._enu[2]))/2)
            # Using the same version of refraction coefficient as SNAP, hence multiply
            # by 2
            corr=-refcoef*subtend*2.0
            brng=brng+corr;
        brng=math.degrees(brng)
        if not ddxyz:
            return brng

        # Calculate the differential in ENU system
        # (wrt self coordinates)
        diff=np.array([denu[2]*denu[0]/hordist,denu[2]*denu[1]/hordist,-hordist])
        diff /= np.vdot(diff,diff)
        # Convert to XYZ
        diff = self._enu.T.dot(diff)
        # Convert to differential of degrees
        diff *= np.degrees(1)
        return brng, -diff, diff

    def heightDifferenceTo( self, trgtstn=None, instofst=None, trgtofst=None, refcoef=None, ddxyz=False, offsettype=None ):
        '''
        Calculate height difference between trgtstn point and this and optionally
        its differentials wrt X,Y,Z coordinates of base point and target
        point.  
        '''
        # Note: some lack of rigour here in whether we are considering offsets in terms
        # of geometric or gravitational vertical. Strictly should be based on gravitational, 
        # but this will be inconsistent if the geoid height is not also updated to reflect
        # the deflection of vertical...
        #
        # Generally differences will be inconsequential...

        hgtdiff=((trgtstn._llh[2]-trgtstn._geoidhgt)
                 -(self._llh[2]-self._geoidhgt))
        if instofst is not None or trgtofst is not None:
            if offsettype is None:
                offsettype=Station.OFFSET_H
            if offsettype==Station.OFFSET_H:
                hgtdiff += (trgtofst or 0.0) - (instofst or 0.0)
            else:
                instofst=instofst if instofst is not None else  [0.0,0.0,0.0]
                trgtofst=trgtofst if trgtofst is not None else  [0.0,0.0,0.0]
                instofst=np.array(instofst)
                trgtofst=np.array(trgtofst)
                if offsettype == Station.OFFSET_XYZ:
                    instofst = instofst.dot(self._genu.T)
                    trgtofst = trgtofst.dot(self._genu.T)
                elif offsettype == Station.OFFSET_ENU:
                    instofst = instofst.dot(self._enu).dot(self._genu.T)
                    trgtofst = trgtofst.dot(self._enu).dot(self._genu.T)
                elif offsettype != Station.OFFSET_GENU:
                    raise RuntimeError("Invalid offsettype in Station.heightDifferentTo")
                hgtdiff += trgtofst[2] - instofst[2]

        if not ddxyz:
            return hgtdiff
        # Technically differentials should be based on genu rather than enu, but
        # not changing in geoid height as point is shifted, so this is mathematically
        # correct with current formulation.
        return hgtdiff, -self._enu[2], trgtstn._enu[2].copy()

    # Note - trgtstn and trgtofst included to allow common function call with
    # other observation types
    def calcXYZ( self, trgtstn=None, instofst=None, trgtofst=None, refcoef=None, ddxyz=False, offsettype=None ):
        assert trgtstn is None
        xyz=self.xyz(instofst, offsettype=offsettype)
        if not ddxyz:
            return xyz
        return xyz, np.identity(3), None
