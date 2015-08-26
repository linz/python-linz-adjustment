
# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from collections import namedtuple
from .Station import Station

class ObservationValue( object ):

    def __init__( self, inststn, trgtstn=None, value=0.0, stderr=None, insthgt=0.0, trgthgt=0.0 ):
        self.inststn=inststn
        self.trgtstn=trgtstn
        self.value=value
        self.stderr=stderr
        self.insthgt=insthgt
        self.trgthgt=trgthgt

class Observation( object ):

    '''
    Class representing an observation

    Each observation is defined by 
       an observation type (obstype)
       a observation date/time
       an array of observation values with elements
         from station code
         from station height
         to station (may be None if not appropriate to type)
         to station height
         value (single value or array)
         stderr (error of observation)
       a covariance matrix 

    The covariance, if defined, overrides the stderr.
    '''

    ObservationType=namedtuple('ObservationType','code name nvalue calcobs')

    ObservationTypes={
        'HA': ObservationType('HA','Horizontal angle',1,Station.azimuthTo),
        'AZ': ObservationType('AZ','Azimuth',1,Station.geodeticAzimuthTo),
        'SD': ObservationType('SD','Slope distance',1,Station.distanceTo),
        'ZD': ObservationType('ZD','Zenith distance',1,Station.zenithDistanceTo),
        'LV': ObservationType('LV','Height difference',1,Station.heightDifferenceTo),
        'GX': ObservationType('GX','XYZ coordinate',3,Station.calcXYZ),
        }

    def __init__( self, obstype, obsdate=None, obsvalue=None, covariance=None ):
        if obstype not in Observation.ObservationTypes:
            raise ValueError("Invalid observation type: "+obstype)
        self.obstype=Observation.ObservationTypes[obstype]
        self.obsdate=obsdate
        self.obsvalues=[obsvalue] if isinstance(obsvalue,ObservationValue) else (obsvalue or [])
        self.covariance=covariance

    def addObservation( self, obs ):
        self.obsvalues.append(obs)

    def setCovariance( self, covar ):
        self.covariance=covar

