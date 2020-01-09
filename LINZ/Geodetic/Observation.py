from collections import namedtuple
from .Station import Station


class ObservationValue(object):
    """
    Represents a single observation value within a set. Has members

      inststn
      trgtstn
      value
      stderr
      insthgt
      trgthgt
      attributes

    Note that stderr may be overridden by the covariance of the set
    """

    def __init__(
        self,
        inststn,
        trgtstn=None,
        value=0.0,
        stderr=None,
        insthgt=0.0,
        trgthgt=0.0,
        attributes={},
    ):
        self.inststn = inststn
        self.trgtstn = trgtstn
        self.value = value
        self.stderr = stderr
        self.insthgt = insthgt
        self.trgthgt = trgthgt
        self.attributes = attributes

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if isinstance(self.value, float):
            strobs = "{0:.4f}".format(self.value)
        else:
            strobs = "[" + ",".join(["{0:.5f}".format(x) for x in self.value]) + "]"

        if self.stderr is None:
            strerr = "None"
        elif isinstance(self.stderr, float):
            strerr = "{0:.4f}".format(self.stderr)
        else:
            strerr = "[" + ",".join(["{0:.5f}".format(x) for x in self.stderr]) + "]"

        strvalue = "{0},{1},{2:.4f},{3:.4f},{4},{5}".format(
            self.inststn, self.trgtstn, self.insthgt, self.trgthgt, strobs, strerr
        )
        delim = ","
        if self.attributes:
            attr = self.attributes
            strvalue = (
                strvalue
                + ","
                + " ".join(str(k) + "=" + str(attr[k]) for k in sorted(attr.keys()))
            )
        return strvalue


class Observation(object):

    """
    Class representing an observation

    Each observation is defined by 
       obstype   an observation type (obstype)
       obsdate   a observation date/time
       obsvalues an array of observation values with elements
         inststn    from station code
         insthgt    from station height
         trgtstn    to station (may be None if not appropriate to type)
         trgthgt    to station height
         value      observation value (single value or array)
         stderr     standard error of observation
         attributes dictionary of custom attributes
       covariance   the observation matrix 

    The covariance, if defined, overrides the individual stderr.
    """

    ObservationType = namedtuple("ObservationType", "code name nvalue calcobs")

    ObservationTypes = {
        "HA": ObservationType("HA", "Horizontal angle", 1, Station.azimuthTo),
        "AZ": ObservationType("AZ", "Azimuth", 1, Station.geodeticAzimuthTo),
        "SD": ObservationType("SD", "Slope distance", 1, Station.distanceTo),
        "HD": ObservationType(
            "HD", "Horizontal distance", 1, Station.horizontalDistanceTo
        ),
        "ZD": ObservationType("ZD", "Zenith distance", 1, Station.zenithDistanceTo),
        "LV": ObservationType("LV", "Height difference", 1, Station.heightDifferenceTo),
        "GX": ObservationType("GX", "XYZ coordinate", 3, Station.calcXYZ),
        "GB": ObservationType("GB", "XYZ baseline", 3, Station.vectorTo),
    }

    def __init__(self, obstype, obsdate=None, obsvalue=None, covariance=None):
        if obstype not in Observation.ObservationTypes:
            raise ValueError("Invalid observation type: " + obstype)
        self.obstype = Observation.ObservationTypes[obstype]
        self.obsdate = obsdate
        self.obsvalues = (
            [obsvalue] if isinstance(obsvalue, ObservationValue) else (obsvalue or [])
        )
        self.covariance = covariance

    def addObservation(self, obs):
        self.obsvalues.append(obs)
        return self

    def setCovariance(self, covar):
        self.covariance = covar
        return self

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "".join(
            [self.obstype.code + "," + str(o) + "\n" for o in self.obsvalues]
        )
