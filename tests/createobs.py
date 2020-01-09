# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os.path
import numpy as np

sys.path.append(
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "LINZ")
)
from Geodetic.Observation import Observation as Obs, ObservationValue as ObsVal


def HA(
    stfrom, targets, offset=0.0, error=0.001, insthgt=0.0, trgthgt=0.0, obsdate=None
):
    ha = Obs("HA", obsdate=obsdate)
    refha = None
    for stto in targets:
        angle = stfrom.azimuthTo(stto, instofst=insthgt, trgtofst=trgthgt)
        if refha is None:
            refha = angle - offset
        angle = angle - refha
        if angle < 0:
            angle += 360.0
        ha.addObservation(
            ObsVal(
                stfrom.code(),
                stto.code(),
                angle,
                error,
                insthgt=insthgt,
                trgthgt=trgthgt,
            )
        )
    return ha


def AZ(stfrom, stto, error=0.001, insthgt=0.0, trgthgt=0.0, obsdate=None):
    angle = stfrom.azimuthTo(stto, instofst=insthgt, trgtofst=trgthgt)
    return Obs(
        "AZ",
        obsdate=obsdate,
        obsvalue=ObsVal(
            stfrom.code(), stto.code(), angle, error, insthgt=insthgt, trgthgt=trgthgt
        ),
    )


def SD(stfrom, stto, error=0.001, insthgt=0.0, trgthgt=0.0, obsdate=None):
    angle = stfrom.distanceTo(stto, instofst=insthgt, trgtofst=trgthgt)
    return Obs(
        "SD",
        obsdate=obsdate,
        obsvalue=ObsVal(
            stfrom.code(), stto.code(), angle, error, insthgt=insthgt, trgthgt=trgthgt
        ),
    )


def ZD(stfrom, stto, error=0.001, insthgt=0.0, trgthgt=0.0, refcoef=0.0, obsdate=None):
    angle = stfrom.zenithDistanceTo(
        stto, instofst=insthgt, trgtofst=trgthgt, refcoef=refcoef
    )
    return Obs(
        "ZD",
        obsdate=obsdate,
        obsvalue=ObsVal(
            stfrom.code(), stto.code(), angle, error, insthgt=insthgt, trgthgt=trgthgt
        ),
    )


def LV(stfrom, stto, error=0.001, insthgt=0.0, trgthgt=0.0, obsdate=None):
    angle = stfrom.heightDifferenceTo(stto, instofst=insthgt, trgtofst=trgthgt)
    return Obs(
        "LV",
        obsdate=obsdate,
        obsvalue=ObsVal(
            stfrom.code(), stto.code(), angle, error, insthgt=insthgt, trgthgt=trgthgt
        ),
    )


def GB(stfrom, stto, error=0.001, insthgt=0.0, trgthgt=0.0, obsdate=None):
    xyz0 = stfrom.xyz(offset=insthgt)
    xyz2 = stto.xyz(offset=trgthgt)
    dxyz = xyz2 - xyz0
    covar = np.identity(3) * error * error
    obs = Obs(
        "GB",
        obsdate=obsdate,
        obsvalue=ObsVal(
            stfrom.code(), stto.code(), value=dxyz, insthgt=insthgt, trgthgt=trgthgt
        ),
    )
    obs.setCovariance(covar)
    return obs


def GX(stfrom, error=0.001, insthgt=0.0, trgthgt=0.0, obsdate=None):
    xyz = stfrom.xyz(offset=insthgt)
    covar = np.identity(3) * error * error
    obs = Obs(
        "GX",
        obsdate=obsdate,
        obsvalue=ObsVal(stfrom.code(), value=xyz, insthgt=insthgt),
    )
    obs.setCovariance(covar)
    return obs


def Traverse(stations):
    obs = []
    for back, base, fore in zip(stations[:-2], stations[1:-1], stations[2:]):
        obs.append(HA(base, [back, fore]))
        obs.append(SD(base, back))
        obs.append(SD(base, fore))
        obs.append(ZD(base, back))
        obs.append(ZD(base, fore))
    return obs
