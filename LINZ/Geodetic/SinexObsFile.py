"""
SinexObsFile: module to read coordinate observations from SINEX file
"""

import sys
import os.path
from datetime import timedelta
from LINZ.Geodetic.Sinex import Reader as SinexReader, COVAR_FULL

from .Observation import Observation, ObservationValue


def read(sinexFile, useMonument=False):
    snx = SinexReader(sinexFile, velocities=False, covariance=COVAR_FULL)

    ids = {}
    for id, code, soln in snx.solutions():
        key = (id, code) if useMonument else (id, None)
        if key in ids:
            raise RuntimeError(
                "SINEX observation files do not support multiple solutions for"
                + str(id)
            )
        ids[key] = 1

    ids = sorted(ids.keys())
    covariance = None
    obsvalues = {}
    startdate = None
    enddate = None
    for id in ids:
        data = snx.get(ptid=id[0], ptcode=id[1])
        prmid0 = data.prmids[0]
        code = data.monument if useMonument else data.id
        obsvalues[prmid0] = ObservationValue(code, value=data.xyz)
        covariance = data.covariance
        if startdate is None or data.crddate < startdate:
            startdate = data.crddate
        if enddate is None or data.crddate > enddate:
            enddate = data.crddate

    obsdate = startdate + timedelta(seconds=(enddate - startdate).total_seconds() / 2.0)
    observation = Observation("GX", obsdate=obsdate)

    for prmid in sorted(obsvalues.keys()):
        observation.addObservation(obsvalues[prmid])
    observation.setCovariance(covariance)
    yield observation


def main():
    import csv
    import argparse

    parser = argparse.ArgumentParser(description="Test load data from SINEX file")
    parser.add_argument("snx_file", help="Name of SINEX file")
    args = parser.parse_args()

    for obs in read(args.snx_file):
        print(obs)


if __name__ == "__main__":
    main()
