# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from .CsvRecord import Reader
from .Observation import Observation, ObservationValue

def read( csvfile, colnames=None ):
    csvreader=Reader('observation','''
                     fromstn tostn date? 
                     fromhgt:float? tohgt:float?
                     obstype obsset? 
                     value:float error:float
                     ''')

    obs=None
    lastset=None
    for row in csvreader.readCsv(csvfile,colnames):
        value=ObservationValue(row.fromstn,row.tostn,row.value,row.error,
                               row.fromhgt or 0.0, row.tohgt or 0.0)
        if obs is not None and (row.obstype != 'HA' or row.obsset != lastset):
            yield obs
            obs=None
        lastset=row.obsset
        if obs is None:
            obs=Observation(row.obstype)
        obs.addObservation(value)

    if obs is not None:
        yield obs
