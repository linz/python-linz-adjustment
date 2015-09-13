# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from .CsvRecord import Reader
from .Observation import Observation, ObservationValue

def read( csvfile, colnames=None ):
    '''
    Reads observations from a data file.  The data file can have columns:
    
       fromstn   Code of instrument station
       fromhgt   Height of instrument above mark
       tostn     Code of target station
       tohgt     Height of target above mark
       date      Date of observation (yyyy-mm-dd hh:mm:ss)
       obsset    Set id of observations (required for HA obserations to identify the set)
       obstype   Observation type code (eg HA, SD, ZD, LV, AZ ...)
       value     Observation value (metres, decimal degrees)
       error     Observtion error (metres, decimal degrees) 
       xx_value  value for xx type observation 
       xx_error  error for xx type observation

    Each row can contain multiple observations usign the xx_value, xx_error format,
    or one row entered using the obstype, value, error.

    '''
    csvreader=Reader('observation','''
                     fromstn tostn date? 
                     fromhgt:float? tohgt:float?
                     obstype? obsset? 
                     value:float? error:float?
                     ha_value:float? ha_error:float?
                     sd_value:float? sd_error:float?
                     zd_value:float? zd_error:float?
                     az_value:float? az_error:float?
                     lv_value:float? lv_error:float?
                     ''')

    csvfile=csvreader.open(csvfile,colnames)
    fields=csvfile.fieldsDefined()
    valuefuncs=[]
    for collist in (
        'value error obstype',
        'ha_value ha_error obsset',
        'sd_value sd_error',
        'zd_value zd_error',
        'az_value az_error',
        'lv_value lv_error'):
        names=collist.split()
        if names[0] in fields:
            for needed in names[1:]:
                if needed not in fields:
                    raise RuntimeError("Cannot have column {0} without {1} in {2}"
                                       .format(names[0],needed,csvfile))
            if names[0][2] == '_':
                obstype="'"+names[0][:2].upper()+"'"
            else:
                obstype="x.obstype"
            funcstr="lambda x: ("+obstype+", x."+names[0]+", x."+names[1]+")"
            valuefuncs.append(eval(funcstr))

    haobs=None
    lastfrom=None
    lastset=None

    for row in csvfile.records():
        # Send pending HA observations
        if haobs is not None and (row.fromstn != lastfrom or row.obsset != lastset):
            yield(haobs)
            haobs=None
        lastfrom=row.fromstn
        lastset=row.obsset

        for f in valuefuncs:
            obstype,obsvalue,obserror=f(row)
            if obsvalue is None:
                continue
            value=ObservationValue(row.fromstn,row.tostn,obsvalue,obserror,
                                   row.fromhgt or 0.0, row.tohgt or 0.0)
            if obstype == 'HA':
                if haobs is None:
                    haobs=Observation('HA')
                haobs.addObservation(value)
            else:
                yield Observation(obstype).addObservation(value)
                continue

    if haobs is not None:
        yield haobs
