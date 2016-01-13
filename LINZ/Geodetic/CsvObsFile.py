# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import re
import math
from .CsvRecord import Reader
from .Observation import Observation, ObservationValue

approxEarthRadius=6371000.0

def read( csvfile, colnames=None, attributes=None ):
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

    Default column names can be overridden using the colnames parameter.  
    This is a dictionary mapping input column names to required column names.

    Additionally attributes can be specified as a list of names or a string
    separated by non-word characters.
    Each which will be defined as observation value attributes.  
    Observations will also have an attribute source containing file name 
    and line number.  Attributes can be defined with 
    type and optionality in the same way as CsvRecord columns.

    Each row can contain multiple observations usign the xx_value, 
    xx_error format, or one row entered using the obstype, value, error.

    '''

    columns='''
         fromstn tostn date? 
         fromhgt:float? tohgt:float?
         obstype? obsset? 
         value:float? error:float?
         ha_value:float? ha_error:float?
         hd_value:float? hd_error:float?
         sd_value:float? sd_error:float?
         zd_value:float? zd_error:float?
         az_value:float? az_error:float?
         lv_value:float? lv_error:float?
         '''

    attributes=attributes or []
    if isinstance(attributes,basestring):
        attributes=[x for x in re.split(r'\W+',attributes) if len(x) > 0]
    columns=columns+' '+' '.join(attributes)
    attfuncs={}
    for a in attributes:
        a=re.sub(r'[:?].*','',a)
        attfuncs[a]=eval('lambda x: x.'+a)

    csvreader=Reader('observation',columns)
    csvsrc=csvreader.open(csvfile,colnames)
    fields=csvsrc.fieldsDefined()
    valuefuncs=[]
    for collist in (
        'value error obstype',
        'ha_value ha_error',
        'hd_value hd_error',
        'sd_value sd_error',
        'zd_value zd_error',
        'az_value az_error',
        'lv_value lv_error'):
        names=collist.split()
        if names[0] in fields:
            for needed in names[1:]:
                if needed not in fields:
                    raise RuntimeError("Cannot have column {0} without {1} in {2}"
                                       .format(names[0],needed,csvsrc.filename()))
            if names[0][2] == '_':
                obstype="'"+names[0][:2].upper()+"'"
            else:
                obstype="x.obstype"
            funcstr="lambda x: ("+obstype+", x."+names[0]+", x."+names[1]+")"
            valuefuncs.append(eval(funcstr))

    haobs=None
    lastfrom=None
    lastset=None

    recordno=0
    filename=csvsrc.filename()

    for row in csvsrc.records():
        # Send pending HA observations
        if haobs is not None and (row.fromstn != lastfrom or row.obsset != lastset):
            yield(haobs)
            haobs=None
        lastfrom=row.fromstn
        lastset=row.obsset
        recordno += 1
        attributes={'source': filename+':'+str(recordno)}
        for k,f in attfuncs.iteritems():
            attributes[k]=f(row)

        for f in valuefuncs:
            obstype,obsvalue,obserror=f(row)
            if obsvalue is None:
                continue
            value=ObservationValue(row.fromstn,row.tostn,obsvalue,obserror,
                                   row.fromhgt or 0.0, row.tohgt or 0.0,
                                  attributes)
            if obstype == 'HA':
                if haobs is None:
                    haobs=Observation('HA')
                haobs.addObservation(value)
            else:
                yield Observation(obstype).addObservation(value)
                continue

    if haobs is not None:
        yield haobs

def readConvertToHorDist( csvfile, colnames=None, attributes=None ):
    global approxEarthRadius
    savedsd=None
    lastsource=None
    for obs in read( csvfile, colnames, attributes ):
        type=obs.obstype.code
        source=obs.obsvalues[0].attributes.get('source')
        if source != lastsource:
            lastsource=source
            if savedsd is not None:
                yield savedsd
                savedsd=None
        if type == 'SD' and savedsd is None:
            savedsd = obs
        elif type in ('ZD','LV') and savedsd is not None:
            sdist=savedsd.obsvalues[0].value
            hdiff=obs.obsvalues[0].value
            if type == 'ZD':
                zdist=math.radians(hdiff)
                subtend=sdist*math.sin(zdist)/(2*approxEarthRadius)
                hdiff=sdist*math.cos(zdist+subtend)
                hdifferr=math.hypot(math.radians(obs.obsvalues[0].stderr)*sdist,
                                    savedsd.obsvalues[0].stderr*subtend)
                obs.obsvalues[0].value=hdiff
                obs.obsvalues[0].stderr=hdifferr
                obs.obstype=Observation.ObservationTypes['LV']
            yield obs
            hdist=math.sqrt(sdist*sdist-hdiff*hdiff)
            subtend=hdist/(2*approxEarthRadius)
            hdist += hdiff*subtend
            savedsd.obsvalues[0].value=hdist
            savedsd.obstype=Observation.ObservationTypes['HD']
            yield savedsd
            savedsd=None
        else:
            yield obs
    if savedsd is not None:
        yield savedsd
