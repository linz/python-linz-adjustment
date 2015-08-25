#!/usr/bin/python

import sys
from collections import namedtuple
import csv
import re

# Class to read a CSV file into typed records.

class Reader( object ):
    
    '''
    CsvRecords reads records from a CSV file.  This allows for parsing fields
    into data types.  Also allows mapping of column names to field names
    '''

    Field=namedtuple('Field','name parsefunc optional')

    def __init__( self, name='record', fields=None ):
        '''
        name is the name given to the namedtuple return records
        fields are a list of fields as [name, parse function, optional]
        eg
           [name, None, False],
           [value,float,True],

        Alternatively can be supplied as a string
           field:type? field:type?
        where:type specifies the type, float or int currently supported,
        and ? implies optional.
        '''
        self._fields=[]
        self._recordName=name
        if fields is not None:
            self.addFields(fields)

    def addFields( self, fields ):
        '''
        Add one or more field definitions to the list. Fields can 
        be defined as in the constructor.
        '''

        if isinstance( fields, basestring ):
            fielddefs=fields.split()
            fields=[]
            for fielddef in fielddefs:
                m = re.match(r'^(\w+)(?:\:(float|int))?(\?)?$',fielddef)
                if not m:
                    raise ValueError("Invalid value "+fielddef+" for field definition")
                fieldname=m.group(1)
                parsefunc={'': None, 'float':float,'int':int}[m.group(2) or '']
                optional=m.group(3)=='?'
                fields.append((fieldname,parsefunc,optional))
        for f in fields:
            self.addField( *f )

    def addField( self, field, parsefunc=None, optional=None ):
        '''
        Add a field to the record type.  Field defined by field name,
        and optional a parsefunc (defines how to convert string to field
        type) and optional status (boolean)
        '''
        self._fields.append(Reader.Field(field,parsefunc,optional))

    def readCsv( self, filename, colnames=None ):
        '''
        Reads a CSV file and yield a record for each row 
        corresponding the defined fields.  Assumes the first
        record defines column names.

        if colnames is supplied it maps from the field names to the
        column names, ie
          { fieldname: colname, ... }
        '''

        if len(self._fields) < 0:
            raise RuntimeError("Cannot read CSV file without defining fields first")

        recordtype=namedtuple(self._recordName, [f.name for f in self._fields])

        with open(filename,'rb') as csvfh:
            csvf=csv.reader(csvfh)
            header=[f.lower() for f in csvf.next()]
            colnos=[];
            colnames=colnames or {}
            for i,f in enumerate(self._fields):
                fname=colnames(f.name) if f.name in colnames else f.name
                colno=-1
                try:
                    colno=header.index(fname.lower())
                except:
                    if not f.optional:
                        raise RuntimeError(filename+' is missing required field '+fname)
                colnos.append(colno)
            recordno=0
            for row in csvf:
                recordno += 1
                record=[]
                for c,f in zip(colnos,self._fields):
                    if c < 0:
                        record.append(None)
                    else:
                        value=row[c] if len(row) > c else None
                        if f.parsefunc:
                            if value == '':
                                value= None
                            else:
                                try:
                                    value=f.parsefunc(value)
                                except:
                                    raise ValueError("Invalid value \""+row[c]+"\" for "+
                                                     header[c]+" in "+filename+
                                                     " at record "+str(recordno))
                        elif value == '' and f.optional:
                            value=None
                        if value is None and not f.optional:
                            raise ValueError("Missing value for "+
                                header[c]+" in "+filename+
                                " at record "+str(recordno))
                        record.append(value)
                yield recordtype(*record)
