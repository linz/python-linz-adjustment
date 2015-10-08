#!/usr/bin/python

# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

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


    def recordType( self ):
        if len(self._fields) < 0:
            raise RuntimeError("Cannot read CSV file without defining fields first")
        return namedtuple(self._recordName, [f.name for f in self._fields])

    def fields( self ):
        return self._fields

    def open( self, source, colnames={}, filename=None ):
        '''
        Returns an open Reader.File object.  
        '''
        return Reader.File( self, source, colnames, filename )

    def readCsv( self, source, colnames={}, filename=None ):
        '''
        Reads a CSV file and yield a record for each row 
        corresponding the defined fields.  Assumes the first
        record defines column names.

        if colnames is supplied it maps from the input column names to the
        field names, ie
          { colname: fieldname, ... }
        '''

        for record in self.open(source,colnames,filename).records():
            yield record

    class File( object ):
        '''
        Reader.File opens a file using the reader record definition
        in readiness for extracting records.
        '''
    
        def __init__( self, csvrecord, source, colnames=None, filename=None ):
            if isinstance(source,basestring):
                if filename is None:
                    filename=source
                source=open(source,'rb')
            if filename is None:
                filename='input'
            self._filename=filename
            self._csvf=csv.reader(source)
            header=[f.lower() for f in self._csvf.next()]
            colnos=[];
            self._fields=csvrecord.fields()
            headercolumns=header
            if colnames is not None:
                colnames={k.lower():v.lower() for k,v in colnames.iteritems()}
                headercolumns=[colnames.get(h,h) for h in header]

            for i,f in enumerate(self._fields):
                colno=-1
                try:
                    colno=headercolumns.index(f.name)
                except:
                    if not f.optional:
                        raise RuntimeError(filename+' is missing required field '+filename)
                colnos.append(colno)
            self._header=header
            self._colnos=colnos
            self._recordtype=csvrecord.recordType()

        def filename( self ):
            return self._filename

        def fieldsDefined( self ):
            '''
            Returns a list of record fields actually defined in the file
            '''
            return [f.name for c,f in zip(self._colnos,self._fields) if c >= 0]

        def records( self ):
            '''
            Iterator yielding each record from the file in turn
            '''
            recordno=0
            for row in self._csvf:
                recordno += 1
                record=[]
                for c,f in zip(self._colnos,self._fields):
                    if c < 0:
                        record.append(None)
                    else:
                        value=row[c].strip() if len(row) > c else ''
                        if f.parsefunc:
                            if value == '':
                                value= None
                            else:
                                try:
                                    value=f.parsefunc(value)
                                except:
                                    raise ValueError("Invalid value \""+value+"\" for "+
                                                     self._header[c]+" in "+self._filename+
                                                     " at record "+str(recordno))
                        elif value == '' and f.optional:
                            value=None
                        if value is None and not f.optional:
                            raise ValueError("Missing value for "+
                                self._header[c]+" in "+filename+
                                " at record "+str(recordno))
                        record.append(value)
                yield self._recordtype(*record)
