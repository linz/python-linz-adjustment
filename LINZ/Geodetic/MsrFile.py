'''
MsrFile: Module to read observations from MSR files used by dynanet
'''

import sys
import os.path
from Observation import Observation, ObservationValue

typemap={
    'D':'HA',
    'S':'SD',
    'V':'ZD',
    'L':'LV',
    'B':'AZ',
    }
heightedtype={ 'S':1, 'V':1 }
angletype={ 'V':1, 'D':1, 'B':1 }

def read( msrfile, skipErrors=False ):
    dtype=''
    stype=''
    inststn=''
    obs=None
    npending=0
    with open(msrfile) as msrf:
        nline=0
        for l in msrf:
            nline += 1
            location=msrfile+':'+str(nline)
            msrtype=l[0:1]
            if msrtype == '*':
                continue
            if msrtype not in typemap:
                if not skipErrors:
                    raise RuntimeError('Unknown type',msrtype,'at',location)
                continue
            obstype=typemap[msrtype]
            # inst=0
            parts=l[1:].split()
            insthgt=0.0
            trgthgt=0.0
            if msrtype in heightedtype:
                trgthgt=float(parts.pop())
                insthgt=float(parts.pop())
            if msrtype == 'D':
                if obs is None:
                    npending=int(parts[2])
                    parts[2:3]=[]
                elif parts[0] != inststn:
                    parts.insert(0,inststn)
            if msrtype in angletype:
                obsval=float(parts[2])+float(parts[3])/60.0+float(parts[4])/3600.0
                obserr=float(parts[-1])/3600.0
                
            else:
                obsval=float(parts[2])
                obserr=float(parts[-1])
            inststn=parts[0]
            trgtstn=parts[1]
            value=ObservationValue(inststn, trgtstn, obsval, obserr, insthgt, trgthgt )
            if obs is None:
                obs=Observation(obstype)
            obs.addObservation(value)
            if npending <= 0:
                yield obs
                obs=None
            else:
                npending -= 1

def main():
    import csv
    import argparse
    parser=argparse.ArgumentParser(description='Load data from MSR file to CSV')
    parser.add_argument('csv_file',help='Name of output CSV file')
    parser.add_argument('msr_files',nargs='+',help='Name(s) of MSR data files')
    parser.add_argument('-o','--overwrite',action='store_true',help='Overwrite CSV file')
    args=parser.parse_args()

    csv_file=args.csv_file
    msr_files=args.msr_files
    overwrite=args.overwrite

    obsdata='fromstn tostn obsset obstype fromhgt tohgt value error'

    if os.path.exists(csv_file) and not overwrite:
        print "Need -o option to overwrite existing CSV file"
        sys.exit()

    for f in msr_files:
        if not os.path.exists(f):
            print "MSR file "+f+" does not exist"
            sys.exit()

    print "Writing",csv_file
    with open(csv_file,'wb') as csvfh:
        csvf=csv.writer(csvfh)
        csvf.writerow(obsdata.split())

        lastsetid=0
        counts={}
        for msrfile in msr_files:
            for obs in read(msrfile):
                # print obs
                typecode=obs.obstype.code
                if typecode not in counts:
                    counts[typecode]=[0,0]
                counts[typecode][0]+=len(obs.obsvalues)
                counts[typecode][1]+=1

                if len(obs.obsvalues) > 1:
                    lastsetid += 1
                    setid=str(lastsetid)
                else:
                    setid=''
                for value in obs.obsvalues:
                    csvf.writerow((
                        value.inststn,
                        value.trgtstn,
                        setid,
                        typecode,
                        value.insthgt,
                        value.trgthgt,
                        value.value,
                        value.stderr))

    for typecode in sorted(counts.keys()):
        if typecode == 'HA':
            print "{0} {1} observations in {2} sets".format(typecode,*counts[typecode])
        else:
            print "{0} {1} observations".format(typecode,*counts[typecode])


if __name__ == '__main__':
    main()
