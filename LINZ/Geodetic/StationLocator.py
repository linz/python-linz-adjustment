#!/usr/bin/python

# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import math
from collections import namedtuple
import numpy as np

from LINZ.Geodetic.Ellipsoid import GRS80

from .Network import Network
from .Observation import Observation
from .Station import Station as NetworkStation

class StationLocator( object ):
    '''
    Class to approximately locate stations based on observations and add them to 
    a network.
    '''

    ObsDef=namedtuple('ObsDef','tostn obstype obsvalue obsset reverse')

    class AngleSet( object ):

        def __init__( self,id=None ):
            self.id=id
            self.orientation=None
            self.defined=False
            self.stations=set()

        def setOrientation( self, orientation ):
            self.orientation=orientation
            self.defined=True

        def addStation( self, station ):
            if station is not None:
                self.stations.add(station)

    class Station( object ):
        
        def __init__( self, locator, code, networkStation=None ):
            self.locator=locator
            self.code=code
            self.networkStation=None
            self.xyz=None if networkStation is None else networkStation.xyz()
            self.enu_axes=None
            self.obs={}
            self.trialXyz={}

        def addObs( self, tostn, obstype, obsvalue, obsset=None, reverse=False ):
            if tostn not in self.obs:
                self.obs[tostn]={}
            toobs=self.obs[tostn]
            if obstype not in toobs:
                toobs[obstype]=[]
            toobs[obstype].append( StationLocator.ObsDef(
                tostn, obstype, obsvalue, obsset, reverse ))

        def setGxLocation( self ):
            '''
            Use GX observations to define location
            '''
            if self.xyz is None:
                gxobs=self.obs.get(None,{}).get('GX')
                if gxobs is not None:
                    gxyz=np.array((0.0,0.0,0.0))
                    for gx in gxobs:
                        gxyz += gx.obsvalue.value
                    self.xyz=gxyz/len(gxobs)
                    self.locator.write("Fixing station {0} using GNSS coordinate observations\n"
                                  .format(self.code))

        def addTrialXyz( self, source, xyz ):
            self.trialXyz[source]=xyz

        def station( self ):
            return self.networkStation or NetworkStation(self.code,xyz=self.xyz)

        def removeUnfixableStations( self ):
            fixable={}
            for tostn,obs in self.obs.iteritems():
                if tostn is not None:
                    if 'SD' not in obs:
                        continue
                    if 'LV' not in obs and 'ZD' not in obs:
                        continue
                    if 'AZ' not in obs and 'HA' not in obs:
                        continue
                fixable[tostn]=obs
            self.obs=fixable

        def fixStation( self, xyz ):
            if self.xyz is None:
                self.locator.write("Fixing station {0}\n".format(self.code))
            self.xyz=xyz
            lon,lat,hgt=GRS80.geodetic(xyz)
            self.enu_axes=GRS80.enu_axes(lon,lat)
            # Fix any HA observation sets that can be defined
            hafixstations=set()

            for tostn in self.obs:
                if tostn is None or tostn.xyz is None:
                    continue
                obs=self.obs[tostn]
                if 'HA' in obs:
                    for haobs in obs['HA']:
                        if not haobs.obsset.defined:
                            stnfrom=self.station()
                            stnto=tostn.station()
                            if haobs.reverse:
                                stnfrom,stnto=stnto,stnfrom
                            az=stnfrom.azimuthTo(stnto)
                            fromcode=stnfrom.code()
                            self.locator.write("Fixing orientation of HA observations {1} at {0}\n"
                                               .format(fromcode,haobs.obsset.id))
                            referredcodes=[stn.code for stn in haobs.obsset.stations if stn.code != fromcode]
                            self.locator.write("   Provides azimuths to {0}\n".format(', '.join(referredcodes)))

                            haobs.obsset.setOrientation(az-haobs.obsvalue.value)
                            for stn in haobs.obsset.stations:
                                if stn != self and stn.xyz is not None:
                                    hafixstations.add(stn)

            # Now work out coordinates can be calculated
            for tostn in self.obs:
                if tostn is None or tostn.xyz is not None:
                    continue
                # self.locator.write("Attempting to coordinate {0}\n".format(tostn.code))
                obs=self.obs[tostn]

                # Work out an azimuth
                azimuths=[]
                if 'AZ' in obs:
                    azimuths=[o.obsvalue.value + (180.0 if o.reverse else 0) for o in obs['AZ']]
                if 'HA' in obs:
                    for o in obs['HA']:
                        if not o.obsset.defined:
                            continue
                        az=o.obsvalue.value+o.obsset.orientation
                        if o.reverse:
                            az += 180.0
                        azimuths.append(az)
                # If no horizontal angles defined yet, then can't set trial coord yet
                if len(azimuths) == 0:
                    # self.locator.write("  No azimuth data available\n")
                    continue
                azimuths=np.array(azimuths)
                azimuths *= math.radians(1.0)
                az0=azimuths[0]
                azimuths=np.remainder(np.array(azimuths)-az0+math.pi,2*math.pi)+az0-math.pi
                azimuth=np.mean(azimuths)

                if 'LV' in obs:
                    hgtdiff=np.mean([o.obsvalue.value for o in obs['LV']])
                # Only use ZD if don't have levelled height difference
                else:
                    distance=np.mean([o.obsvalue.value for o in obs['SD']])
                    hgtdiffs=[]
                    for o in obs['ZD']:
                        hd=math.cos(math.radians(o.obsvalue.value))*distance
                        hd+=o.obsvalue.insthgt
                        hd-=o.obsvalue.trgthgt
                        if o.reverse:
                            hd=-hd
                        hgtdiffs.append(hd)
                    hgtdiff=np.mean(hgtdiffs)

                hordists=[]
                for o in obs['SD']:
                    vdist=hgtdiff
                    if o.reverse:
                        vdist=-vdist
                    vdist += (o.obsvalue.trgthgt-o.obsvalue.insthgt)
                    hdist=o.obsvalue.value*o.obsvalue.value-vdist*vdist
                    if hdist > 0.0:
                        hordists.append(math.sqrt(hdist))
                    else:
                        hordists.append(0.0)
                hordist=np.mean(hordists)

                denu=np.array((hordist*math.sin(azimuth),hordist*math.cos(azimuth),hgtdiff))
                dxyz=self.enu_axes.T.dot(denu)
                tostn.addTrialXyz(self,self.xyz+dxyz)
                self.locator.write("  Trial coordinate determined for {0}: {1}\n".format(tostn.code,str(xyz)))

            # Finally refix any stations 
            for stn in hafixstations:
                stn.fixStation( stn.xyz )


    def __init__( self, network, observations, write=None ):
        self.network=network
        self.observations=observations
        # Replace write with a null function if not defined
        if write is None:
            write=lambda x: x
        self.write=write
        self.nupdated=0
        self.stations={}
        self.obs={}
        self.unlocated=[]
        self.enu_axes=np.identity(3)
        self.buildObservationIndex()
        self.locateInitialStations()
        while len(self.unlocated) > 0:
            self.tryLocateStation()
        self.updateNetwork()

    def getStation( self, code ):
        if code is None:
            return None
        if code not in self.stations:
            networkStation=self.network.get(code)
            self.stations[code]=StationLocator.Station(self,code,self.network.get(code))
        return self.stations[code]

    def buildObservationIndex( self ):
        setid=0
        counts={}
        for obs in self.observations:
            typecode=obs.obstype.code
            if typecode not in counts:
                counts[typecode]=[0,0]
            set=None
            if typecode == 'HA':
                setid+=1
                set=StationLocator.AngleSet(setid)
            counts[typecode][0] += len(obs.obsvalues)
            counts[typecode][1] += 1
            for obsval in obs.obsvalues:
                inststn=self.getStation(obsval.inststn)
                trgtstn=self.getStation(obsval.trgtstn)
                inststn.addObs( trgtstn, typecode, obsval, obsset=set )
                if obsval.trgtstn is not None:
                    trgtstn.addObs( inststn, typecode, obsval, obsset=set, reverse=True )
                if set is not None:
                    set.addStation(inststn)
                    set.addStation(trgtstn)
        for s in self.stations.values():
            s.removeUnfixableStations()
        self.write("Station location observation summary:\n")
        for typecode in sorted(counts.keys()):
            if typecode == 'HA':
                self.write("  {1} {0} observations in {2} sets\n".format(typecode,*counts[typecode]))
            else:
                self.write("  {1} {0} observations\n".format(typecode,counts[typecode][0]))


    def locateInitialStations( self ):
        xyz=np.array((0.0,0.0,0.0))
        nxyz=0
        for s in self.stations.values():
            s.setGxLocation()
            if s.xyz is not None:
                xyz += s.xyz
                nxyz += 1
            else:
                self.unlocated.append(s)
        if nxyz == 0:
            raise RuntimeError('No station data or observations to locate network')
        xyz /= nxyz

        # Define a local topocentric vector
        lon,lat,hgt=GRS80.geodetic(xyz)
        self.enu_axes=GRS80.enu_axes(lon,lat)

        # Fix initial stations to orient horizontal angles
        for s in self.stations.values():
            if s.xyz is not None:
                s.fixStation( s.xyz )

    def tryLocateStation( self ):
        maxtrialcoords=0
        maxtrialstn=None
        for s in self.unlocated:
            ntrial=len(s.trialXyz)
            if ntrial > maxtrialcoords:
                maxtrialcoords=ntrial
                maxtrialstn=s
        if not maxtrialstn:
            raise RuntimeError("Cannot locate {0} remaining stations"
                               .format(len(self.unlocated)))
        xyz=np.array([0.0,0.0,0.0])
        for trialxyz in maxtrialstn.trialXyz.values():
            xyz += trialxyz
        xyz /= maxtrialcoords
        maxtrialstn.fixStation(xyz)
        self.write("   {0} determined from {1} trial coords as {2}\n"
                   .format(maxtrialstn.code,maxtrialcoords,str(xyz)))
        self.unlocated.remove(maxtrialstn)

    def updateNetwork( self ):
        nstations=0
        for s in self.stations.values():
            if s.networkStation is None:
                s.networkStation=NetworkStation(s.code,xyz=s.xyz)
                self.network.addStation( s.networkStation )
                nstations += 1
        self.nupdated += nstations

def locateStations( network, observations, write=None ):
    locator=StationLocator(network,observations,write)
    return locator.nupdated

def main():
    import sys
    import re
    import argparse
    '''
    Program to test StationLocator
    '''
    parser=argparse.ArgumentParser(description='Calculate approximate station coordinates')
    parser.add_argument('coord_file',help="Input coordinate file")
    parser.add_argument('data_files',nargs='*',help="Data files")
    parser.add_argument('-n','--no-input-coords',action='store_true',help='Input coordinate file does not exist')
    parser.add_argument('-g','--geodetic',action='store_true',help='Output geodetic coordinate (ie as lon lat hgt)')
    args=parser.parse_args()

    net=Network()
    if not args.no_input_coords:
        net.loadCsv(args.coord_file)

    observations=[]
    for obsfile in args.data_files:
        if obsfile.lower().endswith('.msr'):
            import MsrFile
            reader=MsrFile.read
        elif re.search(r'\.snx(\.gz)?$',obsfile,re.I):
            import SinexObsFile
            reader=SinexObsFile.read
        else:
            import CsvObsFile
            reader=CsvObsFile.read
        observations.extend(reader(obsfile))

    StationLocator(net,observations,sys.stdout.write)
    net.writeCsv(args.coord_file+'.new',geodetic=args.geodetic)

if __name__ == "__main__":
    main()


