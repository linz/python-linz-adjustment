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

MeanEarthRadius=6371000.0

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
            self.instStation=None
            self.defined=False
            self.stations=set()

        def setOrientation( self, orientation ):
            self.orientation=orientation
            self.defined=True

        def addStation( self, station, instStation=False ):
            if station is not None:
                self.stations.add(station)
                if instStation:
                    self.instStation=station


    class Station( object ):

        # Station class.
        #
        # Holds observations and coordinates relating to a station.
        # code is the code for the station
        # locator is the parent StationLocator object
        # networkStation is the station definition from the network
        # xyz is the set xyz coordinate
        # trialXyz is a dictionary of trial coordinates keyed on the station from
        #   which they are calculated
        # enu_axes is the XYZ-ENU rotation matrix
        # obs is a list of observations.  Obs is a dictionary keyed as
        #   ->tostation->obstype->[ObsDef]
        
        def __init__( self, locator, code, networkStation=None ):
            self.locator=locator
            self.code=code
            self.networkStation=networkStation
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

        # Not right ...  removes stations which are useful for orienting HA
        # observations...
        #
        #def removeUnfixableStations( self ):
        #    fixable={}
        #    for tostn,obs in self.obs.iteritems():
        #        if tostn is not None:
        #            if 'SD' not in obs:
        #                continue
        #            if 'LV' not in obs and 'ZD' not in obs:
        #                continue
        #            if 'AZ' not in obs and 'HA' not in obs:
        #                continue
        #        fixable[tostn]=obs
        #    self.obs=fixable

        def setXyz( self, xyz ):
            self.xyz=np.array(xyz)
            lon,lat,hgt=GRS80.geodetic(xyz)
            self.enu_axes=GRS80.enu_axes(lon,lat)

        def fixStation( self, xyz ):
            if self.xyz is None or True:
                self.locator.write("Fixing station {0}\n".format(self.code))
            self.setXyz( xyz )
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
                            if haobs.reverse:
                                hafixstations.add(tostn)

            # Now work out coordinates can be calculated
            for tostn in self.obs:
                if tostn is None or tostn.xyz is not None:
                    continue
                # self.locator.write("Attempting to coordinate {0}\n".format(tostn.code))
                obs=self.obs[tostn]
                if 'SD' not in obs:
                    continue

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

                llh0=GRS80.geodetic(self.xyz)
                hgtdiff=0.0
                if 'LV' in obs:
                    hgtdiff=np.mean([o.obsvalue.value for o in obs['LV']])
                # Only use ZD if don't have levelled height difference
                elif 'ZD' in obs:
                    distance=np.mean([o.obsvalue.value for o in obs['SD']])
                    hgtdiffs=[]
                    for o in obs['ZD']:
                        angle=math.radians(o.obsvalue.value)
                        sd=math.sin(angle)*distance
                        corr=sd/(2.0*(llh0[2]+MeanEarthRadius))
                        angle -= corr
                        hd=math.cos(angle)*distance
                        hd+=o.obsvalue.insthgt
                        hd-=o.obsvalue.trgthgt
                        if o.reverse:
                            hd=-hd
                        hgtdiffs.append(hd)
                    hgtdiff=np.mean(hgtdiffs)
                if hgtdiff is None:
                    continue

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
                txyz=self.xyz+dxyz
                llh1=GRS80.geodetic(txyz)
                llh1[2]=llh0[2]+hgtdiff
                txyz=GRS80.xyz(llh1)
                tostn.addTrialXyz(self,txyz)
                self.locator.write("  Trial coordinate determined for {0}: {1}\n".format(tostn.code,str(txyz)))

            # Finally refix any stations 
            for stn in hafixstations:
                stn.fixStation( stn.xyz )


    def __init__( self, network, observations, write=None ):
        self.network=network
        self.observations=observations
        # Replace write with a null function if not defined
        if write is None:
            write=lambda x: None
        self.write=write
        self.write("Attempting to locate missing stations\n")
        self.nupdated=0
        self.stations={}
        self.sets=[]
        self.obs={}
        self.unlocated=[]
        self.enu_axes=np.identity(3)
        self.buildObservationIndex()
        self.locateInitialStations()
        nunlocated=len( self.unlocated )
        if len(self.unlocated) > 0:
            while self.tryLocateStation() or self.tryFixAngle():
                # Should not be true - check to avoid infinite loop...
                if len(self.unlocated) >= nunlocated:
                    raise RuntimeError("tryLocateStation or tryFixAngle failed to fix stations")
                nunlocated=len( self.unlocated )
                if nunlocated == 0:
                    break
            if len(self.unlocated) > 0:
                self.write("\nCannot locate {0} remaining stations:\n".format(len(self.unlocated)));
                for code in sorted([s.code for s in self.unlocated]):
                    self.write("  {0}\n".format(code))
                #raise RuntimeError("Cannot locate {0} remaining stations"
                #                   .format(len(self.unlocated)))
        self.updateNetwork()

    def getStation( self, code ):
        if code is None:
            return None
        if code not in self.stations:
            networkStation=self.network.get(code)
            self.stations[code]=StationLocator.Station(self,code,networkStation)
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
                self.sets.append(set)
            counts[typecode][0] += len(obs.obsvalues)
            counts[typecode][1] += 1
            for obsval in obs.obsvalues:
                inststn=self.getStation(obsval.inststn)
                trgtstn=self.getStation(obsval.trgtstn)
                inststn.addObs( trgtstn, typecode, obsval, obsset=set )
                if trgtstn is not None:
                    trgtstn.addObs( inststn, typecode, obsval, obsset=set, reverse=True )
                if set is not None:
                    set.addStation(inststn,True)
                    set.addStation(trgtstn)
        #for s in self.stations.values():
        #    s.removeUnfixableStations()
        self.write("Station location observation summary:\n")
        for typecode in sorted(counts.keys()):
            if typecode == 'HA':
                self.write("  {1} {0} observations in {2} sets\n".format(typecode,*counts[typecode]))
            else:
                self.write("  {1} {0} observations\n".format(typecode,counts[typecode][0]))

    def resetUnlocatedStations( self ):
        self.unlocated=[s for s in self.stations.values() if s.xyz is None]

    def locateInitialStations( self ):
        xyz=np.array((0.0,0.0,0.0))
        nxyz=0
        for s in self.stations.values():
            s.setGxLocation()
            if s.xyz is not None:
                s.fixStation(s.xyz)
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
            return False
        xyz=np.array([0.0,0.0,0.0])
        for trialxyz in maxtrialstn.trialXyz.values():
            xyz += trialxyz
        xyz /= maxtrialcoords
        self.write("Determined {0} from {1} trial coords as {2}\n"
                   .format(maxtrialstn.code,maxtrialcoords,str(xyz)))
        maxtrialstn.fixStation(xyz)
        self.unlocated.remove(maxtrialstn)
        return True

    def getFixedStations( self ):
        fixed={}
        for s in self.stations.values():
            if s.xyz is not None:
                fixed[s]=s.xyz
        return fixed

    def clearFixedStations( self ):
        # Clear all fixed stations and trial coordinates
        for obsset in self.sets:
            obsset.defined=False
        for s in self.stations.values():
            s.xyz=None
            s.trialXyz={}

    def restoreFixedStations( self, fixed ):
        for s in self.stations.values():
            if s in fixed:
                s.fixStation(fixed[s])
            else:
                s.xyz=None
        self.resetUnlocatedStations()

    def unfixedAngles( self ):
        return [s for s in self.sets if not s.defined]

    def tryFixAngle( self ):
        tried=[]
        unfixedAngles=self.unfixedAngles()
        fixedStations=self.getFixedStations()
        success=False
        write=self.write
        self.write=lambda x: None
        write("\nTrying fixing one angle\n")
        try:
            while len(unfixedAngles) > 0:
                maxmatched=0
                matchedset=None
                fixedstation=None
                for haset in unfixedAngles:
                    hafixed=[s for s in haset.stations if s in fixedStations]
                    countfixed=len(hafixed)
                    if countfixed > maxmatched:
                        fixedstation=hafixed[0]
                        maxmatched=countfixed
                        matchedset=haset
                if matchedset is None:
                    break
                # Try fixing one angle
                write("\nFixing set at {0} connected to {1} stations\n"
                           .format(matchedset.instStation.code,maxmatched))
                self.clearFixedStations()
                matchedset.setOrientation(0.0)
                fixedstation.fixStation(fixedStations[fixedstation])
                self.unlocated=[s for s in self.stations.values() if s != fixedstation]
                nunlocated=len(self.unlocated)
                nunlocated0=nunlocated
                while self.tryLocateStation():
                    if len(self.unlocated) >= nunlocated:
                        break
                    nunlocated=len(self.unlocated)
                    if nunlocated == 0:
                        break
                newfixed=self.getFixedStations()
                commonStations=[s for s in newfixed if s in fixedStations]
                write("Fixed angle connects {0} stations - {1} have known coordinates\n"
                           .format(len(newfixed),len(commonStations)))
                # If only one common station then cannot orient new stations...
                # Remove unfixed stations and try again.
                if len(commonStations) < 2:
                    len0=len(unfixedAngles)
                    stillUnfixed=set()
                    for hasset in unfixedAngles:
                        if not haset.defined:
                            stillUnfixed.add(haset)
                    unfixedAngles=stillUnfixed
                    len1=len(unfixedAngles)
                    if len1 >= len0:
                        raise RuntimeError("Failed to remove unfixedAngles in tryFixAngle")
                    write("Failed to use fixed angle\n")
                    continue
                # Find rotation to apply to new stations...
                xyz0=fixedstation.xyz
                enu=fixedstation.enu_axes
                angleref=None
                sumangle=0.0
                sumds=0.0
                for s in commonStations:
                    if s == fixedstation:
                        continue
                    denu0=enu.dot(fixedStations[s]-xyz0)
                    denu1=enu.dot(newfixed[s]-xyz0)
                    angle0=math.atan2(denu0[0],denu0[1])
                    angle1=math.atan2(denu1[0],denu1[1])
                    angdif=angle1-angle0
                    if angleref is None:
                        angleref=angdif
                    angdif-=angleref
                    while angdif < -math.pi:
                        angdif += math.pi*2
                    while angdif > math.pi:
                        angdif -= math.pi*2
                    angdif+=angleref
                    ds=math.hypot(denu1[0],denu1[1])
                    sumangle+=angdif*ds
                    sumds+=ds
                angdif=sumangle/sumds
                # Now restore fixed stations, reset the fixed angle
                # and calculate stations
                self.clearFixedStations()
                matchedset.setOrientation(math.degrees(-angdif))
                self.restoreFixedStations(fixedStations)
                self.write=write
                nunlocated=len(self.unlocated)
                nunlocated0=nunlocated
                success=False
                while self.tryLocateStation():
                    if len(self.unlocated) >= nunlocated:
                        break
                    nunlocated=len(self.unlocated)
                    success=True
                    if nunlocated == 0:
                        break
                nunlocated=len(self.unlocated)
                if nunlocated >= nunlocated0:
                    raise RuntimeError("Failed to fix a  station in tryFixAngle")
                self.write("Successfully fixed {0} stations using fixed angle\n"
                           .format(nunlocated0-nunlocated))
                break

            if not success:
                self.clearFixedStations()
                self.restoreFixedStations(fixedStations)
                write("Restored {0} original fixed stations\n".format(len(fixedStations)))
        finally:
            self.write=write
        return success

    def updateNetwork( self ):
        nstations=0
        for s in self.stations.values():
            if s in self.unlocated:
                continue
            if s.networkStation is None:
                self.write("   Adding network station {0}\n".format(s.code))
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


