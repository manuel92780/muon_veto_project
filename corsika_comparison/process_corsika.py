#!/usr/bin/env python
import os, sys, argparse, time

#import icecube stuff
from I3Tray import I3Tray
from icecube import icetray, dataio, dataclasses, phys_services, simclasses, VHESelfVeto, MuonGun, weighting
from icecube.weighting import weighting, get_weighted_primary
from icecube.weighting.fluxes import GaisserH4a
from icecube.weighting.fluxes import Hoerandel5
from icecube.weighting import CORSIKAWeightCalculator
from icecube.weighting import fluxes
from icecube.icetray import I3Units
from icecube.MuonGun import Cylinder
from icecube.simprod import segments

start_time = time.asctime()
print 'Started:', start_time

def muonstodet(frame, surface):
    detmu = MuonGun.muons_at_surface(frame, surface)
    for i in range(len(detmu)):
        frame['EnteringMuon_'+str(i)] = detmu[i]

def primarytodet(frame, surface):
    primary = frame['MCPrimary']
    intersections = surface.intersection(primary.pos, primary.dir)

def printer(frame):
    print frame
    
def header(frame):
    frame['I3EventHeader'] = dataclasses.I3EventHeader()

def getweights(frame, generator, flux):
    get_weighted_primary(frame)
    energy = frame['MCPrimary'].energy
    ptype = frame['MCPrimary'].type
    weight = (flux(energy,ptype)/generator(energy,ptype))        
    frame['Weight'] = dataclasses.I3Double(weight)

def main():
    
    #arguement parser
    parser = argparse.ArgumentParser(description='Generates muons in IceCube with varying multiplicity')
    parser.add_argument('inp', nargs='+')
    parser.add_argument('--nseed', default=1, type=int,
                        help='seed for randomization')
    parser.add_argument('--out', default='corsika', help='Output file')
    parser.add_argument('--runnum', default=1, type=int,
                        help='Run number for this sim')
    parser.add_argument('--did', default=20209, type=int,
                        help='Simpord Number for this sim')
    #detector args
    parser.add_argument('--gcd', default=os.path.expandvars('$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'),
                        type=str,
                        help='gcd file')
    args = parser.parse_args()

    #geometry parameters
    gcdFile = args.gcd
    surface = Cylinder(1600*I3Units.m, 800*I3Units.m, dataclasses.I3Position(0, 0, 0))
    surface_det = MuonGun.ExtrudedPolygon.from_file(gcdFile)
    
    #setup reconstruction parameters
    pulses    = 'InIcePulses'
    HLCpulses = 'HLCInIcePulses'

    #setup I3Tray
    tray = I3Tray()
    tray.context['I3RandomService'] =  phys_services.I3GSLRandomService(seed = args.nseed)
    #tray.context['I3RandomService'] =  phys_services.I3SPRNGRandomService(seed = args.nseed)
    infiles = args.inp

    tray.Add('I3Reader', Filenamelist=[args.gcd]+infiles)
     
    #weighting
    tray.Add("I3NullSplitter")
    generator = weighting.from_simprod(args.did)
    generator *= 50
    flux = GaisserH4a()
    tray.Add(getweights, generator=generator, flux=flux)

    #find intersecting muons 
    tray.Add(muonstodet, surface=surface_det)
    tray.Add(primarytodet, surface=surface_det)

    #do an N > 0 cut
    def ncut(frame):
        if frame.Has('EnteringMuon_0'):
            return True
        else:
            return False
    tray.AddModule(ncut, 'ncut')
    
    #write everything to file
    tray.AddModule('I3Writer', 'writer',
                   Streams = [icetray.I3Frame.Physics, 
                              #icetray.I3Frame.DAQ
                          ],
                   filename=args.out+"_"+str(args.nseed)+".i3.gz")
    tray.Execute()
    tray.Finish()
    end_time = time.asctime()
    print 'Ended:', end_time

if __name__ == '__main__':
    main()
