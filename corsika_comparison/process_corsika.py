#!/usr/bin/env python
import os, sys, numpy, argparse, time

#import icecube stuff
from I3Tray import I3Tray
from icecube import icetray, dataio, dataclasses, phys_services, simclasses
from icecube.simprod import segments
from icecube import VHESelfVeto, DomTools, trigger_splitter, MuonGun
from icecube.icetray import I3Units
from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw
from icecube.MuonGun.segments import GenerateBundles
from icecube.simprod import segments

start_time = time.asctime()
print 'Started:', start_time

def todet(frame, surface):
    detmu = MuonGun.muons_at_surface(frame, surface)
    #print "multi: " + str(len(detmu))
    for i in range(len(detmu)):
        frame['EnteringMuon_'+str(i)] = detmu[i]
    
def header(frame):
    frame['I3EventHeader'] = dataclasses.I3EventHeader()

@icetray.traysegment
def DetectorSim(tray, name,
                RandomService = None,
                RunID = None,
                KeepMCHits = False,
                KeepPropagatedMCTree = False,
                KeepMCPulses = False,
                SkipNoiseGenerator = False,
                LowMem = False,
                InputPESeriesMapName = "I3MCPESeriesMap",
                BeaconLaunches = True):
    from I3Tray import I3Units

    from icecube import DOMLauncher
    from icecube import topsimulator

    if RunID is None:
        icetray.logging.log_fatal("You *must* set a RunID in production.")

    if not RandomService:
        icetray.logging.log_fatal("You *must* set a RandomService name.")

    MCPESeriesMapNames = [
        InputPESeriesMapName,
        "BackgroundI3MCPESeriesMap",
        "SignalI3MCPEs"
        ]
    MCPulseSeriesMapNames = [
        "I3MCPulseSeriesMap",
        "I3MCPulseSeriesMapParticleIDMap"
        ]
    MCTreeNames = [
        "I3MCTree",
        "BackgroundI3MCTree"
        ]
    MCPMTResponseMapNames = []

    if not SkipNoiseGenerator:
        InputPESeriesMapName_withoutNoise = InputPESeriesMapName + "WithoutNoise"
        tray.Add("Rename", "RenamePESeriesMap",
                 Keys=[InputPESeriesMapName, InputPESeriesMapName_withoutNoise])
        MCPESeriesMapNames.append(InputPESeriesMapName_withoutNoise)

        from icecube import vuvuzela

        tray.AddSegment(vuvuzela.AddNoise, name+"_vuvuzela",
           OutputName = InputPESeriesMapName,
            InputName = InputPESeriesMapName_withoutNoise,
            StartTime = -10.*I3Units.microsecond,
            EndTime = 10.*I3Units.microsecond,
            RandomServiceName = RandomService,
            )

    tray.AddSegment(DOMLauncher.DetectorResponse, "DetectorResponse",
        pmt_config = {'Input':InputPESeriesMapName,
                      'Output':"I3MCPulseSeriesMap",
                      'MergeHits':True,
                      'LowMem':LowMem,
                      'RandomServiceName' : RandomService},
        dom_config = {'Input':'I3MCPulseSeriesMap',
                      'Output':"I3DOMLaunchSeriesMap",
                      'UseTabulatedPT':True,
                      'RandomServiceName' : RandomService,
                      'BeaconLaunches':BeaconLaunches})
    tray.Add(header, streams=[icetray.I3Frame.DAQ])

    from icecube import WaveCalibrator
    tray.AddModule("I3WaveCalibrator", "calibrate",
                   Waveforms='CalibratedWaveforms',
                   CorrectDroop=True,
                   )
    
    tray.AddModule("I3PMTSaturationFlagger", "find_saturation",
                   Waveforms='CalibratedWaveforms',
                   Output="PMTSaturation",
                   )
    
    icetray.load('wavedeform', False)
    tray.AddModule("I3Wavedeform", "deform",
                   Waveforms='CalibratedWaveforms',
                   Output='InIcePulses'
                   )
    
    tray.AddModule("Delete", name+"_cleanup",
                   Keys = ["MCTimeIncEventID",
                           "MCPMTResponseMap",
                           ])
    
    if not KeepMCPulses:
        tray.AddModule("Delete", name+"_cleanup_2",
                       Keys = MCPulseSeriesMapNames + MCPMTResponseMapNames)
        
    if not KeepMCHits:
        tray.AddModule("Delete", name+"_cleanup_I3MCHits_2",
                       Keys = MCPESeriesMapNames)
        
    if not KeepPropagatedMCTree: # Always keep original tree
        tray.AddModule("Delete", name+"_cleanup_I3MCTree_3",
                       Keys = MCTreeNames)
        

def main():
    
    #arguement parser
    parser = argparse.ArgumentParser(description='Generates muons in IceCube with varying multiplicity')
    parser.add_argument('inp', nargs='+')
    parser.add_argument('--nseed', default=1, type=int,
                        help='seed for randomization')
    parser.add_argument('--out', default='corsika', help='Output file')
    parser.add_argument('--runnum', default=1, type=int,
                        help='Run number for this sim')
    #detector args
    parser.add_argument('--gcd', default=os.path.expandvars('$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'),
                        type=str,
                        help='gcd file')
    args = parser.parse_args()

    #geometry parameters
    gcdFile = args.gcd
    surface = Cylinder(1600*I3Units.m, 800*I3Units.m, dataclasses.I3Position(0, 0, 0))
    #surface_det = MuonGun.ExtrudedPolygon.from_file(gcdFile)
    
    #setup reconstruction parameters
    pulses    = 'InIcePulses'
    HLCpulses = 'HLCInIcePulses'

    #setup I3Tray
    tray = I3Tray()
    tray.context['I3RandomService'] =  phys_services.I3GSLRandomService(seed = args.nseed)
    #tray.context['I3RandomService'] =  phys_services.I3SPRNGRandomService(seed = args.nseed)
    infiles = args.inp

    tray.Add('I3Reader', Filenamelist=[args.gcd]+infiles)
     
    #find intersecting muons
    tray.Add("I3NullSplitter")
    tray.Add(todet, surface=surface)

    #write everything to file
    tray.AddModule('I3Writer', 'writer',
                   Streams = [icetray.I3Frame.Physics, 
                              icetray.I3Frame.DAQ],
                   filename=args.out+"_"+str(args.nseed)+".i3.gz")
    tray.Execute()
    tray.Finish()
    end_time = time.asctime()
    print 'Ended:', end_time

if __name__ == '__main__':
    main()
