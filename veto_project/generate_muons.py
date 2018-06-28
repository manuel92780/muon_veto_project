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
    if detmu:
        frame['EntryMuon'] = detmu[0]
        
def header(frame):
    frame['I3EventHeader'] = dataclasses.I3EventHeader()

def printer(frame):
    smuon = frame['EntryMuon']
    if frame.Has('VHESelfVeto_3'):
        print smuon.energy, smuon.dir.zenith, smuon.pos.z, int(frame['VHESelfVeto_3'].value)
    else:
        print smuon.energy, smuon.dir.zenith, smuon.pos.z, 0

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
    tray.Add("I3NullSplitter") #literally converts Q to P frame
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
    parser.add_argument('--nseed', default=1, type=int,
                        help='seed for randomization')
    #muongun args
    parser.add_argument('--model', default='Hoerandel5_atmod12_SIBYLL', type=str)
    parser.add_argument('--multiplicity', default=3, type=int,
                        help='Maximum muon bundle multiplcity')
    parser.add_argument('--emin', default=1e1, type=float,
                        help='Muon min energy (GeV)')
    parser.add_argument('--emax', default=1e5, type=float,
                        help='Muon max energy (GeV)')
    parser.add_argument('--nevents', default=10, type=int,
                        help='Number of events')
    parser.add_argument('--out', default='muongun.i3.gz', help='Output file')
    parser.add_argument('--runnum', default=1, type=int,
                        help='Run number for this sim')
    #detector args
    parser.add_argument('--gcd', default=os.path.expandvars('$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'),
                        type=str,
                        help='gcd file')
    parser.add_argument('--no-hybrid',  action="store_false", default=False,
                        dest='hybrid', help='do not perform a hybrid simulation (i.e. use clsim only)')
    parser.add_argument('--use-gpu',  action='store_true', default=False,
                        help='simulate using GPUs instead of CPU cores')
    args = parser.parse_args()

    #setup muongun parameters
    gcdFile = args.gcd
    model = load_model(args.model)
    model.flux.max_multiplicity = args.multiplicity
    #surface = Cylinder(1600*I3Units.m, 800*I3Units.m, dataclasses.I3Position(0, 0, 0))
    surface = MuonGun.ExtrudedPolygon.from_file(gcdFile)
    #spectrum = OffsetPowerLaw(2, 0, args.emin, args.emax)
    spectrum = OffsetPowerLaw(2, 1*I3Units.TeV, args.emin, args.emax)
    generator = StaticSurfaceInjector(surface, model.flux, spectrum, model.radius)
    
    #setup reconstruction parameters
    icetray.load('VHESelfVeto')
    pulses = 'InIcePulses'
    
    #setup I3Tray
    tray = I3Tray()
    tray.context['I3RandomService'] =  phys_services.I3GSLRandomService(seed = args.nseed)
    #generate events
    tray.AddSegment(GenerateBundles, 'MuonGenerator', Generator=generator, NEvents=args.nevents, GCDFile=gcdFile)
    tray.AddModule('I3MuonGun::WeightCalculatorModule', 'MuonWeight', Model=model,
                   Generator=generator)
    #propagate particles
    tray.Add(segments.PropagateMuons, 'PropagateMuons',
             RandomService='I3RandomService', 
             InputMCTreeName="I3MCTree", 
             OutputMCTreeName="I3MCTree")
    tray.Add(segments.PropagatePhotons, 'PropagatePhotons',
             RandomService='I3RandomService',
             HybridMode=args.hybrid, 
             MaxParallelEvents=100,
             UseAllCPUCores=True,
             UseGPUs=args.use_gpu)
    #detector stuff
    tray.Add(DetectorSim, "DetectorSim",
             RandomService='I3RandomService',                                                                              
             RunID=args.runnum,
             KeepPropagatedMCTree=True,
             KeepMCHits=True,
             KeepMCPulses=True,
             SkipNoiseGenerator=True,
             InputPESeriesMapName="I3MCPESeriesMap")

    #now do the veto
    from icecube.filterscripts import filter_globals
    icetray.load("filterscripts",False)
    icetray.load("cscd-llh",False)

    tray.Add(todet, surface=surface)
    tray.AddModule('HomogenizedQTot', 'qtot_total',
                   Pulses=pulses,
                   Output='HomogenizedQTot',
                   If= lambda frame: frame.Has('EntryMuon'))
    tray.AddModule("VHESelfVeto", 'selfveto_3',
                   VetoThreshold=3,
                   VertexThreshold=3,
                   pulses=pulses,
                   OutputBool="VHESelfVeto_3", OutputVertexPos="VHESelfVetoVertexPos_3", OutputVertexTime="VHESelfVetoVertexTime_3",
                   If = lambda frame: frame.Has('EntryMuon'))
    tray.AddModule("VHESelfVeto", 'selfveto_5',
                   VetoThreshold=5,
                   VertexThreshold=5,
                   pulses=pulses,
                   OutputBool="VHESelfVeto_5", OutputVertexPos="VHESelfVetoVertexPos_5", OutputVertexTime="VHESelfVetoVertexTime_5",
                   If = lambda frame: frame.Has('EntryMuon'))
    tray.AddModule("VHESelfVeto", 'selfveto_10',
                   VetoThreshold=10,
                   VertexThreshold=10,
                   pulses=pulses,
                   OutputBool="VHESelfVeto_10", OutputVertexPos="VHESelfVetoVertexPos_10", OutputVertexTime="VHESelfVetoVertexTime_10",
                   If = lambda frame: frame.Has('EntryMuon'))
    #tray.Add(printer, If = lambda frame:frame.Has('EntryMuon'))    

    #write everything to file
    tray.AddModule('I3Writer', 'writer',
                   Streams = [icetray.I3Frame.Physics, 
                              icetray.I3Frame.DAQ],
                   filename=args.out)
    tray.Execute()
    tray.Finish()

if __name__ == '__main__':
    main()
