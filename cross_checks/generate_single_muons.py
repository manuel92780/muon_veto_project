#!/usr/bin/env python

#import sys defs
import sys
import argparse
from os.path import expandvars

#import icecube defs
from icecube import icetray, dataclasses, dataio, phys_services
from I3Tray import I3Tray
from icecube.icetray import I3Units
from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw
from icecube.MuonGun.segments import GenerateBundles

#propagation 
def make_propagators():
        from icecube.sim_services import I3ParticleTypePropagatorServiceMap
        from icecube.PROPOSAL import I3PropagatorServicePROPOSAL
        from icecube.cmc import I3CascadeMCService
        propagators = I3ParticleTypePropagatorServiceMap()
        muprop = I3PropagatorServicePROPOSAL(type=dataclasses.I3Particle.MuMinus, cylinderHeight=1200, cylinderRadius=700)
        cprop = I3CascadeMCService(phys_services.I3GSLRandomService(1)) # dummy RNG                                              
        for pt in 'MuMinus', 'MuPlus':
                propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = muprop
        for pt in 'DeltaE', 'Brems', 'PairProd', 'NuclInt', 'Hadrons', 'EMinus', 'EPlus':
                propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = cprop
        return propagators
#start main script
def main():
	parser = argparse.ArgumentParser(description='Generates muons in IceCube with varying multiplicity')
	parser.add_argument('--nseed', default=1, type=int,
			    help='seed for randomization')
        parser.add_argument('--model', default='Hoerandel5_atmod12_SIBYLL', type=str)
	parser.add_argument('--gcd', default='auto', type=str)
	parser.add_argument('--ice', default='SpiceLea',
			    help='Either Spice1, SpiceMie or SpiceLea')
	parser.add_argument('--include-gcd', default=False, action='store_true',
			    help='Include GCD in output')
	parser.add_argument('--no-hybrid',  action="store_false", default=False,
			    dest='hybrid', help='do not perform a hybrid simulation (i.e. use clsim only)')
	parser.add_argument('--emin', default=1e3, type=float,
			    help='Muon min energy (GeV)')
	parser.add_argument('--emax', default=1e6, type=float,
			    help='Muon max energy (GeV)')
	parser.add_argument('--nevents', default=1e3, type=int,
			    help='Number of events')
	parser.add_argument('--out', default='muongun.i3.gz', help='Output file')
	parser.add_argument('--runnum', default=1, type=int,
			    help='Run number for this sim')
	parser.add_argument('--use-gpu',  action='store_true', default=False,
			    help='simulate using GPUs instead of CPU cores')
	args = parser.parse_args()
    
	
	if args.gcd == 'auto':
		gcdFilesDefault = {'IC86': os.path.expandvars('$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'),
				   'IC79': os.path.expandvars('$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC79.55380_L2a.i3.gz'), }
		gcdFile = gcdFilesDefault[args.detector]
	else:
		gcdFile = args.gcd

	tray = I3Tray()
	tray.context['I3RandomService'] =  phys_services.I3GSLRandomService(
		seed = args.nseed)
	gcd = expandvars('$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz')
	
	model = load_model(args.model)
	surface = Cylinder(1600*I3Units.m, 800*I3Units.m, dataclasses.I3Position(0, 0, 0))
	spectrum = OffsetPowerLaw(2, 0*I3Units.TeV, args.emin*I3Units.GeV, args.emax*I3Units.GeV)
	generator = StaticSurfaceInjector(surface, model.flux, spectrum, model.radius)

	tray.AddSegment(GenerateBundles, 'MuonGenerator', Generator=generator, NEvents=args.nevents, GCDFile=gcd)
	tray.AddModule('I3MuonGun::WeightCalculatorModule', 'MuonWeight', Model=model,
		       Generator=generator)

	icetray.logging.set_level_for_unit("I3PropagatorService", "TRACE")
	icetray.logging.set_level_for_unit("I3PropagatorMMC", "TRACE")
	icetray.logging.set_level_for_unit("I3PropagatorServiceMMC", "TRACE")

	tray.Add(segments.PropagateMuons, 'PropagateMuons', RandomService='I3RandomService')
	tray.Add(segments.PropagatePhotons, "PropagatePhotons",
		 InputMCTree='I3MCTree',
		 MaxParallelEvents=1,
		 KeepIndividualMaps=False,
		 IceModel=args.ice,
		 UnshadowedFraction=0.9,
		 HybridMode=args.hybrid,
		 IgnoreMuons=False,
		 UseAllCPUCores = True,
		 UseGPUs=args.use_gpu)
	tray.Add(segments.DetectorSim, "DetectorSim",
		 RandomService='I3RandomService',
		 RunID=args.runnum,
		 GCDFile=gcdFile,
		 InputPESeriesMapName="I3MCPESeriesMap",
		 KeepMCHits=True,
		 KeepPropagatedMCTree=True,
		 SkipNoiseGenerator=False)
        #tray.AddModule('I3PropagatorModule', 'propagator', PropagatorServices=make_propagators(),
	#	       RandomService='I3RandomService', RNGStateName="RNGState")
	tray.AddModule('I3Writer', 'writer',
		       Streams=list(map(icetray.I3Frame.Stream, "SQP")),
		       filename=args.model+"_"+args.out)

	
	tray.Execute()

if __name__ == '__main__':
    main()
