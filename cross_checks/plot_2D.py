#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description = "Plot muons.")
parser.add_argument('-d','--debug', default=False, dest = 'BUG')
parser.add_argument('-m','--model', default='Hoerandel5_atmod12_SIBYLL', dest = 'MOD')
parser.add_argument('-k','--kinematic', default='azimuth', dest = 'KIN')
parser.add_argument('-mult','--muonmultiplicity', default=2, type=int, dest = 'MULT')
parser.add_argument('-n','--muonnumber', default=2, type=int, dest = 'NUM')
args = parser.parse_args()

from icecube import icetray, dataclasses, dataio
from icecube import phys_services, simclasses, MuonGun
from I3Tray import I3Tray
from icecube.icetray import I3Units
from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw, BundleConfiguration, BundleEntry

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

#first check to make sure numbers are right
if(args.NUM > args.MULT):
    print "checking for muon #" + str(args.NUM) + " in events with exactly " + str(args.MULT) + " muons, try again"
    exit(0);

#import files
data_name = args.MOD+'_muongun.i3.gz'
infile    = dataio.I3File(data_name)
print "loaded your MC"

#initialize values
m_zen_1  = []; m_zen_1.append([]); m_zen_1.append([])  
m_zen_2  = []; m_zen_2.append([]); m_zen_2.append([])
m_zen_3  = []; m_zen_3.append([]); m_zen_3.append([])
m_a_1  = []; m_a_1.append([]); m_a_1.append([])
m_a_2  = []; m_a_2.append([]); m_a_2.append([])
m_a_3  = []; m_a_3.append([]); m_a_3.append([])
m_r_1  = []; m_r_1.append([]); m_r_1.append([])
m_r_2  = []; m_r_2.append([]); m_r_2.append([])
m_r_3  = []; m_r_3.append([]); m_r_3.append([])
m_e_1  = []; m_e_1.append([]); m_e_1.append([])
m_e_2  = []; m_e_2.append([]); m_e_2.append([])
m_e_3  = []; m_e_3.append([]); m_e_3.append([])
#initialize the loop
event_count   = 1
for frame in infile:
    if(event_count%50000 == 0): 
        print "Event: "+ str(event_count)
        if(args.BUG): break
    if "I3MCTree" in frame: #new frame, reset everything
        event_count += 1
        number_muons = 0; shower_pos = 0.0;
        muon_zen_1 = 0.0; muon_zen_2 = 0.0; muon_zen_3 = 0.0; #zenith
        muon_a_1   = 0.0; muon_a_2   = 0.0; muon_a_3   = 0.0; #azimuthal
        muon_r_1   = 0.0; muon_r_2   = 0.0; muon_r_3   = 0.0; #radius
        muon_e_1   = 0.0; muon_e_2   = 0.0; muon_e_3   = 0.0; #energy
        muon_w_1   = 0.0; muon_w_2   = 0.0; muon_w_3   = 0.0; #weight
        for particle in frame["I3MCTree"]:
            #first load the shower position, since this is recorded at t = 0
            if particle.type == dataclasses.I3Particle.ParticleType.unknown:
                shower_pos = particle.pos;
            #then find all muons in bundle
            if particle.type == dataclasses.I3Particle.ParticleType.MuMinus:
                number_muons += 1 #count the muons
                energy = np.log10(particle.energy / I3Units.GeV) #sort by log10 energy
                if  (energy > muon_e_1):
                    muon_zen_3 = muon_zen_2; muon_zen_2 = muon_zen_1; muon_zen_1 = np.cos(particle.dir.zenith);
                    muon_a_3 = muon_a_2;     muon_a_2 = muon_a_1;     muon_a_1 = particle.dir.azimuth;
                    muon_r_3 = muon_r_2;     muon_r_2 = muon_r_1;     muon_r_1 = phys_services.I3Calculator.closest_approach_distance(particle,shower_pos);
                    muon_e_3 = muon_e_2;     muon_e_2 = muon_e_1;     muon_e_1 = energy;
                    muon_w_3 = muon_w_2;     muon_w_2 = muon_w_1;     muon_w_1 = frame["MuonWeight"].value;
                elif(energy > muon_e_2):
                    muon_zen_3 = muon_zen_2; muon_zen_2 = np.cos(particle.dir.zenith);
                    muon_a_3 = muon_a_2; muon_a_2 = particle.dir.azimuth;
                    muon_r_3 = muon_r_2; muon_r_2 = phys_services.I3Calculator.closest_approach_distance(particle,shower_pos);
                    muon_e_3 = muon_e_2; muon_e_2 = energy;
                    muon_w_3 = muon_w_2; muon_w_2 = frame["MuonWeight"].value;
                elif(energy > muon_e_3):
                    muon_zen_3 = np.cos(particle.dir.zenith);
                    muon_a_3 = particle.dir.azimuth; 
                    muon_r_3 = phys_services.I3Calculator.closest_approach_distance(particle,shower_pos);
                    muon_e_3 = energy;
                    muon_w_3 = frame["MuonWeight"].value;
        if  (number_muons == 3 and args.MULT == 3):
            m_zen_1[0].append(muon_zen_1); m_zen_2[0].append(muon_zen_2); m_zen_3[0].append(muon_zen_3)
            m_a_1[0].append(muon_a_1);     m_a_2[0].append(muon_a_2);     m_a_3[0].append(muon_a_3)
            m_r_1[0].append(muon_r_1);     m_r_2[0].append(muon_r_2);     m_r_3[0].append(muon_r_3)
            m_e_1[0].append(muon_e_1);     m_e_2[0].append(muon_e_2);     m_e_3[0].append(muon_e_3)
            m_e_1[1].append(muon_w_1);     m_e_2[1].append(muon_w_2);     m_e_3[1].append(muon_w_3)
        if  (number_muons == 2 and (args.MULT == 2)):
            m_zen_1[0].append(muon_zen_1); m_zen_2[0].append(muon_zen_2);
            m_a_1[0].append(muon_a_1); m_a_2[0].append(muon_a_2);
            m_r_1[0].append(muon_r_1); m_r_2[0].append(muon_r_2);
            m_e_1[0].append(muon_e_1); m_e_2[0].append(muon_e_2);
            m_e_1[1].append(muon_w_1); m_e_2[1].append(muon_w_2);
        if  (number_muons == 1 and args.MULT == 1):
            m_zen_1[0].append(muon_zen_1);
            m_a_1[0].append(muon_a_1);
            m_r_1[0].append(muon_r_1);
            m_e_1[0].append(muon_e_1);
            m_e_1[1].append(muon_w_1);
print "Processed " +str(event_count)+ " Hoerandel5 events generated by you!"

#make 2D plot
xmin = 3.0; xmax = 6.0; xbins = 50; ymin = 0.0; ymax = 0.0; ybins = 0.0;
m1_vals = 0; m2_vals = 0; m_weights = 0;
znorm = matplotlib.colors.NoNorm()
if(args.KIN == "zenith"):
    if(args.NUM == 1): m1_vals = m_e_1[0]; m2_vals = m_zen_1[0]; m_weights=m_e_1[1]
    if(args.NUM == 2): m1_vals = m_e_2[0]; m2_vals = m_zen_2[0]; m_weights=m_e_2[1]
    if(args.NUM == 3): m1_vals = m_e_3[0]; m2_vals = m_zen_3[0]; m_weights=m_e_3[1]
    ymin = 0.0; ymax = 1.0; ybins = 50;
    znorm = matplotlib.colors.LogNorm()
    plt.ylabel('cos (Zenith)')
if(args.KIN == "azimuth"):
    if(args.NUM == 1): m1_vals = m_e_1[0]; m2_vals = m_a_1[0]; m_weights=m_e_1[1]
    if(args.NUM == 2): m1_vals = m_e_2[0]; m2_vals = m_a_2[0]; m_weights=m_e_2[1]
    if(args.NUM == 3): m1_vals = m_e_3[0]; m2_vals = m_a_3[0]; m_weights=m_e_3[1]
    ymin = 0.0; ymax = 6.5; ybins = 50;
    plt.ylabel('Azimuthal [rad]')
if(args.KIN == "radius"):
    if(args.NUM == 1): m1_vals = m_e_1[0]; m2_vals = m_r_1[0]; m_weights=m_e_1[1]
    if(args.NUM == 2): m1_vals = m_e_2[0]; m2_vals = m_r_2[0]; m_weights=m_e_2[1]
    if(args.NUM == 3): m1_vals = m_e_3[0]; m2_vals = m_r_3[0]; m_weights=m_e_3[1]
    ymin = 0.0; ymax = 100; ybins = 50;
    #znorm = matplotlib.colors.LogNorm()
    plt.ylabel('Lateral Distance [m]')

plt.hist2d(m1_vals, m2_vals, bins=[xbins,ybins], norm=znorm, weights=m_weights, range=[[xmin, xmax], [ymin,ymax]])

#make labels
cbar = plt.colorbar()
cbar.set_label('Rate per Bin [Hz]', labelpad=25, rotation=270)
plt.xlabel('log10 (Muon Energy / [GeV])')
plt.title('Muon #'+str(args.NUM)+' in Events with exactly ' +str(args.MULT)+ " Muon/s")


plt.savefig("energy_vs_"+args.KIN+"_bundle_mult_"+str(args.MULT)+"_"+str(args.NUM)+"_muon_"+args.MOD+".pdf")
print "Made: " + "energy_vs_"+args.KIN+"_bundle_mult_"+str(args.MULT)+"_"+str(args.NUM)+"_muon_"+args.MOD+".pdf" 
