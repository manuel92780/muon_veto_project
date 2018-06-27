#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description = "Plot veto efficiency")
parser.add_argument('-d','--debug', default=False, dest = 'BUG')
parser.add_argument('-k','--kinematic', default='energy', dest = 'KIN')
parser.add_argument('-n','--numbermuons', default='1', type=int, dest = 'NUM')
args = parser.parse_args()

from icecube import icetray, dataclasses, dataio
from icecube import phys_services, simclasses, MuonGun
from I3Tray import I3Tray
from icecube.icetray import I3Units
from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw, BundleConfiguration, BundleEntry

from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def take_ratios(numerator, denominator):
    ratio = list()
    for i in range(len(numerator)):
        if (denominator[i] == 0): n_ratio = 0
        else: n_ratio = numerator[i]/ denominator[i]
        ratio.append(n_ratio)
    return ratio

#load a long list of files
file_name = 'file_list.txt'
file_list = open(file_name).read().splitlines()
infile    = dataio.I3FrameSequence(file_list)
normalization  = len(file_list)
print "loaded your MC in " + str(normalization) +" files!"

#initialize values
muon_1 = [];
muon_1.append([]);#[0] energy
muon_1.append([]);#[1] zenith
muon_1.append([]);#[2] all weights
muon_1.append([]);#[3] passes PE >3 veto
muon_1.append([]);#[4] passes PE >5 veto
muon_1.append([]);#[5] passes PE >10 veto 

#initialize the loop 
event_count = 1
for frame in infile:
    if(event_count%50000 == 0): 
        print "Event: "+ str(event_count); 
        if(args.BUG): break;
    if ("I3MCTree" in frame) and ("EntryMuon" in frame): #count all generated events
        event_count += 1; 
        number_muons = 0;
        energy_1 = 0; zenith_1 = 0; weight_1 = 0;
        veto3_1  = 0; veto5_1  = 0; veto10_1 = 0;
        for particle in frame["I3MCTree"]:
            if particle.type is not dataclasses.I3Particle.ParticleType.MuMinus: continue
            number_muons += 1  #count the muons
            energy = particle.energy  #sort by energy
            if(energy > energy_1):
                energy_1 = particle.energy / I3Units.GeV  #take normalized energy 
                zenith_1 = np.cos(particle.dir.zenith)              #take cosine of zenith
                weight_1 = frame["MuonWeight"].value/normalization #normalized by total files                
                if("VHESelfVeto_3" in frame):  veto3_1  = (int(frame["VHESelfVeto_3"].value) ) * weight_1;
                if("VHESelfVeto_5" in frame):  veto5_1  = (int(frame["VHESelfVeto_5"].value) ) * weight_1;
                if("VHESelfVeto_10" in frame): veto10_1 = (int(frame["VHESelfVeto_10"].value)) * weight_1;
        if(args.NUM == 1 and number_muons == 1):
            muon_1[0].append(energy_1); muon_1[1].append(zenith_1); muon_1[2].append(weight_1);
            muon_1[3].append(veto3_1); muon_1[4].append(veto5_1); muon_1[5].append(veto10_1);

print "Processed " +str(event_count)+ " events generated by you!"

xbins= 50
xmin = 0; xmax = 0;
x_values = 0
if(args.KIN == 'energy'): 
    x_values = muon_1[0]
    xmin = 1e1; xmax = 1e5; xbins = np.logspace(np.log10(xmin),np.log10(xmax), xbins)
if(args.KIN == 'zenith'): 
    x_values = muon_1[1]
    xmin = 0; xmax = 1;
plt.xscale('log')
weight_vals  = plt.hist(x_values, weights=muon_1[2], bins=xbins, range=(xmin,xmax), log="True", histtype='step', label='all events')
veto3_vals   = plt.hist(x_values, weights=muon_1[3], bins=xbins, range=(xmin,xmax), log="True", histtype='step', label='PE > 3')
veto5_vals   = plt.hist(x_values, weights=muon_1[4], bins=xbins, range=(xmin,xmax), log="True", histtype='step', label='PE > 5')
veto10_vals  = plt.hist(x_values, weights=muon_1[5], bins=xbins, range=(xmin,xmax), log="True", histtype='step', label='PE > 10')
acceptance3  = take_ratios(veto3_vals[0],  weight_vals[0])
acceptance5  = take_ratios(veto5_vals[0],  weight_vals[0])
acceptance10 = take_ratios(veto10_vals[0], weight_vals[0])
xvals   = weight_vals[1][:-1]
plt.clf()
#plt.show(); exit(0)
if(args.KIN == 'energy'):
    plt.semilogx(xvals, acceptance3,  label='PE > 3')
    plt.semilogx(xvals, acceptance5,  label='PE > 5')
    #plt.semilogx(xvals, acceptance10, label='PE > 10')
    plt.xlabel("E$_{\mu}$ [GeV]")
    plt.legend(loc='upper left')
if(args.KIN == 'zenith'):
    plt.plot(xvals, acceptance3,  label='PE > 3')
    plt.plot(xvals, acceptance5,  label='PE > 5')
    plt.plot(xvals, acceptance10, label='PE > 10')
    plt.xlabel("cos (theta)")
    plt.legend(loc='lower right')
plt.ylabel('P$_{light}$')
#plt.show()
plt.savefig(args.KIN+"_dist_bundle_"+str(args.NUM)+".pdf")
print "Made file: " + args.KIN+"_dist_bundle_"+str(args.NUM)+".pdf"
