#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description = "Plot muons.")
parser.add_argument('infiles', nargs='+')
parser.add_argument('-d','--debug', default=False, dest = 'BUG')
parser.add_argument('-m','--model', default='GaisserH4a_atmod12_SIBYLL', dest = 'MOD')
parser.add_argument('-k','--kinematic', default='energy', dest = 'KIN')
parser.add_argument('-n','--numbermuons', default='1', type=int, dest = 'NUM')
args = parser.parse_args()

from icecube import icetray, dataclasses, dataio
from icecube import phys_services, simclasses, MuonGun
from I3Tray import I3Tray
from icecube.icetray import I3Units
from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw, BundleConfiguration, BundleEntry

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

#import files
file_name = args.infiles
infile    = dataio.I3FrameSequence(file_name)
normalization  = len(file_name)
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
    if(event_count%10000 == 0): 
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
                    muon_w_3 = muon_w_2;     muon_w_2 = muon_w_1;     muon_w_1 = frame["Weight"].value/normalization;
                elif(energy > muon_e_2):
                    muon_zen_3 = muon_zen_2; muon_zen_2 = np.cos(particle.dir.zenith);
                    muon_a_3 = muon_a_2; muon_a_2 = particle.dir.azimuth;
                    muon_r_3 = muon_r_2; muon_r_2 = phys_services.I3Calculator.closest_approach_distance(particle,shower_pos);
                    muon_e_3 = muon_e_2; muon_e_2 = energy;
                    muon_w_3 = muon_w_2; muon_w_2 = frame["Weight"].value/normalization;
                elif(energy > muon_e_3):
                    muon_zen_3 = np.cos(particle.dir.zenith);
                    muon_a_3 = particle.dir.azimuth; 
                    muon_r_3 = phys_services.I3Calculator.closest_approach_distance(particle,shower_pos);
                    muon_e_3 = energy;
                    muon_w_3 = frame["Weight"].value/normalization;
        if  (number_muons == 3 and args.NUM == 3):
            m_zen_1[0].append(muon_zen_1); m_zen_2[0].append(muon_zen_2); m_zen_3[0].append(muon_zen_3)
            m_a_1[0].append(muon_a_1);     m_a_2[0].append(muon_a_2);     m_a_3[0].append(muon_a_3)
            m_r_1[0].append(muon_r_1);     m_r_2[0].append(muon_r_2);     m_r_3[0].append(muon_r_3)
            m_e_1[0].append(muon_e_1);     m_e_2[0].append(muon_e_2);     m_e_3[0].append(muon_e_3)
            m_e_1[1].append(muon_w_1);     m_e_2[1].append(muon_w_2);     m_e_3[1].append(muon_w_3)
        if  (number_muons == 2 and (args.NUM == 2)):
            m_zen_1[0].append(muon_zen_1); m_zen_2[0].append(muon_zen_2);
            m_a_1[0].append(muon_a_1); m_a_2[0].append(muon_a_2);
            m_r_1[0].append(muon_r_1); m_r_2[0].append(muon_r_2);
            m_e_1[0].append(muon_e_1); m_e_2[0].append(muon_e_2);
            m_e_1[1].append(muon_w_1); m_e_2[1].append(muon_w_2);
        if  (number_muons == 1 and args.NUM == 1):
            m_zen_1[0].append(muon_zen_1);
            m_a_1[0].append(muon_a_1);
            m_r_1[0].append(muon_r_1);
            m_e_1[0].append(muon_e_1);
            m_e_1[1].append(muon_w_1);
print "Processed " +str(event_count)+ " Hoerandel5 events generated by you!"

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[3, 1]})
#make top plot 
dolog=False
xbins= 50
xmin = 0.0; xmax = 0.0;
n1_vals = 0; n2_vals = 0; n3_vals = 0;
m1_vals = 0; m2_vals = 0; m3_vals = 0;
if(args.KIN == "zenith"):
    xmin = 0; xmax = 1; dolog = True;
    m1_vals = m_zen_1[0]; m2_vals = m_zen_2[0]; m3_vals = m_zen_3[0];
if(args.KIN == "azimuthal"):
    xmin = 0; xmax = 6.5; dolog = True;
    m1_vals = m_a_1[0];   m2_vals = m_a_2[0];   m3_vals = m_a_3[0];
if(args.KIN == "radius"):
    xmin = 0; xmax = 50; dolog = False;
    m1_vals = m_r_1[0];   m2_vals = m_r_2[0];   m3_vals = m_r_3[0];
if(args.KIN == "radius2"):
    xmin = 0; xmax = 500; dolog = True;
    m1_vals = np.square(m_r_1[0]);   m2_vals = np.square(m_r_2[0]);   m3_vals = np.square(m_r_3[0]);
if(args.KIN == "energy"):
    xmin = 5; xmax = 9; dolog = True;
    m1_vals = m_e_1[0];   m2_vals = m_e_2[0];   m3_vals = m_e_3[0];

if(args.NUM > 2): n3_vals = ax1.hist(m3_vals, weights=m_e_3[1], bins=xbins, range=(xmin,xmax), log=dolog, histtype='step',
                                     linestyle=('solid'), color=('green'), label='3rd Leading Muon')
if(args.NUM > 1): n2_vals = ax1.hist(m2_vals, weights=m_e_2[1], bins=xbins, range=(xmin,xmax), log=dolog, histtype='step',
                                     linestyle=('solid'), color=('red'), label='2nd Leading Muon')
if(args.NUM > 0): n1_vals = ax1.hist(m1_vals, weights=m_e_1[1], bins=xbins, range=(xmin,xmax), log=dolog, histtype='step',
                                     linestyle=('solid'), color=('blue'), label='Leading Muon')

#make ratio plot
xvals   = n1_vals[1][:-1]
if(args.NUM > 2): 
    ratio_3 = take_ratios(n3_vals[0],n1_vals[0])
    ax2.hist(xvals, bins=xbins, range=(xmin,xmax), weights=ratio_3, log=False, histtype='step',
             linestyle=('solid'), color=('green'))
if(args.NUM > 1): 
    ratio_2 = take_ratios(n2_vals[0],n1_vals[0])
    ax2.hist(xvals, bins=xbins, range=(xmin,xmax), weights=ratio_2, log=False, histtype='step',
             linestyle=('solid'), color=('red'))
if(args.NUM > 0): 
    ratio_1 = take_ratios(n1_vals[0],n1_vals[0])
    ax2.hist(xvals, bins=xbins, range=(xmin,xmax), weights=ratio_1, log=False, histtype='step',
             linestyle=('solid'), color=('blue'))

#make labels
ax1.set_ylabel('Rate per Bin [Hz]')
ax1.grid(linestyle=':', linewidth=0.5)
ax1.set_title('Events with exactly ' +str(args.NUM)+ " Muon/s")
if(args.KIN == "zenith"):
    ax1.legend(loc='lower right')
    ax2.set_xlabel('cos (Zenith)')
if(args.KIN == "azimuthal"):
    ax1.legend(loc='lower left')
    ax2.set_xlabel('Azimuthal [Radian]')
if(args.KIN == "radius"):
    ax1.legend(loc='upper right')
    ax2.set_xlabel('Lateral Distance [m]')
if(args.KIN == "radius2"):
    ax1.legend(loc='upper right')
    ax2.set_xlabel('Lateral Distance Squared [m$^{2}$]')
if(args.KIN == "energy"):
    ax1.legend(loc='upper right')
    ax2.set_xlabel('log10 (Muon Energy / [GeV])')

ax2.set_ylabel('Dist / Leading')
ax2.yaxis.set_ticks(np.arange(0, 1.8, 0.3))
ax2.set_ylim([0,1.8])
ax2.grid(linestyle=':', linewidth=1)

plt.savefig(args.KIN+"_dist_bundle_"+str(args.NUM)+"_"+args.MOD+".pdf")
