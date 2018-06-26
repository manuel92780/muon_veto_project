#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description = "Plot muons.")
parser.add_argument('-o','--outfile', default='tester', dest = 'OUTFILE')
parser.add_argument('-k','--kinematic', default='radial', dest = 'KIN')
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

def dist_cartesian(mx, my, mz, sx, sy, sz):
    length2 = ((mx-sx)*(mx-sx)) + ((my-sy)*(my-sy)) + ((mz-sz)*(mz-sz))
    return np.sqrt(length2)

#import files
hoe_data_name = 'Hoerandel5_atmod12_SIBYLL_muongun.i3.gz'
hoe_infile    = dataio.I3File(hoe_data_name)
gai_data_name = 'GaisserH4a_atmod12_SIBYLL_muongun.i3.gz'
gai_infile    = dataio.I3File(gai_data_name)
print "loaded your MC"

#run the hoe
h_muon_e_1  = list(); h_muon_ze_1  = list(); h_muon_r_1  = list(); h_muon_a_1  = list(); h_muon_w_1  = list()
h_muon_e_2  = list(); h_muon_ze_2  = list(); h_muon_r_2  = list(); h_muon_a_2  = list(); h_muon_w_2  = list()
h_muon_e_3  = list(); h_muon_ze_3  = list(); h_muon_r_3  = list(); h_muon_a_3  = list(); h_muon_w_3  = list()
event_count   = 1
for frame in hoe_infile:
    if(event_count%10000 == 0): print "Event: "+ str(event_count)
    if "I3MCTree" in frame: #new frame, reset everything
        event_count += 1
        number_muons = 0
        muon_e_1 = 0.0; muon_e_2 = 0.0; muon_e_3 = 0.0; #energy
        muon_ze_1 = 0.0; muon_ze_2 = 0.0; muon_ze_3 = 0.0; #zenith
        muon_r_1 = 0.0; muon_r_2 = 0.0; muon_r_3 = 0.0; #radial distance
        muon_a_1 = 0.0; muon_a_2 = 0.0; muon_a_3 = 0.0; #azimuth
        muon_w_1 = 0.0; muon_w_2 = 0.0; muon_w_3 = 0.0; #weight
        #coordinates to find R
        shower_x = 0.0; shower_y = 0.0; shower_z = 0.0;
        muon_x_1 = 0.0; muon_x_2 = 0.0; muon_x_3 = 0.0;
        muon_y_1 = 0.0; muon_y_2 = 0.0; muon_y_3 = 0.0;
        muon_z_1 = 0.0; muon_z_2 = 0.0; muon_z_3 = 0.0;
        for particle in frame["I3MCTree"]: 
            #first find the shower axis
            if particle.type == dataclasses.I3Particle.ParticleType.unknown: 
                shower_x = particle.pos.x; shower_y = particle.pos.y; shower_z = particle.pos.z;
            #then find all muons in bundle
            if particle.type == dataclasses.I3Particle.ParticleType.MuMinus:
                number_muons += 1 #count the muons
                energy = particle.energy / I3Units.GeV #sort by energy
                if  (energy > muon_e_1):
                    muon_e_3 = muon_e_2; muon_e_2 = muon_e_1; muon_e_1 = energy;
                    muon_ze_3 = muon_ze_2; muon_ze_2 = muon_ze_1; muon_ze_1 = particle.dir.zenith;
                    muon_a_3 = muon_a_2; muon_a_2 = muon_a_1; muon_a_1 = particle.dir.azimuth;
                    muon_w_3 = muon_w_2; muon_w_2 = muon_w_1; muon_w_1 = frame["MuonWeight"].value;
                    muon_x_3 = muon_x_2; muon_x_2 = muon_x_1; muon_x_1 = particle.pos.x;
                    muon_y_3 = muon_y_2; muon_y_2 = muon_y_1; muon_y_1 = particle.pos.y;
                    muon_z_3 = muon_z_2; muon_z_2 = muon_z_1; muon_z_1 = particle.pos.z;
                elif(energy > muon_e_2):
                    muon_e_3 = muon_e_2; muon_e_2 = energy;
                    muon_ze_3 = muon_ze_2; muon_ze_2 = particle.dir.zenith;
                    muon_a_3 = muon_a_2; muon_a_2 = particle.dir.azimuth;
                    muon_w_3 = muon_w_2; muon_w_2 = frame["MuonWeight"].value;
                    muon_x_3 = muon_x_2; muon_x_2 = particle.pos.x;
                    muon_y_3 = muon_y_2; muon_y_2 = particle.pos.y;
                    muon_z_3 = muon_z_2; muon_z_2 = particle.pos.z;
                elif(energy > muon_e_3):
                    muon_e_3 = energy;
                    muon_ze_3 = particle.dir.zenith;
                    muon_a_3 = particle.dir.azimuth;
                    muon_w_3 = frame["MuonWeight"].value;
                    muon_x_3 = particle.pos.x;
                    muon_y_3 = particle.pos.y;
                    muon_z_3 = particle.pos.z;
        muon_r_1 = np.sin(muon_ze_1) * dist_cartesian(muon_x_1, muon_y_1, muon_z_1, shower_x, shower_y, shower_z)
        muon_r_2 = np.sin(muon_ze_2) * dist_cartesian(muon_x_2, muon_y_2, muon_z_2, shower_x, shower_y, shower_z)
        muon_r_3 = np.sin(muon_ze_3) * dist_cartesian(muon_x_3, muon_y_3, muon_z_3, shower_x, shower_y, shower_z)
        if  (number_muons > 2):
            h_muon_e_1.append(muon_e_1); h_muon_e_2.append(muon_e_2); h_muon_e_3.append(muon_e_3)
            h_muon_ze_1.append(muon_ze_1); h_muon_ze_2.append(muon_ze_2); h_muon_ze_3.append(muon_ze_3)
            h_muon_r_1.append(muon_r_1); h_muon_r_2.append(muon_r_2); h_muon_r_3.append(muon_r_3)
            h_muon_a_1.append(muon_a_1); h_muon_a_2.append(muon_a_2); h_muon_a_3.append(muon_a_3)
            h_muon_w_1.append(muon_w_1); h_muon_w_2.append(muon_w_2); h_muon_w_3.append(muon_w_3)
        elif(number_muons > 1):
            h_muon_e_1.append(muon_e_1); h_muon_e_2.append(muon_e_2); 
            h_muon_ze_1.append(muon_ze_1); h_muon_ze_2.append(muon_ze_2);
            h_muon_r_1.append(muon_r_1); h_muon_r_2.append(muon_r_2);
            h_muon_a_1.append(muon_a_1); h_muon_a_2.append(muon_a_2);
            h_muon_w_1.append(muon_w_1); h_muon_w_2.append(muon_w_2); 
        elif(number_muons > 0):
            h_muon_e_1.append(muon_e_1); 
            h_muon_ze_1.append(muon_ze_1);
            h_muon_r_1.append(muon_r_1);
            h_muon_a_1.append(muon_a_1);
            h_muon_w_1.append(muon_w_1); 
print "Processed " +str(event_count)+ " Hoerandel5 events generated by you!"

#run gai
g_muon_e_1  = list(); g_muon_ze_1  = list(); g_muon_r_1  = list(); g_muon_a_1  = list(); g_muon_w_1  = list()
g_muon_e_2  = list(); g_muon_ze_2  = list(); g_muon_r_2  = list(); g_muon_a_2  = list(); g_muon_w_2  = list()
g_muon_e_3  = list(); g_muon_ze_3  = list(); g_muon_r_3  = list(); g_muon_a_3  = list(); g_muon_w_3  = list()
event_count   = 1
for frame in gai_infile:
    if(event_count%10000 == 0): print "Event: "+ str(event_count)
    if "I3MCTree" in frame: #new frame, reset everything
        event_count += 1
        number_muons = 0
        muon_e_1 = 0.0; muon_e_2 = 0.0; muon_e_3 = 0.0; #energy
        muon_ze_1 = 0.0; muon_ze_2 = 0.0; muon_ze_3 = 0.0; #zenith
        muon_r_1 = 0.0; muon_r_2 = 0.0; muon_r_3 = 0.0; #radial distance
        muon_a_1 = 0.0; muon_a_2 = 0.0; muon_a_3 = 0.0; #azimuth
        muon_w_1 = 0.0; muon_w_2 = 0.0; muon_w_3 = 0.0; #weight
        #coordinates to find R
        shower_x = 0.0; shower_y = 0.0; shower_z = 0.0;
        muon_x_1 = 0.0; muon_x_2 = 0.0; muon_x_3 = 0.0;
        muon_y_1 = 0.0; muon_y_2 = 0.0; muon_y_3 = 0.0;
        muon_z_1 = 0.0; muon_z_2 = 0.0; muon_z_3 = 0.0;
        for particle in frame["I3MCTree"]:
            #first find the shower axis
            if particle.type == dataclasses.I3Particle.ParticleType.unknown:
                shower_x = particle.pos.x; shower_y = particle.pos.y; shower_z = particle.pos.z;
            #then find all muons in bundle
            if particle.type == dataclasses.I3Particle.ParticleType.MuMinus:
                number_muons += 1 #count the muons
                energy = particle.energy / I3Units.GeV
                if  (energy > muon_e_1):
                    muon_e_3 = muon_e_2; muon_e_2 = muon_e_1; muon_e_1 = energy;
                    muon_ze_3 = muon_ze_2; muon_ze_2 = muon_ze_1; muon_ze_1 = particle.dir.zenith;
                    muon_a_3 = muon_a_2; muon_a_2 = muon_a_1; muon_a_1 = particle.dir.azimuth;
                    muon_w_3 = muon_w_2; muon_w_2 = muon_w_1; muon_w_1 = frame["MuonWeight"].value;
                    muon_x_3 = muon_x_2; muon_x_2 = muon_x_1; muon_x_1 = particle.pos.x;
                    muon_y_3 = muon_y_2; muon_y_2 = muon_y_1; muon_y_1 = particle.pos.y;
                    muon_z_3 = muon_z_2; muon_z_2 = muon_z_1; muon_z_1 = particle.pos.z;
                elif(energy > muon_e_2):
                    muon_e_3 = muon_e_2; muon_e_2 = energy;
                    muon_ze_3 = muon_ze_2; muon_ze_2 = particle.dir.zenith;
                    muon_a_3 = muon_a_2; muon_a_2 = particle.dir.azimuth;
                    muon_w_3 = muon_w_2; muon_w_2 = frame["MuonWeight"].value;
                    muon_x_3 = muon_x_2; muon_x_2 = particle.pos.x;
                    muon_y_3 = muon_y_2; muon_y_2 = particle.pos.y;
                    muon_z_3 = muon_z_2; muon_z_2 = particle.pos.z;
                elif(energy > muon_e_3):
                    muon_e_3 = energy;
                    muon_ze_3 = particle.dir.zenith;
                    muon_a_3 = particle.dir.azimuth;
                    muon_w_3 = frame["MuonWeight"].value;
                    muon_x_3 = particle.pos.x;
                    muon_y_3 = particle.pos.y;
                    muon_z_3 = particle.pos.z;
        #print '--------------'
        #print shower_x; print shower_y; print shower_z;
        #print muon_x_1; print muon_y_1; print muon_z_1;
        muon_r_1 = np.sin(muon_ze_1) * dist_cartesian(muon_x_1, muon_y_1, muon_z_1, shower_x, shower_y, shower_z);
        muon_r_2 = np.sin(muon_ze_2) * dist_cartesian(muon_x_2, muon_y_2, muon_z_2, shower_x, shower_y, shower_z)
        muon_r_3 = np.sin(muon_ze_3) * dist_cartesian(muon_x_3, muon_y_3, muon_z_3, shower_x, shower_y, shower_z)
        if  (number_muons > 2):
            g_muon_e_1.append(muon_e_1); g_muon_e_2.append(muon_e_2); g_muon_e_3.append(muon_e_3)
            g_muon_ze_1.append(muon_ze_1); g_muon_ze_2.append(muon_ze_2); g_muon_ze_3.append(muon_ze_3)
            g_muon_r_1.append(muon_r_1); g_muon_r_2.append(muon_r_2); g_muon_r_3.append(muon_r_3)
            g_muon_a_1.append(muon_a_1); g_muon_a_2.append(muon_a_2); g_muon_a_3.append(muon_a_3)
            g_muon_w_1.append(muon_w_1); g_muon_w_2.append(muon_w_2); g_muon_w_3.append(muon_w_3)
        elif(number_muons > 1):
            g_muon_e_1.append(muon_e_1); g_muon_e_2.append(muon_e_2);
            g_muon_ze_1.append(muon_ze_1); g_muon_ze_2.append(muon_ze_2);
            g_muon_r_1.append(muon_r_1); g_muon_r_2.append(muon_r_2);
            g_muon_a_1.append(muon_a_1); g_muon_a_2.append(muon_a_2);
            g_muon_w_1.append(muon_w_1); g_muon_w_2.append(muon_w_2);
        elif(number_muons > 0):
            g_muon_e_1.append(muon_e_1);
            g_muon_ze_1.append(muon_ze_1);
            g_muon_r_1.append(muon_r_1);
            g_muon_a_1.append(muon_a_1);
            g_muon_w_1.append(muon_w_1);
print "Processed " +str(event_count)+ " GaisserH4a events generated by you!"

h_muon_e_1 = np.log10(h_muon_e_1); h_muon_e_2 = np.log10(h_muon_e_2); h_muon_e_3 = np.log10(h_muon_e_3)
g_muon_e_1 = np.log10(g_muon_e_1); g_muon_e_2 = np.log10(g_muon_e_2); g_muon_e_3 = np.log10(g_muon_e_3)


#plot top histo
xbins= 50
h_1_vals = 0; h_2_vals = 0; h_3_vals = 0
g_1_vals = 0; g_2_vals = 0; g_3_vals = 0
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[3, 1]})
if(args.KIN == "energy"):
    h_1_vals = ax1.hist(h_muon_e_1, weights=h_muon_w_1, bins=xbins, log=True, histtype='step', 
                        linestyle=('solid'), color=('blue'), label='Hoerandel5')
    g_1_vals = ax1.hist(g_muon_e_1, weights=g_muon_w_1, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('blue'), label='GaisserH4a')
    h_2_vals = ax1.hist(h_muon_e_2, weights=h_muon_w_2, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('red'))
    g_2_vals = ax1.hist(g_muon_e_2, weights=g_muon_w_2, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('red'))
    h_3_vals = ax1.hist(h_muon_e_3, weights=h_muon_w_3, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('green'))
    g_3_vals = ax1.hist(g_muon_e_3, weights=g_muon_w_3, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('green'))
if(args.KIN == "zenith"):
    h_1_vals = ax1.hist(np.cos(h_muon_ze_1), weights=h_muon_w_1, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('blue'), label='Hoerandel5')
    g_1_vals = ax1.hist(np.cos(g_muon_ze_1), weights=g_muon_w_1, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('blue'), label='GaisserH4a')
    h_2_vals = ax1.hist(np.cos(h_muon_ze_2), weights=h_muon_w_2, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('red'))
    g_2_vals = ax1.hist(np.cos(g_muon_ze_2), weights=g_muon_w_2, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('red'))
    h_3_vals = ax1.hist(np.cos(h_muon_ze_3), weights=h_muon_w_3, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('green'))
    g_3_vals = ax1.hist(np.cos(g_muon_ze_3), weights=g_muon_w_3, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('green'))
if(args.KIN == "azimuthal"):
    h_1_vals = ax1.hist(h_muon_a_1, weights=h_muon_w_1, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('blue'), label='Hoerandel5')
    g_1_vals = ax1.hist(g_muon_a_1, weights=g_muon_w_1, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('blue'), label='GaisserH4a')
    h_2_vals = ax1.hist(h_muon_a_2, weights=h_muon_w_2, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('red'))
    g_2_vals = ax1.hist(g_muon_a_2, weights=g_muon_w_2, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('red'))
    h_3_vals = ax1.hist(h_muon_a_3, weights=h_muon_w_3, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('green'))
    g_3_vals = ax1.hist(g_muon_a_3, weights=g_muon_w_3, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('green'))
if(args.KIN == "radial"):
    h_1_vals = ax1.hist(h_muon_r_1, weights=h_muon_w_1, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('blue'), label='Hoerandel5')
    g_1_vals = ax1.hist(g_muon_r_1, weights=g_muon_w_1, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('blue'), label='GaisserH4a')
    h_2_vals = ax1.hist(h_muon_r_2, weights=h_muon_w_2, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('red'))
    g_2_vals = ax1.hist(g_muon_r_2, weights=g_muon_w_2, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('red'))
    h_3_vals = ax1.hist(h_muon_r_3, weights=h_muon_w_3, bins=xbins, log=True, histtype='step',
                        linestyle=('solid'), color=('green'))
    g_3_vals = ax1.hist(g_muon_r_3, weights=g_muon_w_3, bins=xbins, log=True, histtype='step',
                        linestyle=('dashed'), color=('green'))
ax1.set_ylabel('Rate per Bin [Hz]')
ax1.grid(linestyle=':', linewidth=1)
if(args.KIN == "energy"):    ax1.legend(loc='upper right')
if(args.KIN == "zenith"):    ax1.legend(loc='lower right')
if(args.KIN == "azimuthal"): ax1.legend(loc='lower left')
if(args.KIN == "radial"):    ax1.legend(loc='upper right')

#plot ratios
xvals = h_1_vals[1][:-1]
ratio_1 = take_ratios(h_1_vals[0],g_1_vals[0])
ratio_2 = take_ratios(h_2_vals[0],g_2_vals[0])
ratio_3 = take_ratios(h_3_vals[0],g_3_vals[0])
ax2.hist(xvals, bins=xbins, weights=ratio_1, log=False, histtype='step',
         linestyle=('solid'), color=('blue'))
ax2.hist(xvals, bins=xbins, weights=ratio_2, log=False, histtype='step',
         linestyle=('solid'), color=('red'))
ax2.hist(xvals, bins=xbins, weights=ratio_3, log=False, histtype='step',
         linestyle=('solid'), color=('green'))
if(args.KIN == "energy"):    ax2.set_xlabel('log10 (Muon Energy [GeV])')
if(args.KIN == "zenith"):    ax2.set_xlabel('cos(Zenith)')
if(args.KIN == "azimuthal"): ax2.set_xlabel('Azimuthal')
if(args.KIN == "radial"):    ax2.set_xlabel('Radial Distance (m)')

ax2.set_ylabel('Hoe / Gai')
ax2.yaxis.set_ticks(np.arange(0, 1.8, 0.3))
ax2.set_ylim([0,1.8])
ax2.grid(linestyle=':', linewidth=1)
plt.savefig(args.KIN+"_dist_"+args.OUTFILE+".pdf")
