#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description = "Plot veto efficiency")
parser.add_argument('-d','--debug', default=False, dest = 'BUG')
parser.add_argument('-p','--PEthreshold', default="3", dest = 'PE')
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
from scipy.optimize import curve_fit

def sigmoid(x, k0, x0, c):
     y = (1-c) / (1 + np.exp(-k0*(np.log10(x)-x0)))+c
     return y

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
depths = [-550,-365,-216,-80,80,216,365,550]
for i in range(len(depths)-1):
     muon_1.append([]);#[i*3+0] energy
     muon_1.append([]);#[i*3+1] weights
     muon_1.append([]);#[i*3+2] passes veto
print muon_1

#initialize the loop 
event_count = 1
for frame in infile:
     if(event_count%50000 == 0): 
          print "Event: "+ str(event_count); 
          if(args.BUG): break;
     if "EnteringMuon_0" in frame: #count all generated events
          event_count += 1;
          energy_1 = frame["EnteringMuon_0"].energy / I3Units.GeV  #take normalized energy 
          depth_1  = frame["EnteringMuon_0"].pos.z                 #find depth
          weight_1 = frame["MuonWeight"].value/normalization #normalized by total files
          if("VHESelfVeto_"+args.PE+"Clean" in frame): 
               veto_1 = (int(frame["VHESelfVeto_"+args.PE+"Clean"].value) ) * weight_1;
          else:
               veto_1 = 0
          for i in range(len(depths)-1):
               if(depth_1 > depths[i]) and (depth_1 < depths[i+1]):
                    muon_1[i*3+0].append(energy_1) 
                    muon_1[i*3+1].append(weight_1) 
                    muon_1[i*3+2].append(veto_1)
print "Processed " +str(event_count)+ " events generated by you!"

xmin = 1e1; xmax = 1e5;
xbins= np.logspace(np.log10(xmin),np.log10(xmax), 50)
plt.xscale('log')
xvals = []; yvals = [];
for i in range(len(depths)-1):
     weight_vals = plt.hist(muon_1[i*3+0], weights=muon_1[i*3+1], bins=xbins, range=(xmin,xmax), log="True", histtype='step', label='all events')
     veto_vals   = plt.hist(muon_1[i*3+0], weights=muon_1[i*3+2], bins=xbins, range=(xmin,xmax), log="True", histtype='step', label='passes veto')
     acceptance  = take_ratios(veto_vals[0], weight_vals[0])
     xvals.append(weight_vals[1][:-1]);
     yvals.append(acceptance)

#do fit
fit_x = []; fit_y = [];
for i in range(len(depths)-1):
     popt, pcov  = curve_fit(sigmoid, xvals[i], yvals[i], maxfev=5000)
     x = np.logspace(np.log10(xmin),np.log10(xmax), len(xbins)*1)
     y  = sigmoid(x, *popt)
     fit_x.append(x); fit_y.append(y);

#do plotting
plt.clf()
#plt.show(); exit(0)
color_map = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
for i in range(len(depths)-1):
     plt.semilogx(xvals[i], yvals[i],  label='Depth: ['+str(depths[i])+","+str(depths[i+1])+ "m]",  ls='-', color=color_map[i])
     plt.semilogx(fit_x[i],fit_y[i],  ls='--', color=color_map[i])

plt.xlabel("E$_{\mu}$ [GeV]")
plt.ylabel('P$_{light}$')
plt.legend(loc='lower right')
plt.ylim(0.0, 1.1)
plt.yticks(np.arange(0, 1.1, 0.1))
plt.grid(linestyle=':', linewidth=0.5)
#plt.show()
plt.savefig("energy_vs_Plight_binned_depth_PE_cut_"+args.PE+".pdf")
print "Made file: " + "energy_vs_Plight_binned_depth_PE_cut_"+args.PE+".pdf"
