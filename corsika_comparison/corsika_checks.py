import time, os
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(description = "Plot corsika")
parser.add_argument('-v','--version', default='23c', type=str, dest = 'VER')
parser.add_argument('--emin', default=1e3, dest = 'EMIN')
parser.add_argument('--emax', default=1e7, dest = 'EMAX')
args = parser.parse_args()

plot_dir="/home/msilva/public_html/corsika_muons_Sibyll"+args.VER+"/2d_plots/"
if not os.path.exists(plot_dir):
       print "path doesn't exist. trying to make"
       os.makedirs(plot_dir)

start_time = time.asctime()
print 'Started:', start_time


file_corsika = 'corsika_'+args.VER+'_array.npy'
print "loading mc: " + file_corsika
array_corsika = np.load(file_corsika)
cor_eventinfo = array_corsika[0]
cor_muon1info = array_corsika[1]
cor_muon2info = array_corsika[2]
cor_weights=np.array(cor_eventinfo[0]); cor_multi=np.array(cor_eventinfo[1]); cor_totE=np.array(cor_eventinfo[2]);
cor_errors=np.array(cor_eventinfo[3]); cor_rad12 = np.array(cor_eventinfo[4]);
cor_type = np.array(cor_eventinfo[5]); cor_energy = np.array(cor_eventinfo[6]);
cor_energy1=np.array(cor_muon1info[0]); cor_zenith1=np.array(np.cos(cor_muon1info[1]));
cor_azimuth1=np.array(cor_muon1info[2]); cor_zpos1=np.array(cor_muon1info[5]); 
cor_energy2=np.array(cor_muon2info[0]); cor_zenith2=np.array(np.cos(cor_muon2info[1]));
cor_azimuth2=np.array(cor_muon2info[2]); cor_zpos2=np.array(cor_muon2info[5]); 

#plot primary energy vs multiplicity
mbins = range(31);
ebins=30; emin=args.EMIN; emax = args.EMAX;
logxbins = np.logspace(np.log10(emin),np.log10(emax), ebins)
cor_2Dhist, xbins, ybins = np.histogram2d(cor_energy, cor_multi, bins=[logxbins, mbins], weights=cor_weights)
x_bincenters = 0.5*(xbins[1:]+xbins[:-1])
y_bincenters = 0.5*(ybins[1:]+ybins[:-1])
fig, ax = plt.subplots()
im = ax.pcolormesh(xbins, ybins, cor_2Dhist.T, norm=LogNorm(vmin=1e-3, vmax=cor_2Dhist.max()) )
ax.set_xscale('log')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Event Rate [Hz]', labelpad=25, rotation=270)
ax.set_xlabel('Primary Energy [GeV]')
ax.set_ylabel('Muon Mutliplicity')
plt.savefig(plot_dir+"prim_energy_vs_muon_multi.pdf")
print "Made file: " + plot_dir+"prim_energy_vs_muon_multi.pdf"

#plot primary energy vs total muon energy
cor_2Dhist, xbins, ybins = np.histogram2d(cor_energy, cor_totE, bins=[logxbins, logxbins], weights=cor_weights)
x_bincenters = 0.5*(xbins[1:]+xbins[:-1])
y_bincenters = 0.5*(ybins[1:]+ybins[:-1])
fig, ax = plt.subplots()
im = ax.pcolormesh(xbins, ybins, cor_2Dhist.T, norm=LogNorm(vmin=1e-3, vmax=cor_2Dhist.max()) )
ax.set_xscale('log')
ax.set_yscale('log')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Event Rate [Hz]', labelpad=25, rotation=270)
ax.set_xlabel('Primary Energy [GeV]')
ax.set_ylabel('Total Muon Energy [GeV]')
plt.savefig(plot_dir+"prim_energy_vs_total_muon_energy.pdf")
print "Made file: " + plot_dir+"prim_energy_vs_muon_multi.pdf"
