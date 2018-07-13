#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(description = "Plot corsika vs muongun comparisons")
parser.add_argument('-n','--numbermuons', default='1', type=int, dest = 'NUM')


def prepare_ratio_plot(top, bottom):
    top.grid(linestyle=':', linewidth=0.5)
    bottom.grid(linestyle=':', linewidth=0.5)
    bottom.set_ylim(0.0, 2.0)
    bottom.yaxis.set_ticks(np.arange(0, 2.0, 0.2))

    return top, bottom
#load MC data
file_corsika = 'corsika_array.npy'
array_corsika = np.load(file_corsika)
cor_eventinfo = array_corsika[0]
cor_muon1info = array_corsika[1]
cor_muon2info = array_corsika[2]
file_muongun = 'muongun_array.npy'
array_muongun = np.load(file_muongun)
mug_eventinfo = array_muongun[0]
mug_muon1info = array_muongun[1]
mug_muon2info = array_muongun[2]

#ratio plots on ax1 and ax2
#muon multiplicity first
xbins=50; xmin=0; xmax = 50;
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
cor_multi_1, cor_edges_1   = np.histogram(cor_eventinfo[1], weights=cor_eventinfo[0], density=True, bins=xbins, range=(xmin,xmax) )
mug_multi_1, mug_edges_1   = np.histogram(mug_eventinfo[1], weights=mug_eventinfo[0], density=True, bins=xbins, range=(xmin,xmax) )
ratio_multi  = np.divide(cor_multi_1, mug_multi_1, out=np.zeros_like(mug_multi_1), where=mug_multi_1!=0)
ax1.hist(cor_edges_1[:-1], weights=cor_multi_1, label=r'Corsika',
         bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
ax1.hist(mug_edges_1[:-1], weights=mug_multi_1,  label=r'MuonGun',
         bins=xbins, histtype='step', linestyle=('solid'), color=('red'))
ax2.hist(mug_edges_1[:-1], weights=ratio_multi,
         bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
ax1, ax2 = prepare_ratio_plot(ax1, ax2)
ax1.legend(loc='upper right')
ax1.set_ylabel('Events per Second (Normed)')
ax2.set_ylabel('Corsika/MuonGun')
ax2.set_xlabel('Event Multiplicity')
plt.savefig("corsika_div_muongun_multiplicity.pdf")
plt.clf()

#energy second, N > 0 
xbins=50; xmin=1e1; xmax = 1e6;
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
logxbins = np.logspace(np.log10(xmin),np.log10(xmax), xbins)
cor_energy_1, cor_edges_1   = np.histogram(cor_muon1info[0], weights=cor_eventinfo[0], density=True, bins=logxbins, range=(xmin,xmax) )
mug_energy_1, mug_edges_1   = np.histogram(mug_muon1info[0], weights=mug_eventinfo[0], density=True, bins=logxbins, range=(xmin,xmax) )
ratio_energy  = np.divide(cor_energy_1, mug_energy_1, out=np.zeros_like(mug_energy_1), where=mug_energy_1!=0)
ax1.semilogx(cor_edges_1[:-1], cor_energy_1, label=r'N$\geq1$, Corsika', ls='-', color='b')
ax1.semilogx(mug_edges_1[:-1], mug_energy_1,  label=r'N$\geq1$, MuonGun', ls='-', color='r')
ax2.semilogx(mug_edges_1[:-1], ratio_energy, color='b')
ax1, ax2 = prepare_ratio_plot(ax1, ax2)
ax1.legend(loc='upper right')
ax1.set_ylabel('Events per Second (Normed)')
ax2.set_ylabel('Corsika/MuonGun')
ax2.set_xlabel('Leading Muon Energy [GeV]')
plt.savefig("corsika_div_muongun_energy_n_geq_1.pdf")
plt.clf()

#cos(zenith) next
xbins=50; xmin=0; xmax = 1;
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
cor_zenith_1, cor_edges_1   = np.histogram(np.cos(cor_muon1info[1]), weights=cor_eventinfo[0], density=True, bins=xbins, range=(xmin,xmax) )
mug_zenith_1, mug_edges_1   = np.histogram(np.cos(mug_muon1info[1]), weights=mug_eventinfo[0], density=True, bins=xbins, range=(xmin,xmax) )
ratio_zenith  = np.divide(cor_zenith_1, mug_zenith_1, out=np.zeros_like(mug_zenith_1), where=mug_zenith_1!=0)
ax1.hist(cor_edges_1[:-1], weights=cor_zenith_1, label=r'N$\geq1$, Corsika',
         bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
ax1.hist(mug_edges_1[:-1], weights=mug_zenith_1,  label=r'N$\geq1$, MuonGun',
         bins=xbins, histtype='step', linestyle=('solid'), color=('red'))
ax2.hist(mug_edges_1[:-1], weights=ratio_zenith,
         bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
ax1, ax2 = prepare_ratio_plot(ax1, ax2)
ax1.legend(loc='upper right')
ax1.set_ylabel('Events per Second (Normed)')
ax2.set_ylabel('Corsika/MuonGun')
ax2.set_xlabel('Leading Muon cos(Zenith)')
plt.savefig("corsika_div_muongun_zenith_n_geq_1.pdf")
#plt.show()
plt.clf()
