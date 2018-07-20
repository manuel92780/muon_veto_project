#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(description = "Plot corsika vs muongun comparisons")
parser.add_argument('-n','--numbermuons', default='1', type=int, dest = 'NUM')

plot_dir="/home/msilva/public_html/corsika_vs_muongun/"

def prepare_ratio_plot(top, bottom):
    top.grid(linestyle=':', linewidth=0.5)
    bottom.grid(linestyle=':', linewidth=0.5)
    bottom.set_ylim(0.0, 2.0)
    bottom.yaxis.set_ticks(np.arange(0, 2.0, 0.2))

    return top, bottom
#load corsika MC data
file_corsika = 'corsika_array.npy'
array_corsika = np.load(file_corsika)
cor_eventinfo = array_corsika[0]
cor_muon1info = array_corsika[1]
cor_muon2info = array_corsika[2]
cor_weights=np.array(cor_eventinfo[0]); cor_multi=np.array(cor_eventinfo[1]);
cor_energy1=np.array(cor_muon1info[0]); cor_zenith1=np.array(np.cos(cor_muon1info[1]));
cor_azimuth1=np.array(cor_muon1info[2]); cor_zpos1=np.array(cor_muon1info[5]); cor_rpos1=np.array(cor_muon1info[6]);
cor_energy2=np.array(cor_muon2info[0]); cor_zenith2=np.array(np.cos(cor_muon2info[1]));
cor_azimuth2=np.array(cor_muon2info[2]); cor_zpos2=np.array(cor_muon2info[5]); cor_rpos2=np.array(cor_muon2info[6]);

#load muongun MC data
file_muongun = 'muongun_array.npy'
array_muongun = np.load(file_muongun)
mug_eventinfo = array_muongun[0]
mug_muon1info = array_muongun[1]
mug_muon2info = array_muongun[2]
mug_weights=np.array(mug_eventinfo[0]); mug_multi=np.array(mug_eventinfo[1]);
mug_energy1=np.array(mug_muon1info[0]); mug_zenith1=np.array(np.cos(mug_muon1info[1]));
mug_azimuth1=np.array(mug_muon1info[2]); mug_zpos1=np.array(mug_muon1info[5]); mug_rpos1=np.array(mug_muon1info[6]);
mug_energy2=np.array(mug_muon2info[0]); mug_zenith2=np.array(np.cos(mug_muon2info[1]));
mug_azimuth2=np.array(mug_muon2info[2]); mug_zpos2=np.array(mug_muon2info[5]); mug_rpos2=np.array(mug_muon2info[6]);

#energy inclusive
xbins=50; xmin=1e2; xmax = 1e5;
logxbins = np.logspace(np.log10(xmin),np.log10(xmax), xbins)
for mbins in range(1,11):
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_temp_location_mbins=np.where(np.logical_and(cor_multi>=mbins, cor_multi<mbins+1))
    cor_location_mbins=np.array(cor_temp_location_mbins[0])
    cor_energy1_mbin=cor_energy1[cor_location_mbins]; cor_weights_mbin=cor_weights[cor_location_mbins]
    mug_temp_location_mbins=np.where(np.logical_and(mug_multi>=mbins, mug_multi<mbins+1))
    mug_location_mbins=np.array(mug_temp_location_mbins[0])
    mug_energy1_mbin=mug_energy1[mug_location_mbins]; mug_weights_mbin=mug_weights[mug_location_mbins]
    cor_energy_1, cor_edges_1   = np.histogram(cor_energy1_mbin, weights=cor_weights_mbin, density=False, bins=logxbins, range=(xmin,xmax) )
    mug_energy_1, mug_edges_1   = np.histogram(mug_energy1_mbin, weights=mug_weights_mbin, density=False, bins=logxbins, range=(xmin,xmax) )
    ratio_energy  = np.divide(cor_energy_1, mug_energy_1, out=np.zeros_like(mug_energy_1), where=mug_energy_1!=0)
    ax1.semilogx(cor_edges_1[:-1], cor_energy_1, label='N='+str(mbins)+', Corsika', ls='-', color='b')
    ax1.semilogx(mug_edges_1[:-1], mug_energy_1,  label='N='+str(mbins)+', MuonGun', ls='-', color='r')
    ax2.semilogx(mug_edges_1[:-1], ratio_energy, color='b')
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    ax1.set_yscale('log')
    ax1.legend(loc='upper right')
    ax1.set_ylabel('Events per Second (Normed)')
    ax2.set_ylabel('Corsika/MuonGun')
    ax2.set_xlabel('Leading Muon Energy [GeV]')
    plt.savefig(plot_dir+"corsika_div_muongun_energy_binned_multi_eq_"+str(mbins)+".pdf")
    plt.clf()
    exit(0);

#now bin by energy
logxbins = np.logspace(np.log10(xmin),np.log10(xmax), 10)
for ebins in range(len(logxbins)-1):
    #first find bin locations
    cor_temp_location_ebins=np.where(np.logical_and(cor_energy1>=logxbins[ebins], cor_energy1<=logxbins[ebins+1]))
    cor_location_ebins = np.array(cor_temp_location_ebins[0])
    cor_multi_ebin=cor_multi[cor_location_ebins]; cor_weights_ebin=cor_weights[cor_location_ebins]
    mug_temp_location_ebins=np.where(np.logical_and(mug_energy1>=logxbins[ebins], mug_energy1<=logxbins[ebins+1]))
    mug_location_ebins = np.array(mug_temp_location_ebins[0])
    mug_multi_ebin=mug_multi[mug_location_ebins]; mug_weights_ebin=mug_weights[mug_location_ebins]
    cor_zenith1_ebin=cor_zenith1[cor_location_ebins];
    mug_zenith1_ebin=mug_zenith1[mug_location_ebins];
    cor_azimuth1_ebin=cor_azimuth1[cor_location_ebins];
    mug_azimuth1_ebin=mug_azimuth1[mug_location_ebins];
    cor_rpos1_ebin=cor_rpos1[cor_location_ebins];
    mug_rpos1_ebin=mug_rpos1[mug_location_ebins];

    xbins=50; xmin=0; xmax = 50;
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_multi_1, cor_edges_1   = np.histogram(cor_multi_ebin, weights=cor_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
    mug_multi_1, mug_edges_1   = np.histogram(mug_multi_ebin, weights=mug_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
    ratio_multi  = np.divide(cor_multi_1, mug_multi_1, out=np.zeros_like(mug_multi_1), where=mug_multi_1!=0)
    ax1.hist(cor_edges_1[:-1], weights=cor_multi_1, label=r'Corsika',
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1.hist(mug_edges_1[:-1], weights=mug_multi_1,  label=r'MuonGun',
             bins=xbins, histtype='step', linestyle=('solid'), color=('red'))
    ax2.hist(mug_edges_1[:-1], weights=ratio_multi,
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    lowE = "%.2f" % (logxbins[ebins]/1000); highE = "%.2f" % (logxbins[ebins+1]/1000);
    ax1.set_title("Binned for range: "+lowE + " < E1 < " + highE + " TeV")
    ax1.set_yscale('log')
    ax1.legend(loc='upper right')
    ax1.set_ylabel('Events per Second (Normed)')
    ax2.set_ylabel('Corsika/MuonGun')
    ax2.set_xlabel('Event Multiplicity')
    plt.savefig(plot_dir+"corsika_div_muongun_multiplicity_binned_energy_"+lowE+"_"+highE+"TeV.pdf")
    plt.clf()
    #zenith now
    xbins=30; xmin=0; xmax = 1;
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_zenith_1, cor_edges_1   = np.histogram(cor_zenith1_ebin, weights=cor_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
    mug_zenith_1, mug_edges_1   = np.histogram(mug_zenith1_ebin, weights=mug_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
    ratio_zenith  = np.divide(cor_zenith_1, mug_zenith_1, out=np.zeros_like(mug_zenith_1), where=mug_zenith_1!=0)
    ax1.hist(cor_edges_1[:-1], weights=cor_zenith_1, label=r'N$\geq1$, Corsika',
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1.hist(mug_edges_1[:-1], weights=mug_zenith_1,  label=r'N$\geq1$, MuonGun',
             bins=xbins, histtype='step', linestyle=('solid'), color=('red'))
    ax2.hist(mug_edges_1[:-1], weights=ratio_zenith,
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    ax1.set_title("Binned for range: "+lowE + " < E1 < " + highE + " TeV")
    ax1.set_yscale('log')
    ax1.legend(loc='upper right')
    ax1.set_title("Binned for range: "+lowE + " < E1 < " + highE + " TeV")
    ax1.set_ylabel('Events per Second (Normed)')
    ax2.set_ylabel('Corsika/MuonGun')
    ax2.set_xlabel('Leading Muon cos(Zenith)')
    plt.savefig(plot_dir+"corsika_div_muongun_zenith_binned_energy_"+lowE+"_"+highE+"TeV..pdf")
    plt.clf()
    #azimuth now
    xbins=30; xmin=0; xmax = 2*np.pi;
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_azimuth_1, cor_edges_1   = np.histogram(cor_azimuth1_ebin, weights=cor_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
    mug_azimuth_1, mug_edges_1   = np.histogram(mug_azimuth1_ebin, weights=mug_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
    ratio_azimuth  = np.divide(cor_azimuth_1, mug_azimuth_1, out=np.zeros_like(mug_azimuth_1), where=mug_azimuth_1!=0)
    ax1.hist(cor_edges_1[:-1], weights=cor_azimuth_1, label=r'N$\geq1$, Corsika',
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1.hist(mug_edges_1[:-1], weights=mug_azimuth_1,  label=r'N$\geq1$, MuonGun',
             bins=xbins, histtype='step', linestyle=('solid'), color=('red'))
    ax2.hist(mug_edges_1[:-1], weights=ratio_azimuth,
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    ax1.legend(loc='upper right')
    ax1.set_title("Binned for range: "+lowE + " < E1 < " + highE + " TeV")
    ax1.set_ylabel('Events per Second (Normed)')
    ax2.set_ylabel('Corsika/MuonGun')
    ax2.set_xlabel('Leading Muon Azimuthal [rad]')
    plt.savefig(plot_dir+"corsika_div_muongun_azimuth_binned_energy_"+lowE+"_"+highE+"TeV..pdf")
    plt.clf()
    #lateral distance now
#    xbins=30; xmin=0; xmax = 50;
#    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
#    print cor_rpos1_ebin
#    print mug_rpos1_ebin
#    cor_rpos_1, cor_edges_1   = np.histogram(cor_rpos1_ebin, weights=cor_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
#    mug_rpos_1, mug_edges_1   = np.histogram(mug_rpos1_ebin, weights=mug_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
#    ratio_rpos  = np.divide(cor_rpos_1, mug_rpos_1, out=np.zeros_like(mug_rpos_1), where=mug_rpos_1!=0)
#    ax1.hist(cor_edges_1[:-1], weights=cor_rpos_1, label=r'N$\geq1$, Corsika',
#             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
#    ax1.hist(mug_edges_1[:-1], weights=mug_rpos_1,  label=r'N$\geq1$, MuonGun',
#             bins=xbins, histtype='step', linestyle=('solid'), color=('red'))
#    ax2.hist(mug_edges_1[:-1], weights=ratio_rpos,
#             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
#    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
#    ax1.set_title("Binned for range: "+lowE + " < E1 < " + highE + " TeV")
#    ax1.legend(loc='upper right')
#    ax1.set_title("Binned for range: "+lowE + " < E1 < " + highE + " TeV")
#    ax1.set_ylabel('Events per Second (Normed)')
#    ax2.set_ylabel('Corsika/MuonGun')
#    ax2.set_xlabel('Leading Muon Lateral Distance [m]')
#    plt.savefig(plot_dir+"corsika_div_muongun_rpos_binned_energy_"+lowE+"_"+highE+"TeV..pdf")
#    plt.clf()
