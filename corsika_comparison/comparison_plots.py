#!/usr/bin/env python
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(description = "Plot corsika vs muongun comparisons")
parser.add_argument('-n','--numbermuons', default='1', type=int, dest = 'NUM')

plot_dir="/home/msilva/public_html/corsika_vs_muongun/"

def ks_w2(data1, data2, wei1, wei2):
    ix1 = np.argsort(data1)
    ix2 = np.argsort(data2)
    data1 = data1[ix1]
    data2 = data2[ix2]
    wei1 = wei1[ix1]
    wei2 = wei2[ix2]
    data = np.concatenate([data1, data2])
    cwei1 = np.hstack([0, np.cumsum(wei1)/sum(wei1)])
    cwei2 = np.hstack([0, np.cumsum(wei2)/sum(wei2)])
    cdf1we = cwei1[[np.searchsorted(data1, data, side='right')]]
    cdf2we = cwei2[[np.searchsorted(data2, data, side='right')]]
    return np.max(np.abs(cdf1we - cdf2we))

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
cor_weights=np.array(cor_eventinfo[0]); cor_multi=np.array(cor_eventinfo[1]); cor_totE=np.array(cor_eventinfo[2]);
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
mug_weights=np.array(mug_eventinfo[0]); mug_multi=np.array(mug_eventinfo[1]); mug_totE=np.array(mug_eventinfo[2]);
mug_energy1=np.array(mug_muon1info[0]); mug_zenith1=np.array(np.cos(mug_muon1info[1]));
mug_azimuth1=np.array(mug_muon1info[2]); mug_zpos1=np.array(mug_muon1info[5]); mug_rpos1=np.array(mug_muon1info[6]);
mug_energy2=np.array(mug_muon2info[0]); mug_zenith2=np.array(np.cos(mug_muon2info[1]));
mug_azimuth2=np.array(mug_muon2info[2]); mug_zpos2=np.array(mug_muon2info[5]); mug_rpos2=np.array(mug_muon2info[6]);

#energy inclusive
xbins=50; emin=3e3; emax = 7e4;
logxbins = np.logspace(np.log10(emin),np.log10(emax), xbins)
print "starting loop through multiplicity bins"
for mbins in range(1,6):
    cor_temp_location_mbins=np.where(np.logical_and(cor_multi>=mbins, cor_multi<mbins+1))
    cor_location_mbins=np.array(cor_temp_location_mbins[0])
    cor_weights_mbin=cor_weights[cor_location_mbins]
    cor_energy1_mbin=cor_energy1[cor_location_mbins]; 
    cor_weights_mbin=cor_weights[cor_location_mbins]
    cor_totE_mbin=cor_totE[cor_location_mbins]
    cor_E1ratio_mbin=np.divide(cor_energy1_mbin, cor_totE_mbin, out=np.zeros_like(cor_totE_mbin), where=cor_totE_mbin!=0)
    mug_temp_location_mbins=np.where(np.logical_and(mug_multi>=mbins, mug_multi<mbins+1))
    mug_location_mbins=np.array(mug_temp_location_mbins[0])
    mug_weights_mbin=mug_weights[mug_location_mbins]
    mug_energy1_mbin=mug_energy1[mug_location_mbins];
    mug_totE_mbin=mug_totE[mug_location_mbins]
    mug_E1ratio_mbin=np.divide(mug_energy1_mbin, mug_totE_mbin, out=np.zeros_like(mug_totE_mbin), where=mug_totE_mbin!=0)

    #total energy
    print "plotting distributions for multiplicity = " + str(mbins)
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_totenergy_1, cor_edges_1   = np.histogram(cor_totE_mbin, weights=cor_weights_mbin, density=False, bins=logxbins, range=(emin,emax) )
    mug_totenergy_1, mug_edges_1   = np.histogram(mug_totE_mbin, weights=mug_weights_mbin, density=False, bins=logxbins, range=(emin,emax) )
    ratio_totenergy  = np.divide(cor_totenergy_1, mug_totenergy_1, out=np.zeros_like(mug_totenergy_1), where=mug_totenergy_1!=0)
    ax1.semilogx(cor_edges_1[:-1], cor_totenergy_1, label='N='+str(mbins)+', Corsika', ls='-', color='b')
    ax1.semilogx(mug_edges_1[:-1], mug_totenergy_1,  label='N='+str(mbins)+', MuonGun', ls='-', color='r')
    ax2.semilogx(mug_edges_1[:-1], ratio_totenergy, color='b')
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    #ks, pval = stats.ks_2samp(cor_totE_mbin, mug_totE_mbin, cor_weights_mbin, mug_weights_mbin)
    ks = ks_w2(cor_totE_mbin, mug_totE_mbin, cor_weights_mbin, mug_weights_mbin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
    ax1.set_yscale('log')
    ax1.legend(loc='upper right')
    ax1.set_ylabel('Events per Second (Normed)')
    ax2.set_ylabel('Corsika/MuonGun')
    ax2.set_xlabel('Sum of Muons Energy [GeV]')
    plt.savefig(plot_dir+"corsika_div_muongun_sum_energy_binned_multi_eq_"+str(mbins)+".pdf")
    plt.clf()
    exit(0)
    #leading energy
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_energy_1, cor_edges_1   = np.histogram(cor_energy1_mbin, weights=cor_weights_mbin, density=True, bins=logxbins, range=(emin,emax) )
    mug_energy_1, mug_edges_1   = np.histogram(mug_energy1_mbin, weights=mug_weights_mbin, density=True, bins=logxbins, range=(emin,emax) )
    ratio_energy  = np.divide(cor_energy_1, mug_energy_1, out=np.zeros_like(mug_energy_1), where=mug_energy_1!=0)
    ax1.semilogx(cor_edges_1[:-1], cor_energy_1, label='N='+str(mbins)+', Corsika', ls='-', color='b')
    ax1.semilogx(mug_edges_1[:-1], mug_energy_1,  label='N='+str(mbins)+', MuonGun', ls='-', color='r')
    ax2.semilogx(mug_edges_1[:-1], ratio_energy, color='b')
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    ks = ks_w2(cor_energy1_mbin, mug_energy1_mbin, cor_weights_mbin, mug_weights_mbin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
    ax1.set_yscale('log')
    ax1.legend(loc='upper right')
    ax1.set_ylabel('Events per Second (Normed)')
    ax2.set_ylabel('Corsika/MuonGun')
    ax2.set_xlabel('Leading Muon Energy [GeV]')
    plt.savefig(plot_dir+"corsika_div_muongun_leading_energy_binned_multi_eq_"+str(mbins)+".pdf")
    plt.clf()

    #leading energy ratios
    xbins=50; xmin=0; xmax = 1;
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_Eratio_1, cor_edges_1   = np.histogram(cor_E1ratio_mbin, weights=cor_weights_mbin, density=True, bins=xbins, range=(xmin,xmax) )
    mug_Eratio_1, mug_edges_1   = np.histogram(mug_E1ratio_mbin, weights=mug_weights_mbin, density=True, bins=xbins, range=(xmin,xmax) )
    ratio_Eratio1  = np.divide(cor_Eratio_1, mug_Eratio_1, out=np.zeros_like(mug_Eratio_1), where=mug_Eratio_1!=0)
    ax1.hist(cor_edges_1[:-1], weights=cor_Eratio_1, label='N='+str(mbins)+', Corsika',
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1.hist(mug_edges_1[:-1], weights=mug_Eratio_1,  label='N='+str(mbins)+', MuonGun',
             bins=xbins, histtype='step', linestyle=('solid'), color=('red'))
    ax2.hist(mug_edges_1[:-1], weights=ratio_Eratio1,
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ks = ks_w2(cor_E1ratio_mbin, mug_E1ratio_mbin, cor_weights_mbin, mug_weights_mbin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
    #ax1.set_yscale('log')
    ax1.legend(loc='upper left')
    ax1.set_ylabel('Events per Second (Normed)')
    ax2.set_ylabel('Corsika/MuonGun')
    ax2.set_xlabel('Leading Muon Energy / All Muons')
    plt.savefig(plot_dir+"corsika_div_muongun_Energy1_div_totE_binned_multi_eq_"+str(mbins)+".pdf")
    plt.clf()


#now bin by energy
xbins=10; emin=3e3; emax = 7e4;
logxbins = np.logspace(np.log10(emin),np.log10(emax), xbins)
for ebins in range(len(logxbins)-1):
    #first find bin locations
    cor_temp_location_ebins=np.where(np.logical_and(cor_totE>=logxbins[ebins], cor_totE<=logxbins[ebins+1]))
    cor_location_ebins = np.array(cor_temp_location_ebins[0])
    cor_multi_ebin=cor_multi[cor_location_ebins]; cor_weights_ebin=cor_weights[cor_location_ebins]
    mug_temp_location_ebins=np.where(np.logical_and(mug_totE>=logxbins[ebins], mug_totE<=logxbins[ebins+1]))
    mug_location_ebins = np.array(mug_temp_location_ebins[0])
    mug_multi_ebin=mug_multi[mug_location_ebins]; mug_weights_ebin=mug_weights[mug_location_ebins]
    cor_zenith1_ebin=cor_zenith1[cor_location_ebins];
    mug_zenith1_ebin=mug_zenith1[mug_location_ebins];
    cor_azimuth1_ebin=cor_azimuth1[cor_location_ebins];
    mug_azimuth1_ebin=mug_azimuth1[mug_location_ebins];
    cor_rpos1_ebin=cor_rpos1[cor_location_ebins];
    mug_rpos1_ebin=mug_rpos1[mug_location_ebins];

    #multiplicity
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
    ks = ks_w2(cor_multi_ebin, mug_multi_ebin, cor_weights_ebin, mug_weights_ebin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
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
    ax1.hist(cor_edges_1[1:], weights=cor_zenith_1, label=r'N$\geq1$, Corsika',
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1.hist(mug_edges_1[1:], weights=mug_zenith_1,  label=r'N$\geq1$, MuonGun',
             bins=xbins, histtype='step', linestyle=('solid'), color=('red'))
    ax2.hist(mug_edges_1[1:], weights=ratio_zenith,
             bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    ks = ks_w2(cor_zenith1_ebin, mug_zenith1_ebin, cor_weights_ebin, mug_weights_ebin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
    ax1.set_title("Binned for range: "+lowE + " < E1 < " + highE + " TeV")
    ax1.set_yscale('log')
    ax1.legend(loc='upper left')
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
    ks = ks_w2(cor_azimuth1_ebin, mug_azimuth1_ebin, cor_weights_ebin, mug_weights_ebin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
    ax1.legend(loc='lower left')
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
