#!/usr/bin/env python
import os
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(description = "Plot corsika vs muongun comparisons")
parser.add_argument('-n','--numbermuons', default='1', type=int, dest = 'NUM')
parser.add_argument('-v','--version', default='23c', type=str, dest = 'VER')
parser.add_argument('--emin', default=1e3, dest = 'EMIN')
parser.add_argument('--emax', default=1e5, dest = 'EMAX')
args = parser.parse_args()

plot_dir="/home/msilva/public_html/corsika21_vs_corsika2x/Sibyll"+args.VER+"/multiplicity_energy/"
if not os.path.exists(plot_dir):
       print "path doesn't exist. trying to make"
       os.makedirs(plot_dir)

#ks test with weights
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

#ratio plots tinkering 
def prepare_ratio_plot(top, bottom):
    top.grid(linestyle=':', linewidth=0.5)
    bottom.grid(linestyle=':', linewidth=0.5)
    bottom.set_ylim(0.0, 2.0)
    bottom.yaxis.set_ticks(np.arange(0, 2.0, 0.2))
    return top, bottom

#propagate errors
def ratio_errors(num, num_err, denom, denom_err):
    yvals  = np.divide(num, denom, out=np.zeros_like(denom), where=denom!=0)
    err_num = np.divide(num_err, num, out=np.zeros_like(num), where=num!=0)
    err_denom = np.divide(denom_err, denom, out=np.zeros_like(denom), where=denom!=0)
    yval_errs = yvals*np.sqrt(np.power(err_num,2) + np.power(err_denom,2))
    num = num/np.sum(num); denom = denom/np.sum(denom);
    yvals_norm = np.divide(num, denom, out=np.zeros_like(denom), where=denom!=0)
    return yvals, yval_errs, yvals_norm


#load corsika MC data
file_corsika = 'corsika_'+args.VER+'_array.npy'
array_corsika = np.load(file_corsika)
cor_eventinfo = array_corsika[0]
cor_muon1info = array_corsika[1]
cor_muon2info = array_corsika[2]
cor_weights=np.array(cor_eventinfo[0]); cor_multi=np.array(cor_eventinfo[1]); cor_totE=np.array(cor_eventinfo[2]);
cor_weights_2=np.array(cor_eventinfo[3]);
cor_energy1=np.array(cor_muon1info[0]); cor_zenith1=np.array(np.cos(cor_muon1info[1]));
cor_azimuth1=np.array(cor_muon1info[2]); cor_zpos1=np.array(cor_muon1info[5]); 
cor_energy2=np.array(cor_muon2info[0]); cor_zenith2=np.array(np.cos(cor_muon2info[1]));
cor_azimuth2=np.array(cor_muon2info[2]); cor_zpos2=np.array(cor_muon2info[5]); 

#load muongun MC data
file_muongun = 'corsika_21_array.npy'
array_muongun = np.load(file_muongun)
mug_eventinfo = array_muongun[0]
mug_muon1info = array_muongun[1]
mug_muon2info = array_muongun[2]
mug_weights=np.array(mug_eventinfo[0]); mug_multi=np.array(mug_eventinfo[1]); mug_totE=np.array(mug_eventinfo[2]);
mug_weights_2=np.array(mug_eventinfo[3]);
mug_energy1=np.array(mug_muon1info[0]); mug_zenith1=np.array(np.cos(mug_muon1info[1]));
mug_azimuth1=np.array(mug_muon1info[2]); mug_zpos1=np.array(mug_muon1info[5]); 
mug_energy2=np.array(mug_muon2info[0]); mug_zenith2=np.array(np.cos(mug_muon2info[1]));
mug_azimuth2=np.array(mug_muon2info[2]); mug_zpos2=np.array(mug_muon2info[5]); 

#bin mutiplicity as function of energy alone
xbins=5; emin=args.EMIN; emax = args.EMAX;
logxbins = np.logspace(np.log10(emin),np.log10(emax), xbins)
print "starting loop through energy bins for multiplicity"
for ebins in range(len(logxbins)-1):
    print "plotting distributions for energy: " + str(logxbins[ebins]/1000.) + " - " + str(logxbins[ebins+1]/1000.) + " TeV"

    #first find bin locations
    cor_temp_location_ebins=np.where(np.logical_and(cor_totE>=logxbins[ebins], cor_totE<=logxbins[ebins+1]))
    cor_location_ebins = np.array(cor_temp_location_ebins[0])
    cor_multi_ebin=cor_multi[cor_location_ebins]; cor_weights_ebin=cor_weights[cor_location_ebins]
    cor_errors_ebin=cor_weights_2[cor_location_ebins]
    mug_temp_location_ebins=np.where(np.logical_and(mug_totE>=logxbins[ebins], mug_totE<=logxbins[ebins+1]))
    mug_location_ebins = np.array(mug_temp_location_ebins[0])
    mug_multi_ebin=mug_multi[mug_location_ebins]; mug_weights_ebin=mug_weights[mug_location_ebins]
    mug_errors_ebin=mug_weights_2[mug_location_ebins]

    #multiplicity
    xbins=50; xmin=0; xmax = 50;
    xbins=np.linspace(xmin,xmax,xbins,dtype='int')
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_multi_1, cor_edges_1   = np.histogram(cor_multi_ebin, weights=cor_weights_ebin, density=False, bins=xbins)
    mug_multi_1, mug_edges_1   = np.histogram(mug_multi_ebin, weights=mug_weights_ebin, density=False, bins=xbins)
    cor_errors_1, _   = np.histogram(cor_multi_ebin, weights=cor_errors_ebin, density=False, bins=xbins)
    mug_errors_1, _   = np.histogram(mug_multi_ebin, weights=mug_errors_ebin, density=False, bins=xbins)
    bincenters = 0.5*(cor_edges_1[1:]+cor_edges_1[:-1])
    ratio_multi, ratio_err, ratio_multi_norm  = ratio_errors(cor_multi_1, np.sqrt(cor_errors_1), mug_multi_1, np.sqrt(mug_errors_1))
    #ax1.hist(bincenters, weights=cor_multi_1, label='Corsika, '+args.VER,
    #         bins=xbins, histtype='step', linestyle=('solid'), color=('blue'))
    ax1.errorbar(bincenters, mug_multi_1/np.sum(mug_multi_1),  yerr = np.sqrt(mug_errors_1), label='Corsika, 2.1',
                 fmt='or', ecolor='r', drawstyle = 'steps-mid-')
    ax1.errorbar(bincenters, cor_multi_1/np.sum(cor_multi_1),  yerr = np.sqrt(cor_errors_1), label='Corsika, '+args.VER,
                 fmt='.b', ecolor='b', drawstyle = 'steps-mid-')
    ax2.errorbar(bincenters, ratio_multi_norm, yerr= ratio_err,
                fmt='ob', ecolor='b', drawstyle = 'steps-mid-')
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    ks = ks_w2(cor_multi_ebin, mug_multi_ebin, cor_weights_ebin, mug_weights_ebin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
    lowE = "%.2f" % (logxbins[ebins]/1000); highE = "%.2f" % (logxbins[ebins+1]/1000);
    ax1.set_title("Binned for range: "+lowE + " < Etot < " + highE + " TeV")
    ax1.set_yscale('log')
    ax1.legend(loc='upper right')
    ax1.set_ylabel('Events per second [normalized]')
    ax2.set_ylabel('Corsika 2.x / 2.1')
    ax2.set_xlabel('Event Multiplicity')
    plt.savefig(plot_dir+"corsika21_div_corsika"+args.VER+"__multiplicity_binned_energy_"+lowE+"_"+highE+"TeV.pdf")
    plt.clf()

#starting a big loop
xbins=30; emin=args.EMIN; emax = args.EMAX;
logxbins = np.logspace(np.log10(emin),np.log10(emax), xbins)
print "starting loop through multiplicity bins"
for mbins in range(1,3):
    cor_temp_location_mbins=np.where(np.logical_and(cor_multi>=mbins, cor_multi<mbins+1))
    cor_location_mbins=np.array(cor_temp_location_mbins[0])

    cor_weights_mbin=cor_weights[cor_location_mbins]
    cor_errors_mbin=cor_weights_2[cor_location_mbins]
    cor_energy1_mbin=cor_energy1[cor_location_mbins]; 
    cor_totE_mbin=cor_totE[cor_location_mbins]
    cor_E1ratio_mbin=np.divide(cor_energy1_mbin, cor_totE_mbin, out=np.zeros_like(cor_totE_mbin), where=cor_totE_mbin!=0)
    cor_zenith1_mbin=cor_zenith1[cor_location_mbins];
    cor_azimuth1_mbin=cor_azimuth1[cor_location_mbins];

    mug_temp_location_mbins=np.where(np.logical_and(mug_multi>=mbins, mug_multi<mbins+1))
    mug_location_mbins=np.array(mug_temp_location_mbins[0])

    mug_weights_mbin=mug_weights[mug_location_mbins]
    mug_errors_mbin=mug_weights_2[mug_location_mbins]
    mug_energy1_mbin=mug_energy1[mug_location_mbins];
    mug_totE_mbin=mug_totE[mug_location_mbins]
    mug_E1ratio_mbin=np.divide(mug_energy1_mbin, mug_totE_mbin, out=np.zeros_like(mug_totE_mbin), where=mug_totE_mbin!=0)
    mug_zenith1_mbin=mug_zenith1[mug_location_mbins];
    mug_azimuth1_mbin=mug_azimuth1[mug_location_mbins];

    #total energy
    print "plotting distributions for multiplicity = " + str(mbins)
    xbins=30; emin=args.EMIN; emax = args.EMAX;
    logxbins = np.logspace(np.log10(emin),np.log10(emax), xbins)
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_totenergy_1, cor_edges_1   = np.histogram(cor_totE_mbin, weights=cor_weights_mbin, density=False, bins=logxbins)
    mug_totenergy_1, mug_edges_1   = np.histogram(mug_totE_mbin, weights=mug_weights_mbin, density=False, bins=logxbins)
    cor_errors_1, _   = np.histogram(cor_totE_mbin, weights=cor_errors_mbin, density=False, bins=logxbins)
    mug_errors_1, _   = np.histogram(mug_totE_mbin, weights=mug_errors_mbin, density=False, bins=logxbins)
    bincenters = 10**(0.5*(np.log10(cor_edges_1[1:])+np.log10(cor_edges_1[:-1])))
    ratio_totenergy, ratio_err, ratio_multi_norm  = ratio_errors(cor_totenergy_1, np.sqrt(cor_errors_1), mug_totenergy_1, np.sqrt(mug_errors_1))
    ax1.errorbar(bincenters, mug_totenergy_1/np.sum(mug_totenergy_1),  yerr = np.sqrt(mug_errors_1), label='N='+str(mbins)+', Corsika, 21',
                 fmt='or', ecolor='r', drawstyle = 'steps-mid-')
    ax1.errorbar(bincenters, cor_totenergy_1/np.sum(cor_totenergy_1),  yerr = np.sqrt(cor_errors_1), label='N='+str(mbins)+', Corsika, '+args.VER,
                 fmt='.b', ecolor='b', drawstyle = 'steps-mid-')
    ax2.errorbar(bincenters, ratio_multi_norm, yerr= ratio_err,
                fmt='ob', ecolor='b', drawstyle = 'steps-mid-')
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    #ks, pval = stats.ks_2samp(cor_totE_mbin, mug_totE_mbin, cor_weights_mbin, mug_weights_mbin)
    ks = ks_w2(cor_totE_mbin, mug_totE_mbin, cor_weights_mbin, mug_weights_mbin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
    ax1.set_xscale("log")
    ax2.set_xscale("log")
    ax1.set_yscale('log')
    ax1.legend(loc='upper right')
    ax1.set_ylabel('Events per second [normalized]')
    ax2.set_ylabel('Corsika 2.x / 2.1')
    ax2.set_xlabel('Sum of Muons Energy [GeV]')
    plt.savefig(plot_dir+"corsika21_div_corsika"+args.VER+"_sum_energy_binned_multi_eq_"+str(mbins)+".pdf")
    plt.clf()

    #leading energy
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_energy_1, cor_edges_1   = np.histogram(cor_energy1_mbin, weights=cor_weights_mbin, density=False, bins=logxbins)
    mug_energy_1, mug_edges_1   = np.histogram(mug_energy1_mbin, weights=mug_weights_mbin, density=False, bins=logxbins)
    cor_errors_1, _   = np.histogram(cor_energy1_mbin, weights=cor_errors_mbin, density=False, bins=logxbins)
    mug_errors_1, _   = np.histogram(mug_energy1_mbin, weights=mug_errors_mbin, density=False, bins=logxbins)
    bincenters = 10**(0.5*(np.log10(cor_edges_1[1:])+np.log10(cor_edges_1[:-1])))
    ratio_energy, ratio_err, ratio_energy_norm  = ratio_errors(cor_energy_1, np.sqrt(cor_errors_1), mug_energy_1, np.sqrt(mug_errors_1))
    ax1.errorbar(bincenters, mug_energy_1/np.sum(mug_energy_1),  yerr = np.sqrt(mug_errors_1), label='N='+str(mbins)+', Corsika, 21',
                 fmt='or', ecolor='r', drawstyle = 'steps-mid-')
    ax1.errorbar(bincenters, cor_energy_1/np.sum(cor_energy_1),  yerr = np.sqrt(cor_errors_1), label='N='+str(mbins)+', Corsika, '+args.VER,
                 fmt='.b', ecolor='b', drawstyle = 'steps-mid-')
    ax2.errorbar(bincenters, ratio_energy_norm, yerr= ratio_err,
                fmt='ob', ecolor='b', drawstyle = 'steps-mid-')
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    ks = ks_w2(cor_energy1_mbin, mug_energy1_mbin, cor_weights_mbin, mug_weights_mbin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
    ax1.set_xscale("log")
    ax2.set_xscale("log")
    ax1.set_yscale('log')
    ax1.legend(loc='lower left')
    ax1.set_ylabel('Events per second [normalized]')
    ax2.set_ylabel('Corsika 2.x / 2.1')
    ax2.set_xlabel('Leading Muon Energy [GeV]')
    plt.savefig(plot_dir+"corsika21_div_corsika"+args.VER+"_leading_energy_binned_multi_eq_"+str(mbins)+".pdf")
    plt.clf()

    #leading energy ratios
    xbins=50; xmin=0; xmax = 1;
    xbins=np.linspace(xmin,xmax,xbins,dtype='float')
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    cor_Eratio_1, cor_edges_1   = np.histogram(cor_E1ratio_mbin, weights=cor_weights_mbin, density=False, bins=xbins)
    mug_Eratio_1, mug_edges_1   = np.histogram(mug_E1ratio_mbin, weights=mug_weights_mbin, density=False, bins=xbins)
    cor_errors_1, _   = np.histogram(cor_E1ratio_mbin, weights=cor_errors_mbin, density=False, bins=xbins)
    mug_errors_1, _   = np.histogram(mug_E1ratio_mbin, weights=mug_errors_mbin, density=False, bins=xbins)
    bincenters = 0.5*(cor_edges_1[1:]+cor_edges_1[:-1])
    ratio_Eratio_1, ratio_err, ratio_Eratio1_norm  = ratio_errors(cor_Eratio_1, np.sqrt(cor_errors_1), mug_Eratio_1, np.sqrt(mug_errors_1))
    ax1.errorbar(bincenters, mug_Eratio_1/np.sum(mug_Eratio_1),  yerr = np.sqrt(mug_errors_1), label='N='+str(mbins)+', Corsika, 21',
                 fmt='or', ecolor='r', drawstyle = 'steps-mid-')
    ax1.errorbar(bincenters, cor_Eratio_1/np.sum(cor_Eratio_1),  yerr = np.sqrt(cor_errors_1), label='N='+str(mbins)+', Corsika, '+args.VER,
                 fmt='.b', ecolor='b', drawstyle = 'steps-mid-')
    ax2.errorbar(bincenters, ratio_Eratio1_norm, yerr= ratio_err,
                fmt='ob', ecolor='b', drawstyle = 'steps-mid-')
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    ks = ks_w2(cor_E1ratio_mbin, mug_E1ratio_mbin, cor_weights_mbin, mug_weights_mbin)
    ax2.set_title('K-S statistic = %0.2f' % ks)
    ax1.legend(loc='upper left')
    ax1.set_ylabel('Events per second [normalized]')
    ax2.set_ylabel('Corsika 2.x / 2.1')
    ax2.set_xlabel('Leading Muon Energy / All Muons')
    plt.savefig(plot_dir+"corsika21_div_corsika"+args.VER+"_Energy1_div_totE_binned_multi_eq_"+str(mbins)+".pdf")
    plt.clf()

    #also bin by energy
    xbins=5; emin=args.EMIN; emax = args.EMAX;
    logxbins = np.logspace(np.log10(emin),np.log10(emax), xbins)
    for ebins in range(len(logxbins)-1):
        lowE = "%.2f" % (logxbins[ebins]/1000); highE = "%.2f" % (logxbins[ebins+1]/1000);
        print "plotting distributions for energy: " + str(logxbins[ebins]/1000.) + " - " + str(logxbins[ebins+1]/1000.) + " TeV"
        #first find bin locations
        cor_temp_location_ebins=np.where(np.logical_and(cor_totE_mbin>=logxbins[ebins], cor_totE_mbin<=logxbins[ebins+1]))
        cor_location_ebins = np.array(cor_temp_location_ebins[0])

        cor_weights_ebin=cor_weights_mbin[cor_location_ebins]
        cor_errors_ebin=cor_errors_mbin[cor_location_ebins]
        cor_zenith1_ebin=cor_zenith1_mbin[cor_location_ebins];
        cor_azimuth1_ebin=cor_azimuth1_mbin[cor_location_ebins];

        mug_temp_location_ebins=np.where(np.logical_and(mug_totE_mbin>=logxbins[ebins], mug_totE_mbin<=logxbins[ebins+1]))
        mug_location_ebins = np.array(mug_temp_location_ebins[0])

        mug_weights_ebin=mug_weights_mbin[mug_location_ebins]
        mug_errors_ebin=mug_errors_mbin[mug_location_ebins]
        mug_zenith1_ebin=mug_zenith1_mbin[mug_location_ebins];
        mug_azimuth1_ebin=mug_azimuth1_mbin[mug_location_ebins];
        
        #zenith now
        xbins=30; xmin=0; xmax = 1;
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
        cor_zenith_1, cor_edges_1   = np.histogram(cor_zenith1_ebin, weights=cor_weights_ebin, density=False, bins=xbins)
        mug_zenith_1, mug_edges_1   = np.histogram(mug_zenith1_ebin, weights=mug_weights_ebin, density=False, bins=xbins)
        cor_errors_1, _   = np.histogram(cor_zenith1_ebin, weights=cor_errors_ebin, density=False, bins=xbins)
        mug_errors_1, _   = np.histogram(mug_zenith1_ebin, weights=mug_errors_ebin, density=False, bins=xbins)
        bincenters = 0.5*(cor_edges_1[1:]+cor_edges_1[:-1])
        ratio_zenith_1, ratio_err, ratio_zenith1_norm  = ratio_errors(cor_zenith_1, np.sqrt(cor_errors_1), mug_zenith_1, np.sqrt(mug_errors_1))
        ax1.errorbar(bincenters, mug_zenith_1/np.sum(mug_zenith_1),  yerr = np.sqrt(mug_errors_1), label='N='+str(mbins)+', Corsika, 21',
                     fmt='or', ecolor='r', drawstyle = 'steps-mid-')
        ax1.errorbar(bincenters, cor_zenith_1/np.sum(cor_zenith_1),  yerr = np.sqrt(cor_errors_1), label='N='+str(mbins)+', Corsika, '+args.VER,
                     fmt='.b', ecolor='b', drawstyle = 'steps-mid-')
        ax2.errorbar(bincenters, ratio_zenith1_norm, yerr= ratio_err,
                     fmt='ob', ecolor='b', drawstyle = 'steps-mid-')
        ax1, ax2 = prepare_ratio_plot(ax1, ax2)
        ks = ks_w2(cor_zenith1_ebin, mug_zenith1_ebin, cor_weights_ebin, mug_weights_ebin)
        ax2.set_title('K-S statistic = %0.2f' % ks)
        ax1.set_title("Binned for range: "+lowE + " < Etot < " + highE + " TeV")
        ax1.set_yscale('log')
        ax1.legend(loc='upper left')
        ax1.set_title("Binned for range: "+lowE + " < Etot < " + highE + " TeV")
        ax1.set_ylabel('Events per second [normalized]')
        ax2.set_ylabel('Corsika 2.x / 2.1')
        ax2.set_xlabel('Leading Muon cos(Zenith)')
        plt.savefig(plot_dir+"corsika21_div_corsika"+args.VER+"_zenith_binned_energy_"+lowE+"_"+highE+"TeV_binned_multi_eq_"+str(mbins)+".pdf")
        plt.clf()

        #azimuth now
        xbins=30; xmin=0; xmax = 2*np.pi;
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
        cor_azimuth_1, cor_edges_1   = np.histogram(cor_azimuth1_ebin, weights=cor_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
        mug_azimuth_1, mug_edges_1   = np.histogram(mug_azimuth1_ebin, weights=mug_weights_ebin, density=True, bins=xbins, range=(xmin,xmax) )
        cor_errors_1, _   = np.histogram(cor_azimuth1_ebin, weights=cor_errors_ebin, density=False, bins=xbins)
        mug_errors_1, _   = np.histogram(mug_azimuth1_ebin, weights=mug_errors_ebin, density=False, bins=xbins)
        bincenters = 0.5*(cor_edges_1[1:]+cor_edges_1[:-1])
        ratio_azimuth_1, ratio_err, ratio_azimuth1_norm  = ratio_errors(cor_azimuth_1, np.sqrt(cor_errors_1), mug_azimuth_1, np.sqrt(mug_errors_1))
        ax1.errorbar(bincenters, mug_azimuth_1/np.sum(mug_azimuth_1),  yerr = np.sqrt(mug_errors_1), label='N='+str(mbins)+', Corsika, 21',
                     fmt='or', ecolor='r', drawstyle = 'steps-mid-')
        ax1.errorbar(bincenters, cor_azimuth_1/np.sum(cor_azimuth_1),  yerr = np.sqrt(cor_errors_1), label='N='+str(mbins)+', Corsika, '+args.VER,
                     fmt='.b', ecolor='b', drawstyle = 'steps-mid-')
        ax2.errorbar(bincenters, ratio_azimuth1_norm, yerr= ratio_err,
                     fmt='ob', ecolor='b', drawstyle = 'steps-mid-')

        ax1, ax2 = prepare_ratio_plot(ax1, ax2)
        ks = ks_w2(cor_azimuth1_ebin, mug_azimuth1_ebin, cor_weights_ebin, mug_weights_ebin)
        ax2.set_title('K-S statistic = %0.2f' % ks)
        ax1.legend(loc='lower left')
        ax1.set_title("Binned for range: "+lowE + " < Etot < " + highE + " TeV")
        ax1.set_ylabel('Events per second [normalized]')
        ax2.set_ylabel('Corsika 2.x / 2.1')
        ax2.set_xlabel('Leading Muon Azimuthal [rad]')
        plt.savefig(plot_dir+"corsika21_div_corsika"+args.VER+"_azimuth_binned_energy_"+lowE+"_"+highE+"TeV_binned_multi_eq_"+str(mbins)+".pdf")
        plt.clf()

