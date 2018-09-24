#!/usr/bin/env python
import time, os
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

plot_dir="/home/msilva/public_html/corsika21_vs_corsika2x/Sibyll"+args.VER+"/depth_zenith/"
if not os.path.exists(plot_dir):
       print "path doesn't exist. trying to make"
       os.makedirs(plot_dir)

start_time = time.asctime()
print 'Started:', start_time

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
print "loading mc"
file_corsika = 'corsika_'+args.VER+'_array.npy'
array_corsika = np.load(file_corsika)
cor_eventinfo = array_corsika[0]
cor_muon1info = array_corsika[1]
cor_muon2info = array_corsika[2]
cor_weights=np.array(cor_eventinfo[0]); cor_multi=np.array(cor_eventinfo[1]); cor_totE=np.array(cor_eventinfo[2]);
cor_errors=np.array(cor_eventinfo[3]);
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
mug_errors=np.array(cor_eventinfo[3]);
mug_energy1=np.array(mug_muon1info[0]); mug_zenith1=np.array(np.cos(mug_muon1info[1]));
mug_azimuth1=np.array(mug_muon1info[2]); mug_zpos1=np.array(mug_muon1info[5]); 
mug_energy2=np.array(mug_muon2info[0]); mug_zenith2=np.array(np.cos(mug_muon2info[1]));
mug_azimuth2=np.array(mug_muon2info[2]); mug_zpos2=np.array(mug_muon2info[5]); 

print "finished loading mc"

#energy bins
ebins=30; emin=args.EMIN; emax = args.EMAX;
logxbins = np.logspace(np.log10(emin),np.log10(emax), ebins)
ebins = np.arange(0,1.05,0.05)

#zenith bins, need irregularly sized bins
zbins = [0, 0.5, 0.75, 1.0]

#depth bins, need irregularly sized bins
dbins = [-550, -150, 200, 550]
#dbins = [-550, -300, -100, 100, 300, 550]

#mutliplicity bin
mbins = range(31)

#loop over cos(zenith) first
print "plotting multiplicity"
for zbin in range(len(zbins)-1):
    lowz = "%.2f" % (zbins[zbin]); highz = "%.2f" % (zbins[zbin+1]);
    print "for cos(zenith) bins: " + lowz + " to " + highz
    cor_zbins=np.where(np.logical_and(cor_zenith1>=zbins[zbin], cor_zenith1<=zbins[zbin+1]))
    cor_zbins = np.array(cor_zbins[0])    
    cor_weights_zbin=cor_weights[cor_zbins]
    cor_errors_zbin=cor_errors[cor_zbins]
    cor_multi_zbin=cor_multi[cor_zbins]
    cor_totE_zbin=cor_totE[cor_zbins]
    cor_energy1_zbin=cor_energy1[cor_zbins]
    cor_depth_zbin=cor_zpos1[cor_zbins]
    cor_energy2_zbin=cor_energy2[cor_zbins]
    print "loaded corsika"

    mug_zbins=np.where(np.logical_and(mug_zenith1>=zbins[zbin], mug_zenith1<=zbins[zbin+1]))
    mug_zbins = np.array(mug_zbins[0])
    mug_weights_zbin=mug_weights[mug_zbins]
    mug_errors_zbin=mug_errors[mug_zbins]
    mug_multi_zbin=mug_multi[mug_zbins]
    mug_totE_zbin=mug_totE[mug_zbins]
    mug_energy1_zbin=mug_energy1[mug_zbins]
    mug_depth_zbin=mug_zpos1[mug_zbins]
    mug_energy2_zbin=mug_energy2[mug_zbins]
    print "loaded muongun"

    #setup the plotting
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    plt.subplots_adjust(wspace=0, hspace=0.05)
    colors = ['b','g','r','c','m','y','k']

    #plot multiplicity first
    cor_2Dhist, xbins, ybins = np.histogram2d(cor_multi_zbin, cor_depth_zbin, bins=[mbins, dbins], weights=cor_weights_zbin)
    mug_2Dhist, xbins, ybins = np.histogram2d(mug_multi_zbin, mug_depth_zbin, bins=[mbins, dbins], weights=mug_weights_zbin)
    cor_errors_1, _, _       = np.histogram2d(cor_multi_zbin, cor_depth_zbin, bins=[mbins, dbins], weights=cor_errors_zbin)
    mug_errors_1, _, _       = np.histogram2d(mug_multi_zbin, mug_depth_zbin, bins=[mbins, dbins], weights=mug_errors_zbin)
    loop_bin = len(cor_2Dhist[0])
    for bin in range(loop_bin):
        print "now loop over depth bin: " + str(dbins[bin]) + ', ' + str(dbins[bin+1]) + " m"
        bincenters = 0.5*(xbins[1:]+xbins[:-1])
        ratio_multi, ratio_err, ratio_multi_norm  = ratio_errors(cor_2Dhist.T[bin], np.sqrt(cor_errors_1.T[bin]), 
                                                                 mug_2Dhist.T[bin], np.sqrt(mug_errors_1.T[bin]))
        ax1.errorbar(bincenters, cor_2Dhist.T[bin]/np.sum(cor_2Dhist.T[bin]),  yerr = np.sqrt(cor_errors_1.T[bin]), 
                     label=str(dbins[bin]) + ', ' + str(dbins[bin+1]) + " m",
                     fmt='o'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
        ax1.errorbar(bincenters, mug_2Dhist.T[bin]/np.sum(mug_2Dhist.T[bin]),  yerr = np.sqrt(mug_errors_1.T[bin]),
                     fmt='v'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
        ax2.errorbar(bincenters, ratio_multi_norm, yerr= ratio_err,
                     fmt='o'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
    #labels etc...
    print "creating labels"
    ax1.set_yscale('log')
    ax1.set_title("Binned for range: "+lowz + " < cos($\Theta_{1}$) < " + highz)
    ax1.legend(loc='upper right', ncol=2, prop={'size': 12})
    ax2.set_xlabel('Bundle Multiplicity')
    ax1.set_ylabel('Events per second [normalized]')
    ax2.set_ylabel('Corsika 2.x / 2.1')
    plt.savefig(plot_dir+"corsika_div_muongun_multiplicity_binned_coszenith_"+lowz+"_"+highz+"_binned_depth.pdf")
    print "made file: " + plot_dir+"corsika_div_muongun_multiplicity_binned_coszenith_"+lowz+"_"+highz+"_binned_depth.pdf"
    plt.close("All")

    #plot lead energy as function of depth now
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    plt.subplots_adjust(wspace=0, hspace=0.05)
    cor_2Dhist, xbins, ybins = np.histogram2d(cor_energy1_zbin, cor_depth_zbin, bins=[logxbins, dbins], weights=cor_weights_zbin)
    mug_2Dhist, xbins, ybins = np.histogram2d(mug_energy1_zbin, mug_depth_zbin, bins=[logxbins, dbins], weights=mug_weights_zbin)
    cor_errors_1, _, _       = np.histogram2d(cor_energy1_zbin, cor_depth_zbin, bins=[logxbins, dbins], weights=cor_errors_zbin)
    mug_errors_1, _, _       = np.histogram2d(mug_energy1_zbin, mug_depth_zbin, bins=[logxbins, dbins], weights=mug_errors_zbin)
    loop_bin = len(cor_2Dhist[0])
    for bin in range(loop_bin):
        print "now loop over depth bin: " + str(dbins[bin]) + ', ' + str(dbins[bin+1]) + " m"
        bincenters = 10**(0.5*(np.log10(xbins[1:])+np.log10(xbins[:-1])) )
        ratio_multi, ratio_err, ratio_multi_norm  = ratio_errors(cor_2Dhist.T[bin], np.sqrt(cor_errors_1.T[bin]),
                                                                 mug_2Dhist.T[bin], np.sqrt(mug_errors_1.T[bin]))
        ax1.errorbar(bincenters, cor_2Dhist.T[bin]/np.sum(cor_2Dhist.T[bin]),  yerr = np.sqrt(cor_errors_1.T[bin]),
                     label=str(dbins[bin]) + ', ' + str(dbins[bin+1]) + " m",
                     fmt='o'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
        ax1.errorbar(bincenters, mug_2Dhist.T[bin]/np.sum(mug_2Dhist.T[bin]),  yerr = np.sqrt(mug_errors_1.T[bin]),
                     fmt='v'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
        ax2.errorbar(bincenters, ratio_multi_norm, yerr= ratio_err,
                     fmt='o'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')

    print "creating labels"
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_title("Binned for range: "+lowz + " < cos($\Theta_{1}$) < " + highz)
    ax1.legend(loc='lower left', ncol=2, prop={'size': 12})
    ax2.set_xlabel('Leading Muon Energy [GeV]')
    ax1.set_ylabel('Events per second [normalized]')
    ax2.set_ylabel('Corsika 2.x / 2.1')
    plt.savefig(plot_dir+"corsika21_div_corsika_leadingenergy_binned_coszenith_"+lowz+"_"+highz+"_binned_depth.pdf")
    print "made file: " + plot_dir+"corsika_div_muongun_leadingenergy_binned_coszenith_"+lowz+"_"+highz+"_binned_depth.pdf"
    plt.close("All")

    #plot total energy as function of depth now
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    plt.subplots_adjust(wspace=0, hspace=0.05)
    cor_2Dhist, xbins, ybins = np.histogram2d(cor_totE_zbin, cor_depth_zbin, bins=[logxbins, dbins], weights=cor_weights_zbin)
    mug_2Dhist, xbins, ybins = np.histogram2d(mug_totE_zbin, mug_depth_zbin, bins=[logxbins, dbins], weights=mug_weights_zbin)
    cor_errors_1, _, _       = np.histogram2d(cor_totE_zbin, cor_depth_zbin, bins=[logxbins, dbins], weights=cor_errors_zbin)
    mug_errors_1, _, _       = np.histogram2d(mug_totE_zbin, mug_depth_zbin, bins=[logxbins, dbins], weights=mug_errors_zbin)
    loop_bin = len(cor_2Dhist[0])
    for bin in range(loop_bin):
        print "now loop over depth bin: " + str(dbins[bin]) + ', ' + str(dbins[bin+1]) + " m"
        bincenters = 10**(0.5*(np.log10(xbins[1:])+np.log10(xbins[:-1])) )
        ratio_multi, ratio_err, ratio_multi_norm  = ratio_errors(cor_2Dhist.T[bin], np.sqrt(cor_errors_1.T[bin]),
                                                                 mug_2Dhist.T[bin], np.sqrt(mug_errors_1.T[bin]))
        ax1.errorbar(bincenters, cor_2Dhist.T[bin]/np.sum(cor_2Dhist.T[bin]),  yerr = np.sqrt(cor_errors_1.T[bin]),
                     label=str(dbins[bin]) + ', ' + str(dbins[bin+1]) + " m",
                     fmt='o'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
        ax1.errorbar(bincenters, mug_2Dhist.T[bin]/np.sum(mug_2Dhist.T[bin]),  yerr = np.sqrt(mug_errors_1.T[bin]),
                     fmt='v'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
        ax2.errorbar(bincenters, ratio_multi_norm, yerr= ratio_err,
                     fmt='o'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
        
    print "creating labels"
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_title("Binned for range: "+lowz + " < cos($\Theta_{1}$) < " + highz)
    ax1.legend(loc='lower left', ncol=2, prop={'size': 12})
    ax2.set_xlabel('Sum of Muons Energy [GeV]')
    ax1.set_ylabel('Events per second [normalized]')
    ax2.set_ylabel('Corsika 2.x / 2.1')
    plt.savefig(plot_dir+"corsika_div_muongun_totalenergy_binned_coszenith_"+lowz+"_"+highz+"_binned_depth.pdf")
    print "made file: " + plot_dir+"corsika_div_muongun_totalenergy_binned_coszenith_"+lowz+"_"+highz+"_binned_depth.pdf"
    plt.close("All")

    #fractional energy
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})
    ax1, ax2 = prepare_ratio_plot(ax1, ax2)
    plt.subplots_adjust(wspace=0, hspace=0.05)
    cor_ratio_energy  = np.divide(cor_energy1_zbin, cor_totE_zbin, out=np.zeros_like(cor_totE_zbin), where=cor_totE_zbin!=0)
    mug_ratio_energy  = np.divide(mug_energy1_zbin, mug_totE_zbin, out=np.zeros_like(mug_totE_zbin), where=mug_totE_zbin!=0)
    cor_2Dhist, xbins, ybins = np.histogram2d(cor_ratio_energy, cor_depth_zbin, bins=[ebins, dbins], weights=cor_weights_zbin)
    mug_2Dhist, xbins, ybins = np.histogram2d(mug_ratio_energy, mug_depth_zbin, bins=[ebins, dbins], weights=mug_weights_zbin)
    cor_errors_1, _, _       = np.histogram2d(cor_ratio_energy, cor_depth_zbin, bins=[ebins, dbins], weights=cor_errors_zbin)
    mug_errors_1, _, _       = np.histogram2d(mug_ratio_energy, mug_depth_zbin, bins=[ebins, dbins], weights=mug_errors_zbin)
    loop_bin = len(cor_2Dhist[0])
    for bin in range(loop_bin):
        print "now loop over depth bin: " + str(dbins[bin]) + ', ' + str(dbins[bin+1]) + " m"
        bincenters = 0.5*(xbins[1:]+xbins[:-1])
        ratio_multi, ratio_err, ratio_multi_norm  = ratio_errors(cor_2Dhist.T[bin], np.sqrt(cor_errors_1.T[bin]),
                                                                 mug_2Dhist.T[bin], np.sqrt(mug_errors_1.T[bin]))
        ax1.errorbar(bincenters, cor_2Dhist.T[bin]/np.sum(cor_2Dhist.T[bin]),  yerr = np.sqrt(cor_errors_1.T[bin]),
                     label=str(dbins[bin]) + ', ' + str(dbins[bin+1]) + " m",
                     fmt='o'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
        ax1.errorbar(bincenters, mug_2Dhist.T[bin]/np.sum(mug_2Dhist.T[bin]),  yerr = np.sqrt(mug_errors_1.T[bin]),
                     fmt='v'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')
        ax2.errorbar(bincenters, ratio_multi_norm, yerr= ratio_err,
                     fmt='o'+colors[bin], ecolor=colors[bin], drawstyle = 'steps-mid-')


    print "creating labels"
    #ax1.set_yscale('log')
    ax1.set_title("Binned for range: "+lowz + " < cos($\Theta_{1}$) < " + highz)
    ax1.legend(loc='upper left', ncol=2, prop={'size': 12})
    ax2.set_xlabel('Leading Energy / Total Energy')
    ax1.set_ylabel('Events per second [normalized]')
    ax2.set_ylabel('Corsika 2.x / 2.1')
    plt.savefig(plot_dir+"corsika_div_muongun_frac_leadenergy_binned_coszenith_"+lowz+"_"+highz+"_binned_depth.pdf")
    print "made file: " + plot_dir+"corsika_div_muongun_frac_leadenergy_binned_coszenith_"+lowz+"_"+highz+"_binned_depth.pdf"
    plt.close("All")



end_time = time.asctime()
print 'Ends:', end_time
