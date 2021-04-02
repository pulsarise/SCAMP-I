#!/bin/python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas
import math

from ScatterMC.data_handling import read_data, read_headerfull
from ScatterMC.model_functions import GxETrain, GxETrain1D


def plot_tausigspectrum(freqGHz, tausec, tauserrsec, sigma, sigma_error, alphaMC, alphaerrMC, ampMC):
    # Convert to ms.
    taus = np.array(tausec)*1000
    tauserr = np.array(tauserrsec)*1000
    sigma = np.array(sigma)*1000
    sigma_error = np.array(sigma_error)*1000
    # Make the plot.
    plt.figure(figsize=(10,6))
    plt.errorbar(freqGHz,taus,yerr=tauserr,fmt='kx', markersize=10.0,capthick=1, capsize=4,linewidth=1.5, label=r'$\tau$')
    plt.errorbar(freqGHz,sigma,yerr=sigma_error,fmt='k.', markersize=10.0,capthick=1, capsize=4,linewidth=1.5, label=r'$\sigma$')
    alpha_ndp = int(-math.floor(math.log10(abs(alphaerr_MC))))
    if alpha_ndp < 0:
        alpha_ndp = 0
    plt.plot(freqGHz,(1000*ampMC*(freqGHz*1000)**(-1*alphaMC)),color='k',linewidth=1.5,label=r'$\alpha = {:.{}f} \pm {:.{}f}$'.format(alphaMC, alpha_ndp, alphaerrMC, alpha_ndp))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Frequency (GHz)',fontsize=22, labelpad=15.0)
    plt.ylabel(r'$\tau$ and $\sigma$ (ms)',fontsize=22)
    plt.legend(fontsize=22)
    plt.grid(True, which="both", ls="-")
    plt.gcf().subplots_adjust(bottom=0.15, top=0.97, right=0.96, left=0.11)


def create_profile_plots(meth, npch, taussec_highsnr, profilexaxis, data_highsnr, model_highsnr, pulseperiod, lmfitstdssec_highsnr, pulsar, freqGHz, proffilename):
    # Set up how many subplots per page.
    dimx, dimy = int(4), int(2)
    numsubplots = dimx*dimy
    numplots = int(np.ceil(npch/numsubplots))

    #"""Plot 1: Pulse profiles and fits"""
    if taussec_highsnr[0] > 1:
        taulabel =  taussec_highsnr
        taulabelerr = lmfitstdssec_highsnr
        taustring = 'sec'
    else:
        taulabel = taussec_highsnr*1000
        taulabelerr = lmfitstdssec_highsnr*1000
        taustring = 'ms'
    for k in range(numplots):
        j = int(numsubplots*k)
        figg = plt.figure(figsize=(int(8),int(12)))
        plots_remaining = int(npch - numsubplots*k)
        for i in range(np.min([int(numsubplots),int(plots_remaining)])):
            # figg.subplots_adjust(left = 0.08, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.15)
            figg.subplots_adjust(left = 0.08, right = 0.98, wspace=0.1,hspace=0.1,bottom=0.15)
            #plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            plt.subplot(dimx,dimy,i+1)
            if np.isnan(model_highsnr[j+i]).all():
                plt.plot(profilexaxis,data_highsnr[j+i]/max(data_highsnr[j+i]),alpha = 0.40, color='k')
                plt.ylim(ymin=-0.5+np.min(data_highsnr[j+i]/max(data_highsnr[j+i])), ymax=1.1*np.max(data_highsnr[j+i]/max(data_highsnr[j+i])))
            else:
                plt.plot(profilexaxis,data_highsnr[j+i]/max(model_highsnr[j+i]),alpha = 0.40, color='k')
                tau_ndp = int(-math.floor(math.log10(abs(taulabelerr[j+i]))))
                if tau_ndp < 0:
                    tau_ndp = 0
                plt.plot(profilexaxis,model_highsnr[j+i]/max(model_highsnr[j+i]),lw = 2.0, label=r'$\tau: {:.{}f} \pm{:.{}f}$ {}'.format(taulabel[j+i], tau_ndp, taulabelerr[j+i], tau_ndp, taustring), color='k') #alpha = 0.85, 
                plt.ylim(ymin=-0.5+np.min(data_highsnr[j+i]/max(model_highsnr[j+i])), ymax=1.1*np.max(data_highsnr[j+i]/max(model_highsnr[j+i])))
                plt.legend(fontsize=14,numpoints=1)
            #plt.ylim(ymin=-0.8, ymax=1.1*np.max(data_highsnr[j+i]/max(model_highsnr[j+i])))
            plt.xlim(min(profilexaxis), max(profilexaxis))
            x1, x2 = plt.xlim()
            y1, y2 = plt.ylim()
            plt.text(0.95*x2,0.9*y2,'{:.2f} GHz'.format(freqGHz[j+i]), fontsize=15, horizontalalignment="right", verticalalignment="top")
            if i == 6 or i == 7:            
                if pulseperiod < 1:
                    #plt.xlim(xmin=0,xmax=pulseperiod*1000)
                    plt.xlabel('Time (ms)',fontsize=15)
                elif pulseperiod >= 1:
                    #plt.xlim(xmin=0,xmax=pulseperiod)
                    plt.xlabel('Time (s)',fontsize=15)
                plt.xticks(fontsize=15)
            else:
                plt.xticks([])
            if i == 0 or i == 2 or i == 4 or i == 6:
                plt.yticks(fontsize=15)
                plt.ylabel('Normalized Intensity',fontsize=15)
            else:
                plt.yticks([])
            plt.tight_layout()
        plt.savefig("{}_{}.png".format(proffilename,k))



if __name__ == '__main__':
    # Define options to the script.
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--paramsfilepath', type=str, help="Information about data and best fit parameters.")
    parser.add_argument('-w', '--writedir', type=str, help="Where to save plots.")
    args = parser.parse_args()
    writedir = args.writedir
    # Load up the parameter file and specify which pulsar using the input index from the command line.
    paramsfilepath = args.paramsfilepath
    psrinfo = pandas.read_csv(paramsfilepath, index_col=0)
    # Get variables from the parameter file.
    tau = psrinfo.TAU.values
    tauerr = psrinfo.TAU_ERROR.values
    mean = psrinfo.MEAN.values
    mean_error = psrinfo.MEAN_ERROR.values
    sigma = psrinfo.SIGMA.values
    sigma_error = psrinfo.SIGMA_ERROR.values
    amplitude = psrinfo.AMPLITUDE.values
    dc = psrinfo.DC.values
    freqMHz = psrinfo.FREQ.values

    freqGHz = freqMHz/1000.
    alpha_MC = psrinfo.ALPHA_MCMC.values[0]
    alphaerr_MC = psrinfo.ALPHA_ERROR_MCMC.values[0]
    amp_MC = psrinfo.AMP_MCMC.values[0]
    method = psrinfo.METHOD.values[0]
    P0 = psrinfo.PERIOD.values[0]
    nchans = psrinfo.shape[0]
    nbins = psrinfo.NBIN.values[0]
    datafilename = psrinfo.DATAFILENAME.values[0]
    datareaddir = psrinfo.DATAREADDIR.values[0]
    samplesfilename = psrinfo.SAMPLESFILENAME.values[0]
    samplesbasename = samplesfilename.split('.h5')[0]
    # Read in the data and the frequencies.
    data_array = np.zeros((nchans, nbins))
    model_array = np.zeros((nchans, nbins))
    for f in range(nchans):
        data_f, _, _ = read_data(datafilename,datareaddir,f,nbins)
        if f == 0:
            peakbin = np.argmax(data_f)
            shift = int(int(nbins/2) -int(peakbin))
        data_array[f,:] = np.roll(data_f, shift)
        # Get the models.
        if np.isnan(tau[f]):
            model_f = np.full([nbins], np.nan)
        else:
            if method == 'iso':
                model_f = GxETrain(0, mean[f]*nbins/P0, sigma[f]*nbins/P0, amplitude[f], tau[f]*nbins/P0, dc[f], nbins)
            elif method == 'onedim':
                model_f = GxETrain1D(0, mean[f]*nbins/P0, sigma[f]*nbins/P0, amplitude[f], tau[f]*nbins/P0, dc[f], nbins)
        model_array[f,:] = model_f
    if P0 < 1:
        profilexaxis = np.linspace(0,P0,nbins)*1000
    elif P0 >= 1:
        profilexaxis = np.linspace(0,P0,nbins)
    psrname = psrinfo.PSRJ.values[0]

    # Do profile plots.
    proffilename = '{}/{}_profiles.png'.format(writedir, samplesbasename)
    create_profile_plots(method, nchans, tau, profilexaxis[:], data_array[:,:], model_array[:,:], P0, tauerr, psrname, freqGHz, proffilename)


    plot_tausigspectrum(freqGHz, tau, tauerr, sigma, sigma_error, alpha_MC, alphaerr_MC, amp_MC)
    plt.savefig('{}/{}_tausig.png'.format(writedir, samplesbasename))

    # Plot DM fit.
    DMdelta = psrinfo.DM_DELTA.values[0]
    DMdeltastd = psrinfo.DM_DELTA_ERROR.values[0]
    DMmodel = mean[~np.isnan(mean)][0] + DMdelta*4148.808*((1/freqMHz)**2-(1/freqMHz[~np.isnan(mean)][0])**2)
    plt.figure()
    plt.errorbar(mean,freqMHz, fmt='x', xerr=mean_error, color='k')
    plt.plot(DMmodel,freqMHz, '-', label=r'DM: ${:.3f} \pm {:.3f}$ $\rm{{cm}}^{{-3}} \rm{{pc}}$'.format(DMdelta,DMdeltastd), color='k')
    plt.xlabel(r'$\mu$ (sec)', fontsize =12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.title(psrname, fontsize=12)
    plt.ylabel(r'$\nu$ (MHz)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
    plt.legend(fontsize = 10, loc='best')
    plt.tight_layout()
    plt.savefig('{}/{}_deltaDM.png'.format(writedir, samplesbasename))        
