#! /usr/bin/env python3

import emcee
import argparse
import h5py
import matplotlib.pyplot as plt
import corner
import sys
import pandas

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input_log', type=str, help="Log information about the alpha MCMC run.")
    parser.add_argument('-w','--writedir', type=str, help="Filepath to save the diagnostic plots.")
    args = parser.parse_args()
    # Get relevant information from info file.
    psrinfo = pandas.read_csv(args.input_log, index_col=0)
    psrname = str(psrinfo.PSRJ.values[0])
    samplesreaddir = str(psrinfo.ALPHASAMPLESREADDIR.values[0])
    samplesfilename = str(psrinfo.ALPHASAMPLESFILENAME.values[0])
    P0 = float(psrinfo.PERIOD.values[0])
    nbins = int(psrinfo.NBIN.values[0])
    # Read in chains and make plots.
    writedir = args.writedir
    samplesbasename = samplesfilename.split('.h5')[0]
    filepath = '{}/{}'.format(samplesreaddir, samplesfilename)
    reader = emcee.backends.HDFBackend(filepath)
    samples = reader.get_chain()
    flatsamples = reader.get_chain(flat=True)
    burnin = int(samples.shape[0]/2)
    flatburnedsamples = reader.get_chain(flat=True, discard=burnin)
    # Chains.
    ndim = samples.shape[2]
    labels=['A', 'alpha']
    fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
    axes[-1].set_xlabel("step number");
    plt.savefig('{}/{}_chains.png'.format(writedir, samplesbasename))
    plt.close()
    # Corner plots.
    labelscorner = ['A', r'$\alpha$']
    fig = corner.corner(flatburnedsamples, labels=labelscorner, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize":16}, label_kwargs={"fontsize":16}, labelpad=1000000)
    plt.savefig('{}/{}_burnin{}_corner.png'.format(writedir, samplesbasename, burnin))
    plt.close()
