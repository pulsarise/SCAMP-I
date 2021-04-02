#! /usr/bin/env python3

import emcee
import argparse
import h5py
import matplotlib.pyplot as plt
import corner
import pandas
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = "14"

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input_log', type=str, help="Log information about the MCMC run.")
    parser.add_argument('-w','--writedir', type=str, help="Filepath to save the diagnostic plots")
    args = parser.parse_args()
    # Get relevant information from info file.
    psrinfo = pandas.read_csv(args.input_log, index_col=0)
    psrname = str(psrinfo.PSRJ.values[0])
    samplesreaddir = str(psrinfo.SAMPLESREADDIR.values[0])
    samplesfilename = str(psrinfo.SAMPLESFILENAME.values[0])
    P0 = float(psrinfo.PERIOD.values[0])
    nbins = int(psrinfo.NBIN.values[0])
    # Read in chains and make plots.
    writedir = args.writedir
    samplesbasename = samplesfilename.split('.h5')[0]
    filepath = '{}/{}'.format(samplesreaddir, samplesfilename)
    with h5py.File(filepath, "r") as f:
        keynames = list(f.keys())
    for key in keynames:
        reader = emcee.backends.HDFBackend(filepath, name=key)
        samples = reader.get_chain()
        flatsamples = reader.get_chain(flat=True)
        burnin = int(samples.shape[0]/2)
        flatburnedsamples = reader.get_chain(flat=True, discard=burnin)
        print("freq: {0}, chain shape: {1}".format(key, samples.shape))
        print("freq: {0}, flat chain shape: {1}".format(key, flatsamples.shape))
        # Chains.
        ndim = samples.shape[2]
        labels=[r'$\sigma$', r'$\mu$', r'$A$', r'$\tau$', 'DC']
        fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.tick_params(which="both", direction="in")
        axes[-1].set_xlabel("Step number");
        fig.subplots_adjust(left=0.1, bottom=0.1, right=0.96, top=0.98, hspace=0)
        plt.savefig('{}/{}_freq{}_burnin{}_chains.png'.format(writedir, samplesbasename, key, burnin))
        plt.close()
        # Corner plots.
        # Rescale the bins to seconds.
        flatburnedsamples[:,0] *= 1000*P0/nbins
        flatburnedsamples[:,1] *= 1000*P0/nbins
        flatburnedsamples[:,3] *= 1000*P0/nbins
        labelscorner = [r'$\sigma$', r'$\mu$', 'A', r'$\tau$', 'DC']
        fig = corner.corner(flatburnedsamples, labels=labelscorner, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize":16}, label_kwargs={"fontsize":16}, labelpad=1000000)
        plt.savefig('{}/{}_freq{}_burnin{}_corner.png'.format(writedir, samplesbasename, key, burnin))
        plt.close()
