#! /usr/bin/env python3

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pandas
import os


from SCAMP-I.data_handling import read_headerfull, read_data
from SCAMP-I.mcmc_functions import tau_fitter_mcmc

if __name__ == '__main__':
    # Define options to the script.
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--filename', type=str, help="Required. Path to config file.")
    parser.add_argument('-o','--inputlog', type=str, help="Filepath where log information about what was input to the run is saved.")
    parser.add_argument('-w', '--writedir', default='.', type=str, help='Location to which chains should be written.')
    parser.add_argument('-t','--runtime',type=int, default=20000, help="Number of steps for the MCMC, default 20,000.")
    parser.add_argument('-nw','--nwalkers',type=int, default=10, help="Number of walkers for the MCMC, default 10.")
    parser.add_argument('-m','--modelfn', default='iso', help="Model function applied. Choose between 'onedim' (extreme anisotropic) and 'iso' (isotropic). Default is 'iso'.")
    parser.add_argument('-n','--nthreads', default=1, help="How many threads to use in running the MCMC. Default is 1.")
    parser.add_argument('--showprogress', default=False, help="Show progress bar for MCMC.", action='store_true')
    args = parser.parse_args()
    # Allocate variable names to the parsed options.
    #readdir = args.readdir
    writedir = args.writedir
    configfilepath = args.filename
    config_df = pandas.read_csv(configfilepath)
    datafilename = config_df.DATAFILENAME.values[0]
    readdir = config_df.DATAREADDIR.values[0]
    if datafilename == None:
        raise RuntimeError('No filename specified. Please use -f option.')
    # Read ascii file header in full.
    pulsar, nch, nbins, nsub, lm_rms, tsub, obsdate = read_headerfull(datafilename, readdir)
    # Continue allocating variables.
    modelfn = args.modelfn
    nthreads = args.nthreads
    nwalkers = args.nwalkers
    runtime = args.runtime
    showprogress = args.showprogress
    # Add MCMC run information to config file and save it as a log.
    config_df["RUNTIME"] = runtime
    config_df["SAMPLESREADDIR"] = writedir
    config_df["METHOD"] = modelfn
    basefilename = os.path.basename(datafilename).split('.')[0]
    samplesfilename = "{}_runtime{}_model_{}.h5".format(basefilename, runtime, modelfn)
    config_df["SAMPLESFILENAME"] = samplesfilename
    config_df.to_csv(args.inputlog)
    # Run the MCMC chains.
    for freqnumber in range(nch):
        data, freqc, freqm = read_data(datafilename,readdir,freqnumber,nbins)
        # Roll the data such that lowest freq channel is centred.
        if freqnumber == 0:
            shift = int(len(data)/2 -int(np.argmax(data)))
        data = np.roll(data,shift)
        tau_fitter_mcmc(data, nbins, freqnumber, runtime, datafilename, writedir, modelfn, nwalkers, nthreads, showprogress)

