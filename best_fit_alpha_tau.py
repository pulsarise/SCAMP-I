#!/bin/python

import argparse
import emcee
import pandas
import numpy as np

if __name__ == '__main__':
    # Define options to the script.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--run_log', type=str, help="All relevant information required (about data and the MCMC run) to get best fit parameters out.")
    args = parser.parse_args()
    psrinfo = pandas.read_csv(args.run_log, index_col=0)
    samplesreaddir = str(psrinfo.ALPHASAMPLESREADDIR.values[0])
    samplesfilename = str(psrinfo.ALPHASAMPLESFILENAME.values[0])
    burnfrac = float(psrinfo.ALPHABURNFRAC.values[0])
    # Read samples out of the h5 file.
    reader = emcee.backends.HDFBackend('{}/{}'.format(samplesreaddir, samplesfilename))
    samples = reader.get_chain()
    burnin = int(samples.shape[0]*burnfrac)
    ndim = samples.shape[2]
    flat_burned_samples = reader.get_chain(discard=burnin, flat=True)
    # Get best fit values.
    for i in range(ndim):
        mcmc = np.percentile(flat_burned_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        if i == 0:
            amp_MCMC = mcmc[1]
            amperr_MCMC = (q[0] + q[1])/2
        elif i == 1:
            alpha_MCMC = mcmc[1]
            alphaerr_MCMC = (q[0] + q[1])/2
    # Assign to results dataframe.
    psrinfo["ALPHA_MCMC"] = alpha_MCMC
    psrinfo["ALPHA_ERROR_MCMC"] = alphaerr_MCMC
    psrinfo["AMP_MCMC"] = amp_MCMC
    psrinfo["AMP_ERROR_MCMC"] = amperr_MCMC
    # Now infer best fit for tau at 1GHz.
    freq1GHz = 1000.
    tau_1GHz_list = []
    for ind in range(len(flat_burned_samples)):
        sample = flat_burned_samples[ind]
        tau_1GHz_list.append(sample[0]*freq1GHz**(-1*sample[1]))
    tau_1GHz_MCMC_values = np.percentile(tau_1GHz_list, [16, 50, 84])
    valdiff = np.diff(tau_1GHz_MCMC_values)
    tau1GHz = tau_1GHz_MCMC_values[1]
    tauerr1GHz = (valdiff[0] + valdiff[1])/2
    psrinfo["TAU_1GHz"] = tau1GHz
    psrinfo["TAU_ERROR_1GHz"] = tauerr1GHz
    # Save.
    psrinfo.to_csv(args.run_log)
