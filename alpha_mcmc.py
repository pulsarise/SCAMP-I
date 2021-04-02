#!/bin/python

import pandas
import numpy as np
import emcee
import argparse
from scipy.optimize import minimize

def powerlaw(x, A, alpha):
    return A*x**(-1*alpha)

def log_prior(theta):
    A, alpha = theta
    if 0. < A < np.inf and 0. < alpha < np.inf:
        return 0.0
    return -np.inf

def log_likelihood(theta, x, y, yerr):
    A, alpha = theta
    model = powerlaw(x, A, alpha)
    return -0.5*np.sum((y-model)**2/yerr**2)

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--parameter_file', type=str, help="All relevant information required to get parameters out of data.")
    parser.add_argument('-o', '--outputfilename', type=str, help='Name of file to which new parameters are written.')
    parser.add_argument('--sampleswritedir', type=str, help="Where to save the chains.")
    parser.add_argument('--resultswritedir', type=str, help="Where to save the results for best fit parameters.")
    parser.add_argument('-t', '--runtime', type=int, help="How many steps for the MCMC run, default 10,000.", default=10000)
    parser.add_argument('-nw','--nwalkers',type=int, default=32, help="Number of walkers for the MCMC, default 32.")
    parser.add_argument('--showprogress', default=False, help="Show progress bar for MCMC.", action='store_true')
    args = parser.parse_args()
    sampleswritedir = args.sampleswritedir
    resultswritedir = args.resultswritedir
    runtime = args.runtime
    outputfilename = args.outputfilename
    paramsfilepath = args.parameter_file
    df = pandas.read_csv(paramsfilepath, index_col=0)
    # Get the useful info from it necessary for calculating alpha.
    psrname = df.PSRJ.values[0]
    tau = df.TAU.values
    tau_error = df.TAU_ERROR.values
    freq = df.FREQ.values
    sigma = df.SIGMA.values
    # Remove nans for fitting purposes.
    tau_error = tau_error[~np.isnan(tau)]
    freq = freq[~np.isnan(tau)]
    tau = tau[~np.isnan(tau)]
    # Calculate alpha from taus and freqs by linalg.
    x = np.log10(freq)
    y = np.log10(tau)
    yerr = tau_error/tau
    A = np.vander(x, 2)
    C = np.diag(yerr * yerr)
    ATA = np.dot(A.T, A / (yerr ** 2)[:, None])
    cov = np.linalg.inv(ATA)
    w = np.linalg.solve(ATA, np.dot(A.T, y / yerr ** 2))
    alpha_linalg = -1*w[0]
    alphaerr_linalg = np.sqrt(cov[0,0])
    amp_linalg = 10**w[1]
    # Now use this to do maximum likelihood.
    nll = lambda *args: -log_likelihood(*args)
    initial = np.array([amp_linalg, alpha_linalg]) + 1 * np.random.randn(2)
    soln = minimize(nll, initial, args=(freq, tau, tau_error))
    amp_maxlik, alpha_maxlik  = soln.x
    print("max likelihood: ", alpha_maxlik)
    # Use this as starting point for MCMC.
    pos = soln.x + 1e-4 * np.random.randn(args.nwalkers, 2)
    nwalkers, ndim = pos.shape
    showprogress = args.showprogress
    # Initialize the walkers
    # Set up the backend
    # Don't forget to clear it before any running in case the file already exists
    filename = "{}_alphachains_runtime{}.h5".format(psrname, runtime)
    print('psrname: {}, nwalkers: {}, runtime: {}'.format(psrname, nwalkers, runtime))
    backend = emcee.backends.HDFBackend("{}/{}".format(sampleswritedir, filename))
    backend.reset(nwalkers, ndim)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(freq, tau, tau_error), backend=backend)
    max_n = runtime
    # We'll track how the average autocorrelation time estimate changes
    index = 0
    autocorr = np.empty(max_n)
    # This will be useful to testing convergence
    old_autocorrtime = np.inf
    # Now we'll sample for up to max_n steps
    for sample in sampler.sample(pos, iterations=max_n, progress=showprogress):
        # Only check convergence every 100 steps
        if sampler.iteration % 100:
            continue
        # Compute the autocorrelation time so far
        # Using tol=0 means that we'll always get an estimate even
        # if it isn't trustworthy
        autocorrtime = sampler.get_autocorr_time(tol=0)
        autocorr[index] = np.mean(autocorrtime)
        index += 1
        # Check convergence
        converged = np.all(autocorrtime * 100 < sampler.iteration)
        converged &= np.all(np.abs(old_autocorrtime - autocorrtime) / autocorrtime < 0.01)
        if converged:
            break
        old_autocorrtime = autocorrtime
    df["ALPHASAMPLESREADDIR"] = sampleswritedir
    df["ALPHASAMPLESFILENAME"] = filename
    # Save.
    df.to_csv("{}/{}".format(resultswritedir,outputfilename))


