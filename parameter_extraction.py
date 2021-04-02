#! /usr/bin/env python3

import numpy as np
import emcee
import argparse
import pandas

from ScatterMC.dmcorrcalc import get_old_new_DM

def get_bestfit_params(flat_burned_samples, P0, nbins):
    # Get the best fit params out from the chains.
    pc = np.percentile(flat_burned_samples[:, 0], [16, 50, 84])
    q = np.diff(pc)
    bestsig = pc[1]
    bestsig_std_n = q[0]
    bestsig_std_p = q[1]
    pc = np.percentile(flat_burned_samples[:, 1], [16, 50, 84])
    q = np.diff(pc)
    bestmu = pc[1]
    bestmu_std_n = q[0]
    bestmu_std_p = q[1]
    pc = np.percentile(flat_burned_samples[:, 2], [16, 50, 84])
    q = np.diff(pc)
    bestA = pc[1]
    bestA_std_n = q[0]
    bestA_std_p = q[1]
    pc = np.percentile(flat_burned_samples[:, 3], [16, 50, 84])
    q = np.diff(pc)
    besttau = pc[1]
    besttau_std_n = q[0]
    besttau_std_p = q[1]
    pc = np.percentile(flat_burned_samples[:, 4], [16, 50, 84])
    q = np.diff(pc)
    bestdc = pc[1]
    bestdc_std_n = q[0]
    bestdc_std_p = q[1]
    # Make arrays.
    best_params = np.array([bestsig, bestmu, bestA, besttau, bestdc])
    positive_error = np.array([bestsig_std_p, bestmu_std_p, bestA_std_p, besttau_std_p, bestdc_std_p])
    negative_error = np.array([bestsig_std_n, bestmu_std_n, bestA_std_n, besttau_std_n, bestdc_std_n])
    # Convert from bins to seconds (only mu, sigma and tau ie the time domain variables).
    bin_to_sec = P0/nbins
    best_params_sec = np.array([bestsig*bin_to_sec, bestmu*bin_to_sec, bestA, besttau*bin_to_sec, bestdc])
    positive_error_sec = np.array([bestsig_std_p*bin_to_sec, bestmu_std_p*bin_to_sec, bestA_std_p, besttau_std_p*bin_to_sec, bestdc_std_p])
    negative_error_sec = np.array([bestsig_std_n*bin_to_sec, bestmu_std_n*bin_to_sec, bestA_std_n, besttau_std_n*bin_to_sec, bestdc_std_n])
    return best_params_sec, positive_error_sec, negative_error_sec, best_params, positive_error, negative_error

def get_samples(readdir, samplesfilename, f, burnin=10000):
    # Read samples out of the h5 file.
    reader = emcee.backends.HDFBackend('{}/{}'.format(readdir, samplesfilename), name='{}'.format(f))
    samples = reader.get_chain()
    flat_burned_samples = reader.get_chain(discard=burnin, flat=True)
    return samples, flat_burned_samples

def get_params_from_samples(samplesfilename, samplesreaddir, passed_list, burnlist, nbins, P0, freqMHz):
    # Set things up for getting the parameters.
    ndim = 5
    nchans = len(freqMHz)
    param_array = np.zeros((nchans, ndim))
    sigma = []
    mean = []
    amplitude = []
    dc = []
    sigma_error = []
    mean_error = []
    amplitude_error = []
    dc_error = []
    tau = []
    tau_error = []
    for f in range(nchans):
        if passed_list[f] == False:
            sigma.append(np.nan)
            mean.append(np.nan)
            amplitude.append(np.nan)
            tau.append(np.nan)
            dc.append(np.nan)
            sigma_error.append(np.nan)
            mean_error.append(np.nan)
            amplitude_error.append(np.nan)
            tau_error.append(np.nan)
            dc_error.append(np.nan)
        else:
            samples, fbs = get_samples(samplesreaddir, samplesfilename, f, burnin=burnlist[f])
            # Parameters.
            params_sec, paramsplus_sec, paramsminus_sec, params_bins, paramsplus_bins, paramsminus_bins = get_bestfit_params(fbs, P0, nbins)
            bestsig, bestmu, bestA, besttau, bestdc = params_bins
            # Get all the parameters.
            sigma.append(params_sec[0])
            mean.append(params_sec[1])
            amplitude.append(params_sec[2])
            tau.append(params_sec[3])
            dc.append(params_sec[4])
            sigma_error.append((paramsplus_sec[0]+paramsminus_sec[0])/2)
            mean_error.append((paramsplus_sec[1]+paramsminus_sec[1])/2)
            amplitude_error.append((paramsplus_sec[2]+paramsminus_sec[2])/2)
            tau_error.append((paramsplus_sec[3]+paramsminus_sec[3])/2)
            dc_error.append((paramsplus_sec[4]+paramsminus_sec[4])/2)

    return mean, mean_error, sigma, sigma_error, amplitude, amplitude_error, dc, dc_error, tau, tau_error


if __name__ == '__main__':
    # Define options to the script.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--run_log', type=str, help="All relevant information required (about data and the MCMC run) to get best fit parameters out.")
    parser.add_argument('-o', '--outputfilename', type=str, help='Name of file to which parameters are written.')
    parser.add_argument('-w', '--writedir', default='.', type=str, help='Location where output should be written.')
    args = parser.parse_args()

    # Read in key info required to extract parameters.
    psrinfo = pandas.read_csv(args.run_log, index_col=0)
    psrname = str(psrinfo.PSRJ.values[0])
    samplesreaddir = str(psrinfo.SAMPLESREADDIR.values[0])
    samplesfilename = str(psrinfo.SAMPLESFILENAME.values[0])
    P0 = float(psrinfo.PERIOD.values[0])
    nbins = int(psrinfo.NBIN.values[0])
    freq = np.array(psrinfo.FREQ.values, dtype=float)
    passed_list = np.array(psrinfo.PASSED.values, dtype=bool)
    burnlist = np.array(psrinfo.BURNFRAC.values, dtype=int)*np.array(psrinfo.RUNTIME.values, dtype=int)

    # Obtain parameters.
    mean, mean_error, sigma, sigma_error, amplitude, amplitude_error, dc, dc_error, tau, tau_error = get_params_from_samples(samplesfilename, samplesreaddir, passed_list, burnlist, nbins, P0, freq)

    # Apply to file.
    psrinfo["MEAN"] = mean
    psrinfo["SIGMA"] = sigma
    psrinfo["AMPLITUDE"] = amplitude
    psrinfo["DC"] = dc
    psrinfo["TAU"] = tau
    psrinfo["MEAN_ERROR"] = mean_error
    psrinfo["SIGMA_ERROR"] = sigma_error
    psrinfo["AMPLITUDE_ERROR"] = amplitude_error
    psrinfo["DC_ERROR"] = dc_error
    psrinfo["TAU_ERROR"] = tau_error

    # Get new DM by fitting means across frequency.
    DMdelta, DMdeltastd, DM_CCval, DM_CCvalstd, DMnew = get_old_new_DM(psrinfo)
    psrinfo['DM_DELTA'] = DMdelta
    psrinfo['DM_DELTA_ERROR'] = DMdeltastd
    psrinfo['DISP_DELTA'] = DM_CCval
    psrinfo['DISP_DELTA_ERROR'] = DM_CCvalstd
    psrinfo['DM_NEW'] = DMnew
    
    psrinfo.to_csv("{}/{}".format(args.writedir,args.outputfilename))
