#!/bin/python

import numpy as np
from lmfit.models import LinearModel

def get_old_new_DM(df_info):
    # Get variables from the parameter file.
    mean = df_info.MEAN.values
    mean_error = df_info.MEAN_ERROR.values
    freqMHz = df_info.FREQ.values
    psrname = df_info.PSRJ.values[0]
    DMorig = df_info.DM_ORIG
    # Fit a DM model to delta mu.
    # Get rid of the nans.
    freqMHz = freqMHz[~np.isnan(mean)]
    mean_error = mean_error[~np.isnan(mean)]
    mean = mean[~np.isnan(mean)]
    delnuarray = np.array([(1/freqMHz[i]**2-1/freqMHz[0]**2) for i in range(len(mean))]) ##in MHz
    delmuarray = np.array([(mean[i] - mean[0]) for i in range(len(mean))]) ##in seconds
    delmu_stdarray = np.array([(mean_error[i] - mean_error[0]) for i in range(len(mean))]) ##in seconds
    if len(delmuarray) < 4:
        print('{}: only {} freqs, not enough to measure DM'.format(psrname, len(delmuarray)))
        return np.nan, np.nan, np.nan, np.nan, np.nan
    else:
        linmod = LinearModel()
        DM_linpars = linmod.guess(delmuarray, x=delnuarray)
        DM_linout  = linmod.fit(delmuarray, DM_linpars, x=delnuarray)
        DM_CCval = DM_linout.best_values['slope']
        DM_CCvalstd = DM_linout.params['slope'].stderr
        DMmodelfit = DM_linout.best_fit ##model gives deltime in seconds (used to shift data)
        DMconstant = 4148.808
        #uncertainty in the constant is 0.003 - only affects the Delta DM value in the 9th decimal
        DMdelta = (DM_CCval/DMconstant)
        DMdeltastd = (DM_CCvalstd/DMconstant)
        # Return the original DM, deltaDM, delta-dispersion-constant (ie the best fit not divided by 4148.808), and the new DM, along with their errors.
        DMnew = DMorig + DMdelta
        return DMdelta, DMdeltastd, DM_CCval, DM_CCvalstd, DMnew

