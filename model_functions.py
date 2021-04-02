#! /usr/bin/env python3

import numpy as np


def find_rms(data,nbins):
    windowsize = 32 
    windows = int(nbins/windowsize)
    rms_loc = np.zeros(windows)
    for i in range(windows):
        start = i*windowsize
        end = start + windowsize
        rms_loc[i] = np.std(data[start:end])
    return np.min(rms_loc)

def makeprofile(nbins = 2**9, ncomps = 1, amps = 1, means = 100, sigmas = 10):
    if ncomps == 1:
        npamps = np.array([amps])
        npmeans = np.array([means])
        npsigmas = np.array([sigmas])
    else:
        npamps = np.array(amps)
        npmeans = np.array(means)
        npsigmas = np.array(sigmas)

    profile = np.zeros(nbins)
    x = np.linspace(1,nbins,nbins)

    for i in range(ncomps):
        profile = profile + \
        npamps[i]*np.exp(-pow(x-npmeans[i],2)/(2*pow(npsigmas[i],2)))
    return x, profile

def pulsetrain(npulses = 10, bins = np.linspace(1,512,512), profile = np.zeros(512)):
    nbins = np.max(bins)
    train = np.zeros(npulses*int(nbins))
    nbins = int(nbins)
    for i in range(npulses):
        startbin = int(i*nbins)
        endbin = int(startbin + nbins)
        train[startbin:endbin] = profile
    return train

def pulsetrain_bins(npulses, numberofbins, profile):
    binsrange = np.linspace(1,numberofbins,numberofbins)
    nbins = np.max(binsrange)
    nbins = int(nbins)
    train = np.zeros(npulses*int(nbins))

    for i in range(npulses):
        startbin = int(i*nbins)
        endbin = int(startbin + nbins)
        train[startbin:endbin] = profile
    return train

def psrscatter(brfunc, profile):
    scattered = np.convolve(profile,brfunc)
    profint = np.sum(profile)
    scint = np.sum(scattered)
    scatterednorm = scattered / scint * profint
    bins = profile.shape[0]
    out = scatterednorm[0:bins]
    return out

def psrscatter_noconserve(brfunc, profile):
    scattered = np.convolve(profile,brfunc)
    bins = profile.shape[0]
    out = scattered[0:bins]
    return out
    
def step(x):
    return 1 * (x >= 0)

def broadfunc(x,tau):
    # 1. Isotropic scattering
    tau = float(tau)
    broadfunc = (1/tau)*np.exp(-x/tau)*step(x)
    return broadfunc   

def broadfunc1D(x,tau):
    # 2. Extremely anisotropic scattering
    broadfunc1 = (1.0/np.sqrt(x*tau*np.pi))*np.exp(-x/tau)
    return broadfunc1

def extractpulse(train, pulsesfromend, binsperpulse):
    if pulsesfromend == 0:
        start = 0
        end = binsperpulse
        zerobpulse = train[start:end]-np.min(train[start:end])
        rectangle = np.min(train[start:end])*binsperpulse
        flux = np.sum(train[start:end]) - rectangle
        return train[start:end], zerobpulse, rectangle, flux

    else:
        start = -pulsesfromend*binsperpulse
        end = start + binsperpulse
        zerobpulse = train[start:end]-np.min(train[start:end])
        rectangle = np.min(train[start:end])*binsperpulse
        flux = np.sum(train[start:end]) - rectangle
        return train[start:end], zerobpulse, rectangle, flux

def GxETrain(x,mu,sigma, A, tau, dc, nbins):
    # This model convolves a pulsetrain with a broadening function
    # It extracts one of the last convolved profiles, subtracts the climbed baseline and then adds noise to it
    mu, sigma, A, tau = float(mu),float(sigma), float(A), float(tau)
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,nbins,nbins)
    npulses = 3
    scat = psrscatter_noconserve(broadfunc(binstau,tau),pulsetrain_bins(npulses, nbins, profile))
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc
    
def GxETrain1D(x,mu, sigma, A, tau1, dc, nbins):
    mu, sigma, A, tau1 = float(mu),float(sigma), float(A), float(tau1)
    bins, profile = makeprofile(nbins = nbins, ncomps = 1, amps = A, means = mu, sigmas = sigma)
    binstau = np.linspace(1,nbins,nbins)
    npulses = 3
    scat = psrscatter(broadfunc1D(binstau,tau1),pulsetrain(npulses, bins, profile))
    climb, observed_nonoise, rec,flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc
    
