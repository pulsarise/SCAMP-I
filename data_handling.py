#! /usr/bin/env python3

import numpy as np
    
### Read ascii files

def read_headerfull(filepath, readdir):
    f = open('{}/{}'.format(readdir,filepath))
    lines = f.readlines()
    header0 = lines[0]
    header1 = lines[1]
    h0_lines = header0.split()
    if h0_lines[0] == '#':
        h0_lines = h0_lines[1:len(h0_lines)]
    else:
        h0_lines = h0_lines    
    file_name = h0_lines[1]
    pulsar_name = h0_lines[3]
    nsub = int(h0_lines[5])
    nch = int(h0_lines[7])
    npol = int(h0_lines[9])
    nbins = int(h0_lines[11])
    rms = float(h0_lines[13])
    h1_lines = header1.split()
    tsub = float(h1_lines[4])
    # Add MJD as a parameter.
    MJD = float(h1_lines[2])
    return pulsar_name, nch, nbins, nsub, rms, tsub, MJD



def read_data(filepath, readdir, profilenumber, nbins):
    d = open('{}/{}'.format(readdir, filepath))
    lines = d.readlines()
    
    profile_start = 2+profilenumber*(nbins+1)
    profile_end = profile_start + nbins
    
    lines_block = lines[profile_start:profile_end]

    if lines[profile_start-1].split()[0] == '#':
        freqc = float(lines[profile_start-1].split()[6])
        bw = float(lines[profile_start-1].split()[8])
        freqm = 10**((np.log10(freqc+ bw/2.)+ np.log10(freqc - bw/2.))/2)
    else:
        freqc = float(lines[profile_start-1].split()[5])
        bw = float(lines[profile_start-1].split()[7])
        freqm = 10**((np.log10(freqc+ bw/2.)+ np.log10(freqc - bw/2.))/2)
    datalist = []
    for i in range(nbins):
        data= float(lines_block[i].split()[3])
        datalist.append(data)

    return np.array(datalist), freqc, freqm

    
