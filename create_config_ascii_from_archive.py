#! usr/bin/env python

import os
import psrchive
import argparse
import subprocess

def archive_to_ascii(archive, outpath, nchan=None, verbose=False):
    if nchan is not None:
        fn, ext = os.path.splitext(archive)
        new_ext = "."+str(nchan)+"ch"
        pam_comand = ['pam', '--setnchn', str(nchan), '-TDp','-e', ext+new_ext, '-u', outpath, archive]
        pam_Popen = subprocess.Popen(pam_comand, shell=False, cwd='.', stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (stdoutdata, stderrdata) = pam_Popen.communicate()
        
        inname = str(archive)+new_ext
    else:
        inname = archive
        
    if verbose == True:
        print("Subprocess output:\n")
        print(stdoutdata+'\n')
    outname = inname+".ascii"
    ascii_file = open(os.path.join(outpath,outname), 'a')
    pdv_comand = ['pdv', '-CKA', inname]
    pdv_Popen = subprocess.Popen(pdv_comand, shell=False, cwd=outpath, stdout=ascii_file)
    (stdoutdata, stderrdata) = pdv_Popen.communicate()
    if verbose == True:
        print "Ascii saved as %s\n" %outname
    
    reduced_archive = inname
    saved_ascii = outname
    
    return reduced_archive, saved_ascii

def get_archive_info(archive, verbose=False):
   """Query archive attributes.
   Input:
       archive: archive filename to be loaded by psrchive
   Output:
       Chosen attrtibutes. Print attributes of the archive if verbose is True.
   """
   arch = psrchive.Archive_load(archive)
   filename = arch.get_filename()
   nbin = arch.get_nbin()
   nchan = arch.get_nchan()
   npol = arch.get_npol()
   nsubint = arch.get_nsubint()
   obs_type = arch.get_type()
   telescope_name = arch.get_telescope()
   source_name = arch.get_source()
   ra = arch.get_coordinates().ra()
   dec = arch.get_coordinates().dec()
   centre_frequency = arch.get_centre_frequency()
   bandwidth = arch.get_bandwidth()
   DM = arch.get_dispersion_measure()
   RM = arch.get_rotation_measure()
   is_dedispersed = arch.get_dedispersed()
   is_faraday_rotated = arch.get_faraday_corrected()
   is_pol_calib = arch.get_poln_calibrated()
   data_units = arch.get_scale()
   data_state = arch.get_state()
   obs_duration = arch.integration_length()
   obs_start = arch.start_time().fracday() + arch.start_time().intday()
   obs_end = arch.end_time().fracday() + arch.end_time().intday()
   receiver_name = arch.get_receiver_name()
   receptor_basis = arch.get_basis()
   backend_name = arch.get_backend_name()
   backend_delay = arch.get_backend_delay()
   # obtain pulse period:
   ephem = arch.get_ephemeris()
   rot_freq = ephem.get_value('F0')
   pulseperiod = 1.0/float(rot_freq)
   # obtain frequencies
   freqs = arch.get_frequencies()

   # low_freq = archive.get_centre_frequency() - archive.get_bandwidth() / 2.0
   # high_freq = archive.get_centre_frequency() + archive.get_bandwidth() / 2.0
   if verbose == True:
       print('Selected metadata for %s' %archive)
       print('------------------------------------')
       print('name             Source name                                %s' % source_name)
       print('file             Name of the file                           %s' % filename)
       print('nbin             Number of pulse phase bins                 %s' % nbin)
       print('nchan            Number of frequency channels               %s' % nchan)
       print('npol             Number of polarizations                    %s' % npol)
       print('nsubint          Number of sub-integrations                 %s' % nsubint)
       print('type             Observation type                           %s' % obs_type)
       print('site             Telescope name                             %s' % telescope_name)     
       print('freq             Centre frequency (MHz)                     %s' % centre_frequency)
       print('bw               Bandwidth (MHz)                            %s' % bandwidth)
       print('P                Pulse period (s)                           %.10f' % pulseperiod)
       print('dm               Dispersion measure (pc/cm^3)               %s' % DM)
       #print('rm               Rotation measure (rad/m^2)                 %s' % RM)
       #print('dmc              Dispersion corrected                       %s' % is_dedispersed)
       #print('rmc              Faraday Rotation corrected                 %s' % is_faraday_rotated)
       #print('polc             Polarization calibrated                    %s' % is_pol_calib)
       #print('scale            Data units                                 %s' % data_units)
       #print('state            Data state                                 %s' % data_state)
       print('length           Observation duration (s)                   %s' % obs_duration)
       print('start            Observation start (MJD)                    %.10f' % obs_start)
       print('end              Observation end (MJD)                      %.10f' % obs_end)
       print('rcvr:name        Receiver name                              %s' % receiver_name)
       #print('rcvr:basis       Basis of receptors                         %s' % receptor_basis)
       #print('be:name          Name of the backend instrument             %s' % backend_name)
       #print('be:delay         Backend propn delay from digi. input.      %s\n' % backend_delay)
       print('\n')
    
   dict={'PSR':source_name, "PERIOD": pulseperiod,"NBIN":nbin,"DM_ORIG":DM,"NCHAN":nchan,"FREQ":freqs}
    
   return dict

def write_config(datadir, archive, outpath, outfile=None, verbose=False):

    if verbose == True:
        verbtag = True
    else:
        verbtag = False
    
    os.chdir(datadir)
    
    header = "PSRJ,DATAREADDIR,DATAFILENAME,PERIOD,NBIN,DM_ORIG,CHAN,FREQ"
    nr_head = len(header.split(","))
    write_ls = ['%s'] * nr_head
    
    d = get_archive_info(archive, verbose=verbtag)
    nchan = d['NCHAN']
    if outfile is not None:
        outfile = outfile
    else:
        outfile = '%s_config_%sch.csv' %(d['PSR'],nchan)

    output = os.path.join(outpath,outfile)
    with open(output,'a') as f:
        f.write('%s\n' %header)
        for i in range(nchan):
            f.write(",".join(write_ls) %(d['PSR'],datadir,archive+'.ascii',d['PERIOD'],d['NBIN'],d['DM_ORIG'],i, d['FREQ'][i]))
            f.write("\n")
        f.close()
    
    if verbtag == True:
        print("Config file saved as %s\n" %outfile)



if __name__ == '__main__':
    # Define options to the script.
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--filename', type=str, help="Required. Provide the filename of the archive data file", required=True)
    parser.add_argument('-d','--dir', type=str, help="Required. Provide the pathway to where the archive data files are", required=True)
    parser.add_argument('-nch','--nchan', type=int, default=4, help="Number of channels across profile to fit scattering values to (default: 4).", required=False)
    parser.add_argument('-w', '--writedir', type=str ,help="Write directory to store channelised archive, corresponding ascii output and config file (default: cwd).", required=False) 
    parser.add_argument('-v', dest='verbose', action="store_true", help='verbose output', required=False)

    args = parser.parse_args()
    verbtag = args.verbose

    if args.writedir is not None:
        writedir = args.writedir
    else:
        writedir = os.getcwd()

    print("Writing to directory %s\n" %writedir ) 

    os.chdir(args.dir)
    
    reduced_ar, ascii_file = archive_to_ascii(args.filename,outpath=writedir,nchan=args.nchan,verbose=verbtag)
    write_config(writedir,reduced_ar,outpath=writedir,verbose=verbtag)







