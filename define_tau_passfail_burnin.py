#! /usr/bin/env python3

import numpy as np
import argparse
import pandas

if __name__ == '__main__':
    # Define options to the script.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input_log', type=str, help="Log information about the MCMC run.")
    parser.add_argument('-o', '--output_log', type=str, help="Output log describing run success/failure and burn-in to be applied")
    parser.add_argument('-w', '--writedir', default='.', type=str, help='Location to whihc output log should be written.')
    args = parser.parse_args()
    # Read in the input log.
    psrinfo = pandas.read_csv(args.input_log, index_col=0)
    # Standard settings.
    psrinfo["PASSED"] = 1
    psrinfo["BURNFRAC"] = 0.5
    # Now have command line input option to change these.
    print("For each of the frequency channels type 0 if it failed.")
    print("If it passed, press enter to move on.")
    for i in psrinfo.index:
        passfail = input("Channel {}, frequency {}:".format(psrinfo.iloc[i].CHAN, psrinfo.iloc[i].FREQ))
        if passfail == "0":
            print("This channel failed")
            psrinfo.at[i, "PASSED"] = 0
        else:
            pass
    print("Now type the appropriate burnin fraction for each channel")
    print("or press enter to keep the default 0.5")
    for i in psrinfo.index:
        burnfrac = input("Channel {}, frequency {}:".format(psrinfo.iloc[i].CHAN, psrinfo.iloc[i].FREQ))
        try:
            bf = float(burnfrac)
            if bf >= 1:
                print("Can't have a burnin fraction longer than the whole chain. Defaulting to 0.5")
                pass
            psrinfo.at[i, "BURNFRAC"] = bf
        except ValueError:
            pass
    # Write to new file.
    psrinfo.to_csv("{}/{}".format(args.writedir,args.output_log))
    
