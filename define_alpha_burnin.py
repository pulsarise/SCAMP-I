#! /usr/bin/env python3

import numpy as np
import argparse
import pandas

if __name__ == '__main__':
    # Define options to the script.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--run_log', type=str, help="Log information about the alpha MCMC run.")
    args = parser.parse_args()
    # Read in the info file.
    psrinfo = pandas.read_csv(args.run_log, index_col=0)
    # Standard settings.
    psrinfo["ALPHABURNFRAC"] = 0.5
    print("Type the desired burnin fraction")
    print("or press enter to keep the default 0.5")
    burnfrac = input("Burnin fraction:")
    try:
        bf = float(burnfrac)
        if bf >= 1:
            print("Can't have a burnin fraction longer than the whole chain. Defaulting to 0.5")
            pass
        psrinfo["ALPHABURNFRAC"] = bf
    except ValueError:
        pass
    # Write to new file.
    psrinfo.to_csv(args.run_log)
    
