#!/bin/env bash

####################################
# POINT TO THE DATA AND CONFIG FILE.
####################################
# Specify base name you want for the logs and output csv files. A standard option would be to use the pulsar's name.
basename="simpsr"
# Give the link to the file containing the basic information required.
configfile="simpsr_config.csv"

################################
# SET OPTIONS FOR THE MODELLING.
################################
# Specify runtime for MCMCs. Defaults 20,000 (for tau) and 10,000 (for alpha) if run from the command line.
taumcruntime=5000
alphamcruntime=5000
# Select model to fit: either "iso" for isotropic, or "onedim" for extremely anisotropic. Default "iso" if not called from command line.
model="iso"
# Choose how many walkers to use for the tau MCMC fit. Default 10 if not called from command line.
tau_nwalkers=20
# Choose how many walkers to use for the alpha MCMC fit. Default 32 if not called from command line.
alpha_nwalkers=20

#########################################################################
# CREATE DIRECTORY STRUCTURE TO HOLD THE OUTPUT OF THE MODELLING PROCESS.
#########################################################################
mkdir -p logs
mkdir -p output/diagnostic_plots
mkdir -p output/chains
mkdir -p output/results/plots

##################
# RUN THE SCRIPTS
##################
echo "Fitting for tau..."
python $MYCODE/SCAMP-I/run_scatter_mcmc.py -f ${configfile} -o logs/${basename}_inputlog.csv -w output/chains -t ${taumcruntime} -nw ${tau_nwalkers} -m ${model} --showprogress
python $MYCODE/SCAMP-I/check_tau_chains.py -i logs/${basename}_inputlog.csv -w output/diagnostic_plots
echo "Now check the chain diagnostic plots saved in the output/diagnostic_plots/ folder and decide whether they have passed or failed and what burnin to apply."
python $MYCODE/SCAMP-I/define_tau_passfail_burnin.py -i logs/${basename}_inputlog.csv -o ${basename}_outputlog.csv -w logs
python $MYCODE/SCAMP-I/parameter_extraction.py -i logs/${basename}_outputlog.csv -o ${basename}_output.csv -w output/results
echo "Fitting for alpha..."
python $MYCODE/SCAMP-I/alpha_mcmc.py -i output/results/${basename}_output.csv -o ${basename}_output.csv --sampleswritedir output/chains --resultswritedir output/results -t ${alphamcruntime} -nw ${alpha_nwalkers} --showprogress
python $MYCODE/SCAMP-I/check_alpha_chains.py -i output/results/${basename}_output.csv -w output/diagnostic_plots
echo "Now check the alpha chain plots and decide what burnin fraction you want."
python $MYCODE/SCAMP-I/define_alpha_burnin.py -i output/results/${basename}_output.csv
python $MYCODE/SCAMP-I/best_fit_alpha_tau.py -i output/results/${basename}_output.csv
python $MYCODE/SCAMP-I/basic_plotting.py -f output/results/${basename}_output.csv -w output/results/plots
echo "Results and plots saved to directories."