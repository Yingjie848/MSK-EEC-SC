#!/bin/bash
# This script runs infSCITE for the EEC14 dataset with specified parameters.

input="EEC14.genotypes.csv"
basename="EEC14"
rep=10              # the desired number of repetitions of the MCMC.
chainLength=10000   # the desired chain length of each MCMC repetition
ad1=0.078           # the estimated rate of missed heterozygous mutations in the sequencing experiment, or allele dropout rate
pct_cells=75        # the percentage of cells to be used in the analysis for bootstrapping

sh run_infSCITE_binary_bootstrap.sh $input $basename $rep $chainLength $ad1 $pct_cells