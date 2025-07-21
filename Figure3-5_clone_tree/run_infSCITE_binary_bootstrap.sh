#!/bin/bash
# This script runs infSCITE in a loop for 100 iterations with different random seeds.

# infSCITE execution file path
infSCITE="infSCITE"

input=$1
basename=$2
rep=$3 # the desired number of repetitions of the MCMC.
chainLength=$4 # the desired chain length of each MCMC repetition
ad1=$5 # the estimated rate of missed heterozygous mutations in the sequencing experiment, or allele dropout rate
pct_cells=$6 # the percentage of cells to be used in the analysis for bootstrapping

sd_alpha=0.1 # Sets standard deviation of prior for false positive error rate (default is 0.1). 
sd_beta=0.1 # Sets standard deviation of prior for false negative error rate (default is 0.1).

outdir=${basename}_pctCells$pct_cells

for i in `seq 1 1 100`
do 
    echo $i

    outdir_MutTree=$outdir/${basename}_binary_r${rep}_l${chainLength}_ad${ad1}_sdAlpha${sd_alpha}_sdBeta${sd_beta}_MutTree_rndSeed/run$i
    mkdir -p $outdir_MutTree

    Rscript prepare_input.R $input $outdir_MutTree $pct_cells

    variants=`cat $outdir_MutTree/summary.txt | cut -d' ' -f1`
    cells=`cat $outdir_MutTree/summary.txt | cut -d' ' -f2`

    # Submit the job to run infSCITE
    ~/bin/sub_now.sh 10:00 10 1 $outdir_MutTree/$basename.err $outdir_MutTree/$basename.out \
            "$infSCITE -i $outdir_MutTree/data.txt \
            -names $outdir_MutTree/data.geneNames \
            -samples $outdir_MutTree/data.sample \
            -n $variants -m $cells -r $rep -l $chainLength -fd 1e-2 -ad $ad1 -e .2 -sd_alpha $sd_alpha -sd_beta $sd_beta -s -p 10000 -g 1 -d -a  \
            -o $outdir_MutTree/$basename
            "

done
