#!/bin/bash

#BSUB -W 2:00                 # for 120 hours of wall clock time
#BSUB -J run_snakemake          # Job title
#BSUB -o run_snakemake.out   # Output file name
#BSUB -e run_snakemake.err   # Error file name

# activate snakemake mamba environment
source $HOME/.bashrc
mamba activate snakemake

# run snakemake for whole pipeline (ending with volcano plot)
snakemake --use-conda --jobs 32 --cluster 'bsub -W 1:00 -n 1 -R "rusage[mem=8G]" -o out.%J.txt -e err.%J.txt' all

# Random comments pertaining to various components of the job submission.
# .%J adds job ID number to output files
# run using bsub < run_snakemake.sh. '<'