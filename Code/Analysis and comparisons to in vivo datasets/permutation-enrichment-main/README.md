# Program-Cluster Gene Enrichment Analysis

This repository contains a snakemake pipeline to perform permutation-based gene enrichment analysis. This pipeline can be used to associate gene programs with cell cluster marker genes from in vivo cell atlas projects.

## Setup

In order to use this repository, you will need to install:
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [conda](https://docs.anaconda.com/) or [mamba](https://mamba.readthedocs.io/en/latest/index.html)

I recommend creating a conda environment with snakemake installed, and then running the pipeline from there. Everything should be taken care of for you!

## Inputs

To run the pipline, you will need a program-gene table. An example can be found in `resources/crispra_program_gene_top50pca.csv`. 

You will also need a cluster-gene table. An example can be found in `resources/cardiac_fibroblast_degs.csv`.

Once you have these two files, place them in the `resources/` folder and then modify the `programs` and `clusters` parameters in the config file. An example config file can be found in `config/config_cardiac_fibro.yaml`. Be sure to also update your config file in the Snakefile (first line).

Then you should be good to run the pipeline! 


## Approach
As a generalization, this procedure can be applied between any two sets of gene-label maps. It will take the genes across all the programs, and resample them to generate a null distribution of gene-program assignments. Then, it will use Jaccard similarity to test if the observed program-cluster overlap is more than would be expected by random chance by computing a permutation testing-based p-value.

## Usage

Once all of this is set up, the pipeline should run as follows:

```snakemake --use-conda all```

The script is written to create 32 separate jobs for parallelization. Therefore, I recommend using the `--jobs` option if running on a cluster. An example job submission script I use is in `run_snakemake.sh` and I submit it using `bsub < run_snakemake.sh`

## Contact
Karthik Guruvayurappan (guruvak@mskcc.org)

