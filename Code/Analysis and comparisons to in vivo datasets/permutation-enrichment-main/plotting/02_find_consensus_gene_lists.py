# This script finds genes that are shared between cluster-gene links in the
# Sankey diagram.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# import results across all atlases
results_dir = 'results/results_archive/'
atlases = [
    'buechler_atlas',
    'cardiac_fibro',
    'myofibroblast',
    'korsunsky_atlas'
]

# read in permutation results from each atlas
results_df = pd.DataFrame()

for atlas in atlases:
    atlas_results = pd.read_csv(
        results_dir + atlas + '/permutation_results.csv'
    )
    atlas_results['atlas'] = atlas
    results_df = pd.concat([results_df, atlas_results], axis = 0)

# filter for significant results
results_df = results_df[results_df['adj_p'] < 0.1]

# convert intersect_genes to list
results_df['intersect_genes'] = results_df['intersect_genes'].apply(lambda x: x.split(' '))

# filter for programs and clusters of interest
# keep clusters of interest 
keep_clusters = [
    # 'CD34+MFAP5+ C9 Synovium', # universal
    # 'CD34+MFAP5+ C9 SalivaryGland', # universal
    'CD34+MFAP5+ C9 Marginal', # universal
    # 'CD34+MFAP5+ C9 Lung', # universal
    # 'CD34+MFAP5+ C9 Gut', # universal
    'PI16+ (n = 779)', # universal
    'POLCE2.MFAP5 Fib', # universal
    'ground state Fib', # universal
    'c05', # universal
    'c03', # universal
    # 'SPARC+COL3A1+ C4 Synovium', # inflammatory
    # 'SPARC+COL3A1+ C4 SalivaryGland', # inflammatory
    'SPARC+COL3A1+ C4 Marginal', # inflammatory
    # 'SPARC+COL3A1+ C4 Lung', # inflammatory
    # 'SPARC+COL3A1+ C4 Gut', # inflammatory
    'COL3A1 (n = 72)', # inflammatory
    'LRRC15+ (n = 472)', # inflammatory
    'POSTN+ Fib', # inflammatory
    'c04', # inflammatory
    # 'MYH11+ C13 Synovium', # myofibroblast
    # 'MYH11+ C13 SalivaryGland', # myofibroblast
    'MYH11+ C13 Marginal', # myofibroblast
    # 'MYH11+ C13 Lung', # myofibroblast
    # 'MYH11+ C13 Gut', # myofibroblast
    'Myofibroblast', # myofibroblast
    'c01', # myofibroblast
    # 'CXCL10+CCL19+ C11 Synovium', # immune/antigen presentation 
    # 'CXCL10+CCL19+ C11 SalivaryGland', # immune
    'CXCL10+CCL19+ C11 Marginal', # immune
    # 'CXCL10+CCL19+ C11 Lung', # immune
    # 'CXCL10+CCL19+ C11 Gut', # immune
    'CCL19+ (n = 400)', # immune
    'Type I IFN Fib', # immune
    'c19' # immune
]

# keep programs of interest
keep_programs = [
    18, # universal
    31, # universal
    27, # inflammatory
    32, # inflammatory
    4, # myofibroblast
    19 # immune
]

results_df = results_df[results_df['Cluster'].isin(keep_clusters)]
results_df = results_df[results_df['Program'].isin(keep_programs)]

# iterate through remaining programs and count gene frequency across clusters
output_df = pd.DataFrame()
kept_programs = pd.Series(keep_programs)[pd.Series(keep_programs).isin(results_df['Program'].unique())]
for program in kept_programs:
    program_results = results_df[results_df['Program'] == program]
    program_genes = program_results['intersect_genes']
    program_genes = program_genes.sum()
    program_gene_counts = pd.Series(program_genes).value_counts()
    program_gene_counts = pd.DataFrame(program_gene_counts)
    program_gene_counts.reset_index(inplace = True)
    program_gene_counts.columns = ['gene', 'count']
    program_gene_counts['program'] = program
    output_df = pd.concat([output_df, program_gene_counts], axis = 0)


# write to output CSV
output_df = output_df[['program', 'gene', 'count']]
output_df.to_csv(
    'results/consensus_genes.csv',
    index = False
)