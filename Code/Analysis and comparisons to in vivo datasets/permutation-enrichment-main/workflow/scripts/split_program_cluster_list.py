# This script creates a dataframe of program-cluster pairs and divides it into
# 32 batches.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# read in program-gene dataframe
program_genes = pd.read_csv(
    snakemake.input[0]
)

# read in cluster-gene dataframe
cluster_genes = pd.read_csv(
    snakemake.input[1]
)

# get unique clusters and programs
programs = program_genes['program'].unique()
clusters = cluster_genes['cluster'].unique()

# get all combinations and put in a dataframe
programs_clusters = pd.DataFrame()

for program in programs:
    program_df = pd.DataFrame(clusters)
    program_df.columns = ['cluster']
    program_df['program'] = program
    programs_clusters = pd.concat([programs_clusters, program_df], axis = 0)


# reorder columns, split into batches, and write to CSV
programs_clusters = programs_clusters[['program', 'cluster']]
programs_clusters_batches = np.array_split(programs_clusters, 32)
for i, batch in enumerate(programs_clusters_batches):
    batch.to_csv(
        snakemake.output[i],
        index = False
    )

