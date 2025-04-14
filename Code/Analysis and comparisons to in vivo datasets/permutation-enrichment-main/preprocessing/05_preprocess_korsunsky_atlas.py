# This script preprocesses the gene lists from the fibroblast atlas in
# Korusnsky et al.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# load differentially expressed gene lists from Korsunsky et al.
in_vivo_degs = pd.read_excel(
    '/data/norman/southark/external_datasets/fibroblast_atlas_med_2022/' + \
    'SupplementaryData.xlsx',
    sheet_name = 'TableS8'
)

# add tisuse name to cluster
in_vivo_degs['Cluster'] = in_vivo_degs['Cluster'] + ' ' + in_vivo_degs['Tissue']

# create dataframe to hold outputs
cluster_df = pd.DataFrame()

# get unique clusters to iterate
clusters = in_vivo_degs['Cluster'].unique()

# create dataframe of cluster-gene assignments
for cluster in clusters:
    cluster_genes = in_vivo_degs[in_vivo_degs['Cluster'] == cluster]['Feature']
    cluster_genes = pd.DataFrame(cluster_genes)
    cluster_genes.columns = ['gene']
    cluster_genes['cluster'] = cluster
    cluster_df = pd.concat([cluster_df, cluster_genes], axis = 0)


# reorder columns and write to output CSV
cluster_df = cluster_df[['cluster', 'gene']]
cluster_df.to_csv(
    'resources/korsunsky_atlas_degs.csv',
    index = False
)

