# This script creates a cluster-gene mapping for the differentially expressed
# genes from the cardiac fibroblast atlas.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# load differentially expressed gene lists from cardiac fibroblast atlas paper
cardiac_degs = pd.read_excel(
    '/data/norman/karthik/permutation_analysis/' + 
        'data/external/41586_2024_8008_MOESM4_ESM_KGversion.xlsx',
    sheet_name = 'Supplementary Table 8',
    header = 0,
    index_col = 0,
    skiprows = 1
)

# create a dataframe of clusters-genes
cluster_df = pd.DataFrame()

clusters = cardiac_degs['cluster'].unique()

for cluster in clusters:
    cluster_genes = list(cardiac_degs[cardiac_degs['cluster'] == cluster]['gene'])
    cluster_genes = pd.DataFrame(cluster_genes)
    cluster_genes.columns = ['gene']
    cluster_genes['cluster'] = cluster
    cluster_df = pd.concat([cluster_df, cluster_genes], axis = 0)


# write to output CSV
cluster_df = cluster_df[['cluster', 'gene']]
cluster_df.to_csv(
    'resources/cardiac_fibroblast_degs.csv',
    index = False
)

