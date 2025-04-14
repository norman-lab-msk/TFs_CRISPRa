# This script processes the DEG lists from the Buechler et al. fibroblast atlas
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# create Excel file object
atlas_sheet = pd.ExcelFile('/data/norman/southark/external_datasets/fibroblast_atlas_2021/supplemental_tables/41586_2021_3549_MOESM10_ESM.xlsx')
sheet_names = atlas_sheet.sheet_names

# read in atlas DEGs
atlas_degs = pd.DataFrame()

for cluster in sheet_names:
    cluster_degs = pd.read_excel(
        '/data/norman/southark/external_datasets/fibroblast_atlas_2021/supplemental_tables/41586_2021_3549_MOESM10_ESM.xlsx',
        sheet_name = cluster
    )
    cluster_degs = cluster_degs.iloc[:, 0:6]
    cluster_degs.columns = [
        'gene',
        'p_val',
        'avg_log2FC',
        'pct.1',
        'pct.2',
        'p_val_adj'
    ]
    cluster_degs['cluster'] = cluster
    cluster_degs = cluster_degs.sort_values(by = 'avg_log2FC', ascending = False)
    cluster_degs = cluster_degs.head(150)
    atlas_degs = pd.concat([atlas_degs, cluster_degs], axis = 0)

# remove "false" clusters
atlas_degs = atlas_degs[
    ~atlas_degs['cluster'].isin(['Up in COL6A3vLRRC15 (n = 137)', 'Up in LRRC15vCOL6A3 (n = 640) '])
]


# organize into cluster-gene format and write to CSV
atlas_degs = atlas_degs[['cluster', 'gene']]

atlas_degs.to_csv(
    'resources/buechler_atlas_degs.csv',
    index = False
)


