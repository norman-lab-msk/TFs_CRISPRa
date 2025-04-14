# This code preprocesses clusters from the myofibroblast atlas into a
# cluster-gene mapping table.
#
# Author: Karthik Guruvayurappan 

import numpy as np
import pandas as pd

# read in cluster DEGs
cluster_degs = pd.read_excel(
    '/data/norman/southark/external_datasets/cancer_cell_fibro_atlas/1-s2.0-S1535610824003192-mmc3.xlsx',
    sheet_name = 'Table S2A',
    skiprows = 2,
    header = 0
)

# filter for necessary columns
cluster_degs = cluster_degs[['cluster', 'gene']]

# write to output file
cluster_degs.to_csv(
    'resources/myofibroblast_atlas_degs.csv',
    index = False
)

