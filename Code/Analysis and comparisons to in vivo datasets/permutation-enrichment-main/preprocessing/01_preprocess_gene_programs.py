# This script creates a list of gene-program assignments based on the CRISPRa
# fibroblast dataset.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# define path to CRISPRa data
data_path = '/data/norman/southark/external_datasets/fibro_CRISPRa_Tfs/'

# read in gene programs in matrix form (rows are programs, columns are genes)
programs_pca = pd.read_csv(
    f'{data_path}20240331_fibroblast_bulk_comps.csv', 
    index_col=0
)

# create dataframe to hold program-gene mapping
program_df = pd.DataFrame()

for i in range(programs_pca.shape[0]):
    program_genes = programs_pca.iloc[i, :]
    program_genes = program_genes[program_genes > 0]
    program_genes = program_genes.sort_values(ascending = False)
    program_genes = program_genes.head(50).index
    program_genes = pd.DataFrame(program_genes)
    program_genes.columns = ['gene']
    program_genes['program'] = i
    program_df = pd.concat([program_df, program_genes], axis = 0)


# write to output CSV
program_df = program_df[['program', 'gene']]
program_df.to_csv(
    'resources/crispra_program_gene_top50pca.csv',
    index = False
)

