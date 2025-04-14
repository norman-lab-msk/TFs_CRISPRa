# This script runs the permutation-based enrichment testing on programs and
# clusters.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# set random seed
np.random.seed(1)

# read in list of program-cluster combos
program_cluster_list = pd.read_csv(
    snakemake.input[0]
)

# read in program-gene assignments
program_genes = pd.read_csv(
    snakemake.input[1]
)

# read in cluster-gene assignments
cluster_genes = pd.read_csv(
    snakemake.input[2]
)

# get number of programs and clusters
num_programs = program_genes['program'].nunique()
num_clusters = cluster_genes['cluster'].nunique()

# create list of genes across all programs
program_gene_list = []

for i in range(num_programs):
    genes = list(program_genes[program_genes['program'] == i]['gene'])
    program_gene_list.append(genes)

program_gene_list = [gene for program in program_gene_list for gene in program]
program_gene_list = np.array(program_gene_list)

# create arrays to hold outputs
program_names = np.empty(program_cluster_list.shape[0])
cluster_names = np.empty(program_cluster_list.shape[0], dtype = object)
overlap_values = np.empty(program_cluster_list.shape[0], dtype = object)
jaccard_similarities = np.empty(program_cluster_list.shape[0])
permutation_pvalues = np.empty(program_cluster_list.shape[0])
overlap_coefs = np.empty(program_cluster_list.shape[0])
intersect_genes = np.empty(program_cluster_list.shape[0], dtype = object)

# define clusters 
clusters = cluster_genes['cluster'].unique()

# define helper function for Jaccard similarity
def compute_jaccard(set1, set2):
    '''Computes Jaccard similarity between two sets'''
    intersection = len(np.intersect1d(set1, set2))
    union = len(np.union1d(set1, set2))
    jaccard = intersection / union
    return jaccard

# use counter to perform permutation testing
counter = 0

for i in range(program_cluster_list.shape[0]):
    program = program_cluster_list['program'][i]
    cluster = program_cluster_list['cluster'][i]
    genes_from_prog = list(program_genes[program_genes['program'] == program]['gene'])
    genes_from_cluster = list(cluster_genes[cluster_genes['cluster'] == cluster]['gene'])
    program_names[counter] = program
    cluster_names[counter] = cluster
    intersect_genes[counter] = " ".join(np.intersect1d(genes_from_prog, genes_from_cluster))
    jaccard_similarity = compute_jaccard(genes_from_prog, genes_from_cluster)
    enough_iterations = False
    num_permutations = 100
    while (enough_iterations == False):  
        permutation_similarities = [None] * num_permutations
        for j in range(num_permutations):
            shuffled_genes = np.random.choice(
                program_gene_list,
                len(genes_from_prog),
                replace = False
            )
            permutation_similarities[j] = compute_jaccard(shuffled_genes, genes_from_cluster) # jaccard for permutations
        permutation_pvalue = sum(np.array(permutation_similarities) >= jaccard_similarity) / num_permutations
        # implement iterative permutation testing to get better p-value estimates
        if num_permutations == 50000:
            enough_iterations = True
        
        if num_permutations == 25000:
            if permutation_pvalue < 0.001:
                num_permutations = 50000
            else:
                enough_iterations = True

        if num_permutations == 2500:
            if permutation_pvalue < 0.01:
                num_permutations = 25000
            else:
                enough_iterations = True

        if num_permutations == 500:
            if permutation_pvalue < 0.05:
                num_permutations = 2500
            else:
                enough_iterations = True

        if num_permutations == 100:
            if permutation_pvalue < 0.1:
                num_permutations = 500
            else:
                enough_iterations = True

    overlap_values[counter] = len(np.intersect1d(genes_from_prog, genes_from_cluster))
    jaccard_similarities[counter] = jaccard_similarity
    overlap_coefs[counter] = len(np.intersect1d(genes_from_prog, genes_from_cluster)) / len(genes_from_prog)
    permutation_pvalues[counter] = permutation_pvalue
    counter += 1

# write to output dataframe
pvalues_df = pd.concat([
    pd.Series(program_names), 
    pd.Series(cluster_names),
    pd.Series(overlap_values),
    pd.Series(jaccard_similarities),
    pd.Series(permutation_pvalues),
    pd.Series(overlap_coefs),
    pd.Series(intersect_genes)
    ], 
    axis = 1)
pvalues_df.columns = ['Program', 'Cluster', 'overlaps', 'jaccard', 'pval', 'overlap_coefs', 'intersect_genes']
pvalues_df.to_csv(
    snakemake.output[0],
    index = False
)

