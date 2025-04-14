# This script creates a Sankey diagram summarizing all of the atlas cell state-
# program links.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

import plotly.graph_objects as go

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

# keep clusters of interest
keep_clusters = [
    # 'CD34+MFAP5+ C9 Synovium', # universal
    # 'CD34+MFAP5+ C9 SalivaryGland', # universal
    'CD34+MFAP5+ C9 Marginal', # universal
    # 'CD34+MFAP5+ C9 Lung', # universal
    # 'CD34+MFAP5+ C9 Gut', # universal
    'PI16+ (n = 779)', # universal
    'PRG4 Fib', # universal
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
    'c16', # inflammatory
    # 'MYH11+ C13 Synovium', # myofibroblast
    # 'MYH11+ C13 SalivaryGland', # myofibroblast
    'MYH11+ C13 Marginal', # myofibroblast
    # 'MYH11+ C13 Lung', # myofibroblast
    # 'MYH11+ C13 Gut', # myofibroblast
    'Myofibroblast', # myofibroblast
    'c01', # myofibroblast
    'ADAMDEC1+ (n = 778)', # immune/antigen presentation
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

# concat nodes and filter for those with significnat results
nodes = pd.concat([pd.Series(keep_programs), pd.Series(keep_clusters)])
results_df = results_df[results_df['Cluster'].isin(keep_clusters)]
results_df = results_df[results_df['Program'].isin(keep_programs)]
sig_prog_clusters = pd.concat([results_df['Program'], results_df['Cluster']]).unique()
nodes = nodes[nodes.isin(sig_prog_clusters)]
kept_programs = pd.Series(keep_programs)[pd.Series(keep_programs).isin(sig_prog_clusters)]
kept_clusters = pd.Series(keep_clusters)[pd.Series(keep_clusters).isin(sig_prog_clusters)]

# filter for necessary columns
results_df = results_df[['Program', 'Cluster', 'jaccard', 'overlaps', 'overlap_coefs', 'pval']]
results_df['log_p'] = -1 * np.log10(results_df['pval'] + 0.00001)
results_df['significant'] = 1

# map to indices
node_idxs = {label: i for i, label in enumerate(nodes)}
# map programs and clusters to indices
results_df['program_idx'] = results_df['Program'].map(node_idxs)
results_df['cluster_idx'] = results_df['Cluster'].map(node_idxs)

# make x positions
num_programs = len(kept_programs)
num_clusters = len(kept_clusters)
x_positions = ([0.7] * num_programs) + ([0.3] * num_clusters)
y_positions = list(np.linspace(0.1, 0.9, num=num_programs)) + \
    list(np.linspace(0.1, 0.9, num=num_clusters))


# write results to output df
results_df.to_csv(
    'results/sankey_diagram_results.csv',
    index = False
)

# read in results_df
results_df = pd.read_csv(
    'results/sankey_diagram_results.csv'
)

# 1. #864D9E - A deep purplish hue.
# 2. #376ADC - A softened medium blue with reduced saturation.
# 3. #AE5981 - A muted pinkish-purple.
# 4. #5AC2B0 - A teal or turquoise-like color.
# 5. #D99828 - A golden yellow.
# 6. #528FF0 - A light blue shade.
# 7. #CE6DA1 - A pinkish red
# 8. #'#D3D3D3', '#808080' #grays

# plot Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=15,
        line=dict(color="black", width=0.5),
        # label=list(nodes),
        x=x_positions,
        y=y_positions,
        color = (['#528FF0', '#D99828', '#D99828', '#CE6DA1', '#5AC2B0']) + (['#528FF0'] * 5) + (['#D99828'] * 6) + (['#CE6DA1'] * 3) + (['#5AC2B0'] * 5),
        align = 'left'
    ),
    link=dict(
        source=results_df['program_idx'],
        target=results_df['cluster_idx'],
        value=results_df['overlap_coefs']
    )
)])

fig.write_image('results/sankey_diagram.pdf', format = 'pdf')

