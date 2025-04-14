# Perturbseq library for loading and manipulating single-cell experiments
# Copyright (C) 2019  Thomas Norman

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import pandas as pd
import numpy as np
from sklearn import preprocessing as pre
from pandas.api.types import is_numeric_dtype
from six.moves import zip as izip
from time import time
import gc
from tqdm import tqdm_notebook
from scipy.spatial import KDTree
#from pandarallel import pandarallel
from .transformers import PCAReducer
from scipy.stats import trim_mean

def strip_low_expression(pop, threshold=0):
    """Remove genes with low or zero expression to reduce memory usage. Modifies the
    target CellPopulation in place.
    
    Args:
        pop: CellPopulation instance
        threshold: all genes with expression <= threshold will be removed
    """
    retain = pop.genes.query('mean > @threshold').index
    remove = pop.genes.query('mean <= @threshold').index
    if len(remove) == 0:
        print('No genes have expression below threshold.')
        return
    pop.matrix = pop.matrix[retain]
    pop.genes.loc[remove, 'in_matrix'] = False
    # set all numeric properties to nan for genes that have been removed 
    for col in np.setdiff1d(pop.genes.columns, ['gene_name', 'in_matrix']):
        if is_numeric_dtype(pop.genes[col]):
            pop.genes.loc[remove, col] = np.nan
    # force garbage collection
    gc.collect()
                
def equalize_UMI_counts(matrix, median_umi_count=None):
    """Normalize all cells in an expression matrix to a specified UMI count
    """
    reads_per_bc = matrix.sum(axis=1)
    if median_umi_count is None:
        median_reads_per_bc = np.median(reads_per_bc)
    else:
        median_reads_per_bc = median_umi_count
    scaling_factors = median_reads_per_bc / reads_per_bc
    m = matrix.astype(np.float64)
    # Normalize expression within each cell by median total count
    m = m.mul(scaling_factors, axis=0)
    if np.mean(median_reads_per_bc) < 5000:
        print("Scaling with a small number of reads. Are you sure this is what you want?")
    return m

def log_normalize_expression(pop, scale_by_total=True, pseudocount=1, median_umi_count=None):
    """ Normalize expression distribution by log transformation.
    The normalization proceeds by first (optionally) normalizing the UMI counts within each cell 
    to the median UMI count within the population. The expression within the population is then 
    log-normalized: i.e., transformed according to Y = log2(X + 1)
       
    Args:
        pop: population to normalize
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
        pseudocount: offset for 0 values (default: 1)
        
    Returns:
        DataFrame of log-normalized expression data
    """
    matrix = pop.matrix
    
    if (scale_by_total):
        m = equalize_UMI_counts(matrix, median_umi_count=None)
    else:
        m = matrix.astype(np.float64) 

    m = np.log2(m + pseudocount)
    
    return pd.DataFrame(m, columns=m.columns, index=m.index)
    
def z_normalize_expression(pop, scale_by_total=True):
    """ Normalize expression distribution by Z-scoring.
    The normalization proceeds by first normalizing the UMI counts within each cell to the 
    median UMI count within the population. The expression within the population is then 
    Z-normalized: i.e., for each gene the mean is subtracted, and then these values are divided 
    by that gene's standard deviation.
       
    Args:
        pop: population to normalize
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
    
    Returns:
        DataFrame of Z-normalized expression data
    """
    matrix = pop.matrix
    
    if (scale_by_total):
        m = equalize_UMI_counts(matrix)
    else:
        m = matrix.astype(np.float64) 

    # Now center and rescale the expression of each gene to average 0 and std 1
    m_out = pre.scale(m.as_matrix(), axis=0)
    
    return pd.DataFrame(m_out, columns=m.columns, index=m.index)

def normalize_matrix_to_control(matrix, control_matrix, scale_by_total=True, median_umi_count=None):
    """ Normalize expression distribution relative to a population of control (unperturbed) cells.
    The normalization proceeds by first normalizing the UMI counts within each cell to the median
    UMI count within the population. The control population is normalized to the same amount. 
    The expression within the population is then Z-normalized with respect to the control 
    distribution: i.e., for each gene the control mean is subtracted, and then these values are
    divided by the control standard deviation.
        
    Args:
        matrix: gene expression matrix to normalize (output from cellranger)
        control_matrix: gene expression matrix of control population
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
        median_umi_count: If provided, set the median to normalize to. This is useful for example
            when normalizing multiple independent lanes/gemgroups.
    
    Returns:
        DataFrame of normalized expression data
        
    Example:
        >>>control_cells = pop.cells[pop.cells['perturbation']=='DMSO_control'].index.values
        >>>pop.normalized_matrix = normalize_expression_to_control(pop.matrix,
                                                                   pop.matrix.loc[control_cells].copy())

    """
    if isinstance(matrix, pd.SparseDataFrame):
        print('     Densifying matrix...')
        matrix = matrix.to_dense()
        
    if isinstance(control_matrix, pd.SparseDataFrame):
        print('     Densifying control matrix...')
        control_matrix = control_matrix.to_dense()

    if (scale_by_total):
        print("     Determining scale factors...")
        reads_per_bc = matrix.sum(axis=1)
        
        if median_umi_count is None:
            median_reads_per_bc = np.median(reads_per_bc)
        else:
            median_reads_per_bc = median_umi_count
        
        scaling_factors = median_reads_per_bc / reads_per_bc
        
        print("     Normalizing matrix to median")
        m = matrix.astype(np.float64)
        # Normalize expression within each cell by median total count
        m = m.mul(scaling_factors, axis=0)
        if np.mean(median_reads_per_bc) < 5000:
            print("Scaling with a small number of reads. Are you sure this is what you want?")
            
        control_reads_per_bc = control_matrix.sum(axis=1)
        
        print("     Normalizing control matrix to median")
        control_scaling_factors = median_reads_per_bc / control_reads_per_bc
        c_m = control_matrix.astype(np.float64)
        c_m = c_m.mul(control_scaling_factors, axis=0)
    else:
        m = matrix.astype(np.float64)
        c_m = matrix.astype(np.float64)

    control_mean = c_m.mean()
    control_std = c_m.std()
    
    print("     Scaling matrix to control")
    # Now center and rescale the expression of each gene to average 0 and std 1
    m_out = (m - control_mean).div(control_std, axis=1)
    
    print("     Done.")
    return pd.DataFrame(m_out, columns=m.columns, index=m.index)


def trim_std(x, proportiontocut):
    mu = trim_mean(x, proportiontocut)
    return np.sqrt(np.mean((x - mu)**2))

def robust_normalize_matrix_to_control(matrix, control_matrix, scale_by_total=True, median_umi_count=None, proportiontocut=0.02):
    """ Normalize expression distribution relative to a population of control (unperturbed) cells.
    The normalization proceeds by first normalizing the UMI counts within each cell to the median
    UMI count within the population. The control population is normalized to the same amount. 
    The expression within the population is then Z-normalized with respect to the control 
    distribution: i.e., for each gene the control mean is subtracted, and then these values are
    divided by the control standard deviation.
        
    Args:
        matrix: gene expression matrix to normalize (output from cellranger)
        control_matrix: gene expression matrix of control population
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
        median_umi_count: If provided, set the median to normalize to. This is useful for example
            when normalizing multiple independent lanes/gemgroups.
    
    Returns:
        DataFrame of normalized expression data
        
    Example:
        >>>control_cells = pop.cells[pop.cells['perturbation']=='DMSO_control'].index.values
        >>>pop.normalized_matrix = normalize_expression_to_control(pop.matrix,
                                                                   pop.matrix.loc[control_cells].copy())

    """
    if isinstance(matrix, pd.SparseDataFrame):
        print('     Densifying matrix...')
        matrix = matrix.to_dense()
        
    if isinstance(control_matrix, pd.SparseDataFrame):
        print('     Densifying control matrix...')
        control_matrix = control_matrix.to_dense()

    if (scale_by_total):
        print("     Determining scale factors...")
        reads_per_bc = matrix.sum(axis=1)
        
        if median_umi_count is None:
            median_reads_per_bc = np.median(reads_per_bc)
        else:
            median_reads_per_bc = median_umi_count
        
        scaling_factors = median_reads_per_bc / reads_per_bc
        
        print("     Normalizing matrix to median")
        m = matrix.astype(np.float64)
        # Normalize expression within each cell by median total count
        m = m.mul(scaling_factors, axis=0)
        if np.mean(median_reads_per_bc) < 5000:
            print("Scaling with a small number of reads. Are you sure this is what you want?")
            
        control_reads_per_bc = control_matrix.sum(axis=1)
        
        print("     Normalizing control matrix to median")
        control_scaling_factors = median_reads_per_bc / control_reads_per_bc
        c_m = control_matrix.astype(np.float64)
        c_m = c_m.mul(control_scaling_factors, axis=0)
    else:
        m = matrix.astype(np.float64)
        c_m = matrix.astype(np.float64)

    control_mean = trim_mean(c_m, proportiontocut)
    control_std = np.apply_along_axis(lambda x: trim_std(x, proportiontocut), 0, c_m)
    
    print("     Robust scaling matrix to control")
    # Now center and rescale the expression of each gene to average 0 and std 1
    m_out = (m - control_mean).div(control_std, axis=1)
    
    print("     Done.")
    return pd.DataFrame(m_out, columns=m.columns, index=m.index)

def normalize_to_control(pop, control_cells, scale_by_total=True, median_umi_count=None, **kwargs):
    """ Normalize expression distribution relative to a population of control (unperturbed) cells.
    The normalization proceeds by first normalizing the UMI counts within each cell to the median
    UMI count within the population. The control population is normalized to the same amount. 
    The expression within the population is then Z-normalized with respect to the control 
    distribution: i.e., for each gene the control mean is subtracted, and then these values are
    divided by the control standard deviation.
        
    Args:
        pop: CellPopulation to normalize
        control_cells: metadata condition passed to pop.where to identify control cell population
            to normalize with respect to
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
        median_umi_count: If provided, set the median to normalize to. This is useful for example
        **kwargs: all other arguments are passed to pop.groupby, and hence to pop.where. These can
            be used to do more refined slicing of the population if desired.
            when normalizing multiple independent lanes/gemgroups.
    
    Returns:
        DataFrame of normalized expression data
        
    Example:
        >>>pop.normalized_matrix = normalize_expression_to_control(pop, 'perturbed == "control"')
    """
    matrix = pop.where(**kwargs)
    control_matrix = pop.where(cells=control_cells, **kwargs)
    return normalize_matrix_to_control(matrix, control_matrix, scale_by_total=scale_by_total, median_umi_count=median_umi_count)

def normalize_to_gemgroup_control(pop, control_cells, median_umi_count=None, robust=False, scale_by_total=True, **kwargs):
    """Normalizes a multi-lane 10x experiment. Cell within each gemgroup are normalized to the 
    control cells within the same gemgroup.
        
    Args:
        pop: CellPopulation to normalize
        control_cells: metadata condition passed to pop.where to identify control cell population
            to normalize with respect to
        median_umi_count: Value to normalize UMI counts to across lanes. If None (the default)
            then all cells are normalized to the median UMI count of control cells within the
            whole experiment.
        **kwargs: all other arguments are passed to pop.groupby, and hence to pop.where. These can
            be used to do more refined slicing of the population if desired.
    All other arguments are passed to normalize_matrix_to_control
    
    Returns:
        DataFrame of normalized expression data    
    
    Example:
        normalized_matrix = normalize_to_gemgroup_control(pop,
                                                          control_cells='guide_identity == "control"')
        will produce a normalized expression matrix where cells in each gemgroup are Z-normalized 
        with respect to the expression distribution of cells in that lane bearing the guide_identity
        "control".
    """
    
    gemgroup_iterator = izip(pop.groupby('gem_group', **kwargs),
                        pop.groupby('gem_group', cells=control_cells, **kwargs))
    
    gem_group_matrices = dict()

    if median_umi_count is None:
        median_umi_count = pop.where(cells=control_cells, **kwargs).sum(axis=1).median()
    if scale_by_total:
        print('Normalizing all cells to {0} UMI...'.format(median_umi_count))
    
    for (i, gemgroup_pop), (_, gemgroup_control_pop) in gemgroup_iterator:
        print('Processing gem group {0}'.format(i))
        t = time()
        if not robust:
            gem_group_matrices[i] = normalize_matrix_to_control(gemgroup_pop,
                                                                gemgroup_control_pop,
                                                                median_umi_count=median_umi_count,
                                                               scale_by_total=scale_by_total)
        else:
            gem_group_matrices[i] = robust_normalize_matrix_to_control(gemgroup_pop,
                                                                gemgroup_control_pop,
                                                                median_umi_count=median_umi_count,
                                                                      scale_by_total=scale_by_total)
        print(time() - t)
    print('Merging submatrices...')
    return pd.concat(gem_group_matrices.values(), verify_integrity=False).loc[pop.matrix.index]

def normalize_matrix_by_key(pop, key):
    subpop_matrices = dict()

    for name, matrix in tqdm_notebook(pop.groupby(key)):
        print(name)
        subpop_matrices[name] = equalize_UMI_counts(matrix)
    print('Merging submatrices...')
    return pd.concat(subpop_matrices.values()).loc[pop.matrix.index]

def inherit_normalized_matrix(pop, parent_pop):
    """Subset a parent population's normalized expression matrix to the cells within a given population
    """
    if isinstance(parent_pop.normalized_matrix,pd.DataFrame):
        pop.normalized_matrix = parent_pop.normalized_matrix.loc[pop.matrix.index, pop.matrix.columns]
    elif isinstance(parent_pop.normalized_matrix,dict):
        pop.normalized_matrix = {}
        for normalized_matrix_name in parent_pop.normalized_matrix.keys():
            pop.normalized_matrix[normalized_matrix_name] = parent_pop.normalized_matrix[normalized_matrix_name].loc[pop.matrix.index, pop.matrix.columns]
            
def UMI_relative_to_gemgroup_control(pop, control_cells, **kwargs):
    """Calculates "fraction of total UMIs" and "Z-score total UMIs" relative to gemgroup control cells. 
        
    Args:
        pop: CellPopulation to normalize
        control_cells: metadata condition passed to pop.where to identify control cell population
            to normalize with respect to
        **kwargs: all other arguments are passed to pop.groupby, and hence to pop.where. These can
            be used to do more refined slicing of the population if desired.
    
    Returns:
        DataFrame of fraction/z-scored UMI counts    
    
    Example:
        pop.cells[['fraction_gemgroup_UMI','z_gemgroup_UMI']] = UMI_relative_to_gemgroup_control(pop,
                                                            control_cells='guide_identity == "control"')
        will produce a "fraction of total UMIs" matrix where cells in each gemgroup are  
        with respect to the expression distribution of cells in that lane bearing the guide_identity
        "control".
    """
    
    gemgroup_iterator = izip(pop.groupby('gem_group', **kwargs),
                        pop.groupby('gem_group', cells=control_cells, **kwargs))
    
    gem_group_relative_matrices = dict()
    gem_group_z_matrices = dict()
    
    for (i, gemgroup_pop), (_, gemgroup_control_pop) in gemgroup_iterator:
        gem_group_relative_matrices[i] = gemgroup_pop.sum(axis=1) / gemgroup_control_pop.sum(axis=1).mean()
        gem_group_z_matrices[i] = (gemgroup_pop.sum(axis=1) - gemgroup_control_pop.sum(axis=1).mean())/gemgroup_control_pop.sum(axis=1).std()
    print('Merging submatrices...')
    f = pd.concat(gem_group_relative_matrices.values())
    z = pd.concat(gem_group_z_matrices.values())
    return pd.concat([z,f],axis=1).loc[pop.matrix.index]

def normalize_to_near_neighbor_control(pop, control_cells, normalized_name, n_near_neighbors=50, n_components=50, nthreads_pca=50, nthreads_apply=10, **kwargs):
    """Normalizes a perturb-seq experiment to control cells. Input data should first be UMI normalized/z-scored to enable PCA. Cells will be normalized to the 
    nearest neighbor control cells, determined via kd-tree nearest neighbor search.
    
    inspired by mixscape with some modification: https://www.biorxiv.org/content/10.1101/2020.06.28.175596v1
    
    Currently, I first UMI normalize and z-score with respect to non-targeting control (normalize_to_control).
    Then, I perform PCA on control cells and then embed perturbed cells in this space.
    My reasoning for this strategy is that we intend to correct for sources of variation present in the control cells.
        
    Args:
        pop: CellPopulation to normalize
        control_cells: metadata condition passed to pop.where to identify control cell population
            to normalize with respect to
        normalized_name: name of the UMI normalized/z-scored input matrix
        n_near_neighbors: number of nearest neighbors to normalize against
        n_components: number of principal components to use for nearest neighbor search
        nthreads: number of threads for parallel processing
        **kwargs: all other arguments are passed to pop.groupby, and hence to pop.where. These can
            be used to do more refined slicing of the population if desired.
    
    Returns:
        DataFrame of normalized expression data    
    
    Example:
    ##scale data to equal number of UMIs per cell
    ##scale data to equal number of UMIs per cell, then gene expression to mean=0/std=1
    pop.normalized_matrix['nontargeting'] = normalize_to_control(pop, control_cells='guide_target == "non-targeting"')
    pop.normalized_matrix['near_neighbors'] = normalize_to_near_neighbor_control(pop, 
                                                                                control_cells='guide_target == "non-targeting"',
                                                                                normalized_name='nontargeting',
                                                                                n_near_neighbors=50, n_components=50)
    """
    
    control_pop = pop.subpopulation(cells=control_cells, normalized_matrix='inherit', **kwargs)
    perturbed_pop = pop.subpopulation(cells='index not in @control',control=control_pop.cells.index, normalized_matrix='inherit', **kwargs)
    print('Normalizing {0} perturbed cells to {1} control cells...'.format(perturbed_pop.cells.shape[0],control_pop.cells.shape[0]))
    
    ##PCA of scaled data
    transformer = PCAReducer(n_components=n_components)
    control_reduced = control_pop.fit_transform(transformer,
                                                normalized=True,
                                                normalized_name=normalized_name)
    perturbed_reduced = transformer.transform(perturbed_pop.normalized_matrix[normalized_name])
    print('Completed PCA.')
    
    tree = KDTree(control_reduced)
    ##remove self from near neighbor for control cells
    distances, control_near_neighbor_index = tree.query(control_reduced, k=n_near_neighbors+1, workers=nthreads_pca)
    control_near_neighbor_index = pd.DataFrame(control_near_neighbor_index,index=control_pop.cells.index).drop(0,axis=1)
    distances, perturbed_near_neighbor_index = tree.query(perturbed_reduced, k=n_near_neighbors, workers=nthreads_pca)
    perturbed_near_neighbor_index = pd.DataFrame(perturbed_near_neighbor_index,index=perturbed_pop.cells.index)
    control_near_neighbor_index.columns = perturbed_near_neighbor_index.columns
    near_neighbor_index = pd.concat([control_near_neighbor_index,perturbed_near_neighbor_index])
    print('Completed near neighbor search.')

    ##parallel apply to normalize
    control_matrix=control_pop.normalized_matrix[normalized_name]
    pandarallel.initialize(nb_workers=nthreads_apply)
    near_neighbor_normalized = pop.normalized_matrix[normalized_name].parallel_apply(normalize_cell_to_near_neighbor, args=(control_matrix,near_neighbor_index,), axis=1)
    
    return near_neighbor_normalized

def normalize_cell_to_near_neighbor(cell,control_matrix,near_neighbor_index):
    avg_nt = control_matrix.iloc[near_neighbor_index.loc[cell.name]].mean()
    cell_corrected = cell-avg_nt
    return cell_corrected
            