from __future__ import annotations  # Helps with forward references
from typing import (
    Sequence, Generator, Any, Optional, Union, 
    Literal, Tuple
)
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from anndata import AnnData
from scanpy import logging as logg
from scanpy._utils import (
    _check_use_raw,
    AnyRandom,
    is_backed_type,
)
from scanpy.get import _get_obs_rep
from numpy.typing import DTypeLike
import numpy.typing as npt

def score_genes_relative_to_control(
    adata: AnnData,
    gene_list: Union[Sequence[str], pd.Index],
    *,
    cell_pool: Optional[Sequence[int]] = None,
    ctrl_as_ref: bool = True,
    ctrl_size: int = 50,
    gene_pool: Union[Sequence[str], pd.Index, None] = None,
    n_bins: int = 25,
    score_name: str = "score",
    random_state: AnyRandom = 0,
    copy: bool = False,
    use_raw: Optional[bool] = None,
) -> Optional[AnnData]:
    """\
    Score a set of genes relative to a control cell pool.

    The score is the average expression of a set of genes subtracted by the
    average expression of a reference set of genes computed within a specified
    control cell pool. This allows comparison against control cells in perturbation
    experiments.

    Parameters
    ----------
    adata
        The annotated data matrix.
    gene_list
        The list of gene names used for score calculation.
    cell_pool
        Indices or boolean mask specifying the control cells to use for background
        gene expression calculation. If None, all cells are used.
    ctrl_as_ref
        Allow the algorithm to use the control genes as reference.
    ctrl_size
        Number of reference genes to be sampled from each bin.
    gene_pool
        Genes for sampling the reference set. Default is all genes.
    n_bins
        Number of expression level bins for sampling.
    score_name
        Name of the field to be added in .obs.
    random_state
        The random seed for sampling.
    copy
        Copy adata or modify it inplace.
    use_raw
        Whether to use raw attribute of adata. Defaults to True if .raw is present.

    Returns
    -------
    Returns None if copy=False, else returns an AnnData object. Sets the following field:

    adata.obs[score_name] : :class:`numpy.ndarray` (dtype float)
        Scores of each cell relative to control cells.
    """
    start = logg.info(f"computing score {score_name!r} relative to control cells")
    adata = adata.copy() if copy else adata
    use_raw = _check_use_raw(adata, use_raw)
    if is_backed_type(adata.X) and not use_raw:
        raise NotImplementedError(
            f"score_genes_relative_to_control is not implemented for matrices of type {type(adata.X)}"
        )

    if random_state is not None:
        np.random.seed(random_state)

    # Changed helper name
    gene_list, gene_pool, get_subset = _check_relative_score_genes_args(
        adata, gene_list, gene_pool, use_raw=use_raw
    )
    del use_raw, random_state

    control_genes = pd.Index([], dtype="string")
    # Changed helper name
    for r_genes in _relative_score_genes_bins(
        gene_list,
        gene_pool,
        cell_pool=cell_pool,
        ctrl_as_ref=ctrl_as_ref,
        ctrl_size=ctrl_size,
        n_bins=n_bins,
        get_subset=get_subset,
    ):
        control_genes = control_genes.union(r_genes)

    if len(control_genes) == 0:
        msg = "No control genes found in any cut."
        if ctrl_as_ref:
            msg += " Try setting ctrl_as_ref=False."
        raise RuntimeError(msg)

    # Calculate means_list as average of gene_list in each cell (all cells)
    means_list = _relative_nan_means(get_subset(gene_list), axis=1)

    # Calculate means_control as average of control_genes in cell_pool cells
    control_data = get_subset(control_genes)
    if cell_pool is not None:
        control_data = control_data[cell_pool, :]
    # Compute mean per gene in control cells, then take the average of those
    means_control_per_gene = _relative_nan_means(control_data, axis=0)
    means_control = np.nanmean(means_control_per_gene)

    score = means_list - means_control

    adata.obs[score_name] = pd.Series(
        np.array(score).ravel(), index=adata.obs_names, dtype="float64"
    )

    logg.info(
        "    finished",
        time=start,
        deep=(
            "added\n"
            f"    {score_name!r}, score of gene set relative to control cells (adata.obs).\n"
            f"    {len(control_genes)} total control genes are used."
        ),
    )
    return adata if copy else None


def _check_relative_score_genes_args(  # Renamed
    adata: AnnData,
    gene_list: Union[pd.Index, Sequence[str]],
    gene_pool: Union[pd.Index, Sequence[str], None],
    *,
    use_raw: bool,
) -> Tuple[pd.Index, pd.Index, Any]:
    """Helper renamed to avoid conflicts"""
    var_names = adata.raw.var_names if use_raw else adata.var_names
    gene_list = pd.Index([gene_list] if isinstance(gene_list, str) else gene_list)
    genes_to_ignore = gene_list.difference(var_names, sort=False)
    gene_list = gene_list.intersection(var_names)
    if len(genes_to_ignore) > 0:
        logg.warning(f"genes are not in var_names and ignored: {genes_to_ignore}")
    if len(gene_list) == 0:
        raise ValueError("No valid genes were passed for scoring.")

    if gene_pool is None:
        gene_pool = var_names.astype("string")
    else:
        gene_pool = pd.Index(gene_pool, dtype="string").intersection(var_names)
    if len(gene_pool) == 0:
        raise ValueError("No valid genes were passed for reference set.")

    def get_subset(genes: pd.Index):
        x = _get_obs_rep(adata, use_raw=use_raw)
        if len(genes) == len(var_names):
            return x
        idx = var_names.get_indexer(genes)
        return x[:, idx]

    return gene_list, gene_pool, get_subset


def _relative_score_genes_bins(  # Renamed
    gene_list: pd.Index,
    gene_pool: pd.Index,
    *,
    cell_pool: Optional[Sequence[int]] = None,
    ctrl_as_ref: bool,
    ctrl_size: int,
    n_bins: int,
    get_subset: Any,
) -> Generator[pd.Index, None, None]:
    # Get the data for gene_pool genes
    data = get_subset(gene_pool)
    if cell_pool is not None:
        data = data[cell_pool, :]
    # Changed helper name
    obs_avg = pd.Series(_relative_nan_means(data, axis=0), index=gene_pool)
    obs_avg = obs_avg[np.isfinite(obs_avg)]

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    keep_ctrl_in_obs_cut = False if ctrl_as_ref else obs_cut.index.isin(gene_list)

    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes: pd.Index = obs_cut[(obs_cut == cut) & ~keep_ctrl_in_obs_cut].index
        if len(r_genes) == 0:
            logg.warning(f"No control genes for cut={cut}")
            continue
        if ctrl_size < len(r_genes):
            r_genes = r_genes.to_series().sample(ctrl_size, random_state=None).index
        if ctrl_as_ref:
            r_genes = r_genes.difference(gene_list)
        yield r_genes


def _relative_nan_means(  # Renamed
    x, *, axis: Literal[0, 1], dtype: Optional[DTypeLike] = None
) -> npt.NDArray[np.float64]:
    if issparse(x):
        return np.array(_relative_sparse_nanmean(x, axis=axis)).flatten()
    return np.nanmean(x, axis=axis, dtype=dtype)


def _relative_sparse_nanmean(  # Renamed
    X, axis: Literal[0, 1]
) -> npt.NDArray[np.float64]:
    """Renamed to avoid conflicts"""
    if axis not in (0, 1):
        raise ValueError("Axis must be 0 or 1")
    if issparse(X):
        X = X.copy()
        X.data = np.where(np.isnan(X.data), 0, X.data)
        counts = X.getnnz(axis=axis)
        sums = X.sum(axis=axis)
        means = np.array(sums / counts).flatten()
        means[counts == 0] = np.nan
        return means
    else:
        return np.nanmean(X, axis=axis)