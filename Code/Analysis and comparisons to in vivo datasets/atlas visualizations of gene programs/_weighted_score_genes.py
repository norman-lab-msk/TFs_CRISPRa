"""Calculate scores based on the expression of gene lists."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.sparse import issparse

from scanpy._utils import _check_use_raw

from scanpy import logging as logg
from scanpy.get import _get_obs_rep

import scipy

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from anndata import AnnData
    from numpy.typing import DTypeLike, NDArray
    from scipy.sparse import csc_matrix, csr_matrix

    from .._utils import AnyRandom


def _sparse_nanmean(
    X: csr_matrix | csc_matrix, axis: Literal[0, 1]
) -> NDArray[np.float64]:
    """
    np.nanmean equivalent for sparse matrices
    """
    if not issparse(X):
        raise TypeError("X must be a sparse matrix")

    # count the number of nan elements per row/column (dep. on axis)
    Z = X.copy()
    Z.data = np.isnan(Z.data)
    Z.eliminate_zeros()
    n_elements = Z.shape[axis] - Z.sum(axis)

    # set the nans to 0, so that a normal .sum() works
    Y = X.copy()
    Y.data[np.isnan(Y.data)] = 0
    Y.eliminate_zeros()

    # the average
    s = Y.sum(axis, dtype="float64")  # float64 for score_genes function compatibility)
    m = s / n_elements

    return m

def _sparse_nansum(
    X: csr_matrix | csc_matrix, axis: Literal[0, 1]
) -> NDArray[np.float64]:
    """
    np.nanmean equivalent for sparse matrices
    """
    if not issparse(X):
        raise TypeError("X must be a sparse matrix")

    # count the number of nan elements per row/column (dep. on axis)
    Z = X.copy()
    Z.data = np.isnan(Z.data)
    Z.eliminate_zeros()
    n_elements = Z.shape[axis] - Z.sum(axis)

    # set the nans to 0, so that a normal .sum() works
    Y = X.copy()
    Y.data[np.isnan(Y.data)] = 0
    Y.eliminate_zeros()

    # the average
    s = Y.sum(axis, dtype="float64")  # float64 for score_genes function compatibility)

    return s


def weighted_score_genes(
    adata: AnnData,
    gene_list: Sequence[str] | pd.Index[str],
    weights: Sequence[float],
    *,
    ctrl_size: int = 50,
    gene_pool: Sequence[str] | pd.Index[str] | None = None,
    n_bins: int = 25,
    score_name: str = "score",
    random_state: AnyRandom = 0,
    copy: bool = False,
    use_raw: bool | None = None,
    layer: str | None = None,
) -> AnnData | None:
    """\
    Score a set of genes :cite:p:`Satija2015`.

    The score is the average expression of a set of genes subtracted with the
    average expression of a reference set of genes. The reference set is
    randomly sampled from the `gene_pool` for each binned expression value.

    This reproduces the approach in Seurat :cite:p:`Satija2015` and has been implemented
    for Scanpy by Davide Cittaro.

    Parameters
    ----------
    adata
        The annotated data matrix.
    gene_list
        The list of gene names used for score calculation.
    ctrl_size
        Number of reference genes to be sampled from each bin. If `len(gene_list)` is not too
        low, you can set `ctrl_size=len(gene_list)`.
    gene_pool
        Genes for sampling the reference set. Default is all genes.
    n_bins
        Number of expression level bins for sampling.
    score_name
        Name of the field to be added in `.obs`.
    random_state
        The random seed for sampling.
    copy
        Copy `adata` or modify it inplace.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.

        .. versionchanged:: 1.4.5
           Default value changed from `False` to `None`.
    layer
        Key from `adata.layers` whose value will be used to perform tests on.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following field:

    `adata.obs[score_name]` : :class:`numpy.ndarray` (dtype `float`)
        Scores of each cell.

    Examples
    --------
    See this `notebook <https://github.com/scverse/scanpy_usage/tree/master/180209_cell_cycle>`__.
    """
    start = logg.info(f"computing score {score_name!r}")
    adata = adata.copy() if copy else adata
    use_raw = _check_use_raw(adata, use_raw)

    if random_state is not None:
        np.random.seed(random_state)

    var_names = adata.raw.var_names if use_raw else adata.var_names
    gene_list = pd.Index([gene_list] if isinstance(gene_list, str) else gene_list)
    weights = pd.Series([weights] if isinstance(weights, float) else weights)
    genes_to_ignore = gene_list.difference(var_names, sort=False)  # first get missing
    gene_list = gene_list.intersection(var_names)  # then restrict to present
    weights = weights[[gene_list.get_loc(x) for x in gene_list.intersection(var_names)]]
    if len(genes_to_ignore) > 0:
        logg.warning(f"genes are not in var_names and ignored: {genes_to_ignore}")
    if len(gene_list) == 0:
        raise ValueError("No valid genes were passed for scoring.")

    if gene_pool is None:
        gene_pool = pd.Index(var_names, dtype="string")
    else:
        gene_pool = pd.Index(gene_pool, dtype="string").intersection(var_names)
    if len(gene_pool) == 0:
        raise ValueError("No valid genes were passed for reference set.")

    # Trying here to match the Seurat approach in scoring cells.
    # Basically we need to compare genes against random genes in a matched
    # interval of expression.

    def get_subset(genes: pd.Index[str]):
        x = _get_obs_rep(adata, use_raw=use_raw, layer=layer)
        if len(genes) == len(var_names):
            return x
        idx = var_names.get_indexer(genes)
        return x[:, idx]

    # average expression of genes
    obs_avg = pd.Series(_nan_means(get_subset(gene_pool), axis=0), index=gene_pool)
    # Sometimes (and I don’t know how) missing data may be there, with NaNs for missing entries
    obs_avg = obs_avg[np.isfinite(obs_avg)]

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    control_genes = []

    # now pick `ctrl_size` genes from every cut
    for cut in obs_cut.loc[gene_list]:
        r_genes = obs_cut[obs_cut == cut].index
        if ctrl_size < len(r_genes):
            r_genes = r_genes.to_series().sample(ctrl_size).index
        control_genes.append(r_genes.difference(gene_list))
    
    means_control = 0
    for i, ctrl_genes in enumerate(control_genes):
        means_control += _nan_means(get_subset(ctrl_genes), axis=1, dtype="float64") * weights[i]
        
    means_list = _nan_sums(get_subset(gene_list) * scipy.sparse.diags(weights, 0), axis=1, dtype="float64")

    score = (means_list - means_control) / (weights.mean()*len(gene_list))

    adata.obs[score_name] = pd.Series(
        np.array(score).ravel(), index=adata.obs_names, dtype="float64"
    )

    logg.info(
        "    finished",
        time=start,
        deep=(
            "added\n"
            f"    {score_name!r}, score of gene set (adata.obs).\n"
            f"    {len(control_genes)} total control genes are used."
        ),
    )
    return adata if copy else None


def _nan_means(
    x, *, axis: Literal[0, 1], dtype: DTypeLike | None = None
) -> NDArray[np.float64]:
    if issparse(x):
        return np.array(_sparse_nanmean(x, axis=axis)).flatten()
    return np.nanmean(x, axis=axis, dtype=dtype)

def _nan_sums(
    x, *, axis: Literal[0, 1], dtype: DTypeLike | None = None
) -> NDArray[np.float64]:
    if issparse(x):
        return np.array(_sparse_nansum(x, axis=axis)).flatten()
    return np.nansum(x, axis=axis, dtype=dtype)
