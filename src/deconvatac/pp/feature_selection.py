import numpy as np
from anndata import AnnData
from scipy import sparse


def _as_1d_array(values):
    """Convert sparse/dense matrix-like values to a one-dimensional ndarray."""
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def _top_indices(values, n_top_features: int):
    n_top = min(int(n_top_features), values.shape[0])
    if n_top <= 0:
        raise ValueError("n_top_features must be positive.")
    return np.argpartition(values, -n_top)[-n_top:]


def highly_variable_peaks(
    adata: AnnData, cluster_key: str, layer: str = None, scale: float = 1, n_top_features: int = 20000
):
    """
    Selects highly variable features the "var" way. 

    Adapted from: https://github.com/GreenleafLab/ArchR/blob/c61b0645d1482f80dcc24e25fbd915128c1b2500/R/IterativeLSI.R#L1015

    Parameters
    -----------

    adata: AnnData
        AnnData object of the (reference) scATAC data.
    cluster_key: str
        Name of column in adata.obs containing the clusters.
    layer: str 
        Layer of the raw counts. If None, uses .X
    scale: float
        Scale factor, i.e. log2((sums/feature_sums) *scale + 1).
    n_top_features: int
        How many features to select.

        
    Returns
    -------
    
    Saves a boolean indicator of the HVPs in place to adata.var['highly_variable'].
    """
    if cluster_key not in adata.obs:
        raise KeyError(f"cluster_key '{cluster_key}' is missing from adata.obs.")

    # step 1: get matrix of shape (clusters x features) (i.e. sum up features per cluster)
    clusters = adata.obs[cluster_key].unique()
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"layer '{layer}' is missing from adata.layers.")
        matrix = adata.layers[layer]
    else:
        matrix = adata.X

    # number of rows becomes the number of clusters
    group_matrix = [_as_1d_array(matrix[np.asarray(adata.obs[cluster_key] == clus)].sum(axis=0)) for clus in clusters]
    group_matrix = np.vstack(group_matrix).astype(float)

    # step 2: log normalize
    # divide each count by sum(all counts in a row/cluster)
    row_sums = group_matrix.sum(axis=1, keepdims=True)
    group_matrix = np.divide(group_matrix, row_sums, out=np.zeros_like(group_matrix), where=row_sums != 0)
    # scale & log normalize
    group_matrix = np.log2(group_matrix * scale + 1)

    # step 3:
    # compute variance of each feature, variance computed all clusters for each feature
    var = np.var(group_matrix, axis=0)

    # step 4
    # get indices of the top n features with the highest variance
    idx = _top_indices(var, n_top_features)

    # step 5: save HVF to anndata, var stores feature names
    hv_features = [adata.var.index[i] for i in idx]
    # var["highly_variable"] is a boolean vector of length n_features, 
    # True for highly variable features, False otherwise
    adata.var["highly_variable"] = adata.var.index.isin(hv_features)

    return


def highly_accessible_peaks(adata: AnnData, layer: str = None, n_top_features: int = 20000, copy: bool = False):
    """
    Selects the most accessible peaks from the given AnnData object.

    Parameters
    ----------

    adata: AnnData
        Annotated data object containing the peaks.
    layer: str 
        Name of the layer to use for peak accessibility.
    n_top_features: int
        Number of top accessible peaks to select.
    copy: bool
        Whether to copy the AnnData object. 

        
    Returns
    -------

    AnnData object with the highly accessible peaks saved to adata.var['highly_accessible'].
    """
    if copy:
        adata = adata.copy()

    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"layer '{layer}' is missing from adata.layers.")
        matrix = adata.layers[layer]
    else:
        matrix = adata.X

    # binarize the matrix
    binary_matrix = matrix > 0

    # get indices of the top n features with the highest sum
    idx = _top_indices(_as_1d_array(binary_matrix.sum(axis=0)), n_top_features)
    # step 5: save HVF to anndata
    hv_features = adata.var.index.values[idx]

    adata.var["highly_accessible"] = adata.var.index.isin(hv_features)
    if copy:
        return adata
    else:
        return
