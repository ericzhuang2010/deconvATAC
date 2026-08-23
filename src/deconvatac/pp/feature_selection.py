from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from ..data.validators import ordered_feature_sha256


@dataclass(frozen=True)
class ReferencePeakSelectionResult:
    """Deterministic reference-only peak ranking and its audit values.

    ``indices`` are positions on the input ``adata.var`` axis, but all selected
    arrays and ``peak_ids`` are stored in ranked output order.  The candidate
    arrays retain the original ``adata.var`` order so the complete ranking
    inputs can be audited without recomputation.
    """

    peak_ids: tuple[str, ...]
    indices: np.ndarray
    scores: np.ndarray
    nonzero_reference_cells: np.ndarray
    total_reference_counts: np.ndarray
    candidate_peak_ids: tuple[str, ...]
    candidate_scores: np.ndarray
    candidate_nonzero_reference_cells: np.ndarray
    candidate_total_reference_counts: np.ndarray
    eligible_mask: np.ndarray
    cell_types: tuple[str, ...]
    candidate_feature_sha256: str
    selected_feature_sha256: str
    n_top_peaks: int
    min_reference_cells: int
    scale: float

    @property
    def candidate_universe_sha256(self) -> str:
        """Return the ordered candidate-universe hash."""
        return self.candidate_feature_sha256

    @property
    def selected_peak_sha256(self) -> str:
        """Return the ordered selected-peak hash."""
        return self.selected_feature_sha256

    def to_frame(self) -> pd.DataFrame:
        """Return the selected peaks and audit values in ranked output order."""
        return pd.DataFrame(
            {
                "rank": np.arange(1, len(self.peak_ids) + 1, dtype=np.int64),
                "peak_id": self.peak_ids,
                "original_index": self.indices,
                "score": self.scores,
                "nonzero_reference_cells": self.nonzero_reference_cells,
                "total_reference_count": self.total_reference_counts,
            }
        )

    def candidate_frame(self) -> pd.DataFrame:
        """Return audit values for every candidate in input feature order."""
        return pd.DataFrame(
            {
                "peak_id": self.candidate_peak_ids,
                "original_index": np.arange(len(self.candidate_peak_ids), dtype=np.int64),
                "score": self.candidate_scores,
                "nonzero_reference_cells": self.candidate_nonzero_reference_cells,
                "total_reference_count": self.candidate_total_reference_counts,
                "eligible": self.eligible_mask,
            }
        )


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


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer.")
    return int(value)


def _validated_reference_count_matrix(adata: AnnData):
    matrix = adata.X
    if matrix.shape != adata.shape:
        raise ValueError("adata.X must have the same shape as the AnnData object.")

    if sparse.issparse(matrix):
        if not np.issubdtype(matrix.dtype, np.number) or np.issubdtype(matrix.dtype, np.complexfloating):
            raise TypeError("adata.X must contain real numeric cut-site counts.")
        if matrix.data.size:
            values = np.asarray(matrix.data)
            if not np.all(np.isfinite(values)):
                raise ValueError("adata.X contains non-finite cut-site counts.")
            if np.any(values < 0):
                raise ValueError("adata.X contains negative cut-site counts.")
            if np.any(values != np.floor(values)):
                raise ValueError("adata.X contains non-integer cut-site counts.")
        # The PBMC candidate matrix is large.  Reuse an already canonical CSR
        # matrix instead of making a full float64 copy; aggregation below still
        # requests float64 accumulators explicitly.
        if sparse.isspmatrix_csr(matrix) and matrix.has_canonical_format:
            return matrix
        validated = sparse.csr_matrix(matrix, copy=True)
        validated.sum_duplicates()
        return validated

    values = np.asarray(matrix)
    if values.ndim != 2 or values.shape != adata.shape:
        raise ValueError("adata.X must be a two-dimensional cell-by-peak matrix.")
    if not np.issubdtype(values.dtype, np.number) or np.issubdtype(values.dtype, np.complexfloating):
        raise TypeError("adata.X must contain real numeric cut-site counts.")
    if not np.all(np.isfinite(values)):
        raise ValueError("adata.X contains non-finite cut-site counts.")
    if np.any(values < 0):
        raise ValueError("adata.X contains negative cut-site counts.")
    if np.any(values != np.floor(values)):
        raise ValueError("adata.X contains non-integer cut-site counts.")
    return values


def _sum_columns_float64(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return _as_1d_array(matrix.sum(axis=0, dtype=np.float64)).astype(np.float64, copy=False)
    return np.asarray(np.sum(matrix, axis=0, dtype=np.float64)).ravel()


def _canonical_cell_types(labels: Sequence[object], cell_types: Optional[Sequence[str]]) -> tuple[str, ...]:
    label_values = np.asarray(labels, dtype=object)
    if label_values.ndim != 1 or label_values.size == 0:
        raise ValueError("The training-reference object must contain at least one cell.")
    if pd.isna(label_values).any():
        raise ValueError("Training-reference cell-type labels must not be missing.")
    if any(not isinstance(value, str) or not value for value in label_values):
        raise TypeError("Training-reference cell-type labels must be non-empty strings.")

    observed = set(label_values.tolist())
    if cell_types is None:
        requested = observed
    else:
        supplied = tuple(cell_types)
        if any(not isinstance(value, str) or not value for value in supplied):
            raise TypeError("cell_types must contain non-empty strings.")
        if len(set(supplied)) != len(supplied):
            raise ValueError("cell_types must not contain duplicates.")
        requested = set(supplied)
        if requested != observed:
            missing = sorted(requested - observed, key=lambda value: value.encode("utf-8"))
            unexpected = sorted(observed - requested, key=lambda value: value.encode("utf-8"))
            raise ValueError(
                "cell_types must exactly match the training-reference label universe; "
                f"missing from observations={missing}, unexpected observations={unexpected}."
            )

    # The score is invariant to type order.  Canonicalizing it also makes the
    # floating-point reduction order identical when callers supply a different
    # ordering of the same type universe.
    return tuple(sorted(requested, key=lambda value: value.encode("utf-8")))


def select_reference_peaks(
    adata: AnnData,
    cell_type_key: str,
    *,
    cell_types: Optional[Sequence[str]] = None,
    n_top_peaks: int = 5000,
    min_reference_cells: int = 10,
    scale: float = 1e4,
) -> ReferencePeakSelectionResult:
    """Rank peaks using only collapsed training-reference counts and labels.

    The implementation follows ShapeMix benchmark protocol version 1.  For
    each type it computes ``log2(1 + scale * K[c,p] / sum_p K[c,p])``, scores
    peaks with population variance across types, filters on nonzero-cell
    coverage, and applies the complete deterministic tie key.

    Parameters
    ----------
    adata
        Training-reference cells.  Only ``adata.X``, ``adata.obs`` at
        ``cell_type_key``, and ``adata.var_names`` are read.
    cell_type_key
        Observation column containing the cell-type labels.
    cell_types
        Optional declared universe.  When supplied, it must exactly match the
        observed training-reference labels.  Its order cannot affect ranking.
    n_top_peaks
        Exact number of ranked peaks to return.  The function fails if fewer
        eligible peaks are available.
    min_reference_cells
        Minimum number of training-reference cells with a nonzero count.
    scale
        Library normalization scale, fixed to ``10^4`` in protocol version 1.

    Returns
    -------
    ReferencePeakSelectionResult
        Selected identifiers and original indices in ranked order, together
        with selected and candidate-level audit values and ordered hashes.
    """
    n_top_peaks = _positive_integer(n_top_peaks, "n_top_peaks")
    min_reference_cells = _positive_integer(min_reference_cells, "min_reference_cells")
    if isinstance(scale, (bool, np.bool_)) or not np.isscalar(scale):
        raise ValueError("scale must be a positive finite number.")
    scale = float(scale)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be a positive finite number.")
    if not isinstance(cell_type_key, str) or not cell_type_key:
        raise TypeError("cell_type_key must be a non-empty string.")
    if cell_type_key not in adata.obs:
        raise KeyError(f"cell_type_key '{cell_type_key}' is missing from adata.obs.")
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("The training-reference object must contain cells and candidate peaks.")
    if not adata.var_names.is_unique:
        raise ValueError("Candidate peak IDs in adata.var_names must be unique.")

    candidate_peak_ids = tuple(adata.var_names.tolist())
    if any(not isinstance(peak_id, str) or not peak_id for peak_id in candidate_peak_ids):
        raise TypeError("Candidate peak IDs must be non-empty strings.")

    labels = np.asarray(adata.obs[cell_type_key], dtype=object)
    canonical_types = _canonical_cell_types(labels, cell_types)
    matrix = _validated_reference_count_matrix(adata)

    group_counts = np.empty((len(canonical_types), adata.n_vars), dtype=np.float64)
    for type_index, cell_type in enumerate(canonical_types):
        type_mask = labels == cell_type
        group_counts[type_index] = _sum_columns_float64(matrix[type_mask])

    type_totals = group_counts.sum(axis=1)
    zero_types = [canonical_types[index] for index in np.flatnonzero(type_totals == 0)]
    if zero_types:
        raise ValueError(f"Every cell type must have a positive total training count; zero totals: {zero_types}.")
    if not np.all(np.isfinite(type_totals)):
        raise ValueError("Cell-type total training counts are non-finite.")

    log_normalized = np.log2(1.0 + scale * group_counts / type_totals[:, np.newaxis])
    scores = np.var(log_normalized, axis=0, ddof=0)
    coverage = _as_1d_array((matrix > 0).sum(axis=0)).astype(np.int64, copy=False)
    total_counts_float = group_counts.sum(axis=0)
    if np.any(total_counts_float > np.iinfo(np.int64).max):
        raise OverflowError("Peak total training counts exceed signed 64-bit integer range.")
    total_counts = total_counts_float.astype(np.int64)
    eligible_mask = coverage >= min_reference_cells
    eligible_indices = np.flatnonzero(eligible_mask)
    if eligible_indices.size < n_top_peaks:
        raise ValueError(
            f"Only {eligible_indices.size} peaks are nonzero in at least {min_reference_cells} "
            f"training-reference cells; {n_top_peaks} are required."
        )

    ranked_indices = np.asarray(
        sorted(
            eligible_indices.tolist(),
            key=lambda index: (
                -scores[index],
                -coverage[index],
                -total_counts[index],
                candidate_peak_ids[index].encode("utf-8"),
            ),
        )[:n_top_peaks],
        dtype=np.int64,
    )
    selected_peak_ids = tuple(candidate_peak_ids[index] for index in ranked_indices)

    return ReferencePeakSelectionResult(
        peak_ids=selected_peak_ids,
        indices=ranked_indices,
        scores=scores[ranked_indices],
        nonzero_reference_cells=coverage[ranked_indices],
        total_reference_counts=total_counts[ranked_indices],
        candidate_peak_ids=candidate_peak_ids,
        candidate_scores=scores,
        candidate_nonzero_reference_cells=coverage,
        candidate_total_reference_counts=total_counts,
        eligible_mask=eligible_mask,
        cell_types=canonical_types,
        candidate_feature_sha256=ordered_feature_sha256(candidate_peak_ids),
        selected_feature_sha256=ordered_feature_sha256(selected_peak_ids),
        n_top_peaks=n_top_peaks,
        min_reference_cells=min_reference_cells,
        scale=scale,
    )


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
