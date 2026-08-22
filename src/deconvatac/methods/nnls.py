from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import nnls

from deconvatac.data import DeconvolutionInput, DeconvolutionResult, normalize_proportions

from .base import BaseDeconvolver


def _matrix_from_layer(adata, layer: Optional[str]):
    matrix = adata.layers[layer] if layer is not None else adata.X
    if sparse.issparse(matrix):
        return matrix.tocsr().astype(float)
    return np.asarray(matrix, dtype=float)


def _row_to_array(matrix, idx: int) -> np.ndarray:
    row = matrix[idx]
    if sparse.issparse(row):
        return row.toarray().ravel()
    return np.asarray(row).ravel()


def _as_array(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values, dtype=float)


def _cell_type_signatures(reference_x, labels: pd.Series, cell_types: list[str]) -> np.ndarray:
    signatures = []
    for cell_type in cell_types:
        mask = labels == cell_type
        mean = reference_x[mask.values].mean(axis=0)
        signatures.append(_as_array(mean).ravel())
    return np.vstack(signatures)


class NNLSDeconvolver(BaseDeconvolver):
    """Plain non-negative least-squares baseline."""

    method_name = "nnls"

    def run(self, data: DeconvolutionInput) -> DeconvolutionResult:
        layer_ref = self.config.get("layer_ref")
        layer_spatial = self.config.get("layer_spatial")
        zero_policy = self.config.get("zero_policy", "zeros")

        reference_x = _matrix_from_layer(data.reference, layer_ref)
        spatial_x = _matrix_from_layer(data.spatial, layer_spatial)
        labels = data.reference.obs[data.labels_key].astype(str)
        cell_types = list(pd.Index(labels).drop_duplicates())
        signatures = _cell_type_signatures(reference_x, labels, cell_types)

        design = signatures.T
        weights = []
        residuals = []
        for idx in range(data.spatial.n_obs):
            y = _row_to_array(spatial_x, idx)
            coef, residual = nnls(design, y)
            weights.append(coef)
            residuals.append(float(residual))

        abundance = pd.DataFrame(weights, index=data.spatial.obs_names, columns=cell_types)
        proportions = normalize_proportions(abundance, zero_policy=zero_policy)
        diagnostics = {
            "method": self.method_name,
            "model": "non_negative_least_squares",
            "params": {
                "layer_ref": layer_ref,
                "layer_spatial": layer_spatial,
                "zero_policy": zero_policy,
            },
            "n_features": data.reference.n_vars,
            "n_cell_types": len(cell_types),
            "mean_residual": float(np.mean(residuals)) if residuals else None,
        }

        return DeconvolutionResult(
            method=self.method_name,
            dataset_id=data.dataset_id,
            modality=data.modality,
            feature_set=data.feature_set,
            proportions=proportions,
            abundance=abundance,
            diagnostics=diagnostics,
        )
