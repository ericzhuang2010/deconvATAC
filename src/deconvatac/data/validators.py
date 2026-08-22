from __future__ import annotations

import numpy as np
import pandas as pd

from .schemas import DeconvolutionInput


def normalize_proportions(values: pd.DataFrame, zero_policy: str = "zeros") -> pd.DataFrame:
    """Row-normalize non-negative abundance-like values into proportions."""
    if zero_policy not in {"zeros", "uniform"}:
        raise ValueError("zero_policy must be one of {'zeros', 'uniform'}.")

    numeric = values.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    numeric = numeric.clip(lower=0)
    row_sums = numeric.sum(axis=1)

    proportions = numeric.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)
    if zero_policy == "uniform":
        zero_rows = row_sums == 0
        if zero_rows.any() and proportions.shape[1] > 0:
            proportions.loc[zero_rows, :] = 1.0 / proportions.shape[1]

    return proportions


def validate_deconvolution_input(data: DeconvolutionInput) -> None:
    """Validate the shared input contract for deconvolution methods."""
    if data.modality not in {"atac", "rna"}:
        raise ValueError("modality must be 'atac' or 'rna'.")
    if data.reference.n_vars == 0:
        raise ValueError("reference contains no features.")
    if data.spatial.n_vars == 0:
        raise ValueError("spatial contains no features.")
    if list(data.reference.var_names) != list(data.spatial.var_names):
        raise ValueError("reference and spatial features must be aligned in the same order.")
    if data.labels_key not in data.reference.obs:
        raise ValueError(f"labels_key '{data.labels_key}' is missing from reference.obs.")
    if data.spatial_key not in data.spatial.obsm:
        raise ValueError(f"spatial coordinates '{data.spatial_key}' are missing from spatial.obsm.")

    if data.truth is not None:
        missing = data.spatial.obs_names.difference(data.truth.index)
        if len(missing) > 0:
            raise ValueError("truth is missing rows for one or more spatial observations.")
