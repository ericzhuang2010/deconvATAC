from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon

from deconvatac.data import normalize_proportions


def align_proportions(true: pd.DataFrame, predicted: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Align true and predicted proportion matrices by shared spots and cell types."""
    shared_index = true.index.intersection(predicted.index)
    shared_columns = true.columns.intersection(predicted.columns)
    if len(shared_index) == 0:
        raise ValueError("No shared observations between true and predicted proportions.")
    if len(shared_columns) == 0:
        raise ValueError("No shared cell-type columns between true and predicted proportions.")

    true_aligned = true.loc[shared_index, shared_columns]
    predicted_aligned = predicted.loc[shared_index, shared_columns]
    return normalize_proportions(true_aligned), normalize_proportions(predicted_aligned)


def rmse(true: pd.DataFrame, predicted: pd.DataFrame) -> float:
    """Compute RMSE after aligning spots and cell-type columns."""
    true_aligned, predicted_aligned = align_proportions(true, predicted)
    return float(np.sqrt(np.mean((true_aligned.values - predicted_aligned.values) ** 2)))


def jsd(true: pd.DataFrame, predicted: pd.DataFrame) -> float:
    """Compute mean Jensen-Shannon divergence after alignment."""
    true_aligned, predicted_aligned = align_proportions(true, predicted)
    values = jensenshannon(true_aligned.values, predicted_aligned.values, axis=1, base=2)
    return float(np.mean(values[np.isfinite(values)]))
