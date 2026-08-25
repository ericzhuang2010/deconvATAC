from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.prepare_shapemix_pbmc_sensitivity import _design_rows
from scripts.summarize_shapemix_pbmc_sensitivity import (
    METHOD_IDS,
    METRICS,
    control_differences,
    expanded_factor_scores,
    factor_level_scores,
    paired_shape_count,
)


def _synthetic_scores() -> pd.DataFrame:
    records = []
    offsets = {
        "shapemix_length": 0.0,
        "shapemix_count_only": 0.1,
        "nnls": 0.2,
    }
    for dataset_index, row in enumerate(_design_rows()):
        for method in METHOD_IDS:
            for metric_index, metric in enumerate(METRICS):
                records.append(
                    {
                        "dataset_id": row.dataset_id,
                        "factor": row.factor,
                        "level": row.level,
                        "control_level": row.control_level,
                        "mixture_seed": row.mixture_seed,
                        "method_run_id": method,
                        "metric": metric,
                        "value": dataset_index / 1000 + metric_index + offsets[method],
                    }
                )
    return pd.DataFrame.from_records(records)


def test_factor_expansion_reuses_anchor_as_five_controls() -> None:
    expanded = expanded_factor_scores(_synthetic_scores())
    assert expanded.groupby("factor")["level"].nunique().to_dict() == {
        "bins": 3,
        "cells": 4,
        "depth": 4,
        "features": 3,
        "rare_nk": 4,
        "reference_support": 4,
        "subtype": 2,
    }
    assert set(factor_level_scores(expanded)["mixture_seeds"]) == {2}


def test_control_differences_are_zero_at_every_control() -> None:
    differences = control_differences(expanded_factor_scores(_synthetic_scores()))
    controls = differences[differences["level"] == differences["control_level"]]
    assert len(controls) == 7 * 2 * len(METHOD_IDS) * len(METRICS)
    assert (controls["value_minus_control"] == 0).all()


def test_shape_count_effect_is_paired_per_dataset() -> None:
    paired = paired_shape_count(_synthetic_scores())
    assert len(paired) == 40 * len(METRICS)
    assert np.allclose(paired["length_minus_count_only"], -0.1)
    assert paired["length_improved"].all()
