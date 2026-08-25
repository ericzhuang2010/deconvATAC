import numpy as np
import pandas as pd

from scripts.summarize_shapemix_gse194122 import (
    donor_level_summary,
    paired_effects,
)


def test_paired_effects_use_length_minus_count_only():
    rows = []
    for method, value in (
        ("shapemix_length", 0.12),
        ("shapemix_count_only", 0.15),
    ):
        rows.append(
            {
                "dataset_id": "dataset",
                "donor": 1,
                "condition": "equal_celltype",
                "mixture_seed": 101,
                "metric": "rmse_v1",
                "method_run_id": method,
                "value": value,
            }
        )
    paired = paired_effects(pd.DataFrame(rows))
    assert np.isclose(paired.loc[0, "length_minus_count_only"], -0.03)
    assert bool(paired.loc[0, "length_improved"])


def test_donor_summary_uses_ten_donor_units():
    values = np.linspace(-0.05, 0.04, 10)
    donors = pd.DataFrame(
        {
            "donor": np.arange(1, 11),
            "condition": "all",
            "metric": "rmse_v1",
            "mean_length_minus_count_only": values,
            "inner_mixtures": 4,
            "inner_mixture_sd": 0.01,
        }
    )
    summary = donor_level_summary(donors)
    assert summary.loc[0, "donors"] == 10
    assert np.isclose(
        summary.loc[0, "mean_length_minus_count_only"],
        values.mean(),
    )
    assert summary.loc[0, "ci95_lower"] < values.mean()
    assert summary.loc[0, "ci95_upper"] > values.mean()
    assert summary.loc[0, "donors_length_improved"] == 5
