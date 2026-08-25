import numpy as np
import pandas as pd

from scripts.summarize_shapemix_gse129785 import (
    collapse_prediction,
    paired_nominal_effects,
)


CELL_TYPES = [
    "Dendritic cells",
    "Monocytes",
    "B cells",
    "Regulatory T cells",
    "Naive CD4 T cells",
    "Memory CD4 T cells",
    "NK cells",
    "Naive CD8 T cells",
    "Memory CD8 T cells",
]


def test_nominal_collapses_preserve_off_target_mass():
    row = pd.Series(
        [0.05, 0.20, 0.05, 0.10, 0.10, 0.15, 0.05, 0.20, 0.10],
        index=CELL_TYPES,
    )
    names, cd4 = collapse_prediction(row, "cd4_memory_cd8_naive")
    assert names == ("CD4 Memory", "CD8 Naive")
    np.testing.assert_allclose(cd4, [0.15, 0.20, 0.65])

    names, mono_t = collapse_prediction(row, "monocyte_t")
    assert names == ("Monocytes", "T cells")
    np.testing.assert_allclose(mono_t, [0.20, 0.65, 0.15])


def test_paired_nominal_effects_are_length_minus_count_only():
    scores = pd.DataFrame(
        [
            {
                "dataset_id": "dataset",
                "sample": "sample",
                "family": "family",
                "method_id": "shapemix_length",
                "nominal_rmse_v1_descriptive": 0.10,
                "nominal_jsd_v2_descriptive": 0.20,
                "predicted_off_target_mass": 0.30,
                "rare_absolute_error": 0.04,
            },
            {
                "dataset_id": "dataset",
                "sample": "sample",
                "family": "family",
                "method_id": "shapemix_count_only",
                "nominal_rmse_v1_descriptive": 0.15,
                "nominal_jsd_v2_descriptive": 0.25,
                "predicted_off_target_mass": 0.20,
                "rare_absolute_error": 0.01,
            },
        ]
    )
    paired = paired_nominal_effects(scores)
    assert np.isclose(
        paired.loc[0, "nominal_rmse_v1_descriptive_length_minus_count_only"], -0.05
    )
    assert np.isclose(
        paired.loc[0, "predicted_off_target_mass_length_minus_count_only"], 0.10
    )
