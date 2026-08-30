import pandas as pd

from scripts.summarize_shapemix_full_evaluation import (
    diagnostic_effects,
    nominal_effects,
    normalize_resources,
)


SOURCE = {
    "stage": "E",
    "family": "family",
    "evidence_class": "class",
    "analysis_unit": "unit",
}


def test_nominal_effects_remain_separate_by_family():
    table = pd.DataFrame(
        {
            "family": ["a", "a", "b"],
            "nominal_rmse_v1_descriptive_length_minus_count_only": [1.0, 3.0, -2.0],
            "nominal_jsd_v2_descriptive_length_minus_count_only": [2.0, 4.0, -3.0],
            "predicted_off_target_mass_length_minus_count_only": [3.0, 5.0, -4.0],
            "rare_absolute_error_length_minus_count_only": [4.0, 6.0, -5.0],
        }
    )
    rows = nominal_effects(SOURCE, table)
    rmse = [row for row in rows if row["endpoint"] == "nominal_rmse_v1_descriptive"]
    assert [(row["context"], row["estimate"], row["n_units"]) for row in rmse] == [
        ("a", 2.0, 2),
        ("b", -2.0, 1),
    ]


def test_diagnostic_effects_keep_factor_levels_separate():
    table = pd.DataFrame(
        {
            "factor": ["depth", "depth", "depth"],
            "level": ["low", "low", "high"],
            "metric": ["rmse_v1", "rmse_v1", "rmse_v1"],
            "length_minus_count_only": [1.0, 3.0, -1.0],
        }
    )
    rows = diagnostic_effects(SOURCE, table)
    assert [(row["context"], row["estimate"], row["n_units"]) for row in rows] == [
        ("depth=high", -1.0, 1),
        ("depth=low", 2.0, 2),
    ]


def test_resource_normalization_accepts_legacy_column_names():
    table = pd.DataFrame(
        {
            "dataset_id": ["d"],
            "method_run_id": ["shapemix_length"],
            "runtime_seconds": [7.5],
            "peak_memory_bytes": [1024],
            "status": ["success"],
        }
    )
    result = normalize_resources(SOURCE, table)
    assert result.loc[0, "method_id"] == "shapemix_length"
    assert result.loc[0, "wall_runtime_seconds"] == 7.5
    assert result.loc[0, "peak_rss_bytes"] == 1024
