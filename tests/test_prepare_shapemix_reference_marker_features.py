from __future__ import annotations

import numpy as np
import pytest

from scripts.prepare_shapemix_reference_marker_features import (
    marker_support_minimum,
    rank_marker_indices,
    validate_marker_reference_feature_count,
    yaml_builtin,
)


def test_reference_marker_ranker_is_reference_only_and_identifier_deterministic() -> None:
    mean_type = np.asarray([4.0, 2.0, 2.0, 1.0])
    mean_rest = np.asarray([1.0, 1.0, 1.0, 2.0])
    coverage = np.asarray([10, 10, 10, 100])
    totals = np.asarray([40.0, 20.0, 20.0, 100.0])
    names = ["peak-a", "peak-c", "peak-b", "peak-d"]

    assert rank_marker_indices(
        mean_type,
        mean_rest,
        coverage,
        totals,
        names,
        n_markers=2,
    ) == [0, 2]


def test_reference_marker_ranker_fails_when_specific_support_is_insufficient() -> None:
    with pytest.raises(ValueError, match="marker-support gate"):
        rank_marker_indices(
            np.asarray([2.0, 0.5]),
            np.asarray([1.0, 1.0]),
            np.asarray([9, 100]),
            np.asarray([2.0, 2.0]),
            ["a", "b"],
            n_markers=1,
        )


def test_marker_support_minimum_adapts_only_for_small_reference_types() -> None:
    assert marker_support_minimum(22) == 3
    assert marker_support_minimum(99) == 10
    assert marker_support_minimum(1) == 1


def test_yaml_builtin_normalizes_nested_numpy_scalars() -> None:
    value = {"version": np.int64(1), "nested": [np.bool_(True), np.float32(2.5)]}
    assert yaml_builtin(value) == {"version": 1, "nested": [True, 2.5]}
    assert type(yaml_builtin(value)["version"]) is int


def test_reference_marker_ranker_accepts_recorded_adaptive_support() -> None:
    assert rank_marker_indices(
        np.asarray([2.0, 1.5, 0.5]),
        np.asarray([1.0, 1.0, 1.0]),
        np.asarray([3, 2, 20]),
        np.asarray([6.0, 3.0, 10.0]),
        ["a", "b", "c"],
        n_markers=1,
        minimum_nonzero_cells=3,
    ) == [0]


def test_marker_reference_accepts_audited_positive_support_subset() -> None:
    manifest = {
        "counts": {"peaks": 4996},
        "feature_support_gate": {
            "source_selected_intervals": 5000,
            "minimum_reconstructed_reference_total": 1,
            "retained_intervals": 4996,
            "excluded_zero_total_intervals": ["p1", "p2", "p3", "p4"],
            "outcome_data_used": False,
        },
    }
    validate_marker_reference_feature_count(manifest, 4996)
    manifest["feature_support_gate"]["retained_intervals"] = 4995
    with pytest.raises(ValueError, match="feature_support_gate is inconsistent"):
        validate_marker_reference_feature_count(manifest, 4996)
