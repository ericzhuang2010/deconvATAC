from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import jensenshannon

from deconvatac.metrics import (
    PROPORTION_CONTRACT_VERSION,
    PROPORTION_METRICS,
    align_proportions,
    declared_cell_types_from_metadata,
    evaluate_proportion_metric,
    js_distance_v1,
    jsd_v2,
    resolve_proportion_metric,
    rmse_v1,
)


CELL_TYPES = ["A", "B", "C"]


@pytest.fixture
def truth() -> pd.DataFrame:
    return pd.DataFrame(
        [[0.6, 0.3, 0.1], [0.1, 0.7, 0.2]],
        index=["s1", "s2"],
        columns=CELL_TYPES,
    )


@pytest.fixture
def prediction() -> pd.DataFrame:
    return pd.DataFrame(
        [[0.2, 0.5, 0.3], [0.5, 0.4, 0.1]],
        index=["s2", "s1"],
        columns=["C", "A", "B"],
    )


def test_alignment_reorders_exact_spot_set_and_declared_prediction_columns(truth, prediction):
    true_aligned, predicted_aligned = align_proportions(truth, prediction, CELL_TYPES)

    assert true_aligned.index.tolist() == ["s1", "s2"]
    assert predicted_aligned.index.tolist() == ["s1", "s2"]
    assert predicted_aligned.columns.tolist() == CELL_TYPES
    np.testing.assert_allclose(predicted_aligned.loc["s1"], [0.4, 0.1, 0.5])


def test_omitted_declared_prediction_type_is_zero_filled_and_penalized(truth):
    predicted = pd.DataFrame(
        [[0.7, 0.3], [0.2, 0.8]],
        index=truth.index,
        columns=["A", "B"],
    )
    _, aligned = align_proportions(truth, predicted, CELL_TYPES)

    assert aligned["C"].eq(0).all()
    expected = np.sqrt(np.mean((truth.to_numpy() - aligned.to_numpy()) ** 2))
    assert rmse_v1(truth, predicted, CELL_TYPES) == pytest.approx(expected)
    assert expected > np.sqrt(np.mean((truth[["A", "B"]].to_numpy() - predicted.to_numpy()) ** 2))


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda frame: frame.rename(columns={"C": "unknown"}), "outside the declared universe"),
        (lambda frame: frame.drop(index="s1"), "exactly the truth spot set"),
        (
            lambda frame: pd.concat(
                [frame, pd.DataFrame([[0.2, 0.3, 0.5]], index=["extra"], columns=frame.columns)]
            ),
            "exactly the truth spot set",
        ),
    ],
)
def test_unknown_types_or_changed_spot_set_fail(truth, mutation, message):
    with pytest.raises(ValueError, match=message):
        align_proportions(truth, mutation(truth.copy()), CELL_TYPES)


@pytest.mark.parametrize("target", ["truth_index", "prediction_index"])
def test_duplicate_spots_fail(truth, target):
    predicted = truth.copy()
    if target == "truth_index":
        truth.index = ["s1", "s1"]
    else:
        predicted.index = ["s1", "s1"]
    with pytest.raises(ValueError, match="spot names must be unique"):
        align_proportions(truth, predicted, CELL_TYPES)


@pytest.mark.parametrize("target", ["truth_columns", "prediction_columns"])
def test_duplicate_columns_fail(truth, target):
    predicted = truth.copy()
    if target == "truth_columns":
        truth.columns = ["A", "A", "C"]
    else:
        predicted.columns = ["A", "A", "C"]
    with pytest.raises(ValueError, match="cell-type columns must be unique"):
        align_proportions(truth, predicted, CELL_TYPES)


def test_duplicate_declared_type_fails(truth):
    with pytest.raises(ValueError, match="cell_types must not contain duplicates"):
        align_proportions(truth, truth, ["A", "A", "C"])


def test_truth_columns_must_equal_declared_order(truth):
    with pytest.raises(ValueError, match="exactly match the declared cell_types order"):
        align_proportions(truth[["B", "A", "C"]], truth, CELL_TYPES)


@pytest.mark.parametrize(
    ("value", "message"),
    [
        (np.nan, "finite"),
        (np.inf, "finite"),
        (-0.1, "nonnegative"),
    ],
)
@pytest.mark.parametrize("target", ["truth", "prediction"])
def test_nan_infinity_and_negative_values_fail(truth, value, message, target):
    predicted = truth.copy()
    selected = truth if target == "truth" else predicted
    selected.iloc[0, 0] = value
    with pytest.raises(ValueError, match=message):
        align_proportions(truth, predicted, CELL_TYPES)


@pytest.mark.parametrize("target", ["truth", "prediction"])
def test_all_zero_rows_fail(truth, target):
    predicted = truth.copy()
    selected = truth if target == "truth" else predicted
    selected.iloc[0] = 0.0
    with pytest.raises(ValueError, match="all-zero rows"):
        align_proportions(truth, predicted, CELL_TYPES)


@pytest.mark.parametrize("target", ["truth", "prediction"])
def test_rows_are_not_silently_normalized(truth, target):
    predicted = truth.copy()
    selected = truth if target == "truth" else predicted
    selected.iloc[0] *= 0.9
    with pytest.raises(ValueError, match="rows must sum to one"):
        align_proportions(truth, predicted, CELL_TYPES)


def test_frozen_absolute_row_sum_tolerance(truth):
    within = truth.copy()
    within.iloc[0, 0] += 0.5e-6
    align_proportions(truth, within, CELL_TYPES)

    outside = truth.copy()
    outside.iloc[0, 0] += 1.1e-6
    with pytest.raises(ValueError, match="atol=1e-06"):
        align_proportions(truth, outside, CELL_TYPES)


def test_nonnumeric_and_boolean_values_fail(truth):
    text = truth.copy().astype(object)
    text.iloc[0, 0] = "0.6"
    with pytest.raises(TypeError, match="must be numeric"):
        align_proportions(truth, text, CELL_TYPES)

    boolean = pd.DataFrame(
        [[True, False, False], [False, True, False]],
        index=truth.index,
        columns=CELL_TYPES,
    )
    with pytest.raises(TypeError, match="must be numeric"):
        align_proportions(truth, boolean, CELL_TYPES)


def test_versioned_jsd_is_mean_base2_divergence_not_distance(truth, prediction):
    true_aligned, predicted_aligned = align_proportions(truth, prediction, CELL_TYPES)
    distances = jensenshannon(
        true_aligned.to_numpy(),
        predicted_aligned.to_numpy(),
        axis=1,
        base=2,
    )

    assert jsd_v2(truth, prediction, CELL_TYPES) == pytest.approx(np.mean(distances**2))
    assert js_distance_v1(truth, prediction, CELL_TYPES) == pytest.approx(np.mean(distances))
    assert jsd_v2(truth, prediction, CELL_TYPES) != pytest.approx(
        js_distance_v1(truth, prediction, CELL_TYPES)
    )
    assert 0 <= jsd_v2(truth, prediction, CELL_TYPES) <= 1


def test_registry_resolves_legacy_selectors_but_reports_canonical_metadata(truth, prediction):
    assert set(PROPORTION_METRICS) == {"rmse_v1", "jsd_v2", "js_distance_v1"}
    assert resolve_proportion_metric("rmse").metric_id == "rmse_v1"
    assert resolve_proportion_metric("jsd").metric_id == "jsd_v2"

    result = evaluate_proportion_metric("jsd", truth, prediction, CELL_TYPES)
    assert result.metric_id == "jsd_v2"
    assert result.metric_name == "jsd"
    assert result.metric_version == "v2"
    assert result.contract_version == PROPORTION_CONTRACT_VERSION
    assert result.cell_types == tuple(CELL_TYPES)
    assert result.n_spots == 2
    assert result.n_cell_types == 3
    assert np.isfinite(result.value)


def test_declared_universe_is_resolved_without_inference():
    direct = {
        "modality": "atac",
        "proportion_evaluation": {"cell_types": CELL_TYPES},
    }
    assert declared_cell_types_from_metadata(direct) == tuple(CELL_TYPES)

    historical_inputs = {
        "modality": "atac",
        "dataset_config": {
            "modalities": {"atac": {"truth": {"cell_types": CELL_TYPES}}}
        },
    }
    assert declared_cell_types_from_metadata({}, historical_inputs) == tuple(CELL_TYPES)
    assert declared_cell_types_from_metadata(override=["X", "Y"]) == ("X", "Y")

    with pytest.raises(ValueError, match="No declared cell-type universe"):
        declared_cell_types_from_metadata({}, {})
