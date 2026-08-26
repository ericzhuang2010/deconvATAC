from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
from scipy import sparse

from deconvatac.shapemix.config import ShapeMixConfig
from deconvatac.shapemix.map import (
    ABUNDANCE_EPSILON,
    MAPOptimizationError,
    _seeded_restart,
    fit_shapemix_map,
)


SPOT_NAMES = ("spot_b", "spot_a", "spot_c")
FEATURE_NAMES = ("peak_3", "peak_1", "peak_5", "peak_0", "peak_4", "peak_2")
CELL_TYPES = ("type_0", "type_1")


def test_restart_zero_is_exact_nnls_and_later_restarts_are_perturbed() -> None:
    initial = np.asarray([[0.05, 25.0, 4566.895409590894]], dtype=np.float64)
    exact_raw = _seeded_restart(initial, (20260822, 0, 0, 0, 0))
    perturbed_raw = _seeded_restart(initial, (20260822, 0, 0, 0, 1))
    exact = np.exp(exact_raw) + ABUNDANCE_EPSILON
    perturbed = np.exp(perturbed_raw) + ABUNDANCE_EPSILON

    np.testing.assert_allclose(exact, initial, rtol=0.0, atol=5.0e-12)
    assert not np.array_equal(perturbed, initial)
    np.testing.assert_array_equal(
        perturbed_raw,
        _seeded_restart(initial, (20260822, 0, 0, 0, 1)),
    )


def _test_config(use_shape: bool, *, restarts: int = 3) -> ShapeMixConfig:
    return ShapeMixConfig(
        use_shape=use_shape,
        learning_rate=0.05,
        max_steps=1_000,
        patience=25,
        tolerance=5.0e-5,
        restarts=restarts,
        spot_batch_size=2,
        peak_chunk_size=3,
    )


def _shape_identifiable_toy():
    accessibility = np.full((2, 6), 10.0)
    omega = np.empty((2, 6, 3), dtype=float)
    omega[0, :, :] = (0.8, 0.1, 0.1)
    omega[1, :, :] = (0.1, 0.1, 0.8)
    abundance = np.asarray([[8.0, 2.0], [2.0, 8.0], [6.0, 4.0]])
    counts = np.rint(
        np.einsum("sc,cp,cpb->spb", abundance, accessibility, omega)
    ).astype(np.int64)
    return counts, accessibility, omega, abundance / abundance.sum(axis=1, keepdims=True)


def _count_identifiable_toy():
    accessibility = np.asarray(
        [
            [30.0, 5.0, 30.0, 5.0, 25.0, 5.0],
            [5.0, 30.0, 5.0, 30.0, 5.0, 25.0],
        ]
    )
    omega = np.empty((2, 6, 3), dtype=float)
    omega[:, :, :] = (0.2, 0.4, 0.4)
    abundance = np.asarray([[8.0, 2.0], [2.0, 8.0], [6.0, 4.0]])
    counts = np.rint(
        np.einsum("sc,cp,cpb->spb", abundance, accessibility, omega)
    ).astype(np.int64)
    return counts, accessibility, omega, abundance / abundance.sum(axis=1, keepdims=True)


@pytest.fixture(scope="module")
def shape_fits():
    counts, accessibility, omega, truth = _shape_identifiable_toy()
    common = {
        "outer_split_seed": 0,
        "inner_mixture_seed": 0,
        "spot_names": SPOT_NAMES,
        "feature_names": FEATURE_NAMES,
        "cell_types": CELL_TYPES,
    }
    aware = fit_shapemix_map(
        counts,
        accessibility,
        omega,
        100.0,
        config=_test_config(True),
        **common,
    )
    peak_only = fit_shapemix_map(
        counts,
        accessibility,
        omega,
        100.0,
        config=_test_config(False),
        **common,
    )
    permuted = fit_shapemix_map(
        counts,
        accessibility,
        omega[::-1].copy(),
        100.0,
        config=_test_config(True, restarts=2),
        **common,
    )
    return counts, accessibility, omega, truth, aware, peak_only, permuted


def test_shape_identifiable_toy_and_permutation_negative_control(shape_fits) -> None:
    _, _, _, truth, aware, peak_only, permuted = shape_fits
    aware_error = float(np.mean(np.abs(aware.proportions - truth)))
    peak_error = float(np.mean(np.abs(peak_only.proportions - truth)))
    permuted_error = float(np.mean(np.abs(permuted.proportions - truth)))

    assert aware_error < 0.01
    assert peak_error > 0.15
    assert permuted_error > 0.25
    assert permuted_error > 20 * aware_error
    assert aware.shape_log_likelihood < 0
    assert peak_only.shape_log_likelihood == 0.0


def test_deterministic_restarts_sparse_input_and_complete_diagnostics(shape_fits) -> None:
    counts, accessibility, omega, _, aware, _, _ = shape_fits
    layers = tuple(sparse.csr_matrix(counts[:, :, index]) for index in range(counts.shape[2]))
    signatures = SimpleNamespace(
        A=accessibility,
        omega=omega,
        phi_ref=100.0,
        cell_types=CELL_TYPES,
        peak_ids=FEATURE_NAMES,
    )
    repeat = fit_shapemix_map(
        layers,
        signatures=signatures,
        config=_test_config(True),
        outer_split_seed=0,
        inner_mixture_seed=0,
        spot_names=SPOT_NAMES,
    )

    np.testing.assert_array_equal(repeat.abundance, aware.abundance)
    np.testing.assert_array_equal(repeat.proportions, aware.proportions)
    assert repeat.selected_restart == aware.selected_restart
    assert repeat.diagnostics.success
    assert repeat.diagnostics.selected_restart_record is not None
    assert len(repeat.restart_diagnostics) == 3
    assert [record.seed_tuple for record in repeat.restart_diagnostics] == [
        (20260822, 0, 0, 0, restart) for restart in range(3)
    ]
    for record in repeat.restart_diagnostics:
        assert record.history[0].step == 0
        assert len(record.history) == record.steps + 1
        assert record.history[-1].step == record.steps
        assert record.total_log_objective is not None
        assert record.count_log_likelihood is not None
        assert record.abundance_log_prior is not None
    assert repeat.to_diagnostics_dict()["restarts"][0]["history"][0]["step"] == 0
    assert repeat.abundance.flags.writeable is False
    assert repeat.proportions.flags.writeable is False
    with pytest.raises(ValueError):
        repeat.proportions[0, 0] = 0.0


def test_count_identifiable_toy_is_recovered_by_both_nested_arms() -> None:
    counts, accessibility, omega, truth = _count_identifiable_toy()
    results = []
    for use_shape in (True, False):
        results.append(
            fit_shapemix_map(
                counts,
                accessibility,
                omega,
                100.0,
                config=_test_config(use_shape, restarts=2),
                outer_split_seed=0,
                inner_mixture_seed=0,
            )
        )
    for result in results:
        assert np.mean(np.abs(result.proportions - truth)) < 0.01
    # Shared omega makes the conditional term a data-only constant, so both
    # arms must follow exactly the same optimizer path.
    np.testing.assert_array_equal(results[0].proportions, results[1].proportions)
    np.testing.assert_array_equal(results[0].abundance, results[1].abundance)
    assert results[0].selected_restart == results[1].selected_restart
    assert [record.steps for record in results[0].restart_diagnostics] == [
        record.steps for record in results[1].restart_diagnostics
    ]
    assert [
        record.count_log_likelihood for record in results[0].restart_diagnostics
    ] == [
        record.count_log_likelihood for record in results[1].restart_diagnostics
    ]
    assert results[0].shape_log_likelihood < 0.0
    assert results[1].shape_log_likelihood == 0.0


def test_all_zero_spots_use_uniform_positive_nnls_fallback() -> None:
    accessibility = np.asarray([[3.0, 1.0], [1.0, 3.0]])
    omega = np.empty((2, 2, 3), dtype=float)
    omega[:, :, :] = (0.2, 0.3, 0.5)
    counts = np.zeros((2, 2, 3), dtype=np.int64)
    config = ShapeMixConfig(
        use_shape=True,
        learning_rate=0.05,
        max_steps=600,
        patience=20,
        tolerance=5.0e-5,
        restarts=2,
        spot_batch_size=2,
        peak_chunk_size=2,
    )

    result = fit_shapemix_map(
        counts,
        accessibility,
        omega,
        10.0,
        config=config,
        outer_split_seed=0,
        inner_mixture_seed=0,
        spot_names=("zero_2", "zero_1"),
    )

    assert result.nnls_fallback_spots == ("zero_1", "zero_2")
    assert np.isfinite(result.abundance).all()
    assert (result.abundance > 0).all()
    np.testing.assert_allclose(result.proportions, 0.5, atol=1.0e-3)
    np.testing.assert_allclose(result.proportions.sum(axis=1), 1.0, atol=1.0e-12)


def test_spot_and_peak_order_invariance() -> None:
    counts, accessibility, omega, _ = _shape_identifiable_toy()
    config = _test_config(True, restarts=1)
    baseline = fit_shapemix_map(
        counts,
        accessibility,
        omega,
        100.0,
        config=config,
        outer_split_seed=0,
        inner_mixture_seed=0,
        spot_names=SPOT_NAMES,
        feature_names=FEATURE_NAMES,
    )
    spot_permutation = np.asarray([2, 0, 1])
    peak_permutation = np.asarray([4, 1, 5, 0, 3, 2])
    reordered_names = tuple(SPOT_NAMES[index] for index in spot_permutation)
    reordered_features = tuple(FEATURE_NAMES[index] for index in peak_permutation)
    reordered_layers = tuple(
        sparse.csr_matrix(counts[spot_permutation][:, peak_permutation, bin_index])
        for bin_index in range(counts.shape[2])
    )
    reordered = fit_shapemix_map(
        reordered_layers,
        accessibility[:, peak_permutation],
        omega[:, peak_permutation, :],
        100.0,
        config=config,
        outer_split_seed=0,
        inner_mixture_seed=0,
        spot_names=reordered_names,
        feature_names=reordered_features,
    )
    back_to_original = [reordered_names.index(name) for name in SPOT_NAMES]
    np.testing.assert_array_equal(
        baseline.proportions,
        reordered.proportions[back_to_original],
    )


def test_no_converged_restart_is_a_hard_failure() -> None:
    counts, accessibility, omega, _ = _count_identifiable_toy()
    config = ShapeMixConfig(
        use_shape=True,
        max_steps=1,
        patience=2,
        restarts=1,
        spot_batch_size=3,
        peak_chunk_size=6,
    )
    with pytest.raises(MAPOptimizationError, match="No finite converged") as caught:
        fit_shapemix_map(
            counts,
            accessibility,
            omega,
            100.0,
            config=config,
            outer_split_seed=0,
            inner_mixture_seed=0,
        )
    assert caught.value.diagnostics.success is False
    assert caught.value.diagnostics.selected_restart is None
    assert caught.value.diagnostics.restarts[0].stopping_reason == "max_steps"


def test_signature_axis_mismatch_is_rejected_before_fitting() -> None:
    counts, accessibility, omega, _ = _shape_identifiable_toy()
    signatures = SimpleNamespace(
        A=accessibility,
        omega=omega,
        phi_ref=100.0,
        cell_types=CELL_TYPES,
        peak_ids=FEATURE_NAMES,
    )
    with pytest.raises(ValueError, match="cell_types"):
        fit_shapemix_map(
            counts,
            signatures=signatures,
            config=_test_config(True, restarts=1),
            outer_split_seed=0,
            inner_mixture_seed=0,
            cell_types=CELL_TYPES[::-1],
        )


def test_frozen_map_defaults_remain_declared() -> None:
    config = ShapeMixConfig()
    assert config.device == "cpu"
    assert config.dtype == "float32"
    assert config.learning_rate == 0.03
    assert config.max_steps == 2_000
    assert config.patience == 100
    assert config.tolerance == 1.0e-5
    assert config.restarts == 3
    assert config.spot_batch_size == 64
    assert config.peak_chunk_size == 512
    assert config.seed == 0
