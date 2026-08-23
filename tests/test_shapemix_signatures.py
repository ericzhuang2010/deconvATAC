import json
from dataclasses import FrozenInstanceError, replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from deconvatac.data import ordered_feature_sha256
from deconvatac.shapemix.config import (
    ShapeMixConfig,
    validate_nested_ablation_configs,
)
from deconvatac.shapemix.signatures import (
    DEFAULT_BIN_NAMES,
    DEFAULT_LAYER_NAMES,
    assign_dispersion_folds,
    estimate_crossfit_dispersion,
    estimate_reference_signatures,
    estimate_reference_signatures_from_array,
    estimate_reference_signatures_from_layers,
    signature_parameter_sha256,
)


CELL_TYPES = ("A", "B")
PEAK_IDS = ("peak_2", "peak_1")
BARCODES = tuple(f"{cell_type}_{index}" for cell_type in CELL_TYPES for index in range(4))
LABELS = np.asarray([cell_type for cell_type in CELL_TYPES for _ in range(4)])


def toy_shape_counts():
    values = np.zeros((8, 2, 3), dtype=np.int64)
    values[0, 0] = [4, 0, 0]
    values[1, 0] = [0, 2, 0]
    values[2, 0] = [0, 0, 1]
    values[3, 0] = [1, 0, 0]
    # Cell type A has a genuine zero-count type/peak combination at peak_1.
    values[4, 0] = [0, 1, 0]
    values[4, 1] = [0, 3, 0]
    values[5, 1] = [0, 0, 2]
    values[6, 1] = [1, 0, 0]
    values[7, 0] = [1, 0, 0]
    values[7, 1] = [0, 1, 0]
    return values


def estimate_toy(
    values=None,
    *,
    labels=LABELS,
    barcodes=BARCODES,
    cell_types=CELL_TYPES,
    peak_ids=PEAK_IDS,
    outer_split_seed=31,
    config=None,
):
    values = toy_shape_counts() if values is None else values
    return estimate_reference_signatures_from_layers(
        [sparse.csr_matrix(values[:, :, index]) for index in range(values.shape[2])],
        labels,
        barcodes,
        cell_types=cell_types,
        peak_ids=peak_ids,
        outer_split_seed=outer_split_seed,
        config=config,
        layer_names=DEFAULT_LAYER_NAMES,
        bin_names=DEFAULT_BIN_NAMES,
    )


def test_config_is_strict_typed_hashed_and_supports_nested_yaml_shape():
    config = ShapeMixConfig.from_mapping(
        {
            "method": "shapemix",
            "params": {
                "use_shape": False,
                "max_steps": 7,
                "patience": 3,
                "restarts": 1,
                "spot_batch_size": 2,
                "peak_chunk_size": 4,
                "learning_rate": 0.1,
                "tolerance": 1.0e-4,
            },
        }
    )
    assert config.use_shape is False
    assert config.max_steps == 7
    assert len(config.content_sha256) == 64
    assert config == ShapeMixConfig.from_mapping(config.to_dict())
    assert config.content_sha256 == ShapeMixConfig.from_mapping(
        dict(reversed(list(config.to_dict().items())))
    ).content_sha256
    assert config.is_protocol_v1 is False
    with pytest.raises(ValueError, match="max_steps"):
        config.validate_protocol_v1()
    with pytest.raises(FrozenInstanceError):
        config.max_steps = 10


@pytest.mark.parametrize(
    ("mapping", "error", "message"),
    [
        ({"learnng_rate": 0.1}, ValueError, "Unknown.*learnng_rate"),
        ({"total_likelihood": "poisson"}, ValueError, "not supported"),
        ({"signature_shape_concentration": 2.0}, ValueError, "frozen"),
        ({"dispersion_crossfit_folds": 3}, ValueError, "exactly two"),
        ({"use_shape": 1}, TypeError, "boolean"),
        ({"max_steps": 0}, ValueError, "strictly positive"),
        ({"learning_rate": np.nan}, ValueError, "finite"),
        ({"seed": 1}, ValueError, "fixes the method seed"),
    ],
)
def test_config_rejects_unknown_or_nonconforming_values(mapping, error, message):
    with pytest.raises(error, match=message):
        ShapeMixConfig.from_mapping(mapping)


def test_nested_ablation_validation_allows_only_use_shape_difference():
    shape = ShapeMixConfig(use_shape=True)
    count = ShapeMixConfig(use_shape=False)
    validate_nested_ablation_configs(shape, count)
    shape.validate_protocol_v1()
    count.validate_protocol_v1()

    with pytest.raises(ValueError, match="learning_rate"):
        validate_nested_ablation_configs(
            shape, replace(count, learning_rate=0.01)
        )
    with pytest.raises(ValueError, match="ordered as"):
        validate_nested_ablation_configs(count, shape)


def test_full_reference_A_and_hierarchical_omega_match_frozen_formulas():
    values = toy_shape_counts()
    result = estimate_toy(values)
    total = values.sum(axis=2, dtype=np.float64)
    g = total.sum(axis=0) / values.shape[0]
    type_bin = np.stack(
        [values[0:4].sum(axis=0), values[4:8].sum(axis=0)]
    ).astype(np.float64)
    type_total = type_bin.sum(axis=2)
    expected_A = (type_total + 0.5 * g[None, :]) / 4.5

    global_counts = values.sum(axis=(0, 1), dtype=np.float64)
    expected_global = global_counts / global_counts.sum()
    pooled_peak_bin = values.sum(axis=0, dtype=np.float64)
    expected_peak = (pooled_peak_bin + expected_global[None, :]) / (
        pooled_peak_bin.sum(axis=1, keepdims=True) + 1.0
    )
    expected_omega = (type_bin + expected_peak[None, :, :]) / (
        type_total[:, :, None] + 1.0
    )

    assert result.cell_types == CELL_TYPES
    assert result.peak_ids == PEAK_IDS
    assert result.bin_names == DEFAULT_BIN_NAMES
    assert result.layer_names == DEFAULT_LAYER_NAMES
    assert np.allclose(result.A, expected_A, rtol=0, atol=1.0e-15)
    assert np.allclose(result.u_global, expected_global, rtol=0, atol=1.0e-15)
    assert np.allclose(result.u_peak, expected_peak, rtol=0, atol=1.0e-15)
    assert np.allclose(result.omega, expected_omega, rtol=0, atol=1.0e-15)
    assert np.all(result.A > 0)
    assert np.all(result.u_peak > 0)
    assert np.all(result.omega > 0)
    assert np.allclose(result.omega.sum(axis=2), 1.0, rtol=0, atol=1.0e-15)
    # The zero-count A/peak_1 combination is finite and shrunk twice.
    assert type_total[0, 1] == 0
    assert np.all(result.omega[0, 1] > 0)

    dense_result = estimate_reference_signatures_from_array(
        values,
        LABELS,
        BARCODES,
        cell_types=CELL_TYPES,
        peak_ids=PEAK_IDS,
        outer_split_seed=31,
        layer_names=DEFAULT_LAYER_NAMES,
        bin_names=DEFAULT_BIN_NAMES,
    )
    assert dense_result.content_sha256 == result.content_sha256


def _manual_crossfit(values, folds):
    total = values.sum(axis=2, dtype=np.float64)
    numerator = 0.0
    denominator = 0.0
    leaky_numerator = 0.0
    full_g = total.mean(axis=0)
    for target_fold in (0, 1):
        opposite = folds == 1 - target_fold
        target = folds == target_fold
        g_opposite = total[opposite].mean(axis=0)
        for cell_type in CELL_TYPES:
            type_mask = LABELS == cell_type
            opposite_type = opposite & type_mask
            target_type = target & type_mask
            type_sum = total[opposite_type].sum(axis=0)
            expected = (type_sum + 0.5 * g_opposite) / (opposite_type.sum() + 0.5)
            leaky_expected = (type_sum + 0.5 * full_g) / (opposite_type.sum() + 0.5)
            observed = total[target_type]
            numerator += np.sum((observed - expected) ** 2 - expected)
            denominator += observed.shape[0] * np.sum(expected**2)
            leaky_numerator += np.sum(
                (observed - leaky_expected) ** 2 - leaky_expected
            )
    return float(numerator), float(denominator), float(leaky_numerator)


def test_crossfit_dispersion_matches_direct_per_cell_opposite_fold_formula():
    values = toy_shape_counts()
    assignment = assign_dispersion_folds(LABELS, BARCODES, CELL_TYPES, 31)
    diagnostics = estimate_crossfit_dispersion(
        values.sum(axis=2), LABELS, BARCODES, CELL_TYPES, 31
    )
    numerator, denominator, leaky_numerator = _manual_crossfit(
        values, assignment.folds
    )

    assert diagnostics.fold_seed == (20260822, 31, 17)
    assert diagnostics.bit_generator == "PCG64"
    assert diagnostics.fold_counts == (("A", 2, 2), ("B", 2, 2))
    assert diagnostics.fold_membership_sha256 == assignment.fold_membership_sha256
    assert diagnostics.fold_membership_sha256 == (
        "f6f3992ab298da18644f10ea0e3dfb5f0e4cbdbfe5369535c1e6bcad811f9956"
    )
    assert diagnostics.numerator == pytest.approx(numerator, abs=1.0e-12)
    assert diagnostics.denominator == pytest.approx(denominator, abs=1.0e-12)
    assert diagnostics.alpha_ref_raw == pytest.approx(
        numerator / denominator, abs=1.0e-15
    )
    assert diagnostics.alpha_ref == max(diagnostics.alpha_ref_raw, 1.0e-8)
    assert diagnostics.phi_ref == pytest.approx(1.0 / diagnostics.alpha_ref)
    # A pooled target computed from all cells would leak target-fold counts and
    # gives a measurably different scalar on this toy case.
    assert abs(numerator - leaky_numerator) > 1.0e-3


def test_fold_and_signature_results_are_invariant_to_input_row_order():
    values = toy_shape_counts()
    baseline_assignment = assign_dispersion_folds(LABELS, BARCODES, CELL_TYPES, 31)
    baseline = estimate_toy(values)
    permutation = np.asarray([7, 0, 5, 2, 4, 1, 6, 3])
    reordered_assignment = assign_dispersion_folds(
        LABELS[permutation],
        [BARCODES[index] for index in permutation],
        CELL_TYPES,
        31,
    )
    reordered = estimate_toy(
        values[permutation],
        labels=LABELS[permutation],
        barcodes=tuple(BARCODES[index] for index in permutation),
    )

    baseline_by_barcode = dict(zip(BARCODES, baseline_assignment.folds))
    reordered_by_barcode = dict(
        zip((BARCODES[index] for index in permutation), reordered_assignment.folds)
    )
    assert baseline_by_barcode == reordered_by_barcode
    assert baseline_assignment.fold_membership_sha256 == (
        reordered_assignment.fold_membership_sha256
    )
    assert np.array_equal(baseline.A, reordered.A)
    assert np.array_equal(baseline.omega, reordered.omega)
    assert baseline.dispersion == reordered.dispersion
    assert baseline.content_sha256 == reordered.content_sha256


def test_declared_type_peak_and_bin_orders_are_preserved_not_inferred_from_sets():
    values = toy_shape_counts()
    baseline = estimate_toy(values)
    reversed_types = estimate_toy(values, cell_types=("B", "A"))
    reversed_peaks = estimate_toy(
        values[:, ::-1, :], peak_ids=tuple(reversed(PEAK_IDS))
    )
    reversed_bins = estimate_reference_signatures_from_layers(
        [sparse.csr_matrix(values[:, :, index]) for index in (2, 1, 0)],
        LABELS,
        BARCODES,
        cell_types=CELL_TYPES,
        peak_ids=PEAK_IDS,
        outer_split_seed=31,
        layer_names=tuple(reversed(DEFAULT_LAYER_NAMES)),
        bin_names=tuple(reversed(DEFAULT_BIN_NAMES)),
    )

    assert reversed_types.cell_types == ("B", "A")
    assert np.array_equal(reversed_types.A, baseline.A[::-1])
    assert np.array_equal(reversed_types.omega, baseline.omega[::-1])
    assert reversed_peaks.peak_ids == tuple(reversed(PEAK_IDS))
    assert np.array_equal(reversed_peaks.A, baseline.A[:, ::-1])
    assert np.array_equal(reversed_peaks.omega, baseline.omega[:, ::-1, :])
    assert reversed_bins.bin_names == tuple(reversed(DEFAULT_BIN_NAMES))
    assert np.array_equal(reversed_bins.omega, baseline.omega[:, :, ::-1])


def test_signature_hash_is_shared_by_shape_and_count_arms_and_optimizer_overrides():
    shape_config = ShapeMixConfig(use_shape=True)
    count_config = ShapeMixConfig(use_shape=False)
    fast_config = replace(shape_config, max_steps=5, patience=2, restarts=1)
    shape = estimate_toy(config=shape_config)
    count = estimate_toy(config=count_config)
    fast = estimate_toy(config=fast_config)

    assert shape_config.content_sha256 != count_config.content_sha256
    assert signature_parameter_sha256(shape_config) == signature_parameter_sha256(
        count_config
    )
    assert shape.content_sha256 == count.content_sha256 == fast.content_sha256
    assert shape.diagnostics == count.diagnostics == fast.diagnostics
    assert shape.to_metadata() == count.to_metadata() == fast.to_metadata()


def test_signatures_and_diagnostics_are_immutable_and_content_audited():
    result = estimate_toy()
    assert len(result.content_sha256) == 64
    assert len(result.dispersion.fold_membership_sha256) == 64
    assert result.diagnostics.feature_sha256 == ordered_feature_sha256(PEAK_IDS)
    assert result.diagnostics.total_cut_sites == int(toy_shape_counts().sum())
    assert result.diagnostics.config_sha256 == signature_parameter_sha256(result.config)
    assert json.loads(json.dumps(result.to_metadata()))["schema_version"] == 1
    with pytest.raises(ValueError, match="read-only"):
        result.A[0, 0] = 0
    with pytest.raises(ValueError, match="WRITEABLE"):
        result.omega.setflags(write=True)


def _toy_anndata(values=None):
    values = toy_shape_counts() if values is None else values
    layers = {
        layer: sparse.csr_matrix(values[:, :, index])
        for index, layer in enumerate(DEFAULT_LAYER_NAMES)
    }
    reference = ad.AnnData(
        X=sum(
            layers.values(), sparse.csr_matrix(values.shape[:2], dtype=np.int64)
        ).tocsr(),
        obs=pd.DataFrame({"cell_type": LABELS}, index=BARCODES),
        var=pd.DataFrame(index=PEAK_IDS),
    )
    for name, matrix in layers.items():
        reference.layers[name] = matrix
    # Deliberately insert numeric keys out of order, as H5AD groups may return.
    reference.uns["fragment_shape"] = {
        "bins": {
            "2": {"name": "long", "order": 2, "layer": DEFAULT_LAYER_NAMES[2]},
            "0": {"name": "short", "order": 0, "layer": DEFAULT_LAYER_NAMES[0]},
            "1": {"name": "mono", "order": 1, "layer": DEFAULT_LAYER_NAMES[1]},
        }
    }
    return reference


def test_anndata_wrapper_reads_roundtrip_bin_mapping_and_checks_conservation():
    reference = _toy_anndata()
    result = estimate_reference_signatures(
        reference, "cell_type", CELL_TYPES, 31
    )
    assert result.layer_names == DEFAULT_LAYER_NAMES
    assert result.bin_names == DEFAULT_BIN_NAMES

    list_encoded = reference.copy()
    list_encoded.uns["fragment_shape"]["bins"] = [
        {"name": name, "order": order, "layer": layer}
        for order, (name, layer) in enumerate(zip(DEFAULT_BIN_NAMES, DEFAULT_LAYER_NAMES))
    ]
    list_result = estimate_reference_signatures(
        list_encoded, "cell_type", CELL_TYPES, 31
    )
    assert list_result.content_sha256 == result.content_sha256

    with pytest.raises(ValueError, match="Explicit layer_names disagree"):
        estimate_reference_signatures(
            reference,
            "cell_type",
            CELL_TYPES,
            31,
            layer_names=tuple(reversed(DEFAULT_LAYER_NAMES)),
        )
    with pytest.raises(ValueError, match="Explicit bin_names disagree"):
        estimate_reference_signatures(
            reference,
            "cell_type",
            CELL_TYPES,
            31,
            bin_names=tuple(reversed(DEFAULT_BIN_NAMES)),
        )

    broken = reference.copy()
    broken.X[0, 0] += 1
    with pytest.raises(ValueError, match="exact sum"):
        estimate_reference_signatures(broken, "cell_type", CELL_TYPES, 31)


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -1.0, 0.5])
def test_layer_count_contract_rejects_nonfinite_negative_or_fractional_counts(invalid):
    values = toy_shape_counts().astype(np.float64)
    values[0, 0, 0] = invalid
    with pytest.raises(ValueError, match="non-finite|negative|non-integer"):
        estimate_toy(values)


def test_signature_failures_cover_types_support_peaks_barcodes_and_dispersion():
    values = toy_shape_counts()
    one_cell_type = LABELS.copy()
    one_cell_type[1:] = "B"
    with pytest.raises(ValueError, match="at least two"):
        estimate_toy(values, labels=one_cell_type)

    with pytest.raises(ValueError, match="unique"):
        estimate_toy(values, barcodes=(BARCODES[0], *BARCODES[:-1]))

    unsupported_bin = values.copy()
    unsupported_bin[:, :, 2] = 0
    with pytest.raises(ValueError, match="strictly positive global support"):
        estimate_toy(unsupported_bin)

    zero_peak = np.concatenate(
        [values, np.zeros((values.shape[0], 1, values.shape[2]), dtype=np.int64)],
        axis=1,
    )
    with pytest.raises(ValueError, match="zero total reference counts"):
        estimate_toy(zero_peak, peak_ids=("p0", "p1", "empty"))

    all_zero = np.zeros((4, 2), dtype=np.int64)
    with pytest.raises(ValueError, match="denominator"):
        estimate_crossfit_dispersion(
            all_zero,
            ["A", "A", "B", "B"],
            ["a0", "a1", "b0", "b1"],
            CELL_TYPES,
            0,
        )
