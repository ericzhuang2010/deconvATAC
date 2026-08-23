"""Fixed reference signatures and genuinely cross-fitted global dispersion."""

from __future__ import annotations

import hashlib
import json
import math
import struct
from dataclasses import dataclass
from numbers import Integral
from typing import Any, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

from ..data.validators import ordered_feature_sha256
from .config import ShapeMixConfig


DISPERSION_SEED_NAMESPACE = 20260822
DISPERSION_STREAM_ID = 17
DEFAULT_LAYER_NAMES = (
    "fragment_length_lt_100",
    "fragment_length_100_249",
    "fragment_length_ge_250",
)
DEFAULT_BIN_NAMES = ("short", "mono", "long")


def _immutable_float64(values: Any) -> np.ndarray:
    array = np.ascontiguousarray(values, dtype="<f8")
    immutable = np.frombuffer(array.tobytes(order="C"), dtype="<f8").reshape(array.shape)
    return immutable


def _immutable_int8(values: Any) -> np.ndarray:
    array = np.ascontiguousarray(values, dtype=np.int8)
    return np.frombuffer(array.tobytes(order="C"), dtype=np.int8).reshape(array.shape)


def _validate_ordered_names(values: Sequence[str], name: str) -> tuple[str, ...]:
    normalized = tuple(values)
    if not normalized:
        raise ValueError(f"{name} must be a non-empty ordered sequence.")
    if any(not isinstance(value, str) or not value for value in normalized):
        raise TypeError(f"{name} must contain non-empty strings.")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{name} must not contain duplicates.")
    return normalized


def _validate_outer_split_seed(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError("outer_split_seed must be an integer.")
    normalized = int(value)
    if normalized < 0:
        raise ValueError("outer_split_seed must be nonnegative.")
    return normalized


def _validate_labels_and_barcodes(
    labels: Sequence[Any],
    barcodes: Sequence[Any],
    cell_types: Sequence[str],
) -> tuple[np.ndarray, tuple[str, ...], tuple[str, ...]]:
    label_values = np.asarray(labels, dtype=object)
    barcode_values = np.asarray(barcodes, dtype=object)
    if label_values.ndim != 1 or barcode_values.ndim != 1:
        raise ValueError("labels and barcodes must be one-dimensional.")
    if label_values.size == 0 or label_values.size != barcode_values.size:
        raise ValueError("labels and barcodes must have the same nonzero length.")
    if pd.isna(label_values).any() or any(
        not isinstance(value, str) or not value for value in label_values
    ):
        raise ValueError("Reference labels must be non-missing, non-empty strings.")
    if pd.isna(barcode_values).any() or any(
        not isinstance(value, str) or not value for value in barcode_values
    ):
        raise ValueError("Reference barcodes must be non-missing, non-empty strings.")
    barcodes_tuple = tuple(str(value) for value in barcode_values)
    if len(set(barcodes_tuple)) != len(barcodes_tuple):
        raise ValueError("Reference barcodes must be unique.")

    ordered_types = _validate_ordered_names(cell_types, "cell_types")
    observed_types = set(label_values.tolist())
    if observed_types != set(ordered_types):
        missing = sorted(set(ordered_types).difference(observed_types))
        unexpected = sorted(observed_types.difference(ordered_types))
        raise ValueError(
            "cell_types must exactly equal the reference label universe; "
            f"missing observations={missing}, unexpected observations={unexpected}."
        )
    return label_values.astype(str), barcodes_tuple, ordered_types


def _validated_count_layer(
    matrix: Any,
    shape: tuple[int, int],
    name: str,
) -> sparse.csr_matrix:
    if getattr(matrix, "shape", None) != shape:
        raise ValueError(f"{name} must have cell-by-peak shape {shape}.")
    if sparse.issparse(matrix):
        if not np.issubdtype(matrix.dtype, np.number) or np.issubdtype(
            matrix.dtype, np.complexfloating
        ):
            raise TypeError(f"{name} must contain real numeric counts.")
        values = np.asarray(matrix.data)
    else:
        array = np.asarray(matrix)
        if array.ndim != 2 or not np.issubdtype(
            array.dtype, np.number
        ) or np.issubdtype(array.dtype, np.complexfloating):
            raise TypeError(f"{name} must be a real numeric two-dimensional matrix.")
        values = array
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} contains non-finite counts.")
    if np.any(values < 0):
        raise ValueError(f"{name} contains negative counts.")
    if np.any(values != np.floor(values)):
        raise ValueError(f"{name} contains non-integer counts.")

    validated = sparse.csr_matrix(matrix, dtype=np.float64, copy=True)
    validated.sum_duplicates()
    validated.eliminate_zeros()
    validated.sort_indices()
    return validated


def _sum_rows(matrix: sparse.csr_matrix, rows: np.ndarray) -> np.ndarray:
    if rows.size == 0:
        return np.zeros(matrix.shape[1], dtype=np.float64)
    return np.asarray(matrix[rows].sum(axis=0, dtype=np.float64)).ravel()


def _sum_squared_entries(matrix: sparse.csr_matrix, rows: np.ndarray) -> float:
    if rows.size == 0:
        return 0.0
    data = np.asarray(matrix[rows].data, dtype=np.float64)
    return float(np.dot(data, data))


def _update_string(digest: Any, value: str) -> None:
    encoded = value.encode("utf-8")
    digest.update(len(encoded).to_bytes(8, byteorder="big", signed=False))
    digest.update(encoded)


def _fold_membership_hash(
    labels: np.ndarray,
    barcodes: tuple[str, ...],
    cell_types: tuple[str, ...],
    folds: np.ndarray,
) -> str:
    digest = hashlib.sha256()
    _update_string(digest, "shapemix_dispersion_fold_membership_v1")
    barcode_array = np.asarray(barcodes, dtype=object)
    for cell_type in cell_types:
        indices = np.flatnonzero(labels == cell_type)
        ordered = indices[np.argsort(barcode_array[indices], kind="mergesort")]
        for index in ordered:
            _update_string(digest, cell_type)
            _update_string(digest, barcodes[int(index)])
            digest.update(int(folds[index]).to_bytes(1, byteorder="big", signed=False))
    return digest.hexdigest()


@dataclass(frozen=True)
class DispersionFoldAssignment:
    """Two-fold membership aligned to caller rows plus canonical audit data."""

    folds: np.ndarray
    seed_tuple: tuple[int, int, int]
    fold_membership_sha256: str
    fold_counts: tuple[tuple[str, int, int], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "folds", _immutable_int8(self.folds))


def assign_dispersion_folds(
    labels: Sequence[Any],
    barcodes: Sequence[Any],
    cell_types: Sequence[str],
    outer_split_seed: int,
) -> DispersionFoldAssignment:
    """Assign balanced within-type folds with the frozen benchmark RNG stream.

    Barcodes are sorted before each type's permutation.  A single PCG64 stream
    is initialized from ``[20260822, outer_split_seed, 17]`` and consumed in
    declared cell-type order.  Permuted positions alternate between folds.
    """
    label_values, barcode_values, ordered_types = _validate_labels_and_barcodes(
        labels, barcodes, cell_types
    )
    outer_seed = _validate_outer_split_seed(outer_split_seed)
    seed_tuple = (DISPERSION_SEED_NAMESPACE, outer_seed, DISPERSION_STREAM_ID)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed_tuple)))
    folds = np.full(label_values.size, -1, dtype=np.int8)
    barcode_array = np.asarray(barcode_values, dtype=object)
    counts: list[tuple[str, int, int]] = []

    for cell_type in ordered_types:
        indices = np.flatnonzero(label_values == cell_type)
        if indices.size < 2:
            raise ValueError(
                f"Cell type {cell_type!r} has {indices.size} reference cell(s); "
                "two-fold cross-fitting requires at least two."
            )
        ordered = indices[np.argsort(barcode_array[indices], kind="mergesort")]
        permuted = ordered[rng.permutation(indices.size)]
        assigned = np.arange(indices.size, dtype=np.int64) % 2
        folds[permuted] = assigned.astype(np.int8)
        n_fold0 = int(np.sum(assigned == 0))
        n_fold1 = int(np.sum(assigned == 1))
        if n_fold0 == 0 or n_fold1 == 0:
            raise ValueError(f"Cell type {cell_type!r} cannot populate both dispersion folds.")
        counts.append((cell_type, n_fold0, n_fold1))

    if np.any(folds < 0):
        raise AssertionError("Dispersion fold assignment left reference cells unassigned.")
    membership_hash = _fold_membership_hash(
        label_values, barcode_values, ordered_types, folds
    )
    return DispersionFoldAssignment(
        folds=folds,
        seed_tuple=seed_tuple,
        fold_membership_sha256=membership_hash,
        fold_counts=tuple(counts),
    )


@dataclass(frozen=True)
class DispersionDiagnostics:
    """Immutable scalar audit trail for the reference dispersion estimator."""

    fold_seed: tuple[int, int, int]
    bit_generator: str
    numpy_version: str
    fold_membership_sha256: str
    fold_counts: tuple[tuple[str, int, int], ...]
    numerator: float
    denominator: float
    alpha_ref_raw: float
    alpha_ref: float
    phi_ref: float
    alpha_floor: float
    rate_pseudocount: float
    crossfit_mean_min: float
    crossfit_mean_max: float

    def to_dict(self) -> dict[str, Any]:
        """Return JSON/YAML-safe scalar diagnostics."""
        return {
            "fold_seed": list(self.fold_seed),
            "bit_generator": self.bit_generator,
            "numpy_version": self.numpy_version,
            "fold_membership_sha256": self.fold_membership_sha256,
            "fold_counts": {
                cell_type: {"fold_0": fold0, "fold_1": fold1}
                for cell_type, fold0, fold1 in self.fold_counts
            },
            "numerator": self.numerator,
            "denominator": self.denominator,
            "alpha_ref_raw": self.alpha_ref_raw,
            "alpha_ref": self.alpha_ref,
            "phi_ref": self.phi_ref,
            "alpha_floor": self.alpha_floor,
            "rate_pseudocount": self.rate_pseudocount,
            "crossfit_mean_min": self.crossfit_mean_min,
            "crossfit_mean_max": self.crossfit_mean_max,
        }


def estimate_crossfit_dispersion(
    collapsed_counts: Any,
    labels: Sequence[Any],
    barcodes: Sequence[Any],
    cell_types: Sequence[str],
    outer_split_seed: int,
    *,
    rate_pseudocount: float = 0.5,
    alpha_floor: float = 1.0e-8,
) -> DispersionDiagnostics:
    """Estimate the global inverse-dispersion with opposite-fold-only means."""
    label_values, barcode_values, ordered_types = _validate_labels_and_barcodes(
        labels, barcodes, cell_types
    )
    rate_pseudocount = float(rate_pseudocount)
    alpha_floor = float(alpha_floor)
    if not math.isfinite(rate_pseudocount) or rate_pseudocount <= 0:
        raise ValueError("rate_pseudocount must be finite and strictly positive.")
    if not math.isfinite(alpha_floor) or alpha_floor <= 0:
        raise ValueError("alpha_floor must be finite and strictly positive.")
    observed_shape = getattr(collapsed_counts, "shape", None)
    if observed_shape is None or len(observed_shape) != 2 or int(observed_shape[1]) <= 0:
        raise ValueError("collapsed_counts must be a non-empty cell-by-peak matrix.")
    shape = (label_values.size, int(observed_shape[1]))
    collapsed = _validated_count_layer(collapsed_counts, shape, "collapsed_counts")
    assignment = assign_dispersion_folds(
        label_values, barcode_values, ordered_types, outer_split_seed
    )
    folds = assignment.folds

    numerator = 0.0
    denominator = 0.0
    mean_min = math.inf
    mean_max = -math.inf
    for target_fold in (0, 1):
        opposite_rows = np.flatnonzero(folds == 1 - target_fold)
        target_rows_all = np.flatnonzero(folds == target_fold)
        if opposite_rows.size == 0 or target_rows_all.size == 0:
            raise ValueError("Both global dispersion folds must contain reference cells.")

        # Both the pooled target g and the type mean exclude every target-fold
        # cell.  This is the key no-self-information cross-fit invariant.
        pooled_opposite = _sum_rows(collapsed, opposite_rows)
        g_opposite = pooled_opposite / float(opposite_rows.size)
        for cell_type in ordered_types:
            opposite_type = np.flatnonzero(
                (folds == 1 - target_fold) & (label_values == cell_type)
            )
            target_type = np.flatnonzero(
                (folds == target_fold) & (label_values == cell_type)
            )
            if opposite_type.size == 0 or target_type.size == 0:
                raise ValueError(
                    f"Cell type {cell_type!r} does not populate both dispersion folds."
                )
            opposite_sum = _sum_rows(collapsed, opposite_type)
            expected = (
                opposite_sum + float(rate_pseudocount) * g_opposite
            ) / (float(opposite_type.size) + float(rate_pseudocount))
            if not np.all(np.isfinite(expected)) or np.any(expected < 0):
                raise ValueError("Cross-fitted reference means must be finite and nonnegative.")
            mean_min = min(mean_min, float(np.min(expected)))
            mean_max = max(mean_max, float(np.max(expected)))

            expected_square_sum = float(np.dot(expected, expected))
            expected_sum = float(np.sum(expected, dtype=np.float64))
            target_sum = _sum_rows(collapsed, target_type)
            observed_square_sum = _sum_squared_entries(collapsed, target_type)
            cross_term = float(np.dot(target_sum, expected))
            n_target = float(target_type.size)
            numerator += (
                observed_square_sum
                - 2.0 * cross_term
                + n_target * expected_square_sum
                - n_target * expected_sum
            )
            denominator += n_target * expected_square_sum

    if not math.isfinite(numerator):
        raise ValueError("Cross-fitted dispersion numerator is non-finite.")
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("Cross-fitted dispersion denominator must be finite and positive.")
    alpha_raw = numerator / denominator
    if not math.isfinite(alpha_raw):
        raise ValueError("Cross-fitted alpha_ref_raw is non-finite.")
    alpha_ref = max(float(alpha_raw), float(alpha_floor))
    phi_ref = 1.0 / alpha_ref
    if not math.isfinite(alpha_ref) or not math.isfinite(phi_ref) or phi_ref <= 0:
        raise ValueError("Cross-fitted dispersion results must be finite and positive.")

    return DispersionDiagnostics(
        fold_seed=assignment.seed_tuple,
        bit_generator="PCG64",
        numpy_version=np.__version__,
        fold_membership_sha256=assignment.fold_membership_sha256,
        fold_counts=assignment.fold_counts,
        numerator=float(numerator),
        denominator=float(denominator),
        alpha_ref_raw=float(alpha_raw),
        alpha_ref=float(alpha_ref),
        phi_ref=float(phi_ref),
        alpha_floor=float(alpha_floor),
        rate_pseudocount=float(rate_pseudocount),
        crossfit_mean_min=float(mean_min),
        crossfit_mean_max=float(mean_max),
    )


@dataclass(frozen=True)
class SignatureDiagnostics:
    """Immutable dimensions and scalar conservation/probability audit values."""

    n_cells: int
    n_cell_types: int
    n_peaks: int
    n_bins: int
    cell_counts: tuple[tuple[str, int], ...]
    bin_totals: tuple[tuple[str, int], ...]
    total_cut_sites: int
    accessibility_min: float
    accessibility_max: float
    omega_min: float
    omega_max: float
    omega_max_sum_error: float
    feature_sha256: str
    signature_parameter_sha256: str

    @property
    def config_sha256(self) -> str:
        """Backward-friendly alias for the signature-relevant parameter hash."""
        return self.signature_parameter_sha256

    def to_dict(self) -> dict[str, Any]:
        """Return JSON/YAML-safe signature diagnostics."""
        return {
            "n_cells": self.n_cells,
            "n_cell_types": self.n_cell_types,
            "n_peaks": self.n_peaks,
            "n_bins": self.n_bins,
            "cell_counts": dict(self.cell_counts),
            "bin_totals": dict(self.bin_totals),
            "total_cut_sites": self.total_cut_sites,
            "accessibility_min": self.accessibility_min,
            "accessibility_max": self.accessibility_max,
            "omega_min": self.omega_min,
            "omega_max": self.omega_max,
            "omega_max_sum_error": self.omega_max_sum_error,
            "feature_sha256": self.feature_sha256,
            "signature_parameter_sha256": self.signature_parameter_sha256,
        }


def signature_parameter_sha256(config: ShapeMixConfig) -> str:
    """Hash only settings that mathematically determine fixed signatures."""
    relevant = {
        "signature_rate_pseudocount": config.signature_rate_pseudocount,
        "signature_shape_concentration": config.signature_shape_concentration,
        "exposure_mode": config.exposure_mode,
        "dispersion_mode": config.dispersion_mode,
        "dispersion_crossfit_folds": config.dispersion_crossfit_folds,
        "dispersion_alpha_floor": config.dispersion_alpha_floor,
    }
    encoded = json.dumps(
        relevant, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _update_array(digest: Any, name: str, values: np.ndarray) -> None:
    _update_string(digest, name)
    normalized = np.ascontiguousarray(values, dtype="<f8")
    digest.update(len(normalized.shape).to_bytes(1, byteorder="big", signed=False))
    for dimension in normalized.shape:
        digest.update(int(dimension).to_bytes(8, byteorder="big", signed=False))
    digest.update(normalized.tobytes(order="C"))


def _signature_content_hash(
    cell_types: tuple[str, ...],
    peak_ids: tuple[str, ...],
    bin_names: tuple[str, ...],
    layer_names: tuple[str, ...],
    A: np.ndarray,
    omega: np.ndarray,
    u_global: np.ndarray,
    u_peak: np.ndarray,
    dispersion: DispersionDiagnostics,
    config: ShapeMixConfig,
) -> str:
    digest = hashlib.sha256()
    _update_string(digest, "shapemix_reference_signatures_v1")
    for domain, names in (
        ("cell_types", cell_types),
        ("peak_ids", peak_ids),
        ("bin_names", bin_names),
        ("layer_names", layer_names),
    ):
        _update_string(digest, domain)
        for name in names:
            _update_string(digest, name)
    for name, values in (
        ("A", A),
        ("omega", omega),
        ("u_global", u_global),
        ("u_peak", u_peak),
    ):
        _update_array(digest, name, values)
    _update_string(digest, signature_parameter_sha256(config))
    _update_string(digest, dispersion.fold_membership_sha256)
    for seed_component in dispersion.fold_seed:
        digest.update(int(seed_component).to_bytes(8, byteorder="big", signed=False))
    for value in (
        dispersion.numerator,
        dispersion.denominator,
        dispersion.alpha_ref_raw,
        dispersion.alpha_ref,
        dispersion.phi_ref,
    ):
        digest.update(struct.pack("<d", float(value)))
    return digest.hexdigest()


@dataclass(frozen=True)
class ReferenceSignatures:
    """Ordered, immutable fixed signatures consumed by both model arms."""

    A: np.ndarray
    omega: np.ndarray
    u_global: np.ndarray
    u_peak: np.ndarray
    cell_types: tuple[str, ...]
    peak_ids: tuple[str, ...]
    bin_names: tuple[str, ...]
    layer_names: tuple[str, ...]
    dispersion: DispersionDiagnostics
    diagnostics: SignatureDiagnostics
    config: ShapeMixConfig
    content_sha256: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "A", _immutable_float64(self.A))
        object.__setattr__(self, "omega", _immutable_float64(self.omega))
        object.__setattr__(self, "u_global", _immutable_float64(self.u_global))
        object.__setattr__(self, "u_peak", _immutable_float64(self.u_peak))

    @property
    def phi_ref(self) -> float:
        """Return the fixed per-reference-cell inverse-dispersion."""
        return self.dispersion.phi_ref

    def to_metadata(self) -> dict[str, Any]:
        """Return compact metadata without embedding the numeric tensors."""
        return {
            "schema_version": 1,
            "content_sha256": self.content_sha256,
            "cell_types": list(self.cell_types),
            "feature_sha256": self.diagnostics.feature_sha256,
            "bin_names": list(self.bin_names),
            "layer_names": list(self.layer_names),
            "signature_parameter_sha256": signature_parameter_sha256(self.config),
            "signature_diagnostics": self.diagnostics.to_dict(),
            "dispersion_diagnostics": self.dispersion.to_dict(),
        }


def estimate_reference_signatures_from_layers(
    layers: Sequence[Any],
    labels: Sequence[Any],
    barcodes: Sequence[Any],
    *,
    cell_types: Sequence[str],
    peak_ids: Sequence[str],
    outer_split_seed: int,
    config: Optional[ShapeMixConfig] = None,
    layer_names: Optional[Sequence[str]] = None,
    bin_names: Optional[Sequence[str]] = None,
) -> ReferenceSignatures:
    """Estimate ``A``, hierarchical ``omega``, and cross-fitted ``phi_ref``.

    This array-level entry point accepts dense or sparse cell-by-peak layers and
    is the reusable numerical core beneath the AnnData wrapper.
    """
    resolved_config = ShapeMixConfig() if config is None else config
    if not isinstance(resolved_config, ShapeMixConfig):
        raise TypeError("config must be a ShapeMixConfig.")
    label_values, barcode_values, ordered_types = _validate_labels_and_barcodes(
        labels, barcodes, cell_types
    )
    ordered_peaks = _validate_ordered_names(peak_ids, "peak_ids")
    if not isinstance(layers, Sequence) or len(layers) == 0:
        raise ValueError("layers must be a non-empty ordered sequence.")
    if layer_names is None:
        ordered_layers = tuple(f"bin_{index}" for index in range(len(layers)))
    else:
        ordered_layers = _validate_ordered_names(layer_names, "layer_names")
        if len(ordered_layers) != len(layers):
            raise ValueError("layer_names must align one-to-one with layers.")
    if bin_names is None:
        ordered_bins = ordered_layers
    else:
        ordered_bins = _validate_ordered_names(bin_names, "bin_names")
        if len(ordered_bins) != len(layers):
            raise ValueError("bin_names must align one-to-one with layers.")

    shape = (label_values.size, len(ordered_peaks))
    validated_layers = tuple(
        _validated_count_layer(layer, shape, f"layers[{ordered_layers[index]!r}]")
        for index, layer in enumerate(layers)
    )
    collapsed = sparse.csr_matrix(shape, dtype=np.float64)
    for layer in validated_layers:
        collapsed = (collapsed + layer).tocsr()
    collapsed.sum_duplicates()
    collapsed.eliminate_zeros()
    collapsed.sort_indices()
    if collapsed.nnz == 0:
        raise ValueError("Training-reference shape layers contain no cut sites.")

    n_types = len(ordered_types)
    n_peaks = len(ordered_peaks)
    n_bins = len(validated_layers)
    type_bin_counts = np.zeros((n_types, n_peaks, n_bins), dtype=np.float64)
    type_cell_counts: list[tuple[str, int]] = []
    for type_index, cell_type in enumerate(ordered_types):
        rows = np.flatnonzero(label_values == cell_type)
        type_cell_counts.append((cell_type, int(rows.size)))
        for bin_index, layer in enumerate(validated_layers):
            type_bin_counts[type_index, :, bin_index] = _sum_rows(layer, rows)

    pooled_peak_bin = np.sum(type_bin_counts, axis=0, dtype=np.float64)
    pooled_peak_total = np.sum(pooled_peak_bin, axis=1, dtype=np.float64)
    if np.any(pooled_peak_total <= 0):
        first = int(np.flatnonzero(pooled_peak_total <= 0)[0])
        raise ValueError(
            f"Selected peak {ordered_peaks[first]!r} has zero total reference counts."
        )

    global_bin_counts = np.sum(pooled_peak_bin, axis=0, dtype=np.float64)
    global_denominator = float(np.sum(global_bin_counts, dtype=np.float64))
    if not math.isfinite(global_denominator) or global_denominator <= 0:
        raise ValueError("Global shape denominator must be finite and positive.")
    u_global = global_bin_counts / global_denominator
    if not np.all(np.isfinite(u_global)) or np.any(u_global <= 0):
        unsupported = [
            ordered_bins[index] for index in np.flatnonzero(u_global <= 0)
        ]
        raise ValueError(
            "Every declared shape bin must have strictly positive global support; "
            f"unsupported bins={unsupported}."
        )

    alpha_omega = resolved_config.signature_shape_concentration
    u_peak = (pooled_peak_bin + alpha_omega * u_global[None, :]) / (
        pooled_peak_total[:, None] + alpha_omega
    )
    type_peak_total = np.sum(type_bin_counts, axis=2, dtype=np.float64)
    omega = (type_bin_counts + alpha_omega * u_peak[None, :, :]) / (
        type_peak_total[:, :, None] + alpha_omega
    )
    omega_sums = np.sum(omega, axis=2, dtype=np.float64)
    if (
        not np.all(np.isfinite(u_peak))
        or not np.all(np.isfinite(omega))
        or np.any(u_peak <= 0)
        or np.any(omega <= 0)
    ):
        raise ValueError("Hierarchically smoothed shape probabilities must be finite and positive.")
    omega_sum_error = float(np.max(np.abs(omega_sums - 1.0)))
    if omega_sum_error > 1.0e-12:
        raise ValueError("Every omega cell-type/peak row must sum to one.")

    collapsed_type_counts = np.sum(type_bin_counts, axis=2, dtype=np.float64)
    g = np.asarray(collapsed.sum(axis=0, dtype=np.float64)).ravel() / float(
        shape[0]
    )
    alpha_A = resolved_config.signature_rate_pseudocount
    cell_count_array = np.asarray(
        [count for _, count in type_cell_counts], dtype=np.float64
    )
    A = (collapsed_type_counts + alpha_A * g[None, :]) / (
        cell_count_array[:, None] + alpha_A
    )
    if not np.all(np.isfinite(A)) or np.any(A <= 0):
        raise ValueError(
            "Smoothed accessibility signatures must be finite and strictly positive "
            "for every selected peak."
        )

    dispersion = estimate_crossfit_dispersion(
        collapsed,
        label_values,
        barcode_values,
        ordered_types,
        outer_split_seed,
        rate_pseudocount=alpha_A,
        alpha_floor=resolved_config.dispersion_alpha_floor,
    )
    feature_hash = ordered_feature_sha256(ordered_peaks)
    bin_totals = tuple(
        (ordered_bins[index], int(global_bin_counts[index]))
        for index in range(n_bins)
    )
    signature_diagnostics = SignatureDiagnostics(
        n_cells=shape[0],
        n_cell_types=n_types,
        n_peaks=n_peaks,
        n_bins=n_bins,
        cell_counts=tuple(type_cell_counts),
        bin_totals=bin_totals,
        total_cut_sites=int(global_denominator),
        accessibility_min=float(np.min(A)),
        accessibility_max=float(np.max(A)),
        omega_min=float(np.min(omega)),
        omega_max=float(np.max(omega)),
        omega_max_sum_error=omega_sum_error,
        feature_sha256=feature_hash,
        signature_parameter_sha256=signature_parameter_sha256(resolved_config),
    )
    content_hash = _signature_content_hash(
        ordered_types,
        ordered_peaks,
        ordered_bins,
        ordered_layers,
        A,
        omega,
        u_global,
        u_peak,
        dispersion,
        resolved_config,
    )
    return ReferenceSignatures(
        A=A,
        omega=omega,
        u_global=u_global,
        u_peak=u_peak,
        cell_types=ordered_types,
        peak_ids=ordered_peaks,
        bin_names=ordered_bins,
        layer_names=ordered_layers,
        dispersion=dispersion,
        diagnostics=signature_diagnostics,
        config=resolved_config,
        content_sha256=content_hash,
    )


def _matrices_equal(left: Any, right: sparse.csr_matrix) -> bool:
    observed = sparse.csr_matrix(left, dtype=np.float64, copy=True)
    observed.sum_duplicates()
    observed.eliminate_zeros()
    difference = observed - right
    difference.eliminate_zeros()
    return difference.nnz == 0


def estimate_reference_signatures_from_array(
    shape_counts: Any,
    labels: Sequence[Any],
    barcodes: Sequence[Any],
    *,
    cell_types: Sequence[str],
    peak_ids: Sequence[str],
    outer_split_seed: int,
    config: Optional[ShapeMixConfig] = None,
    layer_names: Optional[Sequence[str]] = None,
    bin_names: Optional[Sequence[str]] = None,
) -> ReferenceSignatures:
    """Dense ``cells × peaks × bins`` convenience wrapper for toy/core use."""
    values = np.asarray(shape_counts)
    if values.ndim != 3 or min(values.shape) <= 0:
        raise ValueError("shape_counts must have non-empty cells-by-peaks-by-bins axes.")
    return estimate_reference_signatures_from_layers(
        [values[:, :, index] for index in range(values.shape[2])],
        labels,
        barcodes,
        cell_types=cell_types,
        peak_ids=peak_ids,
        outer_split_seed=outer_split_seed,
        config=config,
        layer_names=layer_names,
        bin_names=bin_names,
    )


def _fragment_shape_axes(reference: Any) -> Optional[tuple[tuple[str, ...], tuple[str, ...]]]:
    """Read ordered bin/layer axes from list or H5AD-safe numeric-key metadata."""
    metadata = getattr(reference, "uns", {}).get("fragment_shape")
    if not isinstance(metadata, Mapping) or "bins" not in metadata:
        return None
    raw_bins = metadata["bins"]
    if isinstance(raw_bins, Mapping):
        keys = list(raw_bins)
        if keys and all(str(key).isdigit() for key in keys):
            keys.sort(key=lambda key: int(str(key)))
        elif keys and all(
            isinstance(raw_bins[key], Mapping) and "order" in raw_bins[key]
            for key in keys
        ):
            keys.sort(key=lambda key: int(raw_bins[key]["order"]))
        records = [raw_bins[key] for key in keys]
    elif isinstance(raw_bins, Sequence) and not isinstance(raw_bins, (str, bytes)):
        records = list(raw_bins)
        if records and all(isinstance(record, Mapping) and "order" in record for record in records):
            records.sort(key=lambda record: int(record["order"]))
    else:
        raise TypeError("fragment_shape.bins must be an ordered mapping or sequence.")
    if not records or any(not isinstance(record, Mapping) for record in records):
        raise ValueError("fragment_shape.bins must contain bin mappings.")
    for expected_order, record in enumerate(records):
        if "order" in record:
            observed_order = record["order"]
            if isinstance(observed_order, bool) or not isinstance(
                observed_order, Integral
            ) or int(observed_order) != expected_order:
                raise ValueError(
                    "fragment_shape bin order fields must be consecutive integers "
                    "starting at zero."
                )
    try:
        raw_layer_names = tuple(record["layer"] for record in records)
        raw_bin_names = tuple(record["name"] for record in records)
    except KeyError as error:
        raise ValueError(f"fragment_shape bin is missing {error.args[0]!r}.") from error
    if any(not isinstance(value, str) or not value for value in raw_layer_names):
        raise TypeError("fragment_shape layer names must be non-empty strings.")
    if any(not isinstance(value, str) or not value for value in raw_bin_names):
        raise TypeError("fragment_shape bin names must be non-empty strings.")
    layer_names = tuple(raw_layer_names)
    bin_names = tuple(raw_bin_names)
    return (
        _validate_ordered_names(layer_names, "fragment_shape layer names"),
        _validate_ordered_names(bin_names, "fragment_shape bin names"),
    )


def estimate_reference_signatures(
    reference: Any,
    cell_type_key: str,
    cell_types: Sequence[str],
    outer_split_seed: int,
    *,
    config: Optional[ShapeMixConfig] = None,
    layer_names: Optional[Sequence[str]] = None,
    bin_names: Optional[Sequence[str]] = None,
) -> ReferenceSignatures:
    """AnnData wrapper for fixed signature and dispersion estimation."""
    if not isinstance(cell_type_key, str) or not cell_type_key:
        raise TypeError("cell_type_key must be a non-empty string.")
    if not hasattr(reference, "obs") or not hasattr(reference, "layers"):
        raise TypeError("reference must be an AnnData-like object.")
    if cell_type_key not in reference.obs:
        raise KeyError(f"Reference observations are missing {cell_type_key!r}.")
    stored_axes = _fragment_shape_axes(reference)
    if layer_names is None:
        ordered_layers = stored_axes[0] if stored_axes is not None else DEFAULT_LAYER_NAMES
    else:
        ordered_layers = _validate_ordered_names(layer_names, "layer_names")
        if stored_axes is not None and tuple(ordered_layers) != stored_axes[0]:
            raise ValueError(
                "Explicit layer_names disagree with the ordered fragment_shape bins."
            )
    if bin_names is None:
        if stored_axes is not None and stored_axes[0] == tuple(ordered_layers):
            ordered_bins = stored_axes[1]
        elif tuple(ordered_layers) == DEFAULT_LAYER_NAMES:
            ordered_bins = DEFAULT_BIN_NAMES
        else:
            ordered_bins = tuple(ordered_layers)
    else:
        ordered_bins = _validate_ordered_names(bin_names, "bin_names")
        if stored_axes is not None and tuple(ordered_bins) != stored_axes[1]:
            raise ValueError(
                "Explicit bin_names disagree with the ordered fragment_shape bins."
            )
    if len(ordered_layers) != len(ordered_bins):
        raise ValueError("layer_names and bin_names must have equal lengths.")
    missing_layers = [name for name in ordered_layers if name not in reference.layers]
    if missing_layers:
        raise ValueError(f"Reference is missing shape layers: {missing_layers}.")

    layer_matrices = tuple(reference.layers[name] for name in ordered_layers)
    shape = tuple(reference.shape)
    validated = tuple(
        _validated_count_layer(matrix, shape, f"reference.layers[{name!r}]")
        for name, matrix in zip(ordered_layers, layer_matrices)
    )
    collapsed = sparse.csr_matrix(shape, dtype=np.float64)
    for layer in validated:
        collapsed = (collapsed + layer).tocsr()
    if reference.X is None or not _matrices_equal(reference.X, collapsed):
        raise ValueError("reference.X must equal the exact sum of ordered shape layers.")

    return estimate_reference_signatures_from_layers(
        validated,
        reference.obs[cell_type_key].to_numpy(),
        reference.obs_names.astype(str).tolist(),
        cell_types=cell_types,
        peak_ids=reference.var_names.astype(str).tolist(),
        outer_split_seed=outer_split_seed,
        config=config,
        layer_names=ordered_layers,
        bin_names=ordered_bins,
    )
