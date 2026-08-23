#!/usr/bin/env python
"""Execute the frozen ShapeMix protocol-v1 negative controls.

The fitting path deliberately has no truth argument.  Truth may be loaded by
the standardized data loader, but it is never inspected or supplied to a fit.
The emitted control proportions can be evaluated later by the shared Step 6
evaluator, after all fitting outputs have been frozen.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
import yaml
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data import DeconvolutionInput, load_deconvolution_input
from deconvatac.shapemix.config import (
    ShapeMixConfig,
    validate_nested_ablation_configs,
)
from deconvatac.shapemix.signatures import estimate_reference_signatures


PROTOCOL_ID = "shapemix_negative_controls_v1"
CONTROL_SCHEMA_VERSION = 1
PERMUTATION_SEED_NAMESPACE = 20260822
PERMUTATION_STREAM_ID = 23
PERMUTATION_BIT_GENERATOR = "PCG64"
PERMUTATION_IDENTITY_POLICY = "first_nonidentity_draw"
DEFAULT_PROPORTION_MAX_ABS_TOLERANCE = 1.0e-6
DEFAULT_ONE_BIN_SHAPE_LOG_LIKELIHOOD_ABS_TOLERANCE = 0.0
DEFAULT_POISSON_LOG_LIKELIHOOD_ABS_TOLERANCE = 1.0e-10
KNOWN_CONTROLS = (
    "homogenized_omega",
    "permuted_omega",
    "one_bin",
    "poisson_factorization",
)
DEFAULT_CONFIG_PATH = ROOT / "configs" / "experiments" / "shapemix_negative_controls.yaml"
BENCHMARK_PROTOCOL_PATH = ROOT / "docs" / "ShapeMix" / "benchmark_protocol.md"
EXECUTED_CODE_PATHS = (
    Path(__file__).resolve(),
    ROOT / "src" / "deconvatac" / "data" / "loaders.py",
    ROOT / "src" / "deconvatac" / "data" / "registry.py",
    ROOT / "src" / "deconvatac" / "data" / "schemas.py",
    ROOT / "src" / "deconvatac" / "data" / "validators.py",
    ROOT / "src" / "deconvatac" / "shapemix" / "config.py",
    ROOT / "src" / "deconvatac" / "shapemix" / "diagnostics.py",
    ROOT / "src" / "deconvatac" / "shapemix" / "likelihood.py",
    ROOT / "src" / "deconvatac" / "shapemix" / "map.py",
    ROOT / "src" / "deconvatac" / "shapemix" / "signatures.py",
)


def _canonical_json(value: Any) -> str:
    return json.dumps(
        _plain_data(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    )


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _plain_data(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _plain_data(item) for key, item in value.items()}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return [_plain_data(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_plain_data(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _update_string(digest: Any, value: str) -> None:
    encoded = value.encode("utf-8")
    digest.update(len(encoded).to_bytes(8, byteorder="big", signed=False))
    digest.update(encoded)


def _array_sha256(values: Any, *, domain: str) -> str:
    array = np.ascontiguousarray(values, dtype="<f8")
    digest = hashlib.sha256()
    _update_string(digest, domain)
    digest.update(array.ndim.to_bytes(8, byteorder="big", signed=False))
    for size in array.shape:
        digest.update(int(size).to_bytes(8, byteorder="big", signed=False))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _axis_array_sha256(
    values: Any,
    *,
    domain: str,
    axes: Sequence[Sequence[str]],
) -> str:
    digest = hashlib.sha256()
    _update_string(digest, domain)
    for axis_index, names in enumerate(axes):
        _update_string(digest, f"axis_{axis_index}")
        for name in names:
            _update_string(digest, str(name))
    digest.update(bytes.fromhex(_array_sha256(values, domain=f"{domain}_numeric")))
    return digest.hexdigest()


def _project_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def _display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT.resolve()))
    except ValueError:
        return str(resolved)


def _load_yaml_mapping(path: Path, description: str) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, Mapping):
        raise TypeError(f"{description} must contain a YAML mapping.")
    return dict(value)


@dataclass(frozen=True)
class NegativeControlConfig:
    """Validated, frozen execution configuration for protocol-v1 controls."""

    source_path: Path
    raw: dict[str, Any]
    datasets: tuple[str, ...]
    modality: str
    feature_set: str
    shape_aware_config_path: Path
    count_only_config_path: Path
    controls: tuple[str, ...]
    proportion_max_abs_tolerance: float
    one_bin_shape_log_likelihood_abs_tolerance: float
    poisson_log_likelihood_abs_tolerance: float
    output_root: Path

    @property
    def source_sha256(self) -> str:
        return _sha256_file(self.source_path)

    @property
    def resolved_sha256(self) -> str:
        return _canonical_sha256(self.raw)


def load_negative_control_config(path: str | Path) -> NegativeControlConfig:
    """Load and strictly validate the frozen negative-control declaration."""
    source_path = _project_path(path).resolve()
    raw = _load_yaml_mapping(source_path, "Negative-control config")
    allowed = {
        "schema_version",
        "protocol_id",
        "description",
        "benchmark_scope",
        "datasets",
        "modality",
        "feature_set",
        "shape_aware_config",
        "count_only_config",
        "controls",
        "permutation",
        "tolerances",
        "output_root",
    }
    unknown = sorted(set(raw).difference(allowed))
    if unknown:
        raise ValueError(f"Unknown negative-control config keys: {', '.join(unknown)}.")
    required = allowed.difference({"description"})
    missing = sorted(required.difference(raw))
    if missing:
        raise ValueError(
            f"Negative-control config is missing required keys: {', '.join(missing)}."
        )
    if raw["schema_version"] != CONTROL_SCHEMA_VERSION:
        raise ValueError("Negative-control config requires schema_version: 1.")
    if raw["protocol_id"] != PROTOCOL_ID:
        raise ValueError(f"Negative-control config requires protocol_id: {PROTOCOL_ID}.")
    if raw["benchmark_scope"] != "development_negative_controls":
        raise ValueError(
            "The checked-in control config must be development_negative_controls; "
            "primary dataset IDs are explicit CLI overrides."
        )

    datasets_value = raw["datasets"]
    if not isinstance(datasets_value, list) or not datasets_value:
        raise TypeError("datasets must be a non-empty ordered YAML list.")
    datasets = tuple(datasets_value)
    if any(not isinstance(value, str) or not value for value in datasets):
        raise TypeError("datasets must contain non-empty strings.")
    if len(set(datasets)) != len(datasets):
        raise ValueError("datasets must not contain duplicates.")
    if raw["modality"] != "atac":
        raise ValueError("ShapeMix negative controls require modality: atac.")
    if not isinstance(raw["feature_set"], str) or not raw["feature_set"]:
        raise TypeError("feature_set must be a non-empty string.")

    controls_value = raw["controls"]
    if not isinstance(controls_value, list) or not controls_value:
        raise TypeError("controls must be a non-empty ordered YAML list.")
    controls = tuple(controls_value)
    if len(set(controls)) != len(controls):
        raise ValueError("controls must not contain duplicates.")
    unsupported = sorted(set(controls).difference(KNOWN_CONTROLS))
    if unsupported:
        raise ValueError(f"Unknown negative controls: {', '.join(unsupported)}.")

    permutation = raw["permutation"]
    expected_permutation = {
        "bit_generator": PERMUTATION_BIT_GENERATOR,
        "seed_namespace": PERMUTATION_SEED_NAMESPACE,
        "stream_id": PERMUTATION_STREAM_ID,
        "seed_sequence": "[20260822, outer_split_seed, 23]",
        "identity_policy": PERMUTATION_IDENTITY_POLICY,
    }
    if permutation != expected_permutation:
        raise ValueError(
            "permutation must exactly match the frozen PCG64 SeedSequence "
            "[20260822, outer_split_seed, 23] first-nonidentity policy."
        )

    tolerances = raw["tolerances"]
    expected_tolerances = {
        "proportion_max_abs": DEFAULT_PROPORTION_MAX_ABS_TOLERANCE,
        "one_bin_shape_log_likelihood_abs": (
            DEFAULT_ONE_BIN_SHAPE_LOG_LIKELIHOOD_ABS_TOLERANCE
        ),
        "poisson_log_likelihood_abs": (
            DEFAULT_POISSON_LOG_LIKELIHOOD_ABS_TOLERANCE
        ),
    }
    if tolerances != expected_tolerances:
        raise ValueError(
            "tolerances must exactly match the frozen protocol-v1 control tolerances."
        )

    shape_path = _project_path(raw["shape_aware_config"]).resolve()
    count_path = _project_path(raw["count_only_config"]).resolve()
    output_root = _project_path(raw["output_root"]).resolve()
    for config_path in (shape_path, count_path):
        if not config_path.is_file():
            raise FileNotFoundError(config_path)

    return NegativeControlConfig(
        source_path=source_path,
        raw=raw,
        datasets=datasets,
        modality=raw["modality"],
        feature_set=raw["feature_set"],
        shape_aware_config_path=shape_path,
        count_only_config_path=count_path,
        controls=controls,
        proportion_max_abs_tolerance=float(tolerances["proportion_max_abs"]),
        one_bin_shape_log_likelihood_abs_tolerance=float(
            tolerances["one_bin_shape_log_likelihood_abs"]
        ),
        poisson_log_likelihood_abs_tolerance=float(
            tolerances["poisson_log_likelihood_abs"]
        ),
        output_root=output_root,
    )


def load_ablation_configs(
    control_config: NegativeControlConfig,
) -> tuple[ShapeMixConfig, ShapeMixConfig]:
    """Load the two frozen method arms and verify the nested ablation."""
    aware_raw = _load_yaml_mapping(
        control_config.shape_aware_config_path, "Shape-aware method config"
    )
    count_raw = _load_yaml_mapping(
        control_config.count_only_config_path, "Count-only method config"
    )
    aware = ShapeMixConfig.from_mapping(aware_raw)
    count = ShapeMixConfig.from_mapping(count_raw)
    aware.validate_protocol_v1()
    count.validate_protocol_v1()
    validate_nested_ablation_configs(aware, count)
    return aware, count


def deterministic_cell_type_permutation(
    n_cell_types: int, outer_split_seed: int
) -> tuple[np.ndarray, tuple[int, int, int], int]:
    """Return the first non-identity permutation from the frozen outer stream."""
    if isinstance(n_cell_types, bool) or not isinstance(
        n_cell_types, (int, np.integer)
    ):
        raise TypeError("n_cell_types must be an integer.")
    if int(n_cell_types) < 2:
        raise ValueError("A permutation control requires at least two cell types.")
    if isinstance(outer_split_seed, bool) or not isinstance(
        outer_split_seed, (int, np.integer)
    ):
        raise TypeError("outer_split_seed must be an integer.")
    outer_seed = int(outer_split_seed)
    if outer_seed < 0:
        raise ValueError("outer_split_seed must be nonnegative.")
    seed_tuple = (
        PERMUTATION_SEED_NAMESPACE,
        outer_seed,
        PERMUTATION_STREAM_ID,
    )
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed_tuple)))
    identity = np.arange(int(n_cell_types), dtype=np.int64)
    for draw_index in range(1, 1025):
        permutation = rng.permutation(int(n_cell_types)).astype(np.int64, copy=False)
        if not np.array_equal(permutation, identity):
            return permutation, seed_tuple, draw_index
    raise RuntimeError("Failed to draw a non-identity cell-type permutation.")


def _dataset_seeds(data: DeconvolutionInput) -> tuple[int, int]:
    dataset_config = data.metadata.get("dataset_config")
    if not isinstance(dataset_config, Mapping):
        raise ValueError("Dataset metadata must contain dataset_config.")
    simulation = dataset_config.get("simulation")
    if not isinstance(simulation, Mapping):
        raise ValueError("Dataset config must contain simulation seeds.")
    seeds: list[int] = []
    for field in ("outer_split_seed", "inner_mixture_seed"):
        value = simulation.get(field)
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise TypeError(f"simulation.{field} must be an integer.")
        value = int(value)
        if value < 0:
            raise ValueError(f"simulation.{field} must be nonnegative.")
        seeds.append(value)
    return seeds[0], seeds[1]


def _collapsed_layers(layers: Sequence[Any]) -> sparse.csr_matrix:
    if not layers:
        raise ValueError("At least one fragment-shape layer is required.")
    first = layers[0]
    collapsed = (
        first.tocsr(copy=True)
        if sparse.issparse(first)
        else sparse.csr_matrix(np.asarray(first))
    )
    for layer in layers[1:]:
        collapsed = collapsed + (
            layer.tocsr() if sparse.issparse(layer) else sparse.csr_matrix(np.asarray(layer))
        )
    collapsed = collapsed.tocsr()
    collapsed.sum_duplicates()
    collapsed.eliminate_zeros()
    collapsed.sort_indices()
    return collapsed


def _fit_summary(fit: Any, spot_names: Sequence[str], cell_types: Sequence[str]) -> dict[str, Any]:
    records = fit.restart_diagnostics
    return {
        "selected_restart": int(fit.selected_restart),
        "count_log_likelihood": float(fit.count_log_likelihood),
        "shape_log_likelihood": float(fit.shape_log_likelihood),
        "abundance_log_prior": float(fit.abundance_log_prior),
        "total_log_objective": float(fit.total_log_objective),
        "selected_steps": int(records[fit.selected_restart].steps),
        "selected_stopping_reason": records[fit.selected_restart].stopping_reason,
        "proportions_sha256": _axis_array_sha256(
            fit.proportions,
            domain="shapemix_negative_control_proportions_v1",
            axes=(tuple(map(str, spot_names)), tuple(map(str, cell_types))),
        ),
    }


def _max_abs_difference(left: Any, right: Any) -> float:
    left_array = np.asarray(left, dtype=np.float64)
    right_array = np.asarray(right, dtype=np.float64)
    if left_array.shape != right_array.shape:
        raise ValueError("Compared proportion matrices must have identical shapes.")
    difference = np.abs(left_array - right_array)
    if not np.isfinite(difference).all():
        raise ValueError("Compared proportion matrices must be finite.")
    return float(np.max(difference, initial=0.0))


def poisson_factorization_golden_check(tolerance: float) -> dict[str, Any]:
    """Evaluate the factorized Poisson identity on a frozen toy tensor."""
    import torch

    from deconvatac.shapemix.likelihood import (
        factorized_poisson_log_likelihood,
        independent_poisson_bin_log_likelihood,
    )

    counts = torch.tensor(
        [
            [[2.0, 0.0, 1.0], [0.0, 3.0, 2.0]],
            [[1.0, 4.0, 0.0], [2.0, 1.0, 3.0]],
        ],
        dtype=torch.float64,
    )
    rates = torch.tensor(
        [
            [[1.7, 0.4, 0.9], [0.3, 2.8, 1.6]],
            [[0.8, 3.2, 0.5], [1.3, 0.7, 2.4]],
        ],
        dtype=torch.float64,
    )
    factorized = float(factorized_poisson_log_likelihood(counts, rates).item())
    independent = float(
        independent_poisson_bin_log_likelihood(counts, rates).item()
    )
    difference = abs(factorized - independent)
    return {
        "status": "pass" if difference <= tolerance else "fail",
        "tolerance_abs": float(tolerance),
        "factorized_log_likelihood": factorized,
        "independent_bin_log_likelihood": independent,
        "absolute_difference": difference,
        "counts_sha256": _array_sha256(counts.numpy(), domain="poisson_golden_counts_v1"),
        "rates_sha256": _array_sha256(rates.numpy(), domain="poisson_golden_rates_v1"),
    }


def run_dataset_controls(
    data: DeconvolutionInput,
    shape_aware_config: ShapeMixConfig,
    count_only_config: ShapeMixConfig,
    controls: Sequence[str],
    *,
    proportion_tolerance: float = DEFAULT_PROPORTION_MAX_ABS_TOLERANCE,
    one_bin_shape_tolerance: float = (
        DEFAULT_ONE_BIN_SHAPE_LOG_LIKELIHOOD_ABS_TOLERANCE
    ),
    poisson_tolerance: float = DEFAULT_POISSON_LOG_LIKELIHOOD_ABS_TOLERANCE,
    fit_function: Callable[..., Any] | None = None,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Run selected controls using observations and reference signatures only."""
    selected = tuple(controls)
    if not selected or len(set(selected)) != len(selected):
        raise ValueError("controls must be a non-empty unique sequence.")
    unsupported = sorted(set(selected).difference(KNOWN_CONTROLS))
    if unsupported:
        raise ValueError(f"Unknown negative controls: {', '.join(unsupported)}.")
    if data.modality != "atac" or data.fragment_shape is None or data.cell_types is None:
        raise ValueError("Negative controls require ShapeMix-compatible ATAC data.")
    validate_nested_ablation_configs(shape_aware_config, count_only_config)
    shape_aware_config.validate_protocol_v1()
    count_only_config.validate_protocol_v1()
    outer_seed, inner_seed = _dataset_seeds(data)

    layer_names = data.fragment_shape.layer_names
    bin_names = tuple(bin_spec.name for bin_spec in data.fragment_shape.bins)
    layers = tuple(data.spatial.layers[name] for name in layer_names)
    spot_names = tuple(data.spatial.obs_names.astype(str))
    peak_names = tuple(data.spatial.var_names.astype(str))
    cell_types = tuple(data.cell_types)
    signatures = estimate_reference_signatures(
        data.reference,
        data.labels_key,
        cell_types,
        outer_seed,
        config=shape_aware_config,
        layer_names=layer_names,
        bin_names=bin_names,
    )
    if fit_function is None:
        from deconvatac.shapemix.map import fit_shapemix_map

        fit_function = fit_shapemix_map

    def fit_variant(
        variant_layers: Sequence[Any], variant_omega: np.ndarray, config: ShapeMixConfig
    ) -> Any:
        # The fit API receives no truth object.  A, phi, axes, observations, and
        # the predeclared seeds are the complete fitting inputs.
        return fit_function(
            tuple(variant_layers),
            A=signatures.A,
            omega=variant_omega,
            phi_ref=signatures.phi_ref,
            config=config,
            outer_split_seed=outer_seed,
            inner_mixture_seed=inner_seed,
            spot_names=spot_names,
            feature_names=peak_names,
            cell_types=cell_types,
        )

    predictions: dict[str, np.ndarray] = {}
    fit_summaries: dict[str, Any] = {}
    control_evidence: dict[str, Any] = {}
    count_fit = None

    def require_count_fit() -> Any:
        nonlocal count_fit
        if count_fit is None:
            count_fit = fit_variant(layers, signatures.omega, count_only_config)
            predictions["count_only"] = np.asarray(count_fit.proportions, dtype=np.float64)
            fit_summaries["count_only"] = _fit_summary(
                count_fit, spot_names, cell_types
            )
        return count_fit

    a_sha256 = _axis_array_sha256(
        signatures.A,
        domain="shapemix_accessibility_signature_v1",
        axes=(cell_types, peak_names),
    )
    omega_sha256 = _axis_array_sha256(
        signatures.omega,
        domain="shapemix_shape_signature_v1",
        axes=(cell_types, peak_names, bin_names),
    )
    pooled_sha256 = _axis_array_sha256(
        signatures.u_peak,
        domain="shapemix_pooled_peak_shape_v1",
        axes=(peak_names, bin_names),
    )

    if "homogenized_omega" in selected:
        baseline = require_count_fit()
        homogeneous = np.broadcast_to(
            signatures.u_peak[None, :, :], signatures.omega.shape
        ).copy()
        homogeneous_fit = fit_variant(layers, homogeneous, shape_aware_config)
        predictions["homogenized_omega"] = np.asarray(
            homogeneous_fit.proportions, dtype=np.float64
        )
        fit_summaries["homogenized_omega"] = _fit_summary(
            homogeneous_fit, spot_names, cell_types
        )
        difference = _max_abs_difference(
            homogeneous_fit.proportions, baseline.proportions
        )
        control_evidence["homogenized_omega"] = {
            "status": "pass" if difference <= proportion_tolerance else "fail",
            "construction": "omega[c,p,:] = u_peak[p,:] for every cell type",
            "A_preserved": True,
            "count_likelihood_preserved": True,
            "homogeneous_across_cell_types": bool(
                np.array_equal(
                    homogeneous,
                    np.broadcast_to(homogeneous[0:1], homogeneous.shape),
                )
            ),
            "pooled_u_peak_sha256": pooled_sha256,
            "control_omega_sha256": _axis_array_sha256(
                homogeneous,
                domain="shapemix_homogenized_shape_signature_v1",
                axes=(cell_types, peak_names, bin_names),
            ),
            "comparison": "homogenized_omega_vs_count_only",
            "proportion_max_abs_difference": difference,
            "proportion_max_abs_tolerance": float(proportion_tolerance),
        }

    if "permuted_omega" in selected:
        baseline = require_count_fit()
        permutation, seed_tuple, draw_index = deterministic_cell_type_permutation(
            len(cell_types), outer_seed
        )
        permuted = np.asarray(signatures.omega[permutation], dtype=np.float64)
        permuted_fit = fit_variant(layers, permuted, shape_aware_config)
        predictions["permuted_omega"] = np.asarray(
            permuted_fit.proportions, dtype=np.float64
        )
        fit_summaries["permuted_omega"] = _fit_summary(
            permuted_fit, spot_names, cell_types
        )
        control_evidence["permuted_omega"] = {
            "status": "diagnostic_complete",
            "acceptance_interpretation": (
                "Score these frozen proportions only if the primary shape direction is "
                "positive; a legitimate gain should disappear or reverse."
            ),
            "A_preserved": True,
            "count_likelihood_preserved": True,
            "bit_generator": PERMUTATION_BIT_GENERATOR,
            "seed_tuple": list(seed_tuple),
            "identity_policy": PERMUTATION_IDENTITY_POLICY,
            "selected_draw_index": draw_index,
            "permutation_indices": permutation.tolist(),
            "target_to_source_cell_type": {
                target: cell_types[int(permutation[index])]
                for index, target in enumerate(cell_types)
            },
            "fixed_points": int(np.sum(permutation == np.arange(len(cell_types)))),
            "is_nonidentity": bool(
                not np.array_equal(permutation, np.arange(len(cell_types)))
            ),
            "is_exact_cell_type_axis_permutation": bool(
                np.array_equal(permuted, signatures.omega[permutation])
            ),
            "control_omega_sha256": _axis_array_sha256(
                permuted,
                domain="shapemix_permuted_shape_signature_v1",
                axes=(cell_types, peak_names, bin_names),
            ),
            "numeric_omega_changed": bool(
                not np.array_equal(permuted, signatures.omega)
            ),
            "proportion_max_abs_difference_from_count_only": _max_abs_difference(
                permuted_fit.proportions, baseline.proportions
            ),
        }

    if "one_bin" in selected:
        collapsed = _collapsed_layers(layers)
        one_bin_omega = np.ones(
            (len(cell_types), len(peak_names), 1), dtype=np.float64
        )
        one_bin_count = fit_variant((collapsed,), one_bin_omega, count_only_config)
        one_bin_shape = fit_variant((collapsed,), one_bin_omega, shape_aware_config)
        predictions["one_bin_count_only"] = np.asarray(
            one_bin_count.proportions, dtype=np.float64
        )
        predictions["one_bin_shape"] = np.asarray(
            one_bin_shape.proportions, dtype=np.float64
        )
        fit_summaries["one_bin_count_only"] = _fit_summary(
            one_bin_count, spot_names, cell_types
        )
        fit_summaries["one_bin_shape"] = _fit_summary(
            one_bin_shape, spot_names, cell_types
        )
        difference = _max_abs_difference(
            one_bin_shape.proportions, one_bin_count.proportions
        )
        shape_abs = abs(float(one_bin_shape.shape_log_likelihood))
        passed = (
            difference <= proportion_tolerance
            and shape_abs <= one_bin_shape_tolerance
        )
        control_evidence["one_bin"] = {
            "status": "pass" if passed else "fail",
            "construction": "all declared fragment bins collapsed to one bin",
            "A_preserved": True,
            "count_likelihood_preserved": True,
            "collapsed_total_counts_sha256": _array_sha256(
                collapsed.toarray(), domain="shapemix_collapsed_total_counts_v1"
            ),
            "one_bin_omega_sha256": _axis_array_sha256(
                one_bin_omega,
                domain="shapemix_one_bin_shape_signature_v1",
                axes=(cell_types, peak_names, ("all_fragments",)),
            ),
            "shape_log_likelihood": float(one_bin_shape.shape_log_likelihood),
            "shape_log_likelihood_absolute": shape_abs,
            "shape_log_likelihood_abs_tolerance": float(one_bin_shape_tolerance),
            "comparison": "one_bin_shape_vs_one_bin_count_only",
            "proportion_max_abs_difference": difference,
            "proportion_max_abs_tolerance": float(proportion_tolerance),
        }

    if "poisson_factorization" in selected:
        control_evidence["poisson_factorization"] = poisson_factorization_golden_check(
            poisson_tolerance
        )

    acceptance_statuses = [
        evidence["status"]
        for name, evidence in control_evidence.items()
        if name != "permuted_omega"
    ]
    all_checks_passed = all(status == "pass" for status in acceptance_statuses)
    evidence = {
        "schema_version": CONTROL_SCHEMA_VERSION,
        "protocol_id": PROTOCOL_ID,
        "dataset_id": data.dataset_id,
        "modality": data.modality,
        "feature_set": data.feature_set,
        "selected_controls": list(selected),
        "truth_used_in_fitting": False,
        "all_acceptance_checks_passed": all_checks_passed,
        "axes": {
            "spots": len(spot_names),
            "peaks": len(peak_names),
            "bins": len(bin_names),
            "cell_types": list(cell_types),
        },
        "seeds": {
            "outer_split_seed": outer_seed,
            "inner_mixture_seed": inner_seed,
            "method_seed": shape_aware_config.seed,
            "permutation_seed_tuple": [
                PERMUTATION_SEED_NAMESPACE,
                outer_seed,
                PERMUTATION_STREAM_ID,
            ],
        },
        "fixed_inputs": {
            "reference_signature_sha256": signatures.content_sha256,
            "A_sha256": a_sha256,
            "omega_sha256": omega_sha256,
            "u_peak_sha256": pooled_sha256,
            "phi_ref": float(signatures.phi_ref),
            "total_likelihood": shape_aware_config.total_likelihood,
            "feature_sha256": signatures.diagnostics.feature_sha256,
            "dispersion_fold_membership_sha256": (
                signatures.dispersion.fold_membership_sha256
            ),
        },
        "fit_summaries": fit_summaries,
        "controls": control_evidence,
    }
    return evidence, predictions


def _prediction_frame(
    predictions: Mapping[str, np.ndarray],
    spot_names: Sequence[str],
    cell_types: Sequence[str],
) -> pd.DataFrame:
    frames = []
    for variant, values in predictions.items():
        matrix = np.asarray(values, dtype=np.float64)
        expected = (len(spot_names), len(cell_types))
        if matrix.shape != expected:
            raise ValueError(f"{variant} proportions have shape {matrix.shape}; expected {expected}.")
        frame = pd.DataFrame(matrix, index=pd.Index(spot_names, name="spot_id"), columns=cell_types)
        frame.insert(0, "variant", variant)
        frames.append(frame.reset_index())
    if not frames:
        return pd.DataFrame(columns=["spot_id", "variant", *cell_types])
    return pd.concat(frames, ignore_index=True)


def _git_value(arguments: Sequence[str]) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()


def _software_versions() -> dict[str, str]:
    versions = {"python": platform.python_version()}
    for name, distribution in (
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("scipy", "scipy"),
        ("torch", "torch"),
        ("pyyaml", "PyYAML"),
        ("deconvatac", "deconvATAC"),
    ):
        try:
            versions[name] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "not-installed"
    return versions


def _dataset_config_provenance(data: DeconvolutionInput) -> dict[str, Any]:
    config = data.metadata.get("dataset_config")
    if not isinstance(config, Mapping):
        raise ValueError("Dataset metadata must contain dataset_config.")
    source_value = config.get("config_path")
    if not isinstance(source_value, str) or not source_value:
        raise ValueError("dataset_config.config_path is required for provenance.")
    source_path = Path(source_value).resolve()
    resolved = {key: value for key, value in config.items() if key != "config_path"}
    return {
        "source_path": _display_path(source_path),
        "source_sha256": _sha256_file(source_path),
        "resolved_sha256": _canonical_sha256(resolved),
    }


def write_control_outputs(
    output_dir: Path,
    data: DeconvolutionInput,
    evidence: dict[str, Any],
    predictions: Mapping[str, np.ndarray],
    control_config: NegativeControlConfig,
    shape_aware_config: ShapeMixConfig,
    count_only_config: ShapeMixConfig,
    *,
    overwrite: bool,
) -> Path:
    """Write compact evidence and a non-recursive output hash manifest."""
    output_dir = output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Control output directory already contains files: {output_dir}. "
            "Pass --overwrite to replace the known output files."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    predictions_path = output_dir / "control_proportions.csv"
    evidence_path = output_dir / "control_evidence.yaml"
    manifest_path = output_dir / "output_sha256.yaml"

    spot_names = tuple(data.spatial.obs_names.astype(str))
    cell_types = tuple(data.cell_types or ())
    frame = _prediction_frame(predictions, spot_names, cell_types)
    frame.to_csv(
        predictions_path,
        index=False,
        float_format="%.17g",
        lineterminator="\n",
    )

    script_path = Path(__file__).resolve()
    git_status = _git_value(("status", "--short", "--untracked-files=normal"))
    evidence = {
        **evidence,
        "configuration": {
            "control_config_source_path": _display_path(control_config.source_path),
            "control_config_source_sha256": control_config.source_sha256,
            "control_config_resolved_sha256": control_config.resolved_sha256,
            "shape_aware_config_source_path": _display_path(
                control_config.shape_aware_config_path
            ),
            "shape_aware_config_source_sha256": _sha256_file(
                control_config.shape_aware_config_path
            ),
            "shape_aware_config_resolved_sha256": shape_aware_config.content_sha256,
            "count_only_config_source_path": _display_path(
                control_config.count_only_config_path
            ),
            "count_only_config_source_sha256": _sha256_file(
                control_config.count_only_config_path
            ),
            "count_only_config_resolved_sha256": count_only_config.content_sha256,
        },
        "dataset_config": _dataset_config_provenance(data),
        "code_and_protocol": {
            "git_commit": _git_value(("rev-parse", "HEAD")),
            "git_worktree_dirty": None if git_status is None else bool(git_status),
            "script_path": _display_path(script_path),
            "script_sha256": _sha256_file(script_path),
            "benchmark_protocol_path": _display_path(BENCHMARK_PROTOCOL_PATH),
            "benchmark_protocol_sha256": _sha256_file(BENCHMARK_PROTOCOL_PATH),
            "executed_code_sha256": {
                _display_path(path): _sha256_file(path) for path in EXECUTED_CODE_PATHS
            },
        },
        "software_versions": _software_versions(),
    }
    with evidence_path.open("w") as handle:
        yaml.safe_dump(_plain_data(evidence), handle, sort_keys=False)

    outputs = {}
    for path in (predictions_path, evidence_path):
        outputs[path.name] = {
            "bytes": path.stat().st_size,
            "sha256": _sha256_file(path),
        }
    manifest = {
        "schema_version": CONTROL_SCHEMA_VERSION,
        "protocol_id": PROTOCOL_ID,
        "dataset_id": data.dataset_id,
        "complete": True,
        "all_acceptance_checks_passed": evidence["all_acceptance_checks_passed"],
        "algorithm": "sha256",
        "outputs": outputs,
    }
    with manifest_path.open("w") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False)
    return manifest_path


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument(
        "--dataset-id",
        action="append",
        dest="dataset_ids",
        help=(
            "Override the development dataset; repeat for multiple frozen primary pairs. "
            "Use only after the preregistered primary direction is positive."
        ),
    )
    parser.add_argument("--output-root", help="Override the configured output root.")
    parser.add_argument(
        "--controls",
        nargs="+",
        choices=KNOWN_CONTROLS,
        help="Run a declared subset, e.g. homogenized_omega permuted_omega.",
    )
    parser.add_argument("--registry-path")
    parser.add_argument(
        "--primary-direction-positive",
        action="store_true",
        help=(
            "Required attestation when --dataset-id selects any dataset outside the "
            "checked-in development list."
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    control_config = load_negative_control_config(args.config)
    shape_aware, count_only = load_ablation_configs(control_config)
    dataset_ids = tuple(args.dataset_ids or control_config.datasets)
    if not dataset_ids or len(set(dataset_ids)) != len(dataset_ids):
        raise ValueError("Selected dataset IDs must be non-empty and unique.")
    primary_override = any(
        dataset_id not in control_config.datasets for dataset_id in dataset_ids
    )
    if primary_override and not args.primary_direction_positive:
        raise ValueError(
            "Primary negative controls are gated on a positive preregistered primary "
            "direction; pass --primary-direction-positive to attest that condition."
        )
    controls = tuple(args.controls or control_config.controls)
    output_root = (
        _project_path(args.output_root).resolve()
        if args.output_root
        else control_config.output_root
    )

    all_passed = True
    for dataset_id in dataset_ids:
        data = load_deconvolution_input(
            dataset_id=dataset_id,
            modality=control_config.modality,
            feature_set=control_config.feature_set,
            registry_path=args.registry_path,
            project_root=ROOT,
        )
        evidence, predictions = run_dataset_controls(
            data,
            shape_aware,
            count_only,
            controls,
            proportion_tolerance=control_config.proportion_max_abs_tolerance,
            one_bin_shape_tolerance=(
                control_config.one_bin_shape_log_likelihood_abs_tolerance
            ),
            poisson_tolerance=control_config.poisson_log_likelihood_abs_tolerance,
        )
        evidence["execution_gate"] = {
            "dataset_selection": (
                "primary_cli_override" if primary_override else "checked_in_development"
            ),
            "primary_direction_positive_attested": bool(
                args.primary_direction_positive
            ),
        }
        manifest_path = write_control_outputs(
            output_root / dataset_id,
            data,
            evidence,
            predictions,
            control_config,
            shape_aware,
            count_only,
            overwrite=args.overwrite,
        )
        all_passed = all_passed and evidence["all_acceptance_checks_passed"]
        print(f"{dataset_id}: {manifest_path}")

    return 0 if all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
