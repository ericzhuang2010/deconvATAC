#!/usr/bin/env python
"""Audit fragment-length signal using training-reference cells only."""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data import FragmentShapeSpec
from deconvatac.data.validators import validate_fragment_shape_feature_axis


SEED_NAMESPACE = 20260822
SPLIT_HALF_STREAM = 23


def _as_csr(matrix: Any) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def _sum_axis(matrix: sparse.csr_matrix, axis: int) -> np.ndarray:
    return np.asarray(matrix.sum(axis=axis)).ravel().astype(np.float64, copy=False)


def _entropy_bits(probabilities: np.ndarray, axis: int = -1) -> np.ndarray:
    values = np.asarray(probabilities, dtype=np.float64)
    terms = np.zeros_like(values)
    positive = values > 0
    terms[positive] = -values[positive] * np.log2(values[positive])
    return terms.sum(axis=axis)


def _safe_probabilities(counts: np.ndarray, pseudocount: float = 0.0) -> np.ndarray:
    values = np.asarray(counts, dtype=np.float64)
    if pseudocount < 0 or not math.isfinite(pseudocount):
        raise ValueError("pseudocount must be finite and nonnegative.")
    adjusted = values + pseudocount
    totals = adjusted.sum(axis=-1, keepdims=True)
    return np.divide(adjusted, totals, out=np.zeros_like(adjusted), where=totals > 0)


def _finite_spearman(left: np.ndarray, right: np.ndarray) -> float | None:
    mask = np.isfinite(left) & np.isfinite(right)
    if mask.sum() < 3:
        return None
    left_valid = left[mask]
    right_valid = right[mask]
    if np.ptp(left_valid) == 0 or np.ptp(right_valid) == 0:
        return None
    statistic = spearmanr(left_valid, right_valid).statistic
    return float(statistic) if np.isfinite(statistic) else None


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_reference(
    reference: ad.AnnData,
    labels_key: str,
    cell_types: Sequence[str] | None,
    split_seed: int,
) -> tuple[FragmentShapeSpec, tuple[str, ...], tuple[str, ...]]:
    if labels_key not in reference.obs:
        raise KeyError(f"Training reference is missing obs[{labels_key!r}].")
    preparation = reference.uns.get("shapemix_preparation")
    if not isinstance(preparation, Mapping) or preparation.get("pool") != "reference":
        raise ValueError(
            "Signal audit requires a ShapeMix preparation object explicitly marked as the reference pool."
        )
    if "split_pool" not in reference.obs or set(reference.obs["split_pool"].astype(str)) != {
        "reference"
    }:
        raise ValueError(
            "Signal audit requires every observation to be explicitly marked split_pool='reference'."
        )
    labels = reference.obs[labels_key]
    if labels.isna().any():
        raise ValueError("Training-reference labels contain missing values.")
    if reference.n_obs < 2 or reference.n_vars == 0:
        raise ValueError("Signal audit requires at least two cells and one peak.")
    validate_fragment_shape_feature_axis(reference, "training reference")
    metadata = reference.uns.get("fragment_shape")
    spec = FragmentShapeSpec.from_mapping(metadata)
    if preparation.get("split_sha256") != spec.split_sha256:
        raise ValueError(
            "Reference preparation metadata and fragment-shape metadata disagree on split_sha256."
        )
    if preparation.get("outer_split_seed") != split_seed:
        raise ValueError(
            "Signal-audit split_seed does not match the prepared reference outer_split_seed."
        )
    layer_names = spec.layer_names
    for layer_name in layer_names:
        if layer_name not in reference.layers:
            raise KeyError(f"Training reference is missing layer {layer_name!r}.")
        if not sparse.isspmatrix_csr(reference.layers[layer_name]):
            raise ValueError(f"Training-reference layer {layer_name!r} must be CSR.")

    observed = tuple(pd.unique(labels.astype(str)))
    ordered_types = tuple(observed if cell_types is None else map(str, cell_types))
    if not ordered_types or len(set(ordered_types)) != len(ordered_types):
        raise ValueError("cell_types must be a nonempty ordered unique sequence.")
    if set(observed) != set(ordered_types):
        raise ValueError("cell_types must exactly match the training-reference label universe.")
    return spec, layer_names, ordered_types


def _aggregate_by_type(
    matrices: Sequence[sparse.csr_matrix],
    labels: np.ndarray,
    cell_types: Sequence[str],
) -> np.ndarray:
    counts = np.zeros((len(cell_types), matrices[0].shape[1], len(matrices)), dtype=np.float64)
    for type_index, cell_type in enumerate(cell_types):
        rows = np.flatnonzero(labels == cell_type)
        if rows.size == 0:
            raise ValueError(f"Cell type {cell_type!r} has no training-reference cells.")
        for bin_index, matrix in enumerate(matrices):
            counts[type_index, :, bin_index] = _sum_axis(matrix[rows], axis=0)
    return counts


def _split_half_reproducibility(
    matrices: Sequence[sparse.csr_matrix],
    barcodes: np.ndarray,
    labels: np.ndarray,
    cell_types: Sequence[str],
    split_seed: int,
    pseudocount: float,
) -> dict[str, dict[str, float | int | None]]:
    output: dict[str, dict[str, float | int | None]] = {}
    for type_index, cell_type in enumerate(cell_types, start=1):
        positions = np.flatnonzero(labels == cell_type)
        positions = positions[np.argsort(barcodes[positions], kind="stable")]
        if positions.size < 2:
            raise ValueError(f"Cell type {cell_type!r} needs at least two cells for split-half audit.")
        rng = np.random.Generator(
            np.random.PCG64(
                np.random.SeedSequence(
                    [SEED_NAMESPACE, int(split_seed), SPLIT_HALF_STREAM, type_index]
                )
            )
        )
        permuted = positions[rng.permutation(positions.size)]
        first, second = np.array_split(permuted, 2)

        half_counts = []
        for half in (first, second):
            counts = np.column_stack([_sum_axis(matrix[half], axis=0) for matrix in matrices])
            half_counts.append(counts)
        combined = half_counts[0].sum(axis=1) + half_counts[1].sum(axis=1)
        eligible = combined >= 10
        if eligible.any():
            first_shape = _safe_probabilities(half_counts[0][eligible], pseudocount=pseudocount)
            second_shape = _safe_probabilities(half_counts[1][eligible], pseudocount=pseudocount)
            correlation = _finite_spearman(first_shape.ravel(), second_shape.ravel())
            midpoint = 0.5 * (first_shape + second_shape)
            jsd = 0.5 * (
                _entropy_bits(midpoint)
                - _entropy_bits(first_shape)
                + _entropy_bits(midpoint)
                - _entropy_bits(second_shape)
            )
            mean_jsd = float(np.mean(jsd))
        else:
            correlation = None
            mean_jsd = None
        output[str(cell_type)] = {
            "first_half_cells": int(first.size),
            "second_half_cells": int(second.size),
            "eligible_peaks_min_10_combined_cuts": int(eligible.sum()),
            "shape_spearman": correlation,
            "mean_jsd_bits": mean_jsd,
        }
    return output


def audit_shape_signal(
    reference: ad.AnnData,
    *,
    labels_key: str = "cell_type",
    cell_types: Sequence[str] | None = None,
    split_seed: int = 0,
    pseudocount: float = 0.5,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return peak-level and summary diagnostics from training-reference cells."""
    if isinstance(split_seed, (bool, np.bool_)) or not isinstance(split_seed, (int, np.integer)):
        raise TypeError("split_seed must be an integer.")
    split_seed = int(split_seed)
    if split_seed < 0:
        raise ValueError("split_seed must be nonnegative.")
    spec, layer_names, ordered_types = _validate_reference(
        reference,
        labels_key,
        cell_types,
        split_seed,
    )
    matrices = tuple(_as_csr(reference.layers[name]) for name in layer_names)
    labels = reference.obs[labels_key].astype(str).to_numpy()
    barcodes = reference.obs_names.astype(str).to_numpy()

    per_bin_peak = np.column_stack([_sum_axis(matrix, axis=0) for matrix in matrices])
    peak_totals = per_bin_peak.sum(axis=1)
    peak_probabilities = _safe_probabilities(per_bin_peak)
    peak_entropy = _entropy_bits(peak_probabilities)
    collapsed = _as_csr(reference.X)
    nonzero_cells = np.asarray(collapsed.getnnz(axis=0)).ravel()

    type_counts = _aggregate_by_type(matrices, labels, ordered_types)
    type_totals = type_counts.sum(axis=2)
    type_has_signal = type_totals > 0
    type_probabilities = _safe_probabilities(type_counts, pseudocount=pseudocount)
    mean_type_probability = type_probabilities.mean(axis=0)
    generalized_jsd = _entropy_bits(mean_type_probability) - _entropy_bits(type_probabilities).mean(axis=0)

    peak_data: dict[str, Any] = {
        "peak_id": reference.var_names.astype(str),
        "total_cut_sites": peak_totals.astype(np.int64),
        "nonzero_reference_cells": nonzero_cells.astype(np.int64),
        "shape_entropy_bits": peak_entropy,
        "between_type_generalized_jsd_bits": generalized_jsd,
        "cell_types_with_signal": type_has_signal.sum(axis=0).astype(np.int64),
    }
    for bin_index, layer_name in enumerate(layer_names):
        peak_data[f"count.{layer_name}"] = per_bin_peak[:, bin_index].astype(np.int64)
        peak_data[f"fraction.{layer_name}"] = peak_probabilities[:, bin_index]
    peak_table = pd.DataFrame(peak_data).set_index("peak_id", drop=False)

    per_bin_cell = np.column_stack([_sum_axis(matrix, axis=1) for matrix in matrices])
    cell_totals = per_bin_cell.sum(axis=1)
    cell_probabilities = _safe_probabilities(per_bin_cell)
    positive_cells = cell_totals > 0
    technical: dict[str, Any] = {"positive_count_cells": int(positive_cells.sum())}
    log_depth = np.full(reference.n_obs, np.nan, dtype=np.float64)
    log_depth[positive_cells] = np.log1p(cell_totals[positive_cells])
    for bin_index, layer_name in enumerate(layer_names):
        fractions = np.full(reference.n_obs, np.nan, dtype=np.float64)
        fractions[positive_cells] = cell_probabilities[positive_cells, bin_index]
        technical[f"spearman_log_depth_vs_fraction.{layer_name}"] = _finite_spearman(
            log_depth, fractions
        )

    global_bin_counts = per_bin_peak.sum(axis=0)
    global_bin_probabilities = _safe_probabilities(global_bin_counts[None, :])[0]
    positive_peaks = peak_totals > 0
    summary: dict[str, Any] = {
        "schema_version": 1,
        "scope": "training_reference_only",
        "split_seed": int(split_seed),
        "random_number_generator": {
            "numpy_version": np.__version__,
            "bit_generator": "PCG64",
            "seed_namespace": SEED_NAMESPACE,
            "split_half_stream_tag": SPLIT_HALF_STREAM,
            "seed_streams_by_cell_type": {
                cell_type: [
                    SEED_NAMESPACE,
                    int(split_seed),
                    SPLIT_HALF_STREAM,
                    type_index,
                ]
                for type_index, cell_type in enumerate(ordered_types, start=1)
            },
        },
        "labels_key": labels_key,
        "cell_types": list(ordered_types),
        "cells": int(reference.n_obs),
        "peaks": int(reference.n_vars),
        "positive_peaks": int(positive_peaks.sum()),
        "layer_names": list(layer_names),
        "pseudocount_for_type_comparisons": float(pseudocount),
        "feature_sha256": reference.uns["fragment_shape"]["feature_sha256"],
        "split_sha256": spec.split_sha256,
        "global_bin_counts": {
            layer_name: int(global_bin_counts[index])
            for index, layer_name in enumerate(layer_names)
        },
        "global_bin_fractions": {
            layer_name: float(global_bin_probabilities[index])
            for index, layer_name in enumerate(layer_names)
        },
        "global_shape_entropy_bits": float(_entropy_bits(global_bin_probabilities)),
        "peak_shape_entropy_bits": {
            "mean_positive_peaks": float(peak_entropy[positive_peaks].mean()),
            "median_positive_peaks": float(np.median(peak_entropy[positive_peaks])),
        },
        "between_type_generalized_jsd_bits": {
            "mean_positive_peaks": float(generalized_jsd[positive_peaks].mean()),
            "median_positive_peaks": float(np.median(generalized_jsd[positive_peaks])),
        },
        "split_half_reproducibility": _split_half_reproducibility(
            matrices,
            barcodes,
            labels,
            ordered_types,
            split_seed,
            pseudocount,
        ),
        "technical_confounding": technical,
        "interpretation": (
            "Diagnostics describe one donor and training-reference cells only; they do not "
            "establish donor-level reproducibility or causal nucleosome occupancy."
        ),
    }
    return peak_table, summary


def write_signal_audit(
    peak_table: pd.DataFrame,
    summary: dict[str, Any],
    output_dir: Path,
    *,
    overwrite: bool = False,
) -> tuple[Path, Path]:
    """Write audit outputs without replacing an existing audit by default."""
    output_dir = Path(output_dir)
    peak_path = output_dir / "signal_audit.csv"
    summary_path = output_dir / "signal_audit_summary.yaml"
    if not overwrite and (peak_path.exists() or summary_path.exists()):
        raise FileExistsError(f"Refusing to overwrite signal audit in {output_dir}.")
    output_dir.mkdir(parents=True, exist_ok=True)
    peak_table.to_csv(peak_path, index=False)
    enriched = dict(summary)
    enriched["files"] = {
        "signal_audit.csv": {
            "bytes": peak_path.stat().st_size,
            "sha256": _sha256_file(peak_path),
        }
    }
    with summary_path.open("w") as handle:
        yaml.safe_dump(enriched, handle, sort_keys=False)
    return peak_path, summary_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="Training-reference shape H5AD.")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--labels-key", default="cell_type")
    parser.add_argument("--split-seed", type=int, required=True)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    output_dir = args.output_dir or args.reference.parent
    reference = ad.read_h5ad(args.reference)
    cell_types = reference.uns.get("declared_cell_types")
    peak_table, summary = audit_shape_signal(
        reference,
        labels_key=args.labels_key,
        cell_types=cell_types,
        split_seed=args.split_seed,
        pseudocount=args.pseudocount,
    )
    peak_path, summary_path = write_signal_audit(
        peak_table,
        summary,
        output_dir,
        overwrite=args.overwrite,
    )
    print(f"wrote {peak_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
