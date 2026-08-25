#!/usr/bin/env python
"""Generate ShapeMix PBMC pseudo-spots from held-out single cells.

The primary simulator is deliberately layer-wise: source cells are aggregated
independently in every parent-fragment-length layer and ``X`` is reconstructed
as their exact sum.  All random streams are local PCG64 generators derived
from documented ``SeedSequence`` tuples; importing this module never changes
NumPy's process-global random state.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data import (  # noqa: E402
    FragmentShapeSpec,
    ordered_feature_sha256,
    validate_fragment_shape_spec,
)
from deconvatac.data.validators import (  # noqa: E402
    _validate_shape_object,
    validate_fragment_shape_feature_axis,
)
from scripts.shapemix_provenance import (  # noqa: E402
    code_provenance,
    fragment_shape_declaration,
    matrix_summary,
    software_versions,
)


SEED_NAMESPACE = 20260822
CONDITION_INDICES = {"equal_celltype": 0, "observed_abundance": 1}
PBMC_CELL_TYPE_COUNTS: tuple[tuple[str, int], ...] = (
    ("CD14 Mono", 2551),
    ("CD4 Naive", 1382),
    ("CD8 Naive", 1353),
    ("CD4 TCM", 1113),
    ("CD16 Mono", 442),
    ("NK", 403),
    ("CD8 TEM_1", 322),
    ("CD8 TEM_2", 315),
    ("Intermediate B", 300),
    ("Memory B", 298),
    ("CD4 TEM", 286),
    ("cDC", 180),
    ("Treg", 157),
    ("gdT", 143),
    ("MAIT", 130),
    ("Naive B", 125),
)
PBMC_CELL_TYPES = tuple(cell_type for cell_type, _ in PBMC_CELL_TYPE_COUNTS)
PRIMARY_OUTER_SPLIT_SEEDS = (1103, 2203, 3301, 4409, 5501)
PRIMARY_INNER_MIXTURE_SEEDS = (101, 211)
SMOKE_OUTER_SPLIT_SEED = 0
SMOKE_INNER_MIXTURE_SEED = 0
SMOKE_NUM_SPOTS = 32
SMOKE_MIN_FEATURES = 100
SMOKE_MAX_FEATURES = 500


@dataclass(frozen=True)
class ShapeMixtureSimulation:
    """In-memory result of one deterministic pseudo-spot simulation."""

    spatial: ad.AnnData
    truth: pd.DataFrame
    source_cells_by_spot: tuple[dict[str, Any], ...]
    cell_types: tuple[str, ...]
    sampling_probabilities: dict[str, float]
    condition: str
    outer_split_seed: int
    inner_mixture_seed: int
    seed_streams: dict[str, tuple[int, ...]]
    depth_retain_probability: Optional[float]


def _require_nonnegative_int(value: Any, name: str, *, positive: bool = False) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer.")
    value = int(value)
    minimum = 1 if positive else 0
    if value < minimum:
        qualifier = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be {qualifier}.")
    return value


def _rng(seed_tuple: tuple[int, ...]) -> np.random.Generator:
    """Construct the protocol's required local random generator."""
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed_tuple)))


def _canonical_csr(matrix: Any) -> sparse.csr_matrix:
    result = matrix.tocsr(copy=True) if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    result.sum_duplicates()
    result.eliminate_zeros()
    result.sort_indices()
    return result.astype(np.int64, copy=False)


def _count_matrix_is_valid(matrix: Any) -> bool:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    if not np.issubdtype(np.asarray(values).dtype, np.number):
        return False
    return bool(
        np.isfinite(values).all()
        and np.greater_equal(values, 0).all()
        and np.equal(values, np.floor(values)).all()
    )


def _matrix_equal(left: Any, right: sparse.csr_matrix) -> bool:
    observed = _canonical_csr(left)
    expected = _canonical_csr(right)
    difference = observed - expected
    difference.eliminate_zeros()
    return observed.shape == expected.shape and difference.nnz == 0


def _grid_shape(num_spots: int, grid_shape: Optional[Sequence[int]]) -> tuple[int, int]:
    if grid_shape is not None:
        if len(grid_shape) != 2:
            raise ValueError("grid_shape must contain exactly two dimensions.")
        rows = _require_nonnegative_int(grid_shape[0], "grid_shape[0]", positive=True)
        columns = _require_nonnegative_int(grid_shape[1], "grid_shape[1]", positive=True)
        if rows * columns != num_spots:
            raise ValueError("grid_shape product must equal num_spots.")
        return rows, columns

    # Use the most nearly square exact factorization.  This produces the frozen
    # 32 x 32 primary grid and a useful 4 x 8 grid for the 32-spot smoke path.
    rows = int(np.sqrt(num_spots))
    while num_spots % rows != 0:
        rows -= 1
    return rows, num_spots // rows


def _ordered_probabilities(
    cell_types: tuple[str, ...],
    probabilities: Mapping[str, float],
) -> np.ndarray:
    if set(probabilities) != set(cell_types):
        raise ValueError("Sampling probabilities must provide exactly the declared cell types.")
    values = np.asarray([probabilities[cell_type] for cell_type in cell_types], dtype=np.float64)
    if not np.isfinite(values).all() or (values < 0).any() or values.sum() <= 0:
        raise ValueError("Sampling probabilities must be finite, nonnegative, and have positive sum.")
    return values / values.sum()


def condition_probabilities(
    condition: str,
    cell_types: Sequence[str],
    observed_counts: Optional[Mapping[str, int]] = None,
) -> dict[str, float]:
    """Return ordered probabilities for one frozen benchmark condition."""
    ordered_types = tuple(str(cell_type) for cell_type in cell_types)
    if not ordered_types or len(set(ordered_types)) != len(ordered_types):
        raise ValueError("cell_types must be a non-empty, duplicate-free ordered sequence.")
    if condition == "equal_celltype":
        value = 1.0 / len(ordered_types)
        return {cell_type: value for cell_type in ordered_types}
    if condition != "observed_abundance":
        raise ValueError(f"Unsupported ShapeMix condition: {condition!r}.")

    counts = dict(PBMC_CELL_TYPE_COUNTS) if observed_counts is None else dict(observed_counts)
    if set(counts) != set(ordered_types):
        raise ValueError("Observed counts must provide exactly the declared cell types.")
    values = np.asarray([counts[cell_type] for cell_type in ordered_types], dtype=np.float64)
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("Observed cell-type counts must be finite and positive.")
    values /= values.sum()
    return {cell_type: float(value) for cell_type, value in zip(ordered_types, values)}


def _validate_heldout_cells(
    heldout_cells: ad.AnnData,
    cell_types: tuple[str, ...],
    labels_key: str,
    reference_barcodes: Sequence[str],
) -> tuple[FragmentShapeSpec, tuple[str, ...], np.ndarray]:
    if heldout_cells.n_obs == 0 or heldout_cells.n_vars == 0:
        raise ValueError("Held-out ShapeMix cells must have non-empty observation and feature axes.")
    if not heldout_cells.obs_names.is_unique or not heldout_cells.var_names.is_unique:
        raise ValueError("Held-out observation and feature names must be unique.")
    if labels_key not in heldout_cells.obs:
        raise KeyError(f"Held-out cells are missing obs[{labels_key!r}].")
    labels = heldout_cells.obs[labels_key]
    if labels.isna().any() or labels.astype(str).str.len().eq(0).any():
        raise ValueError(f"Held-out obs[{labels_key!r}] contains missing or empty labels.")
    labels = labels.astype(str)
    if set(labels) != set(cell_types):
        raise ValueError("Held-out labels must exactly match the declared cell-type universe.")

    heldout_names = tuple(str(name) for name in heldout_cells.obs_names)
    if len(set(heldout_names)) != len(heldout_names):
        raise ValueError("Held-out barcodes must remain unique after conversion to strings.")
    overlap = set(heldout_names).intersection(str(name) for name in reference_barcodes)
    if overlap:
        raise ValueError("Reference and held-out barcode pools must be disjoint.")

    metadata = heldout_cells.uns.get("fragment_shape")
    if not isinstance(metadata, Mapping):
        raise ValueError("Held-out cells are missing fragment_shape metadata.")
    spec = FragmentShapeSpec.from_mapping(metadata)
    validate_fragment_shape_spec(spec)
    _validate_shape_object(heldout_cells, spec, "heldout")
    if spec.split_sha256 is None:
        raise ValueError("Held-out fragment_shape metadata must declare split_sha256.")

    layer_sum = sparse.csr_matrix(heldout_cells.shape, dtype=np.int64)
    for layer_name in spec.layer_names:
        if layer_name not in heldout_cells.layers:
            raise ValueError(f"Held-out cells are missing fragment-shape layer {layer_name!r}.")
        layer = heldout_cells.layers[layer_name]
        if not sparse.isspmatrix_csr(layer):
            raise ValueError(f"Held-out layer {layer_name!r} must be CSR.")
        if not _count_matrix_is_valid(layer):
            raise ValueError(f"Held-out layer {layer_name!r} must contain nonnegative integer counts.")
        layer_sum += layer.astype(np.int64)
    if heldout_cells.X is None or not _count_matrix_is_valid(heldout_cells.X):
        raise ValueError("Held-out X must contain nonnegative integer counts.")
    if not _matrix_equal(heldout_cells.X, layer_sum):
        raise ValueError("Held-out X must equal the exact sum of fragment-shape layers.")

    canonical_names = tuple(sorted(heldout_names))
    positions = heldout_cells.obs_names.astype(str).get_indexer(canonical_names)
    if (positions < 0).any():  # Defensive: conversion collisions were checked above.
        raise RuntimeError("Could not restore canonical held-out barcode order.")
    return spec, canonical_names, positions


def _sample_sources(
    sampled_types: np.ndarray,
    pool_positions: Mapping[str, np.ndarray],
    source_rng: np.random.Generator,
) -> np.ndarray:
    """Sample cell positions, avoiding within-spot reuse whenever possible."""
    result = np.empty(sampled_types.size, dtype=np.int64)
    for cell_type in dict.fromkeys(sampled_types.tolist()):
        target_positions = np.flatnonzero(sampled_types == cell_type)
        pool = pool_positions[str(cell_type)]
        requested = target_positions.size
        if requested <= pool.size:
            sampled = source_rng.choice(pool, size=requested, replace=False)
        else:
            # Exhaust the complete pool once before reuse.  The remainder is
            # uniform with replacement, as required when demand exceeds supply.
            sampled = np.concatenate(
                [source_rng.permutation(pool), source_rng.choice(pool, size=requested - pool.size, replace=True)]
            )
        result[target_positions] = sampled
    return result


def _thin_matrix(
    matrix: sparse.csr_matrix,
    retain_probability: float,
    seed_tuple: tuple[int, ...],
) -> sparse.csr_matrix:
    result = matrix.copy()
    if result.nnz:
        result.data = _rng(seed_tuple).binomial(result.data.astype(np.int64), retain_probability)
        result.eliminate_zeros()
        result.sort_indices()
    return result


def _fragment_metadata_for_spatial(
    source_metadata: Mapping[str, Any],
    spatial: ad.AnnData,
    layer_names: tuple[str, ...],
) -> dict[str, Any]:
    metadata = copy.deepcopy(dict(source_metadata))
    metadata["feature_sha256"] = ordered_feature_sha256(spatial.var_names)
    layer_totals = {layer: int(spatial.layers[layer].sum()) for layer in layer_names}
    metadata["matrix_counters"] = {
        "assigned_cut_sites": int(sum(layer_totals.values())),
        **{f"cut_sites_per_bin.{layer}": total for layer, total in layer_totals.items()},
    }
    return metadata


def subset_shape_cells(
    adata: ad.AnnData,
    *,
    observation_mask: Optional[Sequence[bool]] = None,
    feature_names: Optional[Sequence[str]] = None,
) -> ad.AnnData:
    """Subset a shape object while refreshing only matrix-dependent metadata.

    Source-run preprocessing counters and all source/coordinate/split
    provenance are copied unchanged.  The ordered feature digest and matrix
    totals are the only fields whose meaning depends on the stored matrix.
    """
    metadata = adata.uns.get("fragment_shape")
    if not isinstance(metadata, Mapping):
        raise ValueError("Shape cells are missing fragment_shape metadata.")
    spec = FragmentShapeSpec.from_mapping(metadata)
    validate_fragment_shape_spec(spec)
    validate_fragment_shape_feature_axis(adata, "shape subset source")

    obs_selection: Any = slice(None)
    if observation_mask is not None:
        obs_selection = np.asarray(observation_mask, dtype=bool)
        if obs_selection.shape != (adata.n_obs,):
            raise ValueError("observation_mask must have one value per source observation.")
        if not obs_selection.any():
            raise ValueError("Shape subset cannot remove every observation.")

    var_selection: Any = slice(None)
    if feature_names is not None:
        ordered_features = tuple(str(feature) for feature in feature_names)
        if not ordered_features or len(set(ordered_features)) != len(ordered_features):
            raise ValueError("feature_names must be non-empty and duplicate-free.")
        missing = pd.Index(ordered_features).difference(adata.var_names)
        if len(missing):
            raise ValueError(f"Shape subset features are missing from the source: {missing[0]!r}.")
        var_selection = list(ordered_features)

    selected = adata[obs_selection, var_selection].copy()
    for layer_name in spec.layer_names:
        selected.layers[layer_name] = _canonical_csr(selected.layers[layer_name])
    selected.X = _canonical_csr(selected.X)
    selected.uns["fragment_shape"] = _fragment_metadata_for_spatial(
        metadata,
        selected,
        spec.layer_names,
    )
    return selected


def simulate_shapemix_spots(
    heldout_cells: ad.AnnData,
    *,
    cell_types: Sequence[str],
    sampling_probabilities: Mapping[str, float],
    condition: str,
    outer_split_seed: int,
    inner_mixture_seed: int,
    num_spots: int = 1024,
    mean_cells_per_spot: float = 10.0,
    labels_key: str = "cell_type",
    grid_shape: Optional[Sequence[int]] = None,
    reference_barcodes: Sequence[str] = (),
    depth_retain_probability: Optional[float] = None,
) -> ShapeMixtureSimulation:
    """Aggregate held-out cells into deterministic, layer-conserving spots.

    Depth thinning is deliberately isolated in its own seed stream.  Enabling
    it therefore changes matrices but never cell counts, sampled types, or
    source-cell provenance.  Primary datasets must leave it as ``None``.
    """
    if condition not in CONDITION_INDICES:
        raise ValueError(f"Unsupported ShapeMix condition: {condition!r}.")
    outer_split_seed = _require_nonnegative_int(outer_split_seed, "outer_split_seed")
    inner_mixture_seed = _require_nonnegative_int(inner_mixture_seed, "inner_mixture_seed")
    num_spots = _require_nonnegative_int(num_spots, "num_spots", positive=True)
    if isinstance(mean_cells_per_spot, (bool, np.bool_)):
        raise TypeError("mean_cells_per_spot must be numeric.")
    mean_cells_per_spot = float(mean_cells_per_spot)
    if not np.isfinite(mean_cells_per_spot) or mean_cells_per_spot <= 0:
        raise ValueError("mean_cells_per_spot must be finite and positive.")
    if depth_retain_probability is not None:
        if isinstance(depth_retain_probability, (bool, np.bool_)):
            raise TypeError("depth_retain_probability must be numeric.")
        depth_retain_probability = float(depth_retain_probability)
        if not np.isfinite(depth_retain_probability) or not 0 < depth_retain_probability < 1:
            raise ValueError("Secondary depth thinning requires 0 < retain probability < 1.")

    ordered_types = tuple(str(cell_type) for cell_type in cell_types)
    if not ordered_types or len(set(ordered_types)) != len(ordered_types):
        raise ValueError("cell_types must be a non-empty, duplicate-free ordered sequence.")
    probability_array = _ordered_probabilities(ordered_types, sampling_probabilities)
    shape_spec, sorted_barcodes, sorted_positions = _validate_heldout_cells(
        heldout_cells,
        ordered_types,
        labels_key,
        reference_barcodes,
    )
    rows, columns = _grid_shape(num_spots, grid_shape)

    condition_seed = (
        SEED_NAMESPACE,
        outer_split_seed,
        inner_mixture_seed,
        CONDITION_INDICES[condition],
    )
    seed_streams: dict[str, tuple[int, ...]] = {
        "cell_counts_and_types": condition_seed,
        "source_cells": condition_seed + (1,),
    }
    mixture_rng = _rng(seed_streams["cell_counts_and_types"])
    source_rng = _rng(seed_streams["source_cells"])

    sorted_labels = heldout_cells.obs.iloc[sorted_positions][labels_key].astype(str).to_numpy()
    pool_positions = {
        cell_type: sorted_positions[np.flatnonzero(sorted_labels == cell_type)]
        for cell_type in ordered_types
    }
    if any(pool.size == 0 for pool in pool_positions.values()):
        raise ValueError("Every declared cell type must have at least one held-out source cell.")

    cell_counts = np.maximum(1, mixture_rng.poisson(mean_cells_per_spot, size=num_spots)).astype(np.int64)
    layer_rows: dict[str, list[sparse.csr_matrix]] = {layer: [] for layer in shape_spec.layer_names}
    truth_counts = np.zeros((num_spots, len(ordered_types)), dtype=np.int64)
    provenance: list[dict[str, Any]] = []
    type_index = {cell_type: index for index, cell_type in enumerate(ordered_types)}

    for spot_index, cell_count in enumerate(cell_counts):
        sampled_types = mixture_rng.choice(
            np.asarray(ordered_types, dtype=object),
            size=int(cell_count),
            replace=True,
            p=probability_array,
        ).astype(str)
        source_positions = _sample_sources(sampled_types, pool_positions, source_rng)
        source_barcodes = heldout_cells.obs_names[source_positions].astype(str).tolist()
        for layer_name in shape_spec.layer_names:
            row = heldout_cells.layers[layer_name][source_positions].sum(axis=0)
            layer_rows[layer_name].append(_canonical_csr(row))
        for cell_type in sampled_types:
            truth_counts[spot_index, type_index[cell_type]] += 1
        provenance_row = {
            "spot_id": f"spot_{spot_index:04d}",
            "cell_count": int(cell_count),
            "cell_id": source_barcodes,
            "cell_type": sampled_types.tolist(),
        }
        for column in ("sample_key", "site", "donor", "author_cell_type"):
            if column in heldout_cells.obs:
                provenance_row[column] = (
                    heldout_cells.obs.iloc[source_positions][column].astype(str).tolist()
                )
        provenance.append(provenance_row)

    layers = {
        layer_name: sparse.vstack(rows_for_layer, format="csr", dtype=np.int64)
        for layer_name, rows_for_layer in layer_rows.items()
    }
    if depth_retain_probability is not None:
        for layer_index, layer_name in enumerate(shape_spec.layer_names):
            thinning_seed = condition_seed + (2, layer_index)
            seed_streams[f"depth_thinning.{layer_name}"] = thinning_seed
            layers[layer_name] = _thin_matrix(
                layers[layer_name],
                depth_retain_probability,
                thinning_seed,
            )

    x_matrix = sparse.csr_matrix((num_spots, heldout_cells.n_vars), dtype=np.int64)
    for layer_name in shape_spec.layer_names:
        x_matrix = _canonical_csr(x_matrix + layers[layer_name])
    spot_names = pd.Index([f"spot_{index:04d}" for index in range(num_spots)], name="spot_id")
    obs = pd.DataFrame(index=spot_names)
    obs["cell_count"] = cell_counts
    obs["simulation_condition"] = condition
    obs["outer_split_seed"] = outer_split_seed
    obs["inner_mixture_seed"] = inner_mixture_seed
    spatial = ad.AnnData(X=x_matrix, obs=obs, var=heldout_cells.var.copy())
    for layer_name in shape_spec.layer_names:
        spatial.layers[layer_name] = layers[layer_name]
    spatial.obsm["spatial"] = np.asarray(
        [(spot_index % columns, spot_index // columns) for spot_index in range(num_spots)],
        dtype=np.int64,
    )
    truth = pd.DataFrame(
        truth_counts / cell_counts[:, np.newaxis],
        index=spot_names,
        columns=pd.Index(ordered_types, name="cell_type"),
    )
    spatial.obsm["proportions"] = truth.to_numpy(copy=True)
    spatial.uns["proportion_names"] = np.asarray(ordered_types, dtype=str)
    spatial.uns["fragment_shape"] = _fragment_metadata_for_spatial(
        heldout_cells.uns["fragment_shape"],
        spatial,
        shape_spec.layer_names,
    )
    spatial.uns["simulation"] = {
        "schema_version": 1,
        "condition": condition,
        "seed_namespace": SEED_NAMESPACE,
        "outer_split_seed": outer_split_seed,
        "inner_mixture_seed": inner_mixture_seed,
        "num_spots": num_spots,
        "grid_rows": rows,
        "grid_columns": columns,
        "mean_cells_per_spot": mean_cells_per_spot,
        "depth_thinning_enabled": depth_retain_probability is not None,
        "depth_retain_probability": (
            "none" if depth_retain_probability is None else depth_retain_probability
        ),
    }

    return ShapeMixtureSimulation(
        spatial=spatial,
        truth=truth,
        source_cells_by_spot=tuple(provenance),
        cell_types=ordered_types,
        sampling_probabilities={
            cell_type: float(probability_array[index])
            for index, cell_type in enumerate(ordered_types)
        },
        condition=condition,
        outer_split_seed=outer_split_seed,
        inner_mixture_seed=inner_mixture_seed,
        seed_streams=seed_streams,
        depth_retain_probability=depth_retain_probability,
    )


def _repository_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT.resolve()))
    except ValueError:
        return str(resolved)


def _sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _manifest_output(
    split_manifest: Mapping[str, Any],
    role: str,
) -> Mapping[str, Any]:
    outputs = split_manifest.get("outputs")
    if not isinstance(outputs, list):
        raise ValueError("Split manifest outputs must be a list of hashed artifacts.")
    matches = [
        record
        for record in outputs
        if isinstance(record, Mapping) and record.get("role") == role
    ]
    if len(matches) != 1:
        raise ValueError(f"Split manifest must contain exactly one output with role {role!r}.")
    return matches[0]


def _validate_manifest_file(
    path: Path,
    split_dir: Path,
    record: Mapping[str, Any],
    role: str,
) -> None:
    try:
        expected_path = path.resolve().relative_to(split_dir.resolve()).as_posix()
    except ValueError as error:
        raise ValueError(f"Split artifact {role!r} is outside the split directory.") from error
    if record.get("path") != expected_path:
        raise ValueError(
            f"Split manifest path for {role!r} is {record.get('path')!r}, expected {expected_path!r}."
        )
    if isinstance(record.get("bytes"), (bool, np.bool_)) or record.get("bytes") != path.stat().st_size:
        raise ValueError(f"Split artifact byte count does not match the manifest for {role!r}.")
    expected_sha256 = record.get("sha256")
    if not isinstance(expected_sha256, str) or _sha256_file(path) != expected_sha256:
        raise ValueError(f"Split artifact SHA-256 does not match the manifest for {role!r}.")


def _validate_split_bundle(
    *,
    split_dir: Path,
    split_manifest: Mapping[str, Any],
    reference: ad.AnnData,
    heldout: ad.AnnData,
    reference_path: Path,
    heldout_path: Path,
    split_path: Path,
    selected_peaks_path: Path,
    labels_key: str,
    cell_types: tuple[str, ...],
) -> None:
    """Bind the supplied H5ADs to the manifest and canonical split membership."""
    role_paths = {
        "training_reference_fragment_shapes": reference_path,
        "heldout_source_fragment_shapes": heldout_path,
        "canonical_split_membership": split_path,
        "ranked_selected_peak_ids": selected_peaks_path,
    }
    for role, path in role_paths.items():
        _validate_manifest_file(
            path,
            split_dir,
            _manifest_output(split_manifest, role),
            role,
        )

    split_metadata = split_manifest.get("split")
    if not isinstance(split_metadata, Mapping):
        raise ValueError("Split manifest must declare split metadata.")
    observed_split_sha256 = _sha256_file(split_path)
    if split_metadata.get("sha256") != observed_split_sha256:
        raise ValueError("split.csv SHA-256 does not match split.sha256 in the manifest.")

    membership = pd.read_csv(split_path, dtype=str)
    if list(membership.columns) != ["barcode", "cell_type", "pool"]:
        raise ValueError("split.csv must contain exactly barcode, cell_type, and pool columns.")
    if membership.isna().any().any() or membership["barcode"].duplicated().any():
        raise ValueError("split.csv contains missing values or duplicate barcodes.")
    if set(membership["pool"]) != {"reference", "heldout"}:
        raise ValueError("split.csv must contain exactly reference and heldout pools.")
    if set(membership["cell_type"]) != set(cell_types):
        raise ValueError("split.csv cell types do not match the declared ordered universe.")

    expected_by_pool = {
        "reference": reference,
        "heldout": heldout,
    }
    for pool, adata in expected_by_pool.items():
        rows = membership[membership["pool"] == pool]
        observed_names = tuple(adata.obs_names.astype(str))
        if set(rows["barcode"]) != set(observed_names) or len(rows) != len(observed_names):
            raise ValueError(f"{pool} H5AD observations do not exactly match split.csv membership.")
        if labels_key not in adata.obs:
            raise KeyError(f"{pool} H5AD is missing obs[{labels_key!r}].")
        expected_labels = rows.set_index("barcode")["cell_type"].loc[list(observed_names)]
        if not np.array_equal(
            adata.obs[labels_key].astype(str).to_numpy(),
            expected_labels.to_numpy(),
        ):
            raise ValueError(f"{pool} H5AD labels do not match split.csv.")
        if "split_pool" not in adata.obs or set(adata.obs["split_pool"].astype(str)) != {pool}:
            raise ValueError(f"{pool} H5AD observations are not explicitly marked with their split pool.")
        preparation = adata.uns.get("shapemix_preparation")
        if not isinstance(preparation, Mapping) or preparation.get("pool") != pool:
            raise ValueError(f"{pool} H5AD preparation metadata does not identify its split pool.")
        if preparation.get("split_sha256") != observed_split_sha256:
            raise ValueError(f"{pool} H5AD preparation metadata has a stale split SHA-256.")

    if split_metadata.get("reference_cells") != reference.n_obs:
        raise ValueError("Split manifest reference cell count does not match reference H5AD.")
    if split_metadata.get("heldout_cells") != heldout.n_obs:
        raise ValueError("Split manifest held-out cell count does not match heldout H5AD.")
    counts = split_metadata.get("counts_by_cell_type")
    if not isinstance(counts, Mapping):
        raise ValueError("Split manifest must declare counts_by_cell_type.")
    expected_counts = {
        cell_type: {
            pool: int(
                ((membership["cell_type"] == cell_type) & (membership["pool"] == pool)).sum()
            )
            for pool in ("reference", "heldout")
        }
        for cell_type in cell_types
    }
    if counts != expected_counts:
        raise ValueError("Split manifest counts_by_cell_type does not match split.csv.")


def _canonical_json_line(value: Mapping[str, Any]) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def _plain_data(value: Any) -> Any:
    """Convert NumPy/HDF5-returned scalars and arrays to YAML-safe values."""
    if isinstance(value, Mapping):
        return {str(key): _plain_data(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return [_plain_data(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_plain_data(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _shape_declaration(metadata: Mapping[str, Any]) -> dict[str, Any]:
    spec = FragmentShapeSpec.from_mapping(metadata)
    return {
        "schema_version": spec.schema_version,
        "axis": spec.axis,
        "count_unit": spec.count_unit,
        "read_support_policy": spec.read_support_policy,
        "peak_assignment": spec.peak_assignment,
        "bins": [bin_spec.to_dict() for bin_spec in spec.bins],
    }


def dataset_id_for_simulation(
    condition: str,
    outer_split_seed: int,
    inner_mixture_seed: int,
    depth_retain_probability: Optional[float] = None,
    *,
    smoke: bool = False,
) -> str:
    """Build the collision-free primary or secondary dataset ID."""
    if condition not in CONDITION_INDICES:
        raise ValueError(f"Unsupported ShapeMix condition: {condition!r}.")
    dataset_id = (
        "pbmc_granulocyte_sorted_10k_shapemix_"
        f"{condition}_split_{int(outer_split_seed):03d}_mix_{int(inner_mixture_seed):03d}"
    )
    if depth_retain_probability is not None:
        token = f"{float(depth_retain_probability):.6f}".rstrip("0").rstrip(".").replace(".", "p")
        dataset_id += f"_depth_keep_{token}"
    if smoke:
        dataset_id += "_smoke"
    return dataset_id


def is_primary_simulation(simulation: ShapeMixtureSimulation) -> bool:
    """Return whether every frozen primary-dataset condition is satisfied."""
    metadata = simulation.spatial.uns["simulation"]
    if simulation.condition not in CONDITION_INDICES or set(
        simulation.sampling_probabilities
    ) != set(PBMC_CELL_TYPES):
        return False
    expected_probabilities = condition_probabilities(
        simulation.condition,
        PBMC_CELL_TYPES,
    )
    observed_probability_array = np.asarray(
        [simulation.sampling_probabilities[cell_type] for cell_type in PBMC_CELL_TYPES]
    )
    expected_probability_array = np.asarray(
        [expected_probabilities[cell_type] for cell_type in PBMC_CELL_TYPES]
    )
    return bool(
        simulation.cell_types == PBMC_CELL_TYPES
        and simulation.outer_split_seed in PRIMARY_OUTER_SPLIT_SEEDS
        and simulation.inner_mixture_seed in PRIMARY_INNER_MIXTURE_SEEDS
        and simulation.spatial.n_obs == 1024
        and simulation.spatial.n_vars == 5000
        and int(metadata["grid_rows"]) == 32
        and int(metadata["grid_columns"]) == 32
        and float(metadata["mean_cells_per_spot"]) == 10.0
        and simulation.depth_retain_probability is None
        and np.allclose(
            observed_probability_array,
            expected_probability_array,
            rtol=0.0,
            atol=1.0e-15,
        )
    )


def _output_record(path: Path, temporary_root: Path, final_root: Path) -> dict[str, Any]:
    final_path = final_root / path.relative_to(temporary_root)
    return {
        "path": _repository_path(final_path),
        "bytes": path.stat().st_size,
        "sha256": _sha256_file(path),
    }


def write_simulation_dataset(
    simulation: ShapeMixtureSimulation,
    *,
    output_root: Path,
    dataset_id: str,
    reference_path: Path,
    heldout_path: Path,
    split_manifest_path: Path,
    labels_key: str = "cell_type",
    dataset_reference: Optional[ad.AnnData] = None,
    benchmark_scope: Optional[str] = None,
    source: str = "pbmc_granulocyte_sorted_10k_shapemix_heldout_simulation",
    description: Optional[str] = None,
    scientific_scope: Optional[str] = None,
) -> Path:
    """Atomically write one complete dataset and refuse every overwrite."""
    output_root = Path(output_root)
    final_dir = output_root / dataset_id
    if final_dir.exists():
        raise FileExistsError(f"{final_dir} already exists; ShapeMix outputs are never overwritten.")
    for source_path in (reference_path, heldout_path, split_manifest_path):
        if not Path(source_path).is_file():
            raise FileNotFoundError(source_path)
    output_root.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(tempfile.mkdtemp(prefix=f".{dataset_id}.", dir=output_root))

    try:
        atac_dir = temporary_dir / "atac"
        feature_dir = atac_dir / "features"
        truth_dir = temporary_dir / "truth"
        simulation_dir = temporary_dir / "simulation"
        for directory in (feature_dir, truth_dir, simulation_dir):
            directory.mkdir(parents=True, exist_ok=True)

        spatial_path = atac_dir / "spatial.h5ad"
        derived_reference_path = atac_dir / "reference_cells.h5ad"
        feature_path = feature_dir / "highly_variable.txt"
        truth_path = truth_dir / "proportions.csv"
        provenance_path = simulation_dir / "source_cells_by_spot.jsonl"
        dataset_path = temporary_dir / "dataset.yaml"
        manifest_path = simulation_dir / "manifest.yaml"

        simulation.spatial.write_h5ad(spatial_path, compression="gzip")
        if dataset_reference is not None:
            if not dataset_reference.var_names.equals(simulation.spatial.var_names):
                raise ValueError(
                    "A dataset-local reference must match the simulated ordered feature axis."
                )
            dataset_reference.write_h5ad(derived_reference_path, compression="gzip")
        feature_path.write_text("\n".join(map(str, simulation.spatial.var_names)) + "\n")
        simulation.truth.to_csv(truth_path, index_label="spot_id")
        with provenance_path.open("w", encoding="utf-8", newline="\n") as handle:
            for row in simulation.source_cells_by_spot:
                handle.write(_canonical_json_line(row) + "\n")

        final_spatial = final_dir / "atac" / "spatial.h5ad"
        final_reference = final_dir / "atac" / "reference_cells.h5ad"
        final_feature = final_dir / "atac" / "features" / "highly_variable.txt"
        final_truth = final_dir / "truth" / "proportions.csv"
        final_provenance = final_dir / "simulation" / "source_cells_by_spot.jsonl"
        final_manifest = final_dir / "simulation" / "manifest.yaml"
        fragment_shape = _shape_declaration(simulation.spatial.uns["fragment_shape"])
        configured_reference_path = (
            final_reference if dataset_reference is not None else Path(reference_path)
        )
        frozen_primary = is_primary_simulation(simulation)
        if benchmark_scope is None:
            if simulation.depth_retain_probability is not None:
                benchmark_scope = "secondary_depth_sensitivity"
            elif dataset_id.endswith("_smoke"):
                benchmark_scope = "smoke"
            elif frozen_primary:
                benchmark_scope = "primary"
            else:
                benchmark_scope = "development"
        if benchmark_scope == "primary" and not frozen_primary:
            raise ValueError("A primary label requires every frozen primary simulation condition.")
        primary_dataset = benchmark_scope == "primary" and frozen_primary

        dataset_config = {
            "dataset_id": dataset_id,
            "source": source,
            "description": description or (
                "ShapeMix PBMC pseudo-spots aggregated exclusively from held-out cells; "
                "this one-donor benchmark measures conditional resampling variability, "
                "not donor-level or population-level generalization."
            ),
            "labels_key": labels_key,
            "spatial_key": "spatial",
            "benchmark_scope": benchmark_scope,
            "simulation": {
                "condition": simulation.condition,
                "outer_split_seed": simulation.outer_split_seed,
                "inner_mixture_seed": simulation.inner_mixture_seed,
                "source_cells_by_spot": _repository_path(final_provenance),
                "manifest": _repository_path(final_manifest),
                "primary_dataset": primary_dataset,
                "depth_retain_probability": simulation.depth_retain_probability,
            },
            "modalities": {
                "atac": {
                    "reference": {"path": _repository_path(configured_reference_path)},
                    "spatial": {"path": _repository_path(final_spatial)},
                    "labels_key": labels_key,
                    "spatial_key": "spatial",
                    "truth": {
                        "path": _repository_path(final_truth),
                        "cell_types": list(simulation.cell_types),
                    },
                    "fragment_shape": fragment_shape,
                    "feature_sets": {
                        "all": {"mode": "all"},
                        "highly_variable": {"path": _repository_path(final_feature)},
                    },
                }
            },
        }
        with dataset_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(_plain_data(dataset_config), handle, sort_keys=False)

        output_files = {
            "spatial_h5ad": _output_record(spatial_path, temporary_dir, final_dir),
            "selected_peaks": _output_record(feature_path, temporary_dir, final_dir),
            "truth_proportions": _output_record(truth_path, temporary_dir, final_dir),
            "source_cells_by_spot": _output_record(provenance_path, temporary_dir, final_dir),
            "dataset_yaml": _output_record(dataset_path, temporary_dir, final_dir),
        }
        if dataset_reference is not None:
            output_files["reference_h5ad"] = _output_record(
                derived_reference_path,
                temporary_dir,
                final_dir,
            )
        manifest = {
            "schema_version": 1,
            "dataset_id": dataset_id,
            "status": "complete",
            "benchmark_scope": benchmark_scope,
            "scientific_scope": scientific_scope or (
                "Conditional resampling variability within one PBMC Multiome donor; "
                "not donor-level or biological generalization uncertainty."
            ),
            "code_provenance": code_provenance(),
            "software_versions": software_versions(),
            "inputs": {
                "reference_h5ad": {
                    "path": _repository_path(Path(reference_path)),
                    "bytes": Path(reference_path).stat().st_size,
                    "sha256": _sha256_file(Path(reference_path)),
                },
                "heldout_h5ad": {
                    "path": _repository_path(Path(heldout_path)),
                    "bytes": Path(heldout_path).stat().st_size,
                    "sha256": _sha256_file(Path(heldout_path)),
                },
                "split_manifest": {
                    "path": _repository_path(Path(split_manifest_path)),
                    "bytes": Path(split_manifest_path).stat().st_size,
                    "sha256": _sha256_file(Path(split_manifest_path)),
                },
            },
            "fragment_shape": {
                **fragment_shape_declaration(simulation.spatial.uns["fragment_shape"]),
                "split_sha256": simulation.spatial.uns["fragment_shape"]["split_sha256"],
                "source_sha256": copy.deepcopy(
                    simulation.spatial.uns["fragment_shape"]["source_sha256"]
                ),
                "feature_sha256": simulation.spatial.uns["fragment_shape"]["feature_sha256"],
                "matrix_counters": copy.deepcopy(
                    simulation.spatial.uns["fragment_shape"]["matrix_counters"]
                ),
            },
            "simulation": {
                "condition": simulation.condition,
                "num_spots": simulation.spatial.n_obs,
                "num_features": simulation.spatial.n_vars,
                "grid_shape": [
                    int(simulation.spatial.uns["simulation"]["grid_rows"]),
                    int(simulation.spatial.uns["simulation"]["grid_columns"]),
                ],
                "mean_cells_per_spot": float(
                    simulation.spatial.uns["simulation"]["mean_cells_per_spot"]
                ),
                "sampling_probabilities": simulation.sampling_probabilities,
                "cell_types": list(simulation.cell_types),
                "random_number_generator": {
                    "numpy_version": np.__version__,
                    "bit_generator": "PCG64",
                    "seed_namespace": SEED_NAMESPACE,
                    "outer_split_seed": simulation.outer_split_seed,
                    "inner_mixture_seed": simulation.inner_mixture_seed,
                    "seed_streams": {
                        key: list(value) for key, value in simulation.seed_streams.items()
                    },
                },
                "depth_thinning": {
                    "enabled": simulation.depth_retain_probability is not None,
                    "retain_probability": simulation.depth_retain_probability,
                    "primary_dataset": primary_dataset,
                },
                "source_sampling": {
                    "heldout_only": True,
                    "without_replacement_within_spot_when_possible": True,
                    "reuse_across_spots": True,
                },
                "matrix": matrix_summary(
                    simulation.spatial,
                    FragmentShapeSpec.from_mapping(
                        simulation.spatial.uns["fragment_shape"]
                    ).layer_names,
                ),
            },
            "checks": {
                "x_equals_layer_sum": True,
                "truth_rows_sum_to_one": True,
                "reference_heldout_disjoint": True,
                "ordered_feature_hash_matches": True,
            },
            "outputs": output_files,
        }
        with manifest_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(_plain_data(manifest), handle, sort_keys=False)

        # A single same-filesystem rename publishes the complete directory.
        # Path.rename refuses to replace a populated directory on supported
        # platforms, preserving the no-overwrite contract after the preflight.
        temporary_dir.rename(final_dir)
    except BaseException:
        if temporary_dir.exists():
            shutil.rmtree(temporary_dir)
        raise
    return final_dir / "dataset.yaml"


def _validate_reference_and_heldout(
    reference: ad.AnnData,
    heldout: ad.AnnData,
    cell_types: tuple[str, ...],
    labels_key: str,
) -> None:
    if not reference.var_names.equals(heldout.var_names):
        raise ValueError("Reference and held-out objects must have identical ordered features.")
    if labels_key not in reference.obs:
        raise KeyError(f"Reference cells are missing obs[{labels_key!r}].")
    reference_labels = reference.obs[labels_key]
    if reference_labels.isna().any() or set(reference_labels.astype(str)) != set(cell_types):
        raise ValueError("Reference labels must exactly match the declared cell-type universe.")
    overlap = set(reference.obs_names.astype(str)).intersection(heldout.obs_names.astype(str))
    if overlap:
        raise ValueError("Reference and held-out barcode pools must be disjoint.")

    reference_metadata = reference.uns.get("fragment_shape")
    heldout_metadata = heldout.uns.get("fragment_shape")
    if not isinstance(reference_metadata, Mapping) or not isinstance(heldout_metadata, Mapping):
        raise ValueError("Reference and held-out objects must contain fragment_shape metadata.")
    reference_spec = FragmentShapeSpec.from_mapping(reference_metadata)
    heldout_spec = FragmentShapeSpec.from_mapping(heldout_metadata)
    _validate_shape_object(reference, reference_spec, "reference")
    _validate_shape_object(heldout, heldout_spec, "heldout")
    aligned_fields = (
        "schema_version",
        "axis",
        "count_unit",
        "read_support_policy",
        "peak_assignment",
        "bins",
        "left_cut_offset",
        "right_cut_offset",
        "source_sha256",
        "feature_sha256",
        "split_sha256",
        "coordinate_validation",
        "software_versions",
    )
    for field in aligned_fields:
        if getattr(reference_spec, field) != getattr(heldout_spec, field):
            raise ValueError(f"Reference and held-out fragment_shape metadata disagree for {field!r}.")


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regenerate ShapeMix PBMC pseudo-spots from split-specific held-out cells."
    )
    parser.add_argument("--split-dir", type=Path, required=True)
    parser.add_argument("--outer-split-seed", type=int, required=True)
    parser.add_argument("--inner-mixture-seed", type=int, required=True)
    parser.add_argument(
        "--conditions",
        nargs="+",
        choices=tuple(CONDITION_INDICES),
        default=list(CONDITION_INDICES),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "data" / "processed" / "datasets",
    )
    parser.add_argument("--num-spots", type=int, default=1024)
    parser.add_argument("--mean-cells-per-spot", type=float, default=10.0)
    parser.add_argument("--grid-shape", type=int, nargs=2)
    parser.add_argument("--labels-key", default="cell_type")
    parser.add_argument("--cell-types", nargs="+", default=None)
    parser.add_argument(
        "--num-features",
        type=int,
        default=None,
        help="Use the first N ranked selected peaks and write a dataset-local sliced reference.",
    )
    parser.add_argument(
        "--depth-retain-probability",
        type=float,
        default=None,
        help="Secondary sensitivity analysis only; primary datasets never use depth thinning.",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Mark a development smoke dataset explicitly in its ID and metadata.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> list[Path]:
    args = _parse_args(argv)
    split_dir = args.split_dir.resolve()
    reference_path = split_dir / "reference_cells.h5ad"
    heldout_path = split_dir / "heldout_test_cells.h5ad"
    split_manifest_path = split_dir / "manifest.yaml"
    split_path = split_dir / "split.csv"
    selected_peaks_path = split_dir / "selected_peaks.txt"
    for required_path in (
        reference_path,
        heldout_path,
        split_manifest_path,
        split_path,
        selected_peaks_path,
    ):
        if not required_path.is_file():
            raise FileNotFoundError(required_path)

    reference = ad.read_h5ad(reference_path)
    heldout = ad.read_h5ad(heldout_path)
    cell_types = PBMC_CELL_TYPES if args.cell_types is None else tuple(args.cell_types)
    with split_manifest_path.open() as handle:
        split_manifest = yaml.safe_load(handle) or {}
    if split_manifest.get("cell_types") != list(cell_types):
        raise ValueError("Split manifest cell_types must exactly match the requested ordered universe.")
    if (split_manifest.get("rng") or {}).get("outer_split_seed") != args.outer_split_seed:
        raise ValueError("Split manifest rng.outer_split_seed does not match the CLI seed.")
    manifest_split_hash = (split_manifest.get("split") or {}).get("sha256")
    if not isinstance(manifest_split_hash, str):
        raise ValueError("Split manifest must declare split.sha256.")
    manifest_scope = split_manifest.get("benchmark_scope")
    if not isinstance(manifest_scope, str) or not manifest_scope:
        raise ValueError("Split manifest must declare benchmark_scope.")
    if args.smoke and manifest_scope != "development_smoke":
        raise ValueError(
            "--smoke requires a split manifest with benchmark_scope: development_smoke."
        )
    if not args.smoke and manifest_scope == "development_smoke":
        raise ValueError("A smoke split must be generated with --smoke.")
    if args.smoke:
        if args.outer_split_seed != SMOKE_OUTER_SPLIT_SEED or args.inner_mixture_seed != SMOKE_INNER_MIXTURE_SEED:
            raise ValueError("The smoke dataset requires development outer/inner seeds 0/0.")
        if args.num_spots != SMOKE_NUM_SPOTS:
            raise ValueError("The smoke dataset requires exactly 32 spots.")
        if len(cell_types) not in {2, 3}:
            raise ValueError("The smoke dataset requires exactly two or three declared cell types.")
        if args.depth_retain_probability is not None:
            raise ValueError("Depth thinning is a separate sensitivity analysis, not a smoke dataset setting.")
    _validate_split_bundle(
        split_dir=split_dir,
        split_manifest=split_manifest,
        reference=reference,
        heldout=heldout,
        reference_path=reference_path,
        heldout_path=heldout_path,
        split_path=split_path,
        selected_peaks_path=selected_peaks_path,
        labels_key=args.labels_key,
        cell_types=cell_types,
    )
    selected_peaks = [line.strip() for line in selected_peaks_path.read_text().splitlines() if line.strip()]
    if selected_peaks != heldout.var_names.astype(str).tolist():
        raise ValueError("selected_peaks.txt must exactly match the held-out ordered feature axis.")
    if not reference.var_names.equals(heldout.var_names):
        raise ValueError("Reference and held-out objects must have identical ordered features.")

    dataset_local_reference = False
    if args.cell_types is not None:
        if not cell_types or len(set(cell_types)) != len(cell_types):
            raise ValueError("--cell-types must be non-empty, duplicate-free, and ordered.")
        for role, adata in (("Reference", reference), ("Held-out", heldout)):
            if args.labels_key not in adata.obs:
                raise KeyError(f"{role} cells are missing obs[{args.labels_key!r}].")
            missing_types = set(cell_types).difference(adata.obs[args.labels_key].astype(str))
            if missing_types:
                raise ValueError(f"{role} cells are missing requested types: {sorted(missing_types)}.")
        reference = subset_shape_cells(
            reference,
            observation_mask=reference.obs[args.labels_key].astype(str).isin(cell_types),
        )
        heldout = subset_shape_cells(
            heldout,
            observation_mask=heldout.obs[args.labels_key].astype(str).isin(cell_types),
        )
        dataset_local_reference = True

    if args.num_features is not None:
        num_features = _require_nonnegative_int(args.num_features, "num_features", positive=True)
        if num_features > len(selected_peaks):
            raise ValueError("--num-features cannot exceed the ranked selected-peak count.")
        selected_peaks = selected_peaks[:num_features]
        reference = subset_shape_cells(reference, feature_names=selected_peaks)
        heldout = subset_shape_cells(heldout, feature_names=selected_peaks)
        dataset_local_reference = True

    if args.smoke and not SMOKE_MIN_FEATURES <= heldout.n_vars <= SMOKE_MAX_FEATURES:
        raise ValueError("The smoke dataset requires between 100 and 500 selected peaks.")

    _validate_reference_and_heldout(reference, heldout, cell_types, args.labels_key)
    for role, adata in (("reference", reference), ("held-out", heldout)):
        observed_split_hash = FragmentShapeSpec.from_mapping(
            adata.uns["fragment_shape"]
        ).split_sha256
        if observed_split_hash != manifest_split_hash:
            raise ValueError(
                f"The {role} fragment_shape.split_sha256 does not match split manifest split.sha256."
            )

    observed_counts = dict(PBMC_CELL_TYPE_COUNTS)
    if set(cell_types) != set(PBMC_CELL_TYPES):
        # A smoke-specific universe uses its actual pre-split abundance only;
        # scientific runs use the protocol-frozen counts above.
        combined_labels = pd.concat(
            [reference.obs[args.labels_key], heldout.obs[args.labels_key]]
        ).astype(str)
        observed_counts = {
            cell_type: int((combined_labels == cell_type).sum()) for cell_type in cell_types
        }

    targets = []
    for condition in args.conditions:
        dataset_id = dataset_id_for_simulation(
            condition,
            args.outer_split_seed,
            args.inner_mixture_seed,
            args.depth_retain_probability,
            smoke=args.smoke,
        )
        target = args.output_root / dataset_id
        if target.exists():
            raise FileExistsError(f"{target} already exists; no datasets were generated.")
        targets.append((condition, dataset_id))

    written = []
    reference_barcodes = reference.obs_names.astype(str).tolist()
    for condition, dataset_id in targets:
        probabilities = condition_probabilities(condition, cell_types, observed_counts)
        result = simulate_shapemix_spots(
            heldout,
            cell_types=cell_types,
            sampling_probabilities=probabilities,
            condition=condition,
            outer_split_seed=args.outer_split_seed,
            inner_mixture_seed=args.inner_mixture_seed,
            num_spots=args.num_spots,
            mean_cells_per_spot=args.mean_cells_per_spot,
            labels_key=args.labels_key,
            grid_shape=args.grid_shape,
            reference_barcodes=reference_barcodes,
            depth_retain_probability=args.depth_retain_probability,
        )
        config_path = write_simulation_dataset(
            result,
            output_root=args.output_root,
            dataset_id=dataset_id,
            reference_path=reference_path,
            heldout_path=heldout_path,
            split_manifest_path=split_manifest_path,
            labels_key=args.labels_key,
            dataset_reference=reference if dataset_local_reference else None,
            benchmark_scope=(
                "secondary_depth_sensitivity"
                if args.depth_retain_probability is not None
                else ("development_smoke" if args.smoke else None)
            ),
        )
        print(f"wrote {_repository_path(config_path)}", flush=True)
        written.append(config_path)
    return written


if __name__ == "__main__":
    main()
