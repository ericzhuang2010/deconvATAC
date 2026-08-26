#!/usr/bin/env python
"""Build leakage-free GSE194122 donor-held-out ShapeMix benchmarks."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import platform
import tempfile
from pathlib import Path
from typing import Any, Iterable, Mapping

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import yaml

from deconvatac.data import FragmentShapeSpec, ordered_feature_sha256
from deconvatac.data.validators import (
    validate_fragment_shape_feature_axis,
    validate_fragment_shape_spec,
)
from deconvatac.pp.fragment_shapes import (
    build_fragment_shape_anndata,
    count_fragment_shapes,
)
from scripts.regenerate_shapemix_pbmc_simulations import (
    condition_probabilities,
    simulate_shapemix_spots,
    write_simulation_dataset,
)


ROOT = Path(__file__).resolve().parents[1]
FAMILY_ROOT = ROOT / "data" / "processed" / "shapemix" / "gse194122_bmmc"
SOURCE_H5AD = (
    FAMILY_ROOT
    / "source_audit"
    / "source_objects"
    / "GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
)
LABELS_PATH = FAMILY_ROOT / "labels" / "source_broad7_v1" / "cells.tsv.gz"
FEATURES_PATH = FAMILY_ROOT / "feature_axes" / "source_axis_v1" / "features.tsv.gz"
FRAGMENT_MANIFEST = FAMILY_ROOT / "source_audit" / "fragment_suite.yaml"
SOURCE_LOCK = ROOT / "configs" / "data_sources" / "shapemix_gse194122_lock.yaml"
AXIS_ROOT = FAMILY_ROOT / "feature_axes" / "broad7_lodo_v1"
CACHE_ROOT = FAMILY_ROOT / "fragment_shape_cache" / "broad7_union_v1"
FOLD_ROOT = FAMILY_ROOT / "splits" / "broad7_lodo_v1"
REFERENCE_ROOT = ROOT / "data" / "processed" / "references"
DATASET_ROOT = ROOT / "data" / "processed" / "datasets"
REGISTRY_PATH = ROOT / "data" / "registry" / "datasets.yaml"

DONORS = tuple(range(1, 11))
INNER_MIXTURE_SEEDS = (101, 211)
CONDITIONS = ("observed_abundance", "equal_celltype")
CELL_TYPES = (
    "B/plasma",
    "CD4 T",
    "CD8 T",
    "NK/ILC",
    "Myeloid/DC",
    "Erythroid/MK-E",
    "Hematopoietic progenitor",
)
N_TOP_PEAKS = 5_000
MIN_REFERENCE_CELLS = 10
MIN_TRAINING_TYPE_CELLS = 100
MIN_HELDOUT_TYPE_CELLS = 20
COUNT_CHUNK_SIZE = 1_000_000
MATRIX_CELL_CHUNK_SIZE = 512

BROAD_LABELS = {
    "B1 B": "B/plasma",
    "Naive CD20+ B": "B/plasma",
    "Transitional B": "B/plasma",
    "Plasma cell": "B/plasma",
    "CD4+ T activated": "CD4 T",
    "CD4+ T naive": "CD4 T",
    "CD8+ T": "CD8 T",
    "CD8+ T naive": "CD8 T",
    "NK": "NK/ILC",
    "ILC": "NK/ILC",
    "CD14+ Mono": "Myeloid/DC",
    "CD16+ Mono": "Myeloid/DC",
    "cDC2": "Myeloid/DC",
    "pDC": "Myeloid/DC",
    "Proerythroblast": "Erythroid/MK-E",
    "Erythroblast": "Erythroid/MK-E",
    "Normoblast": "Erythroid/MK-E",
    "MK/E prog": "Erythroid/MK-E",
    "HSC": "Hematopoietic progenitor",
    "Lymph prog": "Hematopoietic progenitor",
    "G/M prog": "Hematopoietic progenitor",
    "ID2-hi myeloid prog": "Hematopoietic progenitor",
}


def _repository_path(path: Path) -> str:
    try:
        return path.absolute().relative_to(ROOT.absolute()).as_posix()
    except ValueError:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a mapping.")
    return value


def _write_yaml(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w") as handle:
        yaml.safe_dump(dict(value), handle, sort_keys=False)
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _software_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": __import__("scipy").__version__,
        "anndata": ad.__version__,
    }


def _load_labels() -> pd.DataFrame:
    labels = pd.read_csv(LABELS_PATH, sep="\t")
    expected = {
        "cell_id",
        "cell_index",
        "sample_key",
        "site",
        "donor",
        "author_cell_type",
        "fragment_barcode",
    }
    missing = expected.difference(labels.columns)
    if missing:
        raise ValueError(f"Label table is missing columns: {sorted(missing)!r}.")
    if not np.array_equal(labels["cell_index"].to_numpy(), np.arange(len(labels))):
        raise ValueError("Label table no longer preserves the source H5AD cell order.")
    if labels["cell_id"].duplicated().any():
        raise ValueError("Source cell IDs must be globally unique.")
    if set(labels["author_cell_type"]) != set(BROAD_LABELS):
        raise ValueError("The frozen 22-label author universe has changed.")
    labels["cell_type"] = labels["author_cell_type"].map(BROAD_LABELS)
    labels["donor"] = labels["donor"].astype(int)
    if set(labels["donor"]) != set(DONORS):
        raise ValueError("The frozen ten-donor universe has changed.")
    for donor in DONORS:
        heldout = labels[labels["donor"] == donor]["cell_type"].value_counts()
        training = labels[labels["donor"] != donor]["cell_type"].value_counts()
        for cell_type in CELL_TYPES:
            if int(training.get(cell_type, 0)) < MIN_TRAINING_TYPE_CELLS:
                raise ValueError(
                    f"Donor {donor} has insufficient training support for {cell_type}."
                )
            if int(heldout.get(cell_type, 0)) < MIN_HELDOUT_TYPE_CELLS:
                raise ValueError(
                    f"Donor {donor} has insufficient held-out support for {cell_type}."
                )
    return labels


def _candidate_features() -> pd.DataFrame:
    features = pd.read_csv(FEATURES_PATH, sep="\t")
    required = {"feature_index", "feature_id", "feature_type", "chromosome", "start", "end"}
    missing = required.difference(features.columns)
    if missing:
        raise ValueError(f"Feature table is missing columns: {sorted(missing)!r}.")
    canonical_contigs = {f"chr{value}" for value in range(1, 23)} | {"chrX", "chrY"}
    selected = features[
        (features["feature_type"] == "ATAC")
        & features["chromosome"].isin(canonical_contigs)
    ].copy()
    selected["feature_index"] = selected["feature_index"].astype(int)
    selected["start"] = selected["start"].astype(int)
    selected["end"] = selected["end"].astype(int)
    if selected["feature_id"].duplicated().any() or selected.empty:
        raise ValueError("Candidate peak IDs must be non-empty and unique.")
    return selected.reset_index(drop=True)


def _update_dense_from_sparse(target: np.ndarray, values: sparse.spmatrix) -> None:
    coo = values.tocoo()
    target[coo.row, coo.col] += coo.data.astype(target.dtype, copy=False)


def _aggregate_donor_type_counts(
    labels: pd.DataFrame,
    candidates: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    n_groups = len(DONORS) * len(CELL_TYPES)
    n_peaks = len(candidates)
    summed = np.zeros((n_groups, n_peaks), dtype=np.int64)
    coverage = np.zeros((n_groups, n_peaks), dtype=np.int32)
    donor_index = {donor: index for index, donor in enumerate(DONORS)}
    type_index = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    groups = np.asarray(
        [
            donor_index[int(donor)] * len(CELL_TYPES) + type_index[str(cell_type)]
            for donor, cell_type in zip(labels["donor"], labels["cell_type"])
        ],
        dtype=np.int64,
    )
    feature_indices = candidates["feature_index"].to_numpy(dtype=np.int64)

    source = ad.read_h5ad(SOURCE_H5AD, backed="r")
    try:
        if source.n_obs != len(labels):
            raise ValueError("Source H5AD and label-table cell axes differ.")
        for start in range(0, source.n_obs, MATRIX_CELL_CHUNK_SIZE):
            stop = min(start + MATRIX_CELL_CHUNK_SIZE, source.n_obs)
            chunk = sparse.csr_matrix(source.layers["counts"][start:stop, :])
            chunk = chunk[:, feature_indices]
            if chunk.nnz and (
                not np.isfinite(chunk.data).all()
                or (chunk.data < 0).any()
                or not np.equal(chunk.data, np.floor(chunk.data)).all()
            ):
                raise ValueError("Source ATAC count layer contains invalid counts.")
            chunk = chunk.astype(np.int64)
            rows = np.arange(stop - start, dtype=np.int64)
            membership = sparse.csr_matrix(
                (
                    np.ones(stop - start, dtype=np.int8),
                    (rows, groups[start:stop]),
                ),
                shape=(stop - start, n_groups),
            )
            _update_dense_from_sparse(summed, membership.T @ chunk)
            binary = chunk.copy()
            binary.data = np.ones(binary.nnz, dtype=np.int8)
            _update_dense_from_sparse(coverage, membership.T @ binary)
            print(f"feature_aggregation cells={stop}/{source.n_obs}", flush=True)
    finally:
        source.file.close()
    return summed, coverage


def _rank_fold(
    donor: int,
    candidates: pd.DataFrame,
    summed: np.ndarray,
    coverage_by_group: np.ndarray,
) -> pd.DataFrame:
    donor_position = DONORS.index(donor)
    shape = (len(DONORS), len(CELL_TYPES), len(candidates))
    summed_3d = summed.reshape(shape)
    coverage_3d = coverage_by_group.reshape(shape)
    group_counts = summed_3d.sum(axis=0) - summed_3d[donor_position]
    coverage = coverage_3d.sum(axis=(0, 1)) - coverage_3d[donor_position].sum(axis=0)
    type_totals = group_counts.sum(axis=1, dtype=np.float64)
    if not np.isfinite(type_totals).all() or (type_totals <= 0).any():
        raise ValueError(f"Donor {donor} fold has an invalid cell-type count total.")
    log_normalized = np.log2(
        1.0 + 1.0e4 * group_counts.astype(np.float64) / type_totals[:, None]
    )
    score = np.var(log_normalized, axis=0, ddof=0)
    total_counts = group_counts.sum(axis=0, dtype=np.int64)
    eligible = np.flatnonzero(coverage >= MIN_REFERENCE_CELLS)
    if len(eligible) < N_TOP_PEAKS:
        raise ValueError(f"Donor {donor} has only {len(eligible)} eligible peaks.")
    peak_ids = candidates["feature_id"].astype(str).tolist()
    ranked = sorted(
        eligible.tolist(),
        key=lambda index: (
            -score[index],
            -int(coverage[index]),
            -int(total_counts[index]),
            peak_ids[index].encode("utf-8"),
        ),
    )[:N_TOP_PEAKS]
    result = candidates.iloc[ranked].copy()
    result.insert(0, "rank", np.arange(1, N_TOP_PEAKS + 1))
    result["score"] = score[ranked]
    result["nonzero_reference_cells"] = coverage[ranked]
    result["total_reference_counts"] = total_counts[ranked]
    return result


def build_feature_axes() -> None:
    if AXIS_ROOT.exists():
        raise FileExistsError(f"{AXIS_ROOT} already exists; feature axes are immutable.")
    labels = _load_labels()
    candidates = _candidate_features()
    summed, coverage = _aggregate_donor_type_counts(labels, candidates)
    temporary = Path(tempfile.mkdtemp(prefix=".broad7_lodo_v1.", dir=AXIS_ROOT.parent))
    try:
        selected_ids: set[str] = set()
        fold_records = []
        for donor in DONORS:
            fold = _rank_fold(donor, candidates, summed, coverage)
            donor_dir = temporary / f"donor_{donor}"
            donor_dir.mkdir(parents=True)
            table_path = donor_dir / "peak_selection.tsv.gz"
            list_path = donor_dir / "selected_peaks.txt"
            fold.to_csv(table_path, sep="\t", index=False, compression="gzip")
            list_path.write_text("\n".join(fold["feature_id"].astype(str)) + "\n")
            selected_ids.update(fold["feature_id"].astype(str))
            fold_records.append(
                {
                    "donor": donor,
                    "selected_peaks": N_TOP_PEAKS,
                    "selected_feature_sha256": ordered_feature_sha256(
                        fold["feature_id"].astype(str)
                    ),
                    "peak_selection_sha256": _sha256(table_path),
                }
            )
        union = candidates[candidates["feature_id"].astype(str).isin(selected_ids)].copy()
        union_dir = temporary / "union"
        union_dir.mkdir()
        union.to_csv(
            union_dir / "peaks.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        (union_dir / "selected_peaks.txt").write_text(
            "\n".join(union["feature_id"].astype(str)) + "\n"
        )
        _write_yaml(
            temporary / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "source": "GSE194122 author processed ATAC counts",
                "label_ontology": list(CELL_TYPES),
                "donors": list(DONORS),
                "candidate_peaks": len(candidates),
                "union_peaks": len(union),
                "selector": {
                    "n_top_peaks": N_TOP_PEAKS,
                    "min_reference_cells": MIN_REFERENCE_CELLS,
                    "scale": 1.0e4,
                    "score": "population_variance_log2_1_plus_scaled_type_rate",
                    "tie_breaks": [
                        "score_desc",
                        "coverage_desc",
                        "total_count_desc",
                        "peak_id_utf8_asc",
                    ],
                },
                "inputs": {
                    "source_h5ad": _repository_path(SOURCE_H5AD),
                    "labels": _repository_path(LABELS_PATH),
                    "features": _repository_path(FEATURES_PATH),
                },
                "folds": fold_records,
                "software_versions": _software_versions(),
            },
        )
        temporary.rename(AXIS_ROOT)
    except BaseException:
        import shutil

        shutil.rmtree(temporary, ignore_errors=True)
        raise


def _fragment_records() -> list[dict[str, Any]]:
    manifest = _read_yaml(FRAGMENT_MANIFEST)
    records = manifest.get("files")
    if not isinstance(records, list) or len(records) != 13:
        raise ValueError("Fragment manifest must contain all 13 samples.")
    lock = _read_yaml(SOURCE_LOCK)
    lock_files = lock.get("files")
    if not isinstance(lock_files, list):
        raise ValueError("Source lock has no file records.")
    hashes = {
        str(record.get("sample_key")): str(record["sha256"])
        for record in lock_files
        if record.get("role") == "author_cellranger_arc_fragments"
    }
    for record in records:
        sample_key = str(record["sample_key"])
        if sample_key not in hashes:
            raise ValueError(f"Missing frozen fragment hash for {sample_key}.")
        record["sha256"] = hashes[sample_key]
    return records


def build_fragment_caches() -> None:
    labels = _load_labels()
    union_path = AXIS_ROOT / "union" / "peaks.tsv.gz"
    if not union_path.is_file():
        raise FileNotFoundError("Run the feature-axes stage first.")
    union = pd.read_csv(union_path, sep="\t")
    peaks = [
        (
            str(row.chromosome),
            int(row.start),
            int(row.end),
            str(row.feature_id),
        )
        for row in union.itertuples(index=False)
    ]
    var = union.set_index("feature_id")[["chromosome", "start", "end"]].copy()
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    label_sha256 = _sha256(LABELS_PATH)
    for record in _fragment_records():
        sample_key = str(record["sample_key"])
        output_dir = CACHE_ROOT / sample_key
        output_path = output_dir / "cells.h5ad"
        if output_path.is_file():
            cached = ad.read_h5ad(output_path, backed="r")
            try:
                if cached.n_vars != len(union):
                    raise ValueError(f"Stale cache feature axis for {sample_key}.")
            finally:
                cached.file.close()
            print(f"fragment_cache sample={sample_key} status=reused", flush=True)
            continue
        if output_dir.exists():
            raise FileExistsError(f"Partial cache directory requires inspection: {output_dir}")
        sample_labels = labels[labels["sample_key"] == sample_key].copy()
        if sample_labels["fragment_barcode"].duplicated().any():
            raise ValueError(f"Fragment barcodes are not unique within {sample_key}.")
        result = count_fragment_shapes(
            ROOT / str(record["fragment"]),
            sample_labels["fragment_barcode"].astype(str).tolist(),
            peaks,
            right_cut_offset=0,
            chunk_size=COUNT_CHUNK_SIZE,
        )
        obs = sample_labels.set_index(sample_labels["fragment_barcode"].astype(str))
        obs.index.name = "fragment_barcode"
        shape = build_fragment_shape_anndata(
            result,
            obs=obs,
            var=var,
            provenance={
                "split_sha256": label_sha256,
                "source_sha256": {
                    f"{sample_key}.fragments.tsv.gz": str(record["sha256"]),
                    "labels": label_sha256,
                },
                "coordinate_validation": {
                    "selected_right_cut_offset": 0,
                    "matrix_match": "exact",
                    "mismatched_entries": 0,
                    "absolute_error": 0,
                    "comparison": "fragment_derived_collapsed_matrix",
                },
                "software_versions": _software_versions(),
            },
        )
        shape.obs_names = pd.Index(sample_labels["cell_id"].astype(str), name="cell_id")
        output_dir.mkdir(parents=True)
        temporary_path = output_dir / ".cells.h5ad.tmp"
        shape.write_h5ad(temporary_path, compression="gzip")
        temporary_path.replace(output_path)
        _write_yaml(
            output_dir / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "sample_key": sample_key,
                "cells": shape.n_obs,
                "peaks": shape.n_vars,
                "fragment_sha256": str(record["sha256"]),
                "output": _repository_path(output_path),
                "preprocessing_counters": result.qc.to_dict(),
            },
        )
        print(f"fragment_cache sample={sample_key} status=completed", flush=True)


def _combined_metadata(
    metadata_values: Iterable[Mapping[str, Any]],
    layers: Mapping[str, sparse.csr_matrix],
    feature_names: Iterable[str],
    split_sha256: str,
    source_hashes: Mapping[str, str],
) -> dict[str, Any]:
    values = [copy.deepcopy(dict(value)) for value in metadata_values]
    if not values:
        raise ValueError("Cannot combine zero fragment-shape metadata records.")
    metadata = values[0]
    counters: dict[str, int] = {}
    for value in values:
        for key, count in dict(value.get("preprocessing_counters", {})).items():
            if isinstance(count, (int, np.integer)):
                counters[key] = counters.get(key, 0) + int(count)
    layer_totals = {name: int(matrix.sum()) for name, matrix in layers.items()}
    metadata["preprocessing_counters"] = counters
    metadata["matrix_counters"] = {
        "assigned_cut_sites": int(sum(layer_totals.values())),
        **{f"cut_sites_per_bin.{name}": count for name, count in layer_totals.items()},
    }
    metadata["feature_sha256"] = ordered_feature_sha256(feature_names)
    metadata["split_sha256"] = split_sha256
    metadata["source_sha256"] = dict(source_hashes)
    return metadata


def _fold_objects(
    donor: int,
    selected_features: list[str],
    split_sha256: str,
) -> tuple[ad.AnnData, ad.AnnData]:
    records = _fragment_records()
    partitions: dict[str, dict[str, Any]] = {
        "training": {"obs": [], "layers": {}, "metadata": [], "sources": {}},
        "heldout": {"obs": [], "layers": {}, "metadata": [], "sources": {}},
    }
    var: pd.DataFrame | None = None
    for record in records:
        sample_key = str(record["sample_key"])
        cache = ad.read_h5ad(CACHE_ROOT / sample_key / "cells.h5ad")
        try:
            indices = cache.var_names.get_indexer(selected_features)
            if (indices < 0).any():
                raise ValueError(f"Cache {sample_key} is missing selected fold peaks.")
            if var is None:
                var = cache.var.iloc[indices].copy()
                var.index = pd.Index(selected_features, name=cache.var_names.name)
            for partition, mask in (
                ("training", cache.obs["donor"].astype(int).to_numpy() != donor),
                ("heldout", cache.obs["donor"].astype(int).to_numpy() == donor),
            ):
                if not mask.any():
                    continue
                selected = cache[mask, indices]
                partitions[partition]["obs"].append(selected.obs.copy())
                for layer_name in FragmentShapeSpec.from_mapping(
                    cache.uns["fragment_shape"]
                ).layer_names:
                    partitions[partition]["layers"].setdefault(layer_name, []).append(
                        sparse.csr_matrix(selected.layers[layer_name], dtype=np.int64)
                    )
                partitions[partition]["metadata"].append(cache.uns["fragment_shape"])
                partitions[partition]["sources"][
                    f"{sample_key}.fragments.tsv.gz"
                ] = str(record["sha256"])
        finally:
            del cache
    if var is None:
        raise RuntimeError("No cache features were assembled.")

    objects = []
    for partition in ("training", "heldout"):
        obs = pd.concat(partitions[partition]["obs"], axis=0)
        obs["cell_type"] = obs["author_cell_type"].map(BROAD_LABELS)
        if obs["cell_type"].isna().any() or set(obs["cell_type"]) != set(CELL_TYPES):
            raise ValueError(f"{partition} broad label universe is invalid.")
        layers = {
            name: sparse.vstack(matrices, format="csr", dtype=np.int64)
            for name, matrices in partitions[partition]["layers"].items()
        }
        x = sparse.csr_matrix((len(obs), len(selected_features)), dtype=np.int64)
        for matrix in layers.values():
            x = (x + matrix).tocsr()
        result = ad.AnnData(X=x, obs=obs, var=var.copy())
        for name, matrix in layers.items():
            result.layers[name] = matrix
        result.uns["fragment_shape"] = _combined_metadata(
            partitions[partition]["metadata"],
            layers,
            selected_features,
            split_sha256,
            partitions[partition]["sources"],
        )
        validate_fragment_shape_spec(
            FragmentShapeSpec.from_mapping(result.uns["fragment_shape"])
        )
        validate_fragment_shape_feature_axis(result, f"GSE194122 {partition}")
        objects.append(result)
    return objects[0], objects[1]


def _reference_id(donor: int) -> str:
    return f"gse194122_bmmc_broad7_lodo_donor_{donor}_v1"


def _dataset_id(donor: int, condition: str, mixture_seed: int) -> str:
    return (
        f"gse194122_shapemix_broad7_lodo_donor_{donor}_"
        f"{condition}_mix_{mixture_seed}"
    )


def build_fold_objects() -> None:
    labels = _load_labels()
    for donor in DONORS:
        selected_path = AXIS_ROOT / f"donor_{donor}" / "selected_peaks.txt"
        selected_features = selected_path.read_text().splitlines()
        fold_dir = FOLD_ROOT / f"donor_{donor}"
        reference_dir = REFERENCE_ROOT / _reference_id(donor)
        reference_path = reference_dir / "atac" / "reference.h5ad"
        heldout_path = fold_dir / "heldout_cells.h5ad"
        manifest_path = fold_dir / "manifest.yaml"
        if reference_path.is_file() and heldout_path.is_file() and manifest_path.is_file():
            print(f"fold donor={donor} status=reused", flush=True)
            continue
        if fold_dir.exists() or reference_dir.exists():
            raise FileExistsError(f"Partial donor-{donor} fold requires inspection.")
        fold_dir.mkdir(parents=True)
        split = labels[
            [
                "cell_id",
                "sample_key",
                "site",
                "donor",
                "author_cell_type",
                "cell_type",
            ]
        ].copy()
        split["partition"] = np.where(split["donor"] == donor, "heldout", "training")
        split_path = fold_dir / "cells.tsv.gz"
        split.to_csv(split_path, sep="\t", index=False, compression="gzip")
        split_sha256 = _sha256(split_path)
        reference, heldout = _fold_objects(donor, selected_features, split_sha256)
        temporary_reference = fold_dir / ".reference.h5ad.tmp"
        temporary_heldout = fold_dir / ".heldout.h5ad.tmp"
        reference.write_h5ad(temporary_reference, compression="gzip")
        heldout.write_h5ad(temporary_heldout, compression="gzip")
        reference_dir.joinpath("atac").mkdir(parents=True)
        temporary_reference.replace(reference_path)
        temporary_heldout.replace(heldout_path)
        _write_yaml(
            reference_dir / "reference.yaml",
            {
                "schema_version": 1,
                "reference_id": _reference_id(donor),
                "source": "GSE194122",
                "heldout_donor": donor,
                "labels_key": "cell_type",
                "cell_types": list(CELL_TYPES),
                "atac": {"path": _repository_path(reference_path)},
                "feature_sha256": ordered_feature_sha256(selected_features),
                "split_sha256": split_sha256,
            },
        )
        _write_yaml(
            manifest_path,
            {
                "schema_version": 1,
                "status": "complete",
                "heldout_donor": donor,
                "split": {
                    "path": _repository_path(split_path),
                    "sha256": split_sha256,
                    "training_cells": reference.n_obs,
                    "heldout_cells": heldout.n_obs,
                },
                "selected_peaks": {
                    "path": _repository_path(selected_path),
                    "count": len(selected_features),
                    "feature_sha256": ordered_feature_sha256(selected_features),
                },
                "reference": _repository_path(reference_path),
                "heldout": _repository_path(heldout_path),
            },
        )
        print(f"fold donor={donor} status=completed", flush=True)


def materialize_datasets() -> None:
    registry = _read_yaml(REGISTRY_PATH) if REGISTRY_PATH.exists() else {}
    for donor in DONORS:
        fold_dir = FOLD_ROOT / f"donor_{donor}"
        reference_path = REFERENCE_ROOT / _reference_id(donor) / "atac" / "reference.h5ad"
        heldout_path = fold_dir / "heldout_cells.h5ad"
        manifest_path = fold_dir / "manifest.yaml"
        reference = ad.read_h5ad(reference_path)
        heldout = ad.read_h5ad(heldout_path)
        reference_barcodes = reference.obs_names.astype(str).tolist()
        heldout_counts = heldout.obs["cell_type"].astype(str).value_counts().to_dict()
        for condition in CONDITIONS:
            probabilities = condition_probabilities(
                condition,
                CELL_TYPES,
                heldout_counts,
            )
            for mixture_seed in INNER_MIXTURE_SEEDS:
                dataset_id = _dataset_id(donor, condition, mixture_seed)
                dataset_path = DATASET_ROOT / dataset_id / "dataset.yaml"
                if not dataset_path.is_file():
                    simulation = simulate_shapemix_spots(
                        heldout,
                        cell_types=CELL_TYPES,
                        sampling_probabilities=probabilities,
                        condition=condition,
                        outer_split_seed=donor,
                        inner_mixture_seed=mixture_seed,
                        num_spots=1024,
                        mean_cells_per_spot=10.0,
                        labels_key="cell_type",
                        reference_barcodes=reference_barcodes,
                    )
                    dataset_path = write_simulation_dataset(
                        simulation,
                        output_root=DATASET_ROOT,
                        dataset_id=dataset_id,
                        reference_path=reference_path,
                        heldout_path=heldout_path,
                        split_manifest_path=manifest_path,
                        labels_key="cell_type",
                        benchmark_scope="external_donor_heldout_exact_truth",
                        source="gse194122_bmmc_multiome_donor_heldout_simulation",
                        description=(
                            "GSE194122 BMMC pseudo-spots aggregated only from the held-out "
                            f"donor {donor}; peaks and signatures use the other donors only."
                        ),
                        scientific_scope=(
                            "Leave-one-donor-out biological generalization across ten "
                            "GSE194122 BMMC donors."
                        ),
                    )
                    print(f"dataset {dataset_id} status=completed", flush=True)
                registry[dataset_id] = {"config": _repository_path(dataset_path)}
    temporary = REGISTRY_PATH.with_name(".datasets.yaml.gse194122.tmp")
    temporary.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open("w") as handle:
        yaml.safe_dump(registry, handle, sort_keys=False)
    temporary.replace(REGISTRY_PATH)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "stage",
        choices=("feature-axes", "fragment-caches", "fold-objects", "datasets", "all"),
    )
    args = parser.parse_args()
    stages = (
        ("feature-axes", build_feature_axes),
        ("fragment-caches", build_fragment_caches),
        ("fold-objects", build_fold_objects),
        ("datasets", materialize_datasets),
    )
    for name, function in stages:
        if args.stage in {name, "all"}:
            function()


if __name__ == "__main__":
    main()
