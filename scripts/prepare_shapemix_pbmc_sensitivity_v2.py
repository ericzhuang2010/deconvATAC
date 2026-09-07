#!/usr/bin/env python
"""Build the support-safe, versioned PBMC ShapeMix sensitivity campaign.

Version 1 selected 5,000 peaks before applying the subtype and reference-cell
subsets. Some of those peaks consequently had zero pooled reference counts in
the restricted references. Version 2 freezes one common 5,000-peak axis that
has positive support in every predeclared reference variant, using training
reference counts only, and recounts fragment shapes once on that axis.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Mapping

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import yaml

from deconvatac.data import FragmentShapeSpec, ordered_feature_sha256
from deconvatac.data.validators import validate_fragment_shape_feature_axis
from deconvatac.pp import select_reference_peaks
from deconvatac.pp.fragment_shapes import (
    FragmentLengthBin,
    FragmentShapeResult,
    build_fragment_shape_anndata,
    count_fragment_shapes,
)
from scripts import prepare_shapemix_pbmc_sensitivity as v1


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_ID = "shapemix_pbmc_stress_v2"
TEMPLATE_PATH = ROOT / "configs/datasets/shapemix_pbmc_stress_v2.yaml"
EXPERIMENT_PATH = ROOT / "configs/experiments/shapemix_pbmc_stress_v2.yaml"
REUSABLE_ROOT = (
    ROOT
    / "data/processed/shapemix/pbmc_granulocyte_sorted_10k"
    / "sensitivity/stress_v2"
)
REFERENCE_ROOT = ROOT / "data/processed/references/shapemix_pbmc_stress_v2"
WORK_ROOT = ROOT / "data/work/preprocessing/shapemix_pbmc_stress_v2"
RANKING_REFERENCE = (
    ROOT
    / "data/processed/references/pbmc_granulocyte_sorted_10k_multiome"
    / "atac/reference.h5ad"
)
SOURCE_MANIFEST = REUSABLE_ROOT / "manifests/source_axis.yaml"
FEATURE_AXIS_ROOT = REUSABLE_ROOT / "feature_axis"
THREE_BINS = (
    FragmentLengthBin("short", 0, 100, "fragment_length_lt_100"),
    FragmentLengthBin("mono", 100, 250, "fragment_length_100_249"),
    FragmentLengthBin("long", 250, None, "fragment_length_ge_250"),
)
_V1_DESIGN_ROWS = v1._design_rows


def _ranked_indices(selection: Any) -> np.ndarray:
    eligible = np.flatnonzero(selection.eligible_mask).tolist()
    ids = selection.candidate_peak_ids
    return np.asarray(
        sorted(
            eligible,
            key=lambda index: (
                -selection.candidate_scores[index],
                -selection.candidate_nonzero_reference_cells[index],
                -selection.candidate_total_reference_counts[index],
                str(ids[index]).encode("utf-8"),
            ),
        ),
        dtype=np.int64,
    )


def _reference_rows(
    reference: ad.AnnData,
    cell_types: tuple[str, ...],
    per_type: int | None = None,
) -> np.ndarray:
    labels = reference.obs["cell_type"].astype(str)
    names: list[str] = []
    for cell_type in cell_types:
        candidates = reference.obs_names[labels == cell_type].astype(str).tolist()
        if per_type is None:
            names.extend(candidates)
            continue
        if len(candidates) < per_type:
            raise ValueError(f"{cell_type} has fewer than {per_type} reference cells.")
        names.extend(
            sorted(
                candidates,
                key=lambda name: (hashlib.sha256(name.encode()).hexdigest(), name),
            )[:per_type]
        )
    rows = reference.obs_names.get_indexer(names)
    if np.any(rows < 0) or len(set(rows.tolist())) != len(rows):
        raise ValueError("Reference-variant barcode selection is invalid.")
    return rows


def _derive_three_bin_result(five: FragmentShapeResult) -> FragmentShapeResult:
    if tuple(five.bins) != v1.FIVE_BINS:
        raise ValueError("Three-bin derivation requires the frozen five-bin axis.")
    layers = {
        THREE_BINS[0].layer: sparse.csr_matrix(
            five.layers[v1.FIVE_BINS[0].layer]
            + five.layers[v1.FIVE_BINS[1].layer]
        ),
        THREE_BINS[1].layer: sparse.csr_matrix(
            five.layers[v1.FIVE_BINS[2].layer]
            + five.layers[v1.FIVE_BINS[3].layer]
        ),
        THREE_BINS[2].layer: sparse.csr_matrix(
            five.layers[v1.FIVE_BINS[4].layer]
        ),
    }
    for matrix in layers.values():
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        matrix.sort_indices()
    qc = copy.deepcopy(five.qc)
    qc.cut_sites_per_bin = {
        name: int(matrix.sum()) for name, matrix in layers.items()
    }
    return FragmentShapeResult(
        barcodes=five.barcodes,
        peaks=five.peaks,
        bins=THREE_BINS,
        layers=layers,
        qc=qc,
        right_cut_offset=five.right_cut_offset,
    )


def _feature_axis(
    ranking_reference: ad.AnnData,
    membership: pd.DataFrame,
) -> tuple[list[str], pd.DataFrame, dict[str, Any]]:
    reference_barcodes = membership.loc[
        membership["pool"] == "reference", "barcode"
    ].astype(str).tolist()
    training = ranking_reference[reference_barcodes, :].copy()
    observed = training.obs["cell_type"].astype(str)
    expected = membership.set_index("barcode").loc[reference_barcodes, "cell_type"].astype(str)
    if not np.array_equal(observed.to_numpy(), expected.to_numpy()):
        raise ValueError("Ranking-reference labels disagree with split membership.")
    selection = select_reference_peaks(
        training,
        "cell_type",
        cell_types=v1.PBMC_CELL_TYPES,
        n_top_peaks=5000,
        min_reference_cells=10,
        scale=1.0e4,
    )
    frozen_v1 = (v1.SPLIT_DIR / "selected_peaks.txt").read_text().splitlines()
    if list(selection.peak_ids) != frozen_v1:
        raise ValueError("Recomputed protocol-v1 ranking differs from the frozen split.")

    variants = {
        "all16_all": np.arange(training.n_obs),
        "broad3_all": _reference_rows(training, v1.BROAD3),
        "broad3_per_type_250": _reference_rows(training, v1.BROAD3, 250),
        "broad3_per_type_100": _reference_rows(training, v1.BROAD3, 100),
        "broad3_per_type_50": _reference_rows(training, v1.BROAD3, 50),
        "cd4_related3_all": _reference_rows(training, v1.CD4_RELATED3),
    }
    support: dict[str, np.ndarray] = {}
    for name, rows in variants.items():
        support[name] = np.asarray(training.X[rows, :].sum(axis=0)).ravel() > 0
    universal = np.logical_and.reduce(list(support.values())) & selection.eligible_mask
    ranked = _ranked_indices(selection)
    selected_indices = ranked[universal[ranked]][:5000]
    if len(selected_indices) != 5000:
        raise ValueError(
            f"Only {len(selected_indices)} eligible peaks have support in every reference variant."
        )
    selected = [str(selection.candidate_peak_ids[index]) for index in selected_indices]
    if len(set(selected)) != 5000:
        raise ValueError("Corrected feature axis is not 5,000 unique peaks.")

    rank = np.zeros(len(selection.candidate_peak_ids), dtype=np.int64)
    rank[ranked] = np.arange(1, len(ranked) + 1)
    audit = pd.DataFrame(
        {
            "peak_id": selection.candidate_peak_ids,
            "original_index": np.arange(len(selection.candidate_peak_ids)),
            "protocol_v1_rank": rank,
            "score": selection.candidate_scores,
            "nonzero_reference_cells": selection.candidate_nonzero_reference_cells,
            "total_reference_count": selection.candidate_total_reference_counts,
            "eligible_all16": selection.eligible_mask,
            **{f"supported_{name}": mask for name, mask in support.items()},
            "universal_reference_support": universal,
        }
    )
    selected_set = set(selected)
    frozen_set = set(frozen_v1)
    audit["selected_v1"] = audit["peak_id"].isin(frozen_set)
    audit["selected_v2"] = audit["peak_id"].isin(selected_set)
    unsupported_v1 = {
        name: int(sum(not mask[index] for index in selection.indices))
        for name, mask in support.items()
    }
    metadata = {
        "candidate_peaks": len(selection.candidate_peak_ids),
        "eligible_all16": int(selection.eligible_mask.sum()),
        "universal_supported_eligible": int(universal.sum()),
        "selected_peaks": len(selected),
        "v1_peaks_replaced": len(selected_set.difference(frozen_set)),
        "v1_unsupported_peaks_by_variant": unsupported_v1,
        "candidate_feature_sha256": selection.candidate_feature_sha256,
        "v1_selected_feature_sha256": selection.selected_feature_sha256,
        "v2_selected_feature_sha256": ordered_feature_sha256(selected),
    }
    return selected, audit, metadata


def _validate_cached_source() -> None:
    manifest = v1._read_yaml(SOURCE_MANIFEST)
    if manifest.get("status") != "complete" or manifest.get("campaign") != CAMPAIGN_ID:
        raise ValueError(f"Invalid PBMC v2 source manifest: {SOURCE_MANIFEST}")
    for relative in (
        "feature_axis/reference_cells.h5ad",
        "feature_axis/heldout_cells.h5ad",
        "fragment_shape_cache/bins_three_full_split_1103.h5ad",
        "fragment_shape_cache/bins_two_full_split_1103.h5ad",
        "fragment_shape_cache/bins_five_full_split_1103.h5ad",
    ):
        if not (REUSABLE_ROOT / relative).is_file():
            raise FileNotFoundError(REUSABLE_ROOT / relative)


def materialize_source_axis() -> None:
    if SOURCE_MANIFEST.is_file():
        _validate_cached_source()
        return
    if REUSABLE_ROOT.exists():
        raise FileExistsError(f"Partial PBMC v2 source tree requires inspection: {REUSABLE_ROOT}")
    for required in (
        TEMPLATE_PATH,
        EXPERIMENT_PATH,
        RANKING_REFERENCE,
        v1.RAW_FRAGMENTS,
        v1.SPLIT_DIR / "split.csv",
        v1.SPLIT_DIR / "reference_cells.h5ad",
    ):
        if not required.is_file():
            raise FileNotFoundError(required)

    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    temporary_parent = Path(tempfile.mkdtemp(prefix="source_axis.", dir=WORK_ROOT))
    staging = temporary_parent / "stress_v2"
    staging.mkdir()
    try:
        membership = pd.read_csv(v1.SPLIT_DIR / "split.csv")
        ranking_reference = ad.read_h5ad(RANKING_REFERENCE)
        missing = pd.Index(membership["barcode"].astype(str)).difference(
            ranking_reference.obs_names
        )
        if len(missing):
            raise ValueError(f"Ranking reference lacks split barcode {missing[0]!r}.")
        selected, audit, axis_metadata = _feature_axis(ranking_reference, membership)
        combined_names = membership["barcode"].astype(str).tolist()
        combined_obs = ranking_reference.obs.loc[combined_names].copy()
        expected_labels = membership.set_index("barcode").loc[combined_names, "cell_type"].astype(str)
        if not np.array_equal(
            combined_obs["cell_type"].astype(str).to_numpy(), expected_labels.to_numpy()
        ):
            raise ValueError("Combined source labels disagree with split membership.")
        combined_obs["split_pool"] = (
            membership.set_index("barcode").loc[combined_names, "pool"].astype(str).to_numpy()
        )
        selected_var = ranking_reference.var.loc[selected].copy()
        peaks = [
            (str(row.chrom), int(row.start), int(row.end), str(name))
            for name, row in selected_var.iterrows()
        ]
        five = count_fragment_shapes(
            v1.RAW_FRAGMENTS,
            combined_names,
            peaks,
            right_cut_offset=0,
            bins=v1.FIVE_BINS,
            chunk_size=1_000_000,
        )
        source_reference = ad.read_h5ad(v1.SPLIT_DIR / "reference_cells.h5ad")
        source_spec = FragmentShapeSpec.from_mapping(
            source_reference.uns["fragment_shape"]
        )
        provenance = {
            "split_sha256": source_spec.split_sha256,
            "source_sha256": copy.deepcopy(source_spec.source_sha256),
            "coordinate_validation": copy.deepcopy(source_spec.coordinate_validation),
            "software_versions": {
                "python": v1.platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": __import__("scipy").__version__,
                "anndata": ad.__version__,
            },
        }
        five_full = build_fragment_shape_anndata(
            five, obs=combined_obs, var=selected_var, provenance=provenance
        )
        official = ranking_reference[combined_names, selected].copy().X
        difference = five_full.X != official
        if difference.nnz:
            raise ValueError(
                f"Corrected fragment recount differs from official counts at {difference.nnz} entries."
            )
        three = _derive_three_bin_result(five)
        three_full = build_fragment_shape_anndata(
            three, obs=combined_obs, var=selected_var, provenance=provenance
        )
        two = v1.derive_two_bin_result(five)
        two_full = build_fragment_shape_anndata(
            two, obs=combined_obs, var=selected_var, provenance=provenance
        )

        cache = staging / "fragment_shape_cache"
        cache.mkdir(parents=True)
        cache_paths = {
            "three": cache / "bins_three_full_split_1103.h5ad",
            "two": cache / "bins_two_full_split_1103.h5ad",
            "five": cache / "bins_five_full_split_1103.h5ad",
        }
        for level, value in (
            ("three", three_full),
            ("two", two_full),
            ("five", five_full),
        ):
            v1._atomic_h5ad(value, cache_paths[level])

        reference_mask = three_full.obs["split_pool"].astype(str).eq("reference")
        heldout_mask = three_full.obs["split_pool"].astype(str).eq("heldout")
        reference = v1.subset_shape_cells(
            three_full, observation_mask=reference_mask
        )
        heldout = v1.subset_shape_cells(three_full, observation_mask=heldout_mask)
        axis_root = staging / "feature_axis"
        axis_root.mkdir(parents=True)
        reference_path = axis_root / "reference_cells.h5ad"
        heldout_path = axis_root / "heldout_cells.h5ad"
        selected_path = axis_root / "selected_peaks.txt"
        audit_path = axis_root / "feature_support_audit.csv"
        v1._atomic_h5ad(reference, reference_path)
        v1._atomic_h5ad(heldout, heldout_path)
        selected_path.write_text("\n".join(selected) + "\n")
        audit.to_csv(audit_path, index=False)

        v1._write_yaml(
            cache / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "campaign": CAMPAIGN_ID,
                "raw_fragments": v1._repository_path(v1.RAW_FRAGMENTS),
                "raw_fragment_sha256": v1._sha256(v1.RAW_FRAGMENTS),
                "stream_count": 1,
                "three_bin_cache": v1._repository_path(
                    REUSABLE_ROOT / "fragment_shape_cache" / cache_paths["three"].name
                ),
                "five_bin_cache": v1._repository_path(
                    REUSABLE_ROOT / "fragment_shape_cache" / cache_paths["five"].name
                ),
                "two_bin_cache": v1._repository_path(
                    REUSABLE_ROOT / "fragment_shape_cache" / cache_paths["two"].name
                ),
                "three_bin_derivation": "exact sums of five-bin layers",
                "two_bin_derivation": "exact sums of five-bin layers",
                "official_matrix_mismatch_entries": 0,
            },
        )
        manifest_path = staging / "manifests/source_axis.yaml"
        v1._write_yaml(
            manifest_path,
            {
                "schema_version": 1,
                "status": "complete",
                "campaign": CAMPAIGN_ID,
                "amends_campaign": "shapemix_pbmc_stress_v1",
                "selection_scope": "training_reference_only",
                "selection_rule": (
                    "protocol-v1 ranking filtered to peaks with positive pooled counts "
                    "in every predeclared reference variant"
                ),
                "one_common_axis_across_all_factors": True,
                "source_split": v1._repository_path(v1.SPLIT_DIR),
                "source_split_manifest_sha256": v1._sha256(v1.SPLIT_DIR / "manifest.yaml"),
                "ranking_reference": v1._repository_path(RANKING_REFERENCE),
                "ranking_reference_sha256": v1._sha256(RANKING_REFERENCE),
                "template": v1._repository_path(TEMPLATE_PATH),
                "template_sha256": v1._sha256(TEMPLATE_PATH),
                "experiment": v1._repository_path(EXPERIMENT_PATH),
                "experiment_sha256": v1._sha256(EXPERIMENT_PATH),
                "feature_axis": axis_metadata,
                "outputs": {
                    "selected_peaks": {
                        "path": v1._repository_path(FEATURE_AXIS_ROOT / selected_path.name),
                        "sha256": v1._sha256(selected_path),
                    },
                    "feature_support_audit": {
                        "path": v1._repository_path(FEATURE_AXIS_ROOT / audit_path.name),
                        "sha256": v1._sha256(audit_path),
                    },
                    "reference_cells": {
                        "path": v1._repository_path(FEATURE_AXIS_ROOT / reference_path.name),
                        "sha256": v1._sha256(reference_path),
                    },
                    "heldout_cells": {
                        "path": v1._repository_path(FEATURE_AXIS_ROOT / heldout_path.name),
                        "sha256": v1._sha256(heldout_path),
                    },
                },
            },
        )
        REUSABLE_ROOT.parent.mkdir(parents=True, exist_ok=True)
        os.replace(staging, REUSABLE_ROOT)
    except BaseException:
        shutil.rmtree(temporary_parent, ignore_errors=True)
        raise
    shutil.rmtree(temporary_parent, ignore_errors=True)
    _validate_cached_source()


def _load_base_v2() -> tuple[ad.AnnData, ad.AnnData]:
    materialize_source_axis()
    reference = ad.read_h5ad(FEATURE_AXIS_ROOT / "reference_cells.h5ad")
    heldout = ad.read_h5ad(FEATURE_AXIS_ROOT / "heldout_cells.h5ad")
    if not reference.var_names.equals(heldout.var_names) or reference.n_vars != 5000:
        raise ValueError("PBMC v2 source objects must share one 5,000-peak axis.")
    for role, value, pool in (
        ("reference", reference, "reference"),
        ("heldout", heldout, "heldout"),
    ):
        validate_fragment_shape_feature_axis(value, f"PBMC v2 {role}")
        if set(value.obs["split_pool"].astype(str)) != {pool}:
            raise ValueError(f"Unexpected split pool in PBMC v2 {role}.")
    return reference, heldout


def _design_rows_v2() -> list[v1.DesignRow]:
    rows = []
    for row in _V1_DESIGN_ROWS():
        rows.append(
            v1.DesignRow(
                dataset_id=row.dataset_id.replace(
                    "pbmc_shapemix_stress_v1", "pbmc_shapemix_stress_v2", 1
                ),
                factor=row.factor,
                level=row.level,
                mixture_seed=row.mixture_seed,
                control_level=row.control_level,
                reference_variant=row.reference_variant,
                cell_types=row.cell_types,
                condition=row.condition,
                mean_cells_per_spot=row.mean_cells_per_spot,
                depth_retain_probability=row.depth_retain_probability,
                rare_nk_fraction=row.rare_nk_fraction,
            )
        )
    return rows


def _write_reference_v2(
    variant: str,
    reference: ad.AnnData,
    metadata: Mapping[str, Any],
) -> Path:
    path = REFERENCE_ROOT / variant / "atac/reference.h5ad"
    descriptor = path.parents[1] / "reference.yaml"
    if path.is_file() and descriptor.is_file():
        existing = v1._read_yaml(descriptor)
        if existing.get("reference_id") != f"{CAMPAIGN_ID}/{variant}":
            raise ValueError(f"Stale reference descriptor: {descriptor}")
        return path
    if path.exists() or descriptor.exists() or path.parents[1].exists():
        raise FileExistsError(f"Partial v2 reference requires inspection: {path.parents[1]}")
    validate_fragment_shape_feature_axis(reference, f"PBMC v2 reference {variant}")
    pooled = np.asarray(reference.X.sum(axis=0)).ravel()
    if np.any(pooled <= 0):
        peak = reference.var_names[int(np.flatnonzero(pooled <= 0)[0])]
        raise ValueError(f"PBMC v2 reference {variant} has unsupported peak {peak!r}.")
    v1._atomic_h5ad(reference, path)
    spec = FragmentShapeSpec.from_mapping(reference.uns["fragment_shape"])
    v1._write_yaml(
        descriptor,
        {
            "schema_version": 1,
            "reference_id": f"{CAMPAIGN_ID}/{variant}",
            "source_dataset_id": "pbmc_granulocyte_sorted_10k",
            "source_split": v1._repository_path(v1.SPLIT_DIR),
            "source_feature_axis": v1._repository_path(SOURCE_MANIFEST),
            "labels_key": "cell_type",
            "cell_types": list(dict.fromkeys(reference.obs["cell_type"].astype(str))),
            "cells": reference.n_obs,
            "peaks": reference.n_vars,
            "bins": len(spec.bins),
            "feature_sha256": ordered_feature_sha256(reference.var_names),
            "all_peaks_have_positive_pooled_reference_counts": True,
            "atac": {"path": v1._repository_path(path)},
            "sensitivity": dict(metadata),
        },
    )
    return path


def _heldout_for_row_v2(
    row: v1.DesignRow,
    heldout: ad.AnnData,
) -> tuple[ad.AnnData, Path]:
    if row.factor == "bins":
        path = REUSABLE_ROOT / "splits" / f"bins_{row.level}" / "heldout_cells.h5ad"
        return ad.read_h5ad(path), path
    peaks = 5000
    if row.factor == "features":
        peaks = int(row.level.removeprefix("peaks_"))
    return (
        v1._subset(heldout, cell_types=row.cell_types, num_features=peaks),
        FEATURE_AXIS_ROOT / "heldout_cells.h5ad",
    )


@contextmanager
def _patched_v1_module() -> Iterator[None]:
    replacements = {
        "TEMPLATE_PATH": TEMPLATE_PATH,
        "EXPERIMENT_PATH": EXPERIMENT_PATH,
        "REUSABLE_ROOT": REUSABLE_ROOT,
        "REFERENCE_ROOT": REFERENCE_ROOT,
        "_load_base": _load_base_v2,
        "_design_rows": _design_rows_v2,
        "_write_reference": _write_reference_v2,
        "_heldout_for_row": _heldout_for_row_v2,
    }
    original = {name: getattr(v1, name) for name in replacements}
    try:
        for name, value in replacements.items():
            setattr(v1, name, value)
        yield
    finally:
        for name, value in original.items():
            setattr(v1, name, value)


def prepare_references() -> None:
    with _patched_v1_module():
        v1.prepare_references()


def _enrich_dataset(dataset_path: Path, row: v1.DesignRow) -> None:
    sensitivity = {
        "campaign": CAMPAIGN_ID,
        "factor": row.factor,
        "level": row.level,
        "control_level": row.control_level,
        "outer_split_seed": v1.OUTER_SPLIT_SEED,
        "mixture_seed": row.mixture_seed,
        "reference_variant": row.reference_variant,
        "mean_cells_per_spot": row.mean_cells_per_spot,
        "depth_retain_probability": row.depth_retain_probability,
        "rare_nk_fraction": row.rare_nk_fraction,
        "observed_nk_fraction": v1.OBSERVED_NK_FRACTION,
        "feature_axis_manifest": v1._repository_path(SOURCE_MANIFEST),
    }
    descriptor = v1._read_yaml(dataset_path)
    if descriptor.get("sensitivity") not in (None, sensitivity):
        raise ValueError(f"Stale v2 sensitivity metadata: {dataset_path}")
    descriptor["sensitivity"] = sensitivity
    v1._write_yaml(dataset_path, descriptor)
    manifest_path = dataset_path.parent / "simulation/manifest.yaml"
    manifest = v1._read_yaml(manifest_path)
    manifest["sensitivity"] = sensitivity
    record = manifest["outputs"]["dataset_yaml"]
    record["bytes"] = dataset_path.stat().st_size
    record["sha256"] = v1._sha256(dataset_path)
    v1._write_yaml(manifest_path, manifest)


def materialize_datasets() -> None:
    with _patched_v1_module():
        v1._validate_template()
        _, base_heldout = _load_base_v2()
        rows = _design_rows_v2()
        expected = v1._read_yaml(EXPERIMENT_PATH)["datasets"]
        if [row.dataset_id for row in rows] != expected:
            raise ValueError("PBMC v2 builder order differs from the frozen experiment.")
        registry = v1._read_yaml(v1.REGISTRY_PATH) if v1.REGISTRY_PATH.exists() else {}
        design_records: list[dict[str, Any]] = []
        for row in rows:
            reference_path = REFERENCE_ROOT / row.reference_variant / "atac/reference.h5ad"
            if not reference_path.is_file():
                raise FileNotFoundError(f"Run the v2 references stage first: {reference_path}")
            reference = ad.read_h5ad(reference_path)
            heldout, heldout_path = _heldout_for_row_v2(row, base_heldout)
            if not reference.var_names.equals(heldout.var_names):
                raise ValueError(f"Reference/held-out axes differ for {row.dataset_id}.")
            dataset_path = v1.DATASET_ROOT / row.dataset_id / "dataset.yaml"
            if not dataset_path.is_file():
                simulation = v1.simulate_shapemix_spots(
                    heldout,
                    cell_types=row.cell_types,
                    sampling_probabilities=v1._probabilities(row),
                    condition=row.condition,
                    outer_split_seed=v1.OUTER_SPLIT_SEED,
                    inner_mixture_seed=row.mixture_seed,
                    num_spots=v1.NUM_SPOTS,
                    mean_cells_per_spot=row.mean_cells_per_spot,
                    labels_key="cell_type",
                    grid_shape=(32, 32),
                    reference_barcodes=reference.obs_names.astype(str).tolist(),
                    depth_retain_probability=row.depth_retain_probability,
                )
                dataset_path = v1.write_simulation_dataset(
                    simulation,
                    output_root=v1.DATASET_ROOT,
                    dataset_id=row.dataset_id,
                    reference_path=reference_path,
                    heldout_path=heldout_path,
                    split_manifest_path=SOURCE_MANIFEST,
                    labels_key="cell_type",
                    benchmark_scope="secondary_one_donor_diagnostic_sensitivity",
                    source="pbmc_granulocyte_sorted_10k_shapemix_stress_v2",
                    description=(
                        "Support-safe one-factor-at-a-time PBMC ShapeMix diagnostic "
                        "sensitivity on split-1103 held-out source cells."
                    ),
                    scientific_scope=(
                        "Conditional resampling within one PBMC Multiome donor; explanatory, "
                        "not donor-level generalization."
                    ),
                )
                _enrich_dataset(dataset_path, row)
                print(f"dataset {row.dataset_id} status=completed", flush=True)
            else:
                descriptor = v1._read_yaml(dataset_path)
                if (descriptor.get("sensitivity") or {}).get("campaign") != CAMPAIGN_ID:
                    raise ValueError(f"Existing dataset is not PBMC v2: {dataset_path}")
                print(f"dataset {row.dataset_id} status=reused", flush=True)
            registry[row.dataset_id] = {"config": v1._repository_path(dataset_path)}
            design_records.append(
                {
                    **row.__dict__,
                    "cell_types": "|".join(row.cell_types),
                    "reference_path": v1._repository_path(reference_path),
                }
            )

        design_path = REUSABLE_ROOT / "manifests/design.csv"
        design_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_design = design_path.with_name(".design.csv.tmp")
        pd.DataFrame.from_records(design_records).to_csv(temporary_design, index=False)
        temporary_design.replace(design_path)
        v1._write_yaml(
            REUSABLE_ROOT / "manifests/materialization.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "campaign": CAMPAIGN_ID,
                "template": v1._repository_path(TEMPLATE_PATH),
                "experiment": v1._repository_path(EXPERIMENT_PATH),
                "source_axis": v1._repository_path(SOURCE_MANIFEST),
                "datasets": len(rows),
                "design": {
                    "path": v1._repository_path(design_path),
                    "sha256": v1._sha256(design_path),
                },
            },
        )
        temporary_registry = v1.REGISTRY_PATH.with_name(".datasets.yaml.pbmc_stress_v2.tmp")
        with temporary_registry.open("w") as handle:
            yaml.safe_dump(registry, handle, sort_keys=False)
        temporary_registry.replace(v1.REGISTRY_PATH)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("source", "references", "datasets", "all"))
    args = parser.parse_args()
    if args.stage in {"source", "all"}:
        materialize_source_axis()
    if args.stage in {"references", "all"}:
        prepare_references()
    if args.stage in {"datasets", "all"}:
        materialize_datasets()


if __name__ == "__main__":
    main()
