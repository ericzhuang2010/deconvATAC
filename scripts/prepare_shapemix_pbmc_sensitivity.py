#!/usr/bin/env python
"""Build the frozen one-factor-at-a-time ShapeMix PBMC sensitivities."""

from __future__ import annotations

import argparse
import copy
import hashlib
import os
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

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
    FragmentLengthBin,
    FragmentShapeResult,
    build_fragment_shape_anndata,
    count_fragment_shapes,
)
from scripts.regenerate_shapemix_pbmc_simulations import (
    PBMC_CELL_TYPE_COUNTS,
    PBMC_CELL_TYPES,
    condition_probabilities,
    simulate_shapemix_spots,
    subset_shape_cells,
    write_simulation_dataset,
)


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE_PATH = ROOT / "configs/datasets/shapemix_pbmc_stress_v1.yaml"
EXPERIMENT_PATH = ROOT / "configs/experiments/shapemix_pbmc_stress_v1.yaml"
SPLIT_DIR = (
    ROOT
    / "data/processed/shapemix/pbmc_granulocyte_sorted_10k/split_1103"
)
RAW_FRAGMENTS = (
    ROOT
    / "data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k"
    / "cellranger_arc_2.0.0"
    / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
)
REUSABLE_ROOT = (
    ROOT
    / "data/processed/shapemix/pbmc_granulocyte_sorted_10k"
    / "sensitivity/stress_v1"
)
REFERENCE_ROOT = (
    ROOT / "data/processed/references/shapemix_pbmc_stress_v1"
)
DATASET_ROOT = ROOT / "data/processed/datasets"
REGISTRY_PATH = ROOT / "data/registry/datasets.yaml"
OUTER_SPLIT_SEED = 1103
MIXTURE_SEEDS = (307, 401)
NUM_SPOTS = 1024
BROAD3 = ("CD14 Mono", "CD4 Naive", "CD8 Naive")
CD4_RELATED3 = ("CD4 Naive", "CD4 TCM", "CD4 TEM")
OBSERVED_NK_FRACTION = dict(PBMC_CELL_TYPE_COUNTS)["NK"] / sum(
    count for _, count in PBMC_CELL_TYPE_COUNTS
)
FIVE_BINS = (
    FragmentLengthBin("very_short", 0, 80, "fragment_length_lt_80"),
    FragmentLengthBin("short", 80, 100, "fragment_length_80_99"),
    FragmentLengthBin("mono_short", 100, 180, "fragment_length_100_179"),
    FragmentLengthBin("mono_long", 180, 250, "fragment_length_180_249"),
    FragmentLengthBin("long", 250, None, "fragment_length_ge_250"),
)
TWO_BINS = (
    FragmentLengthBin("short", 0, 100, "fragment_length_lt_100"),
    FragmentLengthBin("long", 100, None, "fragment_length_ge_100"),
)


@dataclass(frozen=True)
class DesignRow:
    dataset_id: str
    factor: str
    level: str
    mixture_seed: int
    control_level: str
    reference_variant: str
    cell_types: tuple[str, ...]
    condition: str
    mean_cells_per_spot: float = 10.0
    depth_retain_probability: float | None = None
    rare_nk_fraction: float | None = None


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


def _repository_path(path: Path) -> str:
    return path.absolute().relative_to(ROOT.absolute()).as_posix()


def _atomic_h5ad(adata: ad.AnnData, path: Path) -> None:
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    try:
        adata.write_h5ad(temporary, compression="gzip")
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _validate_template() -> None:
    template = _read_yaml(TEMPLATE_PATH)
    experiment = _read_yaml(EXPERIMENT_PATH)
    if template.get("status") != "frozen_before_predictions":
        raise ValueError("PBMC sensitivity template is not frozen.")
    if template.get("materialization", {}).get("unique_datasets") != 40:
        raise ValueError("PBMC sensitivity template must declare 40 datasets.")
    if len(experiment.get("datasets", [])) != 40:
        raise ValueError("PBMC sensitivity experiment must declare 40 datasets.")


def _load_base() -> tuple[ad.AnnData, ad.AnnData]:
    reference = ad.read_h5ad(SPLIT_DIR / "reference_cells.h5ad")
    heldout = ad.read_h5ad(SPLIT_DIR / "heldout_test_cells.h5ad")
    if not reference.var_names.equals(heldout.var_names) or reference.n_vars != 5000:
        raise ValueError("Split 1103 must contain one aligned 5,000-peak axis.")
    for role, value, pool in (
        ("reference", reference, "reference"),
        ("heldout", heldout, "heldout"),
    ):
        validate_fragment_shape_feature_axis(value, f"PBMC sensitivity {role}")
        spec = FragmentShapeSpec.from_mapping(value.uns["fragment_shape"])
        validate_fragment_shape_spec(spec)
        if len(spec.bins) != 3 or set(value.obs["split_pool"].astype(str)) != {pool}:
            raise ValueError(f"Unexpected source contract for {role}.")
        if set(value.obs["cell_type"].astype(str)) != set(PBMC_CELL_TYPES):
            raise ValueError(f"Unexpected source cell-type universe for {role}.")
    return reference, heldout


def _subset(
    adata: ad.AnnData,
    *,
    cell_types: Sequence[str] = PBMC_CELL_TYPES,
    num_features: int = 5000,
) -> ad.AnnData:
    mask = adata.obs["cell_type"].astype(str).isin(tuple(cell_types))
    features = adata.var_names[:num_features].astype(str).tolist()
    return subset_shape_cells(
        adata,
        observation_mask=mask,
        feature_names=features,
    )


def _cap_reference(reference: ad.AnnData, per_type: int | None) -> ad.AnnData:
    if per_type is None:
        return reference.copy()
    keep: set[str] = set()
    labels = reference.obs["cell_type"].astype(str)
    for cell_type in BROAD3:
        names = reference.obs_names[labels == cell_type].astype(str).tolist()
        if len(names) < per_type:
            raise ValueError(f"{cell_type} has fewer than {per_type} reference cells.")
        ranked = sorted(
            names,
            key=lambda name: (hashlib.sha256(name.encode()).hexdigest(), name),
        )
        keep.update(ranked[:per_type])
    return subset_shape_cells(
        reference,
        observation_mask=reference.obs_names.astype(str).isin(keep),
    )


def _reference_path(variant: str) -> Path:
    return REFERENCE_ROOT / variant / "atac/reference.h5ad"


def derive_two_bin_result(five: FragmentShapeResult) -> FragmentShapeResult:
    """Collapse the frozen five-bin recount without rereading raw fragments."""
    if tuple(five.bins) != FIVE_BINS:
        raise ValueError("Two-bin derivation requires the frozen five-bin axis.")
    layers = {
        TWO_BINS[0].layer: sparse.csr_matrix(
            five.layers[FIVE_BINS[0].layer] + five.layers[FIVE_BINS[1].layer]
        ),
        TWO_BINS[1].layer: sparse.csr_matrix(
            five.layers[FIVE_BINS[2].layer]
            + five.layers[FIVE_BINS[3].layer]
            + five.layers[FIVE_BINS[4].layer]
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
        bins=TWO_BINS,
        layers=layers,
        qc=qc,
        right_cut_offset=five.right_cut_offset,
    )


def _write_reference(
    variant: str,
    reference: ad.AnnData,
    metadata: Mapping[str, Any],
) -> Path:
    path = _reference_path(variant)
    descriptor = path.parents[1] / "reference.yaml"
    if path.is_file() and descriptor.is_file():
        existing = _read_yaml(descriptor)
        if existing.get("reference_id") != f"shapemix_pbmc_stress_v1/{variant}":
            raise ValueError(f"Stale reference descriptor: {descriptor}")
        return path
    if path.exists() or descriptor.exists() or path.parents[1].exists():
        raise FileExistsError(f"Partial sensitivity reference requires inspection: {path.parents[1]}")
    validate_fragment_shape_feature_axis(reference, f"sensitivity reference {variant}")
    _atomic_h5ad(reference, path)
    spec = FragmentShapeSpec.from_mapping(reference.uns["fragment_shape"])
    _write_yaml(
        descriptor,
        {
            "schema_version": 1,
            "reference_id": f"shapemix_pbmc_stress_v1/{variant}",
            "source_dataset_id": "pbmc_granulocyte_sorted_10k",
            "source_split": _repository_path(SPLIT_DIR),
            "labels_key": "cell_type",
            "cell_types": list(dict.fromkeys(reference.obs["cell_type"].astype(str))),
            "cells": reference.n_obs,
            "peaks": reference.n_vars,
            "bins": len(spec.bins),
            "feature_sha256": ordered_feature_sha256(reference.var_names),
            "atac": {"path": _repository_path(path)},
            "sensitivity": dict(metadata),
        },
    )
    return path


def _alternate_bin_objects(
    reference: ad.AnnData,
    heldout: ad.AnnData,
) -> dict[str, tuple[ad.AnnData, ad.AnnData, Path]]:
    result: dict[str, tuple[ad.AnnData, ad.AnnData, Path]] = {}
    cache_root = REUSABLE_ROOT / "fragment_shape_cache"
    cached_paths = {
        "five": cache_root / "bins_five_full_split_1103.h5ad",
        "two": cache_root / "bins_two_full_split_1103.h5ad",
    }
    if not all(path.is_file() for path in cached_paths.values()):
        if any(path.exists() for path in cached_paths.values()):
            raise FileExistsError("Partial alternate-bin cache requires inspection.")
        combined_obs = pd.concat([reference.obs, heldout.obs], axis=0)
        combined_names = combined_obs.index.astype(str).tolist()
        peaks = [
            (str(row.chrom), int(row.start), int(row.end), str(name))
            for name, row in reference.var.iterrows()
        ]
        five = count_fragment_shapes(
            RAW_FRAGMENTS,
            combined_names,
            peaks,
            right_cut_offset=0,
            bins=FIVE_BINS,
            chunk_size=1_000_000,
        )
        source_spec = FragmentShapeSpec.from_mapping(reference.uns["fragment_shape"])
        provenance = {
            "split_sha256": source_spec.split_sha256,
            "source_sha256": copy.deepcopy(source_spec.source_sha256),
            "coordinate_validation": copy.deepcopy(source_spec.coordinate_validation),
            "software_versions": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": __import__("scipy").__version__,
                "anndata": ad.__version__,
            },
        }
        five_full = build_fragment_shape_anndata(
            five,
            obs=combined_obs,
            var=reference.var,
            provenance=provenance,
        )
        two = derive_two_bin_result(five)
        two_full = build_fragment_shape_anndata(
            two,
            obs=combined_obs,
            var=reference.var,
            provenance=provenance,
        )
        _atomic_h5ad(five_full, cached_paths["five"])
        _atomic_h5ad(two_full, cached_paths["two"])
        _write_yaml(
            cache_root / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "raw_fragments": _repository_path(RAW_FRAGMENTS),
                "raw_fragment_sha256": _sha256(RAW_FRAGMENTS),
                "stream_count": 1,
                "five_bin_cache": _repository_path(cached_paths["five"]),
                "two_bin_cache": _repository_path(cached_paths["two"]),
                "two_bin_derivation": "exact sums of five-bin layers",
            },
        )

    for level, full_path in cached_paths.items():
        full = ad.read_h5ad(full_path)
        ref = subset_shape_cells(
            full,
            observation_mask=full.obs["split_pool"].astype(str).eq("reference"),
        )
        test = subset_shape_cells(
            full,
            observation_mask=full.obs["split_pool"].astype(str).eq("heldout"),
        )
        heldout_path = REUSABLE_ROOT / "splits" / f"bins_{level}" / "heldout_cells.h5ad"
        if not heldout_path.is_file():
            if heldout_path.exists():
                raise FileExistsError(heldout_path)
            _atomic_h5ad(test, heldout_path)
        result[level] = (ref, test, heldout_path)
    return result


def prepare_references() -> None:
    _validate_template()
    reference, heldout = _load_base()
    variants: list[tuple[str, ad.AnnData, Mapping[str, Any]]] = [
        ("anchor_all16_peaks5000_bins3", reference, {"factor": "anchor"}),
        ("features_all16_peaks1000_bins3", _subset(reference, num_features=1000), {"factor": "features", "peaks": 1000}),
        ("features_all16_peaks2500_bins3", _subset(reference, num_features=2500), {"factor": "features", "peaks": 2500}),
        ("subtype_broad3_peaks5000_bins3", _subset(reference, cell_types=BROAD3), {"factor": "subtype", "level": "broad3"}),
        ("subtype_cd4_related3_peaks5000_bins3", _subset(reference, cell_types=CD4_RELATED3), {"factor": "subtype", "level": "cd4_related3"}),
    ]
    broad_reference = _subset(reference, cell_types=BROAD3)
    for support in (50, 100, 250, None):
        token = "all" if support is None else str(support)
        variants.append(
            (
                f"reference_support_broad3_{token}_peaks5000_bins3",
                _cap_reference(broad_reference, support),
                {"factor": "reference_support", "per_type": support or "all"},
            )
        )
    alternate = _alternate_bin_objects(reference, heldout)
    for level in ("two", "five"):
        variants.append(
            (
                f"bins_{level}_all16_peaks5000",
                alternate[level][0],
                {"factor": "bins", "level": level},
            )
        )
    for variant, value, metadata in variants:
        path = _write_reference(variant, value, metadata)
        print(f"reference {variant} status=ready path={_repository_path(path)}", flush=True)


def _design_rows() -> list[DesignRow]:
    rows: list[DesignRow] = []
    for seed in MIXTURE_SEEDS:
        suffix = f"mix_{seed}"
        rows.append(DesignRow(f"pbmc_shapemix_stress_v1_anchor_standard_{suffix}", "anchor", "standard", seed, "standard", "anchor_all16_peaks5000_bins3", PBMC_CELL_TYPES, "observed_abundance"))
        for token, probability in (("0p25", 0.25), ("0p50", 0.50), ("0p75", 0.75)):
            rows.append(DesignRow(f"pbmc_shapemix_stress_v1_depth_keep_{token}_{suffix}", "depth", f"keep_{token}", seed, "keep_1p00", "anchor_all16_peaks5000_bins3", PBMC_CELL_TYPES, "observed_abundance", depth_retain_probability=probability))
        for cells in (2, 5, 20):
            rows.append(DesignRow(f"pbmc_shapemix_stress_v1_cells_mean_{cells}_{suffix}", "cells", f"mean_{cells}", seed, "mean_10", "anchor_all16_peaks5000_bins3", PBMC_CELL_TYPES, "observed_abundance", mean_cells_per_spot=float(cells)))
        for token, fraction in (("0p001", 0.001), ("0p005", 0.005), ("0p010", 0.010)):
            rows.append(DesignRow(f"pbmc_shapemix_stress_v1_rare_nk_fraction_{token}_{suffix}", "rare_nk", f"fraction_{token}", seed, "observed", "anchor_all16_peaks5000_bins3", PBMC_CELL_TYPES, "observed_abundance", rare_nk_fraction=fraction))
        for peaks in (1000, 2500):
            rows.append(DesignRow(f"pbmc_shapemix_stress_v1_features_peaks_{peaks}_{suffix}", "features", f"peaks_{peaks}", seed, "peaks_5000", f"features_all16_peaks{peaks}_bins3", PBMC_CELL_TYPES, "observed_abundance"))
        rows.extend(
            [
                DesignRow(f"pbmc_shapemix_stress_v1_subtype_broad3_{suffix}", "subtype", "broad3", seed, "broad3", "subtype_broad3_peaks5000_bins3", BROAD3, "equal_celltype"),
                DesignRow(f"pbmc_shapemix_stress_v1_subtype_cd4_related3_{suffix}", "subtype", "cd4_related3", seed, "broad3", "subtype_cd4_related3_peaks5000_bins3", CD4_RELATED3, "equal_celltype"),
            ]
        )
        for support in (50, 100, 250, None):
            token = "all" if support is None else f"per_type_{support}"
            ref_token = "all" if support is None else str(support)
            rows.append(DesignRow(f"pbmc_shapemix_stress_v1_reference_support_{token}_{suffix}", "reference_support", token, seed, "all", f"reference_support_broad3_{ref_token}_peaks5000_bins3", BROAD3, "equal_celltype"))
        for level in ("two", "five"):
            rows.append(DesignRow(f"pbmc_shapemix_stress_v1_bins_{level}_{suffix}", "bins", level, seed, "three", f"bins_{level}_all16_peaks5000", PBMC_CELL_TYPES, "observed_abundance"))
    if len(rows) != 40 or len({row.dataset_id for row in rows}) != 40:
        raise RuntimeError("PBMC sensitivity design must contain 40 unique datasets.")
    return rows


def _probabilities(row: DesignRow) -> dict[str, float]:
    if row.condition == "equal_celltype":
        return condition_probabilities("equal_celltype", row.cell_types)
    observed = condition_probabilities(
        "observed_abundance",
        PBMC_CELL_TYPES,
        dict(PBMC_CELL_TYPE_COUNTS),
    )
    if row.rare_nk_fraction is None:
        return observed
    scale = (1.0 - row.rare_nk_fraction) / (1.0 - observed["NK"])
    probabilities = {
        cell_type: probability * scale
        for cell_type, probability in observed.items()
    }
    probabilities["NK"] = row.rare_nk_fraction
    return probabilities


def _heldout_for_row(row: DesignRow, heldout: ad.AnnData) -> tuple[ad.AnnData, Path]:
    if row.factor == "bins":
        path = REUSABLE_ROOT / "splits" / f"bins_{row.level}" / "heldout_cells.h5ad"
        return ad.read_h5ad(path), path
    peaks = 5000
    if row.factor == "features":
        peaks = int(row.level.removeprefix("peaks_"))
    return _subset(heldout, cell_types=row.cell_types, num_features=peaks), SPLIT_DIR / "heldout_test_cells.h5ad"


def _enrich_dataset(dataset_path: Path, row: DesignRow) -> None:
    sensitivity = {
        "campaign": "shapemix_pbmc_stress_v1",
        "factor": row.factor,
        "level": row.level,
        "control_level": row.control_level,
        "outer_split_seed": OUTER_SPLIT_SEED,
        "mixture_seed": row.mixture_seed,
        "reference_variant": row.reference_variant,
        "mean_cells_per_spot": row.mean_cells_per_spot,
        "depth_retain_probability": row.depth_retain_probability,
        "rare_nk_fraction": row.rare_nk_fraction,
        "observed_nk_fraction": OBSERVED_NK_FRACTION,
    }
    descriptor = _read_yaml(dataset_path)
    if descriptor.get("sensitivity") not in (None, sensitivity):
        raise ValueError(f"Stale sensitivity metadata: {dataset_path}")
    descriptor["sensitivity"] = sensitivity
    _write_yaml(dataset_path, descriptor)
    manifest_path = dataset_path.parent / "simulation/manifest.yaml"
    manifest = _read_yaml(manifest_path)
    manifest["sensitivity"] = sensitivity
    record = manifest["outputs"]["dataset_yaml"]
    record["bytes"] = dataset_path.stat().st_size
    record["sha256"] = _sha256(dataset_path)
    _write_yaml(manifest_path, manifest)


def materialize_datasets() -> None:
    _validate_template()
    _, base_heldout = _load_base()
    rows = _design_rows()
    expected = _read_yaml(EXPERIMENT_PATH)["datasets"]
    if [row.dataset_id for row in rows] != expected:
        raise ValueError("Builder design order differs from the frozen experiment.")
    registry = _read_yaml(REGISTRY_PATH) if REGISTRY_PATH.exists() else {}
    design_records: list[dict[str, Any]] = []
    for row in rows:
        reference_path = _reference_path(row.reference_variant)
        if not reference_path.is_file():
            raise FileNotFoundError(f"Run the references stage first: {reference_path}")
        reference = ad.read_h5ad(reference_path)
        heldout, heldout_path = _heldout_for_row(row, base_heldout)
        if not reference.var_names.equals(heldout.var_names):
            raise ValueError(f"Reference/held-out axes differ for {row.dataset_id}.")
        dataset_path = DATASET_ROOT / row.dataset_id / "dataset.yaml"
        if not dataset_path.is_file():
            simulation = simulate_shapemix_spots(
                heldout,
                cell_types=row.cell_types,
                sampling_probabilities=_probabilities(row),
                condition=row.condition,
                outer_split_seed=OUTER_SPLIT_SEED,
                inner_mixture_seed=row.mixture_seed,
                num_spots=NUM_SPOTS,
                mean_cells_per_spot=row.mean_cells_per_spot,
                labels_key="cell_type",
                grid_shape=(32, 32),
                reference_barcodes=reference.obs_names.astype(str).tolist(),
                depth_retain_probability=row.depth_retain_probability,
            )
            dataset_path = write_simulation_dataset(
                simulation,
                output_root=DATASET_ROOT,
                dataset_id=row.dataset_id,
                reference_path=reference_path,
                heldout_path=heldout_path,
                split_manifest_path=SPLIT_DIR / "manifest.yaml",
                labels_key="cell_type",
                benchmark_scope="secondary_one_donor_diagnostic_sensitivity",
                source="pbmc_granulocyte_sorted_10k_shapemix_stress_v1",
                description=(
                    "Frozen one-factor-at-a-time PBMC ShapeMix diagnostic sensitivity; "
                    "pseudo-spots use only split-1103 held-out source cells."
                ),
                scientific_scope=(
                    "Conditional resampling within one PBMC Multiome donor; explanatory, "
                    "not donor-level generalization."
                ),
            )
            _enrich_dataset(dataset_path, row)
            print(f"dataset {row.dataset_id} status=completed", flush=True)
        else:
            descriptor = _read_yaml(dataset_path)
            if (descriptor.get("sensitivity") or {}).get("factor") != row.factor:
                raise ValueError(f"Existing dataset is not the frozen design: {dataset_path}")
            print(f"dataset {row.dataset_id} status=reused", flush=True)
        registry[row.dataset_id] = {"config": _repository_path(dataset_path)}
        design_records.append(
            {
                **row.__dict__,
                "cell_types": "|".join(row.cell_types),
                "reference_path": _repository_path(reference_path),
            }
        )

    design_path = REUSABLE_ROOT / "manifests/design.csv"
    design_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_design = design_path.with_name(".design.csv.tmp")
    pd.DataFrame.from_records(design_records).to_csv(temporary_design, index=False)
    temporary_design.replace(design_path)
    _write_yaml(
        REUSABLE_ROOT / "manifests/materialization.yaml",
        {
            "schema_version": 1,
            "status": "complete",
            "template": _repository_path(TEMPLATE_PATH),
            "experiment": _repository_path(EXPERIMENT_PATH),
            "datasets": len(rows),
            "design": {
                "path": _repository_path(design_path),
                "sha256": _sha256(design_path),
            },
        },
    )
    temporary_registry = REGISTRY_PATH.with_name(".datasets.yaml.pbmc_stress.tmp")
    with temporary_registry.open("w") as handle:
        yaml.safe_dump(registry, handle, sort_keys=False)
    temporary_registry.replace(REGISTRY_PATH)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("references", "datasets", "all"))
    args = parser.parse_args()
    if args.stage in {"references", "all"}:
        prepare_references()
    if args.stage in {"datasets", "all"}:
        materialize_datasets()


if __name__ == "__main__":
    main()
