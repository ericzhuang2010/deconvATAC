#!/usr/bin/env python3
"""Select outcome-blind marker intervals from standardized ShapeMix references."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

import anndata as ad
import numpy as np
from scipy import sparse
import yaml

from deconvatac.data import FragmentShapeSpec, ordered_feature_sha256
from deconvatac.data.validators import validate_fragment_shape_feature_axis


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = ROOT / "data/processed/shapemix/real_spatial_validation/marker_features_v1"
DEFAULT_REFERENCES = (
    "gse246791_adult_mouse_brain_broad9_v1",
    "gse216371_mouse_embryo_e13_major_types_v1",
    "gse244618_human_hippocampus_donor3_region3_v1",
)
MIN_NONZERO_CELLS = 10
MARKERS_PER_TYPE = 25
ROW_CHUNK = 4_096
RATE_PSEUDOCOUNT = 1.0e-4


def repository_path(path: Path) -> str:
    return path.absolute().relative_to(ROOT.absolute()).as_posix()


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def atomic_yaml(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(descriptor, "w") as handle:
            yaml.safe_dump(dict(value), handle, sort_keys=False)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def rank_marker_indices(
    mean_type: np.ndarray,
    mean_rest: np.ndarray,
    coverage: np.ndarray,
    totals: np.ndarray,
    feature_names: Iterable[str],
    *,
    n_markers: int = MARKERS_PER_TYPE,
) -> list[int]:
    """Rank reference-only marker features with deterministic identifier ties."""
    mean_type = np.asarray(mean_type, dtype=np.float64)
    mean_rest = np.asarray(mean_rest, dtype=np.float64)
    coverage = np.asarray(coverage, dtype=np.int64)
    totals = np.asarray(totals, dtype=np.float64)
    names = np.asarray([str(value) for value in feature_names], dtype=object)
    if not (
        mean_type.ndim == 1
        and mean_type.shape == mean_rest.shape == coverage.shape == totals.shape == names.shape
    ):
        raise ValueError("Marker-ranker inputs must be aligned one-dimensional arrays")
    score = np.log2(
        (mean_type + RATE_PSEUDOCOUNT) / (mean_rest + RATE_PSEUDOCOUNT)
    )
    eligible = np.flatnonzero(
        np.isfinite(score)
        & (coverage >= MIN_NONZERO_CELLS)
        & (mean_type > mean_rest)
        & (totals > 0)
    )
    if len(eligible) < n_markers:
        raise ValueError(
            f"Only {len(eligible)} features pass the marker-support gate; need {n_markers}"
        )
    ordered = sorted(
        eligible.tolist(),
        key=lambda index: (
            -float(score[index]),
            -int(coverage[index]),
            -float(totals[index]),
            str(names[index]),
        ),
    )
    return ordered[:n_markers]


def aggregate_reference(
    reference: ad.AnnData,
    labels_key: str,
    cell_types: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if labels_key not in reference.obs:
        raise KeyError(f"Reference lacks labels_key={labels_key!r}")
    labels = reference.obs[labels_key].astype(str).to_numpy()
    unknown = sorted(set(labels).difference(cell_types))
    if unknown:
        raise ValueError(f"Reference has labels outside the frozen universe: {unknown}")
    type_index = {cell_type: index for index, cell_type in enumerate(cell_types)}
    totals = np.zeros((len(cell_types), reference.n_vars), dtype=np.float64)
    coverage = np.zeros((len(cell_types), reference.n_vars), dtype=np.int64)
    cells = np.zeros(len(cell_types), dtype=np.int64)
    for start in range(0, reference.n_obs, ROW_CHUNK):
        stop = min(start + ROW_CHUNK, reference.n_obs)
        block = sparse.csr_matrix(reference.X[start:stop])
        block.sum_duplicates()
        block.eliminate_zeros()
        block_labels = labels[start:stop]
        for cell_type, index in type_index.items():
            selected = np.flatnonzero(block_labels == cell_type)
            if len(selected) == 0:
                continue
            subset = block[selected]
            totals[index] += np.asarray(subset.sum(axis=0)).ravel()
            coverage[index] += np.asarray(subset.getnnz(axis=0)).ravel()
            cells[index] += len(selected)
        if stop % (ROW_CHUNK * 25) == 0 or stop == reference.n_obs:
            print(
                f"reference_marker_scan cells={stop}/{reference.n_obs}",
                flush=True,
            )
    if (cells == 0).any():
        missing = [cell_types[index] for index in np.flatnonzero(cells == 0)]
        raise ValueError(f"Frozen reference types have no cells: {missing}")
    return totals, coverage, cells


def build_reference_markers(reference_id: str) -> Path:
    output_root = OUTPUT_ROOT / reference_id
    marker_path = output_root / "marker_features.yaml"
    manifest_path = output_root / "manifest.yaml"
    if marker_path.is_file() and manifest_path.is_file():
        print(f"reference_markers reference={reference_id} status=reused", flush=True)
        return marker_path
    if output_root.exists():
        raise FileExistsError(f"Partial immutable reference-marker directory: {output_root}")

    reference_root = ROOT / "data/processed/references" / reference_id
    reference_manifest_path = reference_root / "reference.yaml"
    if not reference_manifest_path.is_file():
        raise FileNotFoundError(reference_manifest_path)
    manifest = read_yaml(reference_manifest_path)
    if str(manifest.get("reference_id")) != reference_id:
        raise ValueError(f"Reference manifest ID mismatch: {reference_id}")
    reference_path = ROOT / str(manifest["modalities"]["atac"]["path"])
    if not reference_path.is_file():
        raise FileNotFoundError(reference_path)
    labels_key = str(manifest["labels_key"])
    cell_types = [str(value) for value in manifest["cell_types"]]
    reference = ad.read_h5ad(reference_path, backed="r")
    try:
        if reference.n_vars != 5_000:
            raise ValueError(f"Expected 5,000 reference features: {reference_id}")
        validate_fragment_shape_feature_axis(reference, f"{reference_id} marker source")
        shape_spec = FragmentShapeSpec.from_mapping(reference.uns["fragment_shape"])
        names = reference.var_names.astype(str).tolist()
        totals, coverage, cells = aggregate_reference(reference, labels_key, cell_types)
    finally:
        reference.file.close()
    all_totals = totals.sum(axis=0)
    all_cells = int(cells.sum())
    markers: dict[str, Any] = {}
    for type_index, cell_type in enumerate(cell_types):
        mean_type = totals[type_index] / int(cells[type_index])
        rest_cells = all_cells - int(cells[type_index])
        if rest_cells <= 0:
            raise ValueError("Marker selection requires at least two reference cell types")
        mean_rest = (all_totals - totals[type_index]) / rest_cells
        selected = rank_marker_indices(
            mean_type,
            mean_rest,
            coverage[type_index],
            totals[type_index],
            names,
        )
        markers[cell_type] = {
            "features": [names[index] for index in selected],
            "log2_rate_fold_change": [
                float(
                    np.log2(
                        (mean_type[index] + RATE_PSEUDOCOUNT)
                        / (mean_rest[index] + RATE_PSEUDOCOUNT)
                    )
                )
                for index in selected
            ],
            "nonzero_reference_cells": [int(coverage[type_index, index]) for index in selected],
        }

    output_root.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{reference_id}.", dir=output_root.parent))
    try:
        marker_value = {
            "schema_version": 1,
            "marker_version": "reference_count_rate_log2fc_v1",
            "reference_id": reference_id,
            "source": "standardized_reference_only",
            "outcome_data_used": False,
            "selection": {
                "markers_per_cell_type": MARKERS_PER_TYPE,
                "minimum_nonzero_reference_cells": MIN_NONZERO_CELLS,
                "score": "log2((mean_type+1e-4)/(mean_all_other_types+1e-4))",
                "eligibility": "mean_type_gt_mean_rest",
                "tie_breaks": [
                    "score_descending",
                    "nonzero_cell_coverage_descending",
                    "type_total_descending",
                    "feature_identifier_ascending",
                ],
            },
            "cell_types": cell_types,
            "reference_cell_counts": {
                cell_type: int(cells[index]) for index, cell_type in enumerate(cell_types)
            },
            "feature_sha256": ordered_feature_sha256(names),
            "fragment_shape": shape_spec.to_dict(omit_none=True),
            "markers": markers,
        }
        atomic_yaml(temporary / "marker_features.yaml", marker_value)
        atomic_yaml(
            temporary / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "completed_at": datetime.now(timezone.utc).isoformat(),
                "reference_id": reference_id,
                "reference_manifest": repository_path(reference_manifest_path),
                "reference_manifest_sha256": digest(reference_manifest_path),
                "reference_h5ad": repository_path(reference_path),
                "reference_h5ad_bytes": reference_path.stat().st_size,
                "marker_features": repository_path(output_root / "marker_features.yaml"),
                "marker_features_sha256": digest(temporary / "marker_features.yaml"),
                "outcome_data_used": False,
            },
        )
        temporary.rename(output_root)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(
        f"reference_markers reference={reference_id} cell_types={len(cell_types)} "
        f"markers_per_type={MARKERS_PER_TYPE} status=completed",
        flush=True,
    )
    return marker_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-id", action="append")
    return parser.parse_args()


def main() -> None:
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run marker selection through scripts/run_shapemix_low_impact.sh")
    args = parse_args()
    for reference_id in args.reference_id or DEFAULT_REFERENCES:
        build_reference_markers(reference_id)


if __name__ == "__main__":
    main()
