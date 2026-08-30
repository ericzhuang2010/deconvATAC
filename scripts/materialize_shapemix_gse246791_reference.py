#!/usr/bin/env python3
"""Materialize the audited GSE246791 adult mouse-brain ShapeMix reference."""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import os
import platform
import re
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
import scipy
import yaml

from deconvatac.data import FragmentShapeSpec, ordered_feature_sha256
from deconvatac.data.validators import (
    validate_fragment_shape_feature_axis,
    validate_fragment_shape_spec,
)
from deconvatac.pp.fragment_shapes import (
    build_fragment_shape_anndata,
    count_fragment_shapes_from_records,
)
from scripts.build_shapemix_gse246791_reference import (
    CELL_TYPES,
    LABEL_VERSION,
    build_labels,
)
from scripts.prepare_shapemix_gse246791 import CONFIG, selected_sources


ROOT = Path(__file__).resolve().parents[1]
FAMILY_ROOT = ROOT / str(CONFIG["processed_directory"])
LABEL_ROOT = FAMILY_ROOT / "labels" / LABEL_VERSION
AXIS_ROOT = FAMILY_ROOT / "feature_axes" / LABEL_VERSION
CACHE_ROOT = FAMILY_ROOT / "fragment_shape_cache" / LABEL_VERSION
AUDIT_PATH = FAMILY_ROOT / "manifests/representative_coordinate_audit.yaml"
GATE_PATH = FAMILY_ROOT / "manifests/adult_reference_preacquisition_gate.yaml"
REFERENCE_ROOT = ROOT / "data/processed/references" / str(CONFIG["standardized_reference_id"])
N_TOP_PEAKS = 5_000
MIN_REFERENCE_CELLS_PER_PEAK = 10
H5_AXIS_CHUNK = 100_000
COUNT_CHUNK_SIZE = 1_000_000
INTERVAL_PATTERN = re.compile(r"^([^:]+):(\d+)-(\d+)$")


def repository_path(path: Path) -> str:
    return path.absolute().relative_to(ROOT.absolute()).as_posix()


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


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


def software_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "h5py": h5py.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }


def decode_strings(values: Iterable[Any]) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def parse_interval(value: str) -> tuple[str, int, int]:
    match = INTERVAL_PATTERN.fullmatch(value)
    if match is None:
        raise ValueError(f"Invalid deposited interval identifier: {value!r}")
    chrom, start_text, end_text = match.groups()
    start, end = int(start_text), int(end_text)
    if start < 0 or end <= start:
        raise ValueError(f"Invalid deposited interval coordinates: {value!r}")
    return chrom, start, end


def hash_h5_axis(dataset: h5py.Dataset) -> str:
    digest = hashlib.sha256()
    for start in range(0, len(dataset), H5_AXIS_CHUNK):
        for value in decode_strings(dataset[start : start + H5_AXIS_CHUNK]):
            encoded = value.encode("utf-8")
            digest.update(len(encoded).to_bytes(8, byteorder="big", signed=False))
            digest.update(encoded)
    return digest.hexdigest()


def h5_axis_values(dataset: h5py.Dataset, indices: np.ndarray) -> dict[int, str]:
    ordered = np.asarray(sorted({int(value) for value in indices}), dtype=np.int64)
    if len(ordered) == 0:
        return {}
    values = decode_strings(dataset[ordered])
    return dict(zip(ordered.tolist(), values, strict=True))


def load_labels() -> pd.DataFrame:
    cells_path = LABEL_ROOT / "cells.tsv.gz"
    if not cells_path.is_file():
        build_labels()
    labels = pd.read_csv(cells_path, sep="\t", dtype={"sample": str, "barcode": str})
    required = {"cell_id", "sample", "gsm", "major_region", "barcode", "cell_type"}
    missing = sorted(required.difference(labels.columns))
    if missing:
        raise ValueError(f"GSE246791 labels lack required columns: {missing}")
    if labels["cell_id"].duplicated().any():
        raise ValueError("GSE246791 retained cell IDs are not unique")
    if tuple(dict.fromkeys(labels["cell_type"].astype(str))) != CELL_TYPES:
        observed = set(labels["cell_type"].astype(str))
        if observed != set(CELL_TYPES):
            raise ValueError(f"GSE246791 label universe changed: {sorted(observed)}")
    return labels


def rank_selected_indices(
    score: np.ndarray,
    coverage: np.ndarray,
    total: np.ndarray,
    *,
    n_top: int,
    id_loader: Callable[[np.ndarray], Mapping[int, str]],
) -> list[int]:
    score = np.asarray(score, dtype=np.float64)
    coverage = np.asarray(coverage, dtype=np.int64)
    total = np.asarray(total, dtype=np.int64)
    if score.shape != coverage.shape or score.shape != total.shape or score.ndim != 1:
        raise ValueError("Ranker arrays must be aligned one-dimensional vectors")
    eligible = np.flatnonzero(
        np.isfinite(score) & (coverage >= MIN_REFERENCE_CELLS_PER_PEAK)
    )
    if len(eligible) < n_top:
        raise ValueError(f"Only {len(eligible)} features pass the support gate")
    numeric_order = np.lexsort(
        (
            eligible,
            -total[eligible],
            -coverage[eligible],
            -score[eligible],
        )
    )
    boundary = int(eligible[numeric_order[n_top - 1]])
    boundary_score = score[boundary]
    boundary_coverage = coverage[boundary]
    boundary_total = total[boundary]
    strict_mask = (
        (score[eligible] > boundary_score)
        | (
            (score[eligible] == boundary_score)
            & (coverage[eligible] > boundary_coverage)
        )
        | (
            (score[eligible] == boundary_score)
            & (coverage[eligible] == boundary_coverage)
            & (total[eligible] > boundary_total)
        )
    )
    strict = eligible[strict_mask]
    tied = eligible[
        (score[eligible] == boundary_score)
        & (coverage[eligible] == boundary_coverage)
        & (total[eligible] == boundary_total)
    ]
    remaining = n_top - len(strict)
    if remaining < 1 or remaining > len(tied):
        raise RuntimeError("Numeric feature-rank boundary is inconsistent")
    tied_ids = id_loader(tied)
    chosen_ties = sorted(tied.tolist(), key=lambda index: tied_ids[index])[:remaining]
    selected = np.concatenate(
        [strict.astype(np.int64), np.asarray(chosen_ties, dtype=np.int64)]
    )
    selected_ids = id_loader(selected)
    return sorted(
        selected.tolist(),
        key=lambda index: (
            -score[index],
            -coverage[index],
            -total[index],
            selected_ids[index],
        ),
    )


def aggregate_deposited_feature_statistics(
    labels: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, str, Path]:
    sources = selected_sources()
    type_index = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    counts: np.ndarray | None = None
    coverage: np.ndarray | None = None
    expected_axis_hash: str | None = None
    axis_source: Path | None = None
    for source in sources:
        if not source.output_path.is_file():
            raise FileNotFoundError(
                f"Run prepare_shapemix_gse246791.py first: {source.output_path}"
            )
        sample_labels = labels[labels["sample"] == source.sample]
        barcode_type = {
            str(barcode): type_index[str(cell_type)]
            for barcode, cell_type in zip(
                sample_labels["barcode"], sample_labels["cell_type"], strict=True
            )
        }
        with h5py.File(source.output_path, "r") as handle:
            x = handle["X"]
            shape = tuple(int(value) for value in x.attrs["shape"])
            barcodes = decode_strings(handle["obs/index"][:])
            if shape[0] != len(barcodes):
                raise ValueError(f"Deposited CSR observation shape changed: {source.output_path}")
            if counts is None:
                counts = np.zeros((len(CELL_TYPES), shape[1]), dtype=np.uint64)
                coverage = np.zeros(shape[1], dtype=np.uint32)
                axis_source = source.output_path
            elif shape[1] != counts.shape[1]:
                raise ValueError("GSE246791 deposited feature counts do not share one axis")
            axis_hash = hash_h5_axis(handle["var/index"])
            if expected_axis_hash is None:
                expected_axis_hash = axis_hash
            elif axis_hash != expected_axis_hash:
                raise ValueError(f"GSE246791 deposited feature axis changed in {source.gsm}")
            indptr = np.asarray(x["indptr"][:], dtype=np.int64)
            if len(indptr) != len(barcodes) + 1 or indptr[0] != 0:
                raise ValueError(f"Invalid deposited CSR pointers: {source.output_path}")
            for row_index, barcode in enumerate(barcodes):
                group = barcode_type.get(barcode)
                if group is None:
                    continue
                start, end = int(indptr[row_index]), int(indptr[row_index + 1])
                indices = np.asarray(x["indices"][start:end], dtype=np.int64)
                values = np.asarray(x["data"][start:end], dtype=np.uint64)
                if len(indices) != len(values):
                    raise ValueError("Deposited CSR indices and data do not align")
                if len(indices) and (
                    np.any(indices[1:] <= indices[:-1])
                    or indices[0] < 0
                    or indices[-1] >= shape[1]
                    or np.any(values <= 0)
                ):
                    raise ValueError(f"Invalid deposited CSR row for {source.gsm}:{barcode}")
                counts[group, indices] += values
                coverage[indices] += 1
        print(
            f"gse246791 feature_statistics gsm={source.gsm} cells={len(sample_labels)} "
            "status=completed",
            flush=True,
        )
    if counts is None or coverage is None or expected_axis_hash is None or axis_source is None:
        raise ValueError("No GSE246791 deposited source objects were aggregated")
    return counts, coverage, expected_axis_hash, axis_source


def build_feature_axis() -> pd.DataFrame:
    selected_path = AXIS_ROOT / "selected_500bp_intervals.tsv.gz"
    manifest_path = AXIS_ROOT / "manifest.yaml"
    if selected_path.is_file() and manifest_path.is_file():
        print("gse246791 feature_axis status=reused", flush=True)
        return pd.read_csv(selected_path, sep="\t")
    if AXIS_ROOT.exists():
        raise FileExistsError(f"Partial immutable feature axis: {AXIS_ROOT}")
    labels = load_labels()
    counts, coverage, axis_hash, axis_source = aggregate_deposited_feature_statistics(labels)
    type_totals = counts.sum(axis=1, dtype=np.float64)
    if np.any(type_totals <= 0) or not np.isfinite(type_totals).all():
        raise ValueError("A GSE246791 broad class has no deposited feature counts")
    score = np.empty(counts.shape[1], dtype=np.float64)
    for start in range(0, counts.shape[1], H5_AXIS_CHUNK):
        stop = min(start + H5_AXIS_CHUNK, counts.shape[1])
        normalized = np.log2(
            1.0 + 1.0e4 * counts[:, start:stop].astype(np.float64) / type_totals[:, None]
        )
        score[start:stop] = np.var(normalized, axis=0, ddof=0)
    total = counts.sum(axis=0, dtype=np.uint64)
    with h5py.File(axis_source, "r") as handle:
        axis_dataset = handle["var/index"]

        def load_ids(indices: np.ndarray) -> Mapping[int, str]:
            return h5_axis_values(axis_dataset, indices)

        selected_indices = rank_selected_indices(
            score,
            coverage,
            total,
            n_top=N_TOP_PEAKS,
            id_loader=load_ids,
        )
        selected_ids = h5_axis_values(
            axis_dataset, np.asarray(selected_indices, dtype=np.int64)
        )
    rows = []
    for rank, index in enumerate(selected_indices, start=1):
        peak_id = selected_ids[index]
        chrom, interval_start, interval_end = parse_interval(peak_id)
        rows.append(
            {
                "rank": rank,
                "source_feature_index": index,
                "peak_id": peak_id,
                "chrom": chrom,
                "start": interval_start,
                "end": interval_end,
                "score": float(score[index]),
                "nonzero_reference_cells": int(coverage[index]),
                "total_reference_counts": int(total[index]),
            }
        )
    selected = pd.DataFrame.from_records(rows)
    AXIS_ROOT.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{LABEL_VERSION}.", dir=AXIS_ROOT.parent))
    try:
        selected.to_csv(
            temporary / "selected_500bp_intervals.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        (temporary / "selected_peaks.txt").write_text(
            "\n".join(selected["peak_id"].astype(str)) + "\n"
        )
        atomic_yaml(
            temporary / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "axis_version": LABEL_VERSION,
                "candidate_intervals": int(counts.shape[1]),
                "eligible_intervals": int(
                    np.count_nonzero(coverage >= MIN_REFERENCE_CELLS_PER_PEAK)
                ),
                "selected_intervals": N_TOP_PEAKS,
                "feature_sha256": ordered_feature_sha256(
                    selected["peak_id"].astype(str)
                ),
                "deposited_axis_sha256": axis_hash,
                "label_ontology": list(CELL_TYPES),
                "selector": {
                    "minimum_nonzero_reference_cells": MIN_REFERENCE_CELLS_PER_PEAK,
                    "score": "population_variance_log2_1_plus_scaled_type_rate",
                    "scale": 1.0e4,
                    "tie_breaks": [
                        "score_desc",
                        "coverage_desc",
                        "total_count_desc",
                        "peak_id_asc",
                    ],
                    "outcome_data_used": False,
                },
                "inputs": {
                    "labels": repository_path(LABEL_ROOT / "cells.tsv.gz"),
                    "source_objects": [
                        repository_path(source.output_path)
                        for source in selected_sources()
                    ],
                },
                "software_versions": software_versions(),
            },
        )
        temporary.rename(AXIS_ROOT)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print("gse246791 feature_axis peaks=5000 status=completed", flush=True)
    return selected


def fragment_path(gsm: str) -> Path:
    return FAMILY_ROOT / "normalized_fragments" / gsm / "fragments.tsv.gz"


def fragment_manifest(gsm: str) -> dict[str, Any]:
    path = fragment_path(gsm).parent / "manifest.yaml"
    record = read_yaml(path)
    output = record.get("output", {})
    if (
        not fragment_path(gsm).is_file()
        or output.get("path") != repository_path(fragment_path(gsm))
        or fragment_path(gsm).stat().st_size != int(output.get("bytes", -1))
        or not record.get("concordance", {}).get("passed", False)
    ):
        raise ValueError(f"Invalid normalized-fragment manifest for {gsm}")
    return record


def count_source(
    gsm: str,
    barcodes: list[str],
    peaks: list[tuple[str, int, int, str]],
    *,
    right_cut_offset: int,
):
    fragment_manifest(gsm)
    with gzip.open(fragment_path(gsm), "rt") as handle:
        return count_fragment_shapes_from_records(
            handle,
            barcodes,
            peaks,
            right_cut_offset=right_cut_offset,
            chunk_size=COUNT_CHUNK_SIZE,
        )


def deposited_selected_matrix(
    h5ad: Path,
    barcodes: list[str],
    source_feature_indices: np.ndarray,
) -> sparse.csr_matrix:
    with h5py.File(h5ad, "r") as handle:
        x = handle["X"]
        source_barcodes = decode_strings(handle["obs/index"][:])
        source_rows = {barcode: index for index, barcode in enumerate(source_barcodes)}
        missing = sorted(set(barcodes).difference(source_rows))
        if missing:
            raise ValueError(f"{len(missing)} retained barcodes are absent from {h5ad}")
        n_features = int(x.attrs["shape"][1])
        selected_map = np.full(n_features, -1, dtype=np.int32)
        selected_map[source_feature_indices] = np.arange(
            len(source_feature_indices), dtype=np.int32
        )
        indptr = np.asarray(x["indptr"][:], dtype=np.int64)
        rows: list[np.ndarray] = []
        columns: list[np.ndarray] = []
        values: list[np.ndarray] = []
        for output_row, barcode in enumerate(barcodes):
            source_row = source_rows[barcode]
            start, end = int(indptr[source_row]), int(indptr[source_row + 1])
            source_columns = np.asarray(x["indices"][start:end], dtype=np.int64)
            mapped = selected_map[source_columns]
            keep = mapped >= 0
            if np.any(keep):
                columns.append(mapped[keep].astype(np.int32))
                values.append(np.asarray(x["data"][start:end], dtype=np.int64)[keep])
                rows.append(np.full(np.count_nonzero(keep), output_row, dtype=np.int32))
    if not rows:
        return sparse.csr_matrix((len(barcodes), len(source_feature_indices)), dtype=np.int64)
    return sparse.coo_matrix(
        (np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
        shape=(len(barcodes), len(source_feature_indices)),
        dtype=np.int64,
    ).tocsr()


def matrix_concordance(
    observed: sparse.csr_matrix,
    deposited: sparse.csr_matrix,
) -> dict[str, Any]:
    observed = sparse.csr_matrix(observed, dtype=np.int64)
    deposited = sparse.csr_matrix(deposited, dtype=np.int64)
    if observed.shape != deposited.shape:
        raise ValueError("Observed and deposited representative matrices do not align")
    difference = (observed - deposited).tocsr()
    l1_error = int(np.abs(difference.data).sum(dtype=np.int64))
    deposited_total = int(deposited.sum())
    observed_total = int(observed.sum())
    cell_observed = np.asarray(observed.sum(axis=1)).ravel()
    cell_deposited = np.asarray(deposited.sum(axis=1)).ravel()
    feature_observed = np.asarray(observed.sum(axis=0)).ravel()
    feature_deposited = np.asarray(deposited.sum(axis=0)).ravel()
    return {
        "observed_total": observed_total,
        "deposited_total": deposited_total,
        "l1_error": l1_error,
        "normalized_l1_error": l1_error / max(deposited_total, 1),
        "cell_total_spearman_r": float(
            spearmanr(cell_deposited, cell_observed).statistic
        ),
        "feature_total_spearman_r": float(
            spearmanr(feature_deposited, feature_observed).statistic
        ),
        "exact_matrix_entries": int(
            observed.shape[0] * observed.shape[1] - difference.count_nonzero()
        ),
        "matrix_entries": int(observed.shape[0] * observed.shape[1]),
    }


def write_preacquisition_gate(audit: Mapping[str, Any]) -> None:
    if not audit.get("passed", False):
        raise ValueError("Cannot write an adult reference gate for a failed audit")
    label_manifest = LABEL_ROOT / "manifest.yaml"
    feature_manifest = AXIS_ROOT / "manifest.yaml"
    required = (label_manifest, feature_manifest, AUDIT_PATH)
    missing = [repository_path(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Adult reference gate inputs are absent: {missing}")
    atomic_yaml(
        GATE_PATH,
        {
            "schema_version": 1,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "label_manifest": repository_path(label_manifest),
            "label_manifest_sha256": file_digest(label_manifest),
            "feature_manifest": repository_path(feature_manifest),
            "feature_manifest_sha256": file_digest(feature_manifest),
            "coordinate_audit": repository_path(AUDIT_PATH),
            "coordinate_audit_sha256": file_digest(AUDIT_PATH),
            "selected_right_cut_offset": int(audit["selected_right_cut_offset"]),
            "outcome_data_used": False,
        },
    )


def build_coordinate_audit() -> int:
    if AUDIT_PATH.is_file():
        record = read_yaml(AUDIT_PATH)
        if not record.get("passed", False):
            raise ValueError("Stored representative coordinate audit is gated")
        write_preacquisition_gate(record)
        print("gse246791 coordinate_audit status=reused", flush=True)
        return int(record["selected_right_cut_offset"])
    selected = build_feature_axis()
    source = next(value for value in selected_sources() if value.gsm == "GSM7877011")
    labels = load_labels()
    sample_labels = labels[labels["gsm"] == source.gsm].sort_values("barcode")
    barcodes = sample_labels["barcode"].astype(str).tolist()
    peaks = list(
        selected[["chrom", "start", "end", "peak_id"]].itertuples(index=False, name=None)
    )
    deposited = deposited_selected_matrix(
        source.output_path,
        barcodes,
        selected["source_feature_index"].to_numpy(dtype=np.int64),
    )
    metrics: dict[int, dict[str, Any]] = {}
    for offset in (0, -1):
        result = count_source(
            source.gsm,
            barcodes,
            peaks,
            right_cut_offset=offset,
        )
        observed = sparse.csr_matrix(deposited.shape, dtype=np.int64)
        for layer in result.layers.values():
            observed = (observed + sparse.csr_matrix(layer, dtype=np.int64)).tocsr()
        metrics[offset] = matrix_concordance(observed, deposited)
        print(
            f"gse246791 coordinate_audit offset={offset} "
            f"normalized_l1={metrics[offset]['normalized_l1_error']:.6f}",
            flush=True,
        )
    selected_offset = min(
        metrics,
        key=lambda offset: (
            metrics[offset]["normalized_l1_error"],
            0 if offset == 0 else 1,
        ),
    )
    thresholds = CONFIG["preprocessing_policy"]["representative_matrix_concordance_gate"]
    selected_metrics = metrics[selected_offset]
    gates = {
        "normalized_l1_error": selected_metrics["normalized_l1_error"]
        <= float(thresholds["normalized_l1_error_max"]),
        "cell_total_spearman_r": selected_metrics["cell_total_spearman_r"]
        >= float(thresholds["cell_total_spearman_r_min"]),
        "feature_total_spearman_r": selected_metrics["feature_total_spearman_r"]
        >= float(thresholds["feature_total_spearman_r_min"]),
        "selected_offset_not_worse": all(
            selected_metrics["normalized_l1_error"]
            <= value["normalized_l1_error"]
            for value in metrics.values()
        ),
    }
    record = {
        "schema_version": 1,
        "status": "complete" if all(gates.values()) else "gated",
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "representative_gsm": source.gsm,
        "selected_right_cut_offset": selected_offset,
        "offset_metrics": {str(key): value for key, value in metrics.items()},
        "thresholds": dict(thresholds),
        "gates": gates,
        "passed": all(gates.values()),
        "inputs": {
            "feature_axis": repository_path(
                AXIS_ROOT / "selected_500bp_intervals.tsv.gz"
            ),
            "normalized_fragments": repository_path(fragment_path(source.gsm)),
            "deposited_h5ad": repository_path(source.output_path),
        },
        "outcome_data_used": False,
    }
    atomic_yaml(AUDIT_PATH, record)
    if not record["passed"]:
        raise ValueError(f"Representative coordinate audit failed: {gates}")
    write_preacquisition_gate(record)
    print(
        f"gse246791 coordinate_audit selected_offset={selected_offset} status=completed",
        flush=True,
    )
    return selected_offset


def combined_metadata(
    values: Iterable[Mapping[str, Any]],
    layers: Mapping[str, sparse.csr_matrix],
    feature_names: Iterable[str],
    label_sha256: str,
    source_hashes: Mapping[str, str],
) -> dict[str, Any]:
    records = [copy.deepcopy(dict(value)) for value in values]
    if not records:
        raise ValueError("Cannot combine zero GSE246791 shape metadata records")
    metadata = records[0]
    counters: dict[str, int] = {}
    for value in records:
        for key, count in dict(value.get("preprocessing_counters", {})).items():
            if isinstance(count, (int, np.integer)):
                counters[key] = counters.get(key, 0) + int(count)
    layer_totals = {name: int(matrix.sum()) for name, matrix in layers.items()}
    metadata["preprocessing_counters"] = counters
    metadata["matrix_counters"] = {
        "assigned_cut_sites": int(sum(layer_totals.values())),
        **{f"cut_sites_per_bin.{name}": value for name, value in layer_totals.items()},
    }
    metadata["feature_sha256"] = ordered_feature_sha256(feature_names)
    metadata["split_sha256"] = label_sha256
    metadata["source_sha256"] = dict(source_hashes)
    return metadata


def build_sample_cache(
    source,
    labels: pd.DataFrame,
    selected: pd.DataFrame,
    *,
    right_cut_offset: int,
) -> Path:
    output_dir = CACHE_ROOT / source.gsm
    output_path = output_dir / "cells.h5ad"
    manifest_path = output_dir / "manifest.yaml"
    if output_path.is_file() and manifest_path.is_file():
        print(f"gse246791 fragment_cache gsm={source.gsm} status=reused", flush=True)
        return output_path
    if output_dir.exists():
        raise FileExistsError(f"Partial immutable GSE246791 cache: {output_dir}")
    sample_labels = labels[labels["gsm"] == source.gsm].sort_values("barcode").copy()
    if sample_labels.empty:
        raise ValueError(f"No retained broad labels for {source.gsm}")
    barcodes = sample_labels["barcode"].astype(str).tolist()
    peaks = list(
        selected[["chrom", "start", "end", "peak_id"]].itertuples(index=False, name=None)
    )
    fragment_record = fragment_manifest(source.gsm)
    result = count_source(
        source.gsm,
        barcodes,
        peaks,
        right_cut_offset=right_cut_offset,
    )
    obs = sample_labels.set_index(sample_labels["barcode"].astype(str))
    obs.index.name = "fragment_barcode"
    var = selected.set_index("peak_id")[["chrom", "start", "end"]].copy()
    shape = build_fragment_shape_anndata(
        result,
        obs=obs,
        var=var,
        provenance={
            "split_sha256": file_digest(LABEL_ROOT / "cells.tsv.gz"),
            "source_sha256": {
                fragment_path(source.gsm).name: str(
                    fragment_record["output"]["sha256"]
                )
            },
            "coordinate_validation": {
                "selected_right_cut_offset": right_cut_offset,
                "matrix_match": "representative_recovered_read_concordance",
                "validation_method": "deposited_matrix_concordance_thresholds",
                "passed": True,
                "audit": repository_path(AUDIT_PATH),
            },
            "software_versions": software_versions(),
        },
    )
    shape.obs_names = pd.Index(sample_labels["cell_id"].astype(str), name="cell_id")
    output_dir.mkdir(parents=True)
    temporary = output_dir / ".cells.h5ad.tmp"
    try:
        shape.write_h5ad(temporary, compression="gzip")
        os.replace(temporary, output_path)
        atomic_yaml(
            manifest_path,
            {
                "schema_version": 1,
                "status": "complete",
                "gsm": source.gsm,
                "sample": source.sample,
                "major_region": source.major_region,
                "cells": shape.n_obs,
                "peaks": shape.n_vars,
                "right_cut_offset": right_cut_offset,
                "source_fragment_sha256": fragment_record["output"]["sha256"],
                "preprocessing_counters": result.qc.to_dict(),
                "output": repository_path(output_path),
            },
        )
    except BaseException:
        temporary.unlink(missing_ok=True)
        shutil.rmtree(output_dir, ignore_errors=True)
        raise
    print(
        f"gse246791 fragment_cache gsm={source.gsm} cells={shape.n_obs} "
        "status=completed",
        flush=True,
    )
    return output_path


def requested_sources(gsms: list[str] | None = None):
    sources = selected_sources()
    if not gsms:
        return sources
    requested = set(gsms)
    result = [source for source in sources if source.gsm in requested]
    missing = sorted(requested.difference(source.gsm for source in result))
    if missing:
        raise ValueError(f"Unknown selected GSM accession(s): {', '.join(missing)}")
    return result


def build_fragment_caches(gsms: list[str] | None = None) -> list[Path]:
    labels = load_labels()
    selected = build_feature_axis()
    right_cut_offset = build_coordinate_audit()
    return [
        build_sample_cache(
            source,
            labels,
            selected,
            right_cut_offset=right_cut_offset,
        )
        for source in requested_sources(gsms)
    ]




def build_reference() -> Path:
    reference_path = REFERENCE_ROOT / "atac/reference.h5ad"
    manifest_path = REFERENCE_ROOT / "reference.yaml"
    if reference_path.is_file() and manifest_path.is_file():
        print("gse246791 reference status=reused", flush=True)
        return reference_path
    if REFERENCE_ROOT.exists():
        raise FileExistsError(f"Partial immutable reference directory: {REFERENCE_ROOT}")
    caches = build_fragment_caches()
    selected = build_feature_axis()
    obs_parts: list[pd.DataFrame] = []
    layer_parts: dict[str, list[sparse.csr_matrix]] = {}
    metadata_parts: list[Mapping[str, Any]] = []
    source_hashes: dict[str, str] = {}
    for source, cache_path in zip(selected_sources(), caches, strict=True):
        cache = ad.read_h5ad(cache_path)
        obs_parts.append(cache.obs.copy())
        spec = FragmentShapeSpec.from_mapping(cache.uns["fragment_shape"])
        for layer_name in spec.layer_names:
            layer_parts.setdefault(layer_name, []).append(
                sparse.csr_matrix(cache.layers[layer_name], dtype=np.int64)
            )
        metadata_parts.append(cache.uns["fragment_shape"])
        source_hashes[source.gsm] = str(
            read_yaml(cache_path.parent / "manifest.yaml")["source_fragment_sha256"]
        )
    obs = pd.concat(obs_parts, axis=0)
    if obs.index.duplicated().any() or set(obs["cell_type"].astype(str)) != set(CELL_TYPES):
        raise ValueError("Assembled GSE246791 reference observation axis is invalid")
    layers = {
        name: sparse.vstack(parts, format="csr", dtype=np.int64)
        for name, parts in layer_parts.items()
    }
    var = selected.set_index("peak_id")[["chrom", "start", "end"]].copy()
    x = sparse.csr_matrix((len(obs), len(var)), dtype=np.int64)
    for matrix in layers.values():
        x = (x + matrix).tocsr()
    reference = ad.AnnData(X=x, obs=obs, var=var)
    for name, matrix in layers.items():
        reference.layers[name] = matrix
    label_sha256 = file_digest(LABEL_ROOT / "cells.tsv.gz")
    reference.uns["fragment_shape"] = combined_metadata(
        metadata_parts,
        layers,
        reference.var_names.astype(str),
        label_sha256,
        source_hashes,
    )
    validate_fragment_shape_spec(
        FragmentShapeSpec.from_mapping(reference.uns["fragment_shape"])
    )
    validate_fragment_shape_feature_axis(reference, "GSE246791 reference")
    REFERENCE_ROOT.parent.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(
        tempfile.mkdtemp(
            prefix=f".{CONFIG['standardized_reference_id']}.",
            dir=REFERENCE_ROOT.parent,
        )
    )
    try:
        (temporary_root / "atac").mkdir()
        temporary_path = temporary_root / "atac/reference.h5ad"
        reference.write_h5ad(temporary_path, compression="gzip")
        atomic_yaml(
            temporary_root / "reference.yaml",
            {
                "schema_version": 1,
                "reference_id": CONFIG["standardized_reference_id"],
                "source_dataset_id": CONFIG["source_dataset_id"],
                "description": (
                    "Twelve-region adult mouse-brain snATAC reference on a "
                    "reference-only 5,000-interval mm10 axis."
                ),
                "labels_key": "cell_type",
                "genome_build": "mm10",
                "cell_types": list(CELL_TYPES),
                "counts": {"cells": reference.n_obs, "peaks": reference.n_vars},
                "modalities": {
                    "atac": {
                        "path": repository_path(
                            REFERENCE_ROOT / "atac/reference.h5ad"
                        ),
                        "feature_type": "500bp_intervals",
                    }
                },
                "provenance": {
                    "source_lock": repository_path(
                        ROOT / "configs/data_sources/shapemix_gse246791_reference_lock.yaml"
                    ),
                    "label_manifest": repository_path(LABEL_ROOT / "manifest.yaml"),
                    "feature_manifest": repository_path(AXIS_ROOT / "manifest.yaml"),
                    "coordinate_audit": repository_path(AUDIT_PATH),
                    "fragment_cache_manifests": [
                        repository_path(path.parent / "manifest.yaml") for path in caches
                    ],
                },
            },
        )
        temporary_root.rename(REFERENCE_ROOT)
    except BaseException:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    print(
        f"gse246791 reference cells={reference.n_obs} peaks={reference.n_vars} "
        "status=completed",
        flush=True,
    )
    return reference_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=(
            "feature-axis",
            "coordinate-audit",
            "fragment-cache",
            "reference",
            "all",
        ),
        default="all",
    )
    parser.add_argument("--gsm", action="append")
    return parser.parse_args()


def main() -> None:
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run this materializer through run_shapemix_low_impact.sh")
    args = parse_args()
    if args.stage in {"feature-axis", "all"}:
        build_feature_axis()
    if args.stage in {"coordinate-audit", "all"}:
        build_coordinate_audit()
    if args.stage == "fragment-cache":
        build_fragment_caches(args.gsm)
    if args.stage in {"reference", "all"}:
        build_reference()


if __name__ == "__main__":
    main()
