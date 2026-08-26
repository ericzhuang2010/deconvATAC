#!/usr/bin/env python
"""Preprocess GSE194122 annotations and author ARC fragments for ShapeMix."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import os
import re
import shutil
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pysam
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/data_sources/shapemix_gse194122.yaml"
CHUNK_BYTES = 8 * 1024 * 1024
SAMPLE_PATTERN = re.compile(r"(?<![a-z0-9])s(\d+)d(\d+)(?!\d)", re.IGNORECASE)
BARCODE_PATTERN = re.compile(
    r"(?<![ACGTN])([ACGTN]{16}(?:-\d+)?)(?![ACGTN])", re.IGNORECASE
)
PEAK_PATTERN = re.compile(r"^(.+?)(?::|-)(\d+)-(\d+)$")


def load_config(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict) or value.get("accession") != "GSE194122":
        raise ValueError(f"Not a GSE194122 source manifest: {path}")
    return value


def project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def processed_root(config: Mapping[str, Any]) -> Path:
    return project_path(str(config["processed_directory"]))


def raw_h5ad_gzip(config: Mapping[str, Any]) -> Path:
    return (
        project_path(str(config["raw_directory"]))
        / "processed_downloads"
        / str(config["processed_object"]["filename"])
    )


def source_h5ad(config: Mapping[str, Any]) -> Path:
    filename = str(config["processed_object"]["filename"])
    if not filename.endswith(".gz"):
        raise ValueError("processed_object filename must end in .gz")
    return processed_root(config) / "source_audit" / "source_objects" / filename.removesuffix(".gz")


def hash_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sequence_sha256(values: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        encoded = str(value).encode()
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


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


def atomic_tsv_gzip(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.unlink(missing_ok=True)
    try:
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                    frame.to_csv(text, sep="\t", index=False, lineterminator="\n")
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def ensure_source_h5ad(config: Mapping[str, Any]) -> Path:
    source = raw_h5ad_gzip(config)
    destination = source_h5ad(config)
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        temporary = destination.with_name(f".{destination.name}.tmp")
        temporary.unlink(missing_ok=True)
        try:
            with gzip.open(source, "rb") as input_handle, temporary.open("wb") as output_handle:
                shutil.copyfileobj(input_handle, output_handle, length=CHUNK_BYTES)
            os.replace(temporary, destination)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    with destination.open("rb") as handle:
        if handle.read(8) != b"\x89HDF\r\n\x1a\n":
            raise ValueError(f"Decompressed object is not HDF5: {destination}")
    return destination


def parse_sample_key(value: object) -> Optional[tuple[str, int, int]]:
    match = SAMPLE_PATTERN.search(str(value))
    if match is None:
        return None
    site, donor = int(match.group(1)), int(match.group(2))
    return f"s{site}d{donor}", site, donor


def parse_barcode(value: object) -> Optional[str]:
    matches = BARCODE_PATTERN.findall(str(value))
    unique = {match.upper() for match in matches}
    if not unique:
        return None
    if len(unique) != 1:
        raise ValueError(f"Ambiguous 10x barcode in cell identifier: {value!r}")
    return unique.pop()


def barcode_sequence(value: object) -> Optional[str]:
    barcode = parse_barcode(value)
    return barcode.split("-", 1)[0] if barcode is not None else None


def choose_column(frame: pd.DataFrame, candidates: tuple[str, ...], role: str) -> str:
    normalized = {str(column).lower(): str(column) for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in normalized:
            return normalized[candidate.lower()]
    raise ValueError(f"Could not identify {role} column; observed={list(frame.columns)}")


def sample_records(config: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    records = {
        f"s{int(sample['site'])}d{int(sample['donor'])}": sample
        for sample in config["atac_samples"]
    }
    if len(records) != 13:
        raise ValueError("ATAC manifest does not contain 13 unique site/donor samples")
    return records


def cell_table(adata: ad.AnnData, config: Mapping[str, Any]) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs.insert(0, "cell_id", adata.obs_names.astype(str))
    label_column = choose_column(
        obs,
        ("cell_type", "celltype", "cell_type_l2", "author_cell_type"),
        "author cell-type",
    )
    search_columns = [
        column
        for column in obs.columns
        if str(column).lower() in {"batch", "sample", "sample_id", "dataset", "site"}
    ]
    manifest = sample_records(config)
    rows = []
    for record in obs.to_dict(orient="records"):
        cell_id = str(record["cell_id"])
        parsed = None
        for value in [cell_id, *(record[column] for column in search_columns)]:
            parsed = parse_sample_key(value)
            if parsed is not None:
                break
        if parsed is None:
            raise ValueError(f"Could not resolve site/donor sample for cell {cell_id!r}")
        sample_key, site, donor = parsed
        if sample_key not in manifest:
            raise ValueError(f"H5AD cell resolves outside ATAC manifest: {sample_key}")
        barcode = parse_barcode(cell_id)
        if barcode is None:
            for column in ("barcode", "cell_barcode"):
                if column in record:
                    barcode = parse_barcode(record[column])
                    if barcode is not None:
                        break
        if barcode is None:
            raise ValueError(f"Could not resolve corrected 10x barcode for {cell_id!r}")
        sample = manifest[sample_key]
        label = str(record[label_column])
        rows.append(
            {
                "cell_index": len(rows),
                "cell_id": cell_id,
                "processed_barcode": barcode,
                "barcode_sequence": barcode.split("-", 1)[0],
                "sample_key": sample_key,
                "gsm": sample["gsm"],
                "run": sample["run"],
                "site": site,
                "donor": donor,
                "donor_id": str(sample["donor_id"]),
                "author_cell_type": label,
                "harmonized_cell_type": label.strip(),
            }
        )
    result = pd.DataFrame.from_records(rows)
    if result["cell_id"].duplicated().any():
        raise ValueError("Processed H5AD cell identifiers are not unique")
    if result.duplicated(["sample_key", "barcode_sequence"]).any():
        raise ValueError("Corrected barcode mapping is not injective within sample")
    if set(result["sample_key"]) != set(manifest):
        raise ValueError("Processed H5AD does not contain all 13 declared ATAC samples")
    return result


def feature_table(adata: ad.AnnData) -> pd.DataFrame:
    var = adata.var.copy()
    var.insert(0, "feature_id", adata.var_names.astype(str))
    type_column = choose_column(
        var,
        ("feature_types", "feature_type", "modality"),
        "feature-type",
    )
    feature_types = var[type_column].astype(str)
    result = pd.DataFrame(
        {
            "feature_index": np.arange(adata.n_vars),
            "feature_id": var["feature_id"].astype(str),
            "feature_type": feature_types,
        }
    )
    chromosomes, starts, ends = [], [], []
    for feature_id, feature_type in zip(result["feature_id"], result["feature_type"]):
        match = PEAK_PATTERN.fullmatch(feature_id)
        is_atac = "atac" in feature_type.lower() or "peak" in feature_type.lower()
        if is_atac and match is None:
            raise ValueError(f"ATAC feature is not chr:start-end or chr-start-end: {feature_id!r}")
        chromosomes.append(match.group(1) if match else "")
        starts.append(int(match.group(2)) if match else pd.NA)
        ends.append(int(match.group(3)) if match else pd.NA)
    result["chromosome"] = chromosomes
    result["start"] = pd.array(starts, dtype="Int64")
    result["end"] = pd.array(ends, dtype="Int64")
    if result["feature_id"].duplicated().any():
        raise ValueError("Processed H5AD feature identifiers are not unique")
    return result


def write_lodo_folds(cells: pd.DataFrame, root: Path) -> None:
    for donor in sorted(cells["donor"].unique()):
        fold = cells[["cell_id", "sample_key", "site", "donor", "harmonized_cell_type"]].copy()
        fold.insert(1, "role", np.where(fold["donor"] == donor, "heldout", "training"))
        atomic_tsv_gzip(fold, root / f"donor_{donor}" / "cells.tsv.gz")


def fragment_manifest(config: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    records = {
        str(record["sample_key"]): record
        for record in config["openproblems_fragments"]["files"]
    }
    if len(records) != 13:
        raise ValueError("Fragment manifest must contain 13 unique sample keys")
    return records


def source_sample_file(config: Mapping[str, Any], gsm: str, name: str) -> Path:
    return (
        project_path(str(config["raw_directory"]))
        / "samples"
        / gsm
        / "source_files"
        / name
    )


def retain_hard_link(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        os.link(source, destination)
    if source.stat().st_size != destination.stat().st_size or not os.path.samefile(
        source, destination
    ):
        raise ValueError(f"Processed retention is not an exact hard link: {destination}")


def validate_fragment_prefix(tabix: pysam.TabixFile, sample_key: str) -> int:
    observed = 0
    for contig in tabix.contigs:
        for line in tabix.fetch(contig):
            fields = line.split("\t")
            if len(fields) != 5:
                raise ValueError(f"{sample_key}: fragment row is not five-column TSV")
            start, end, support = int(fields[1]), int(fields[2]), int(fields[4])
            if start < 0 or end <= start or support < 1 or parse_barcode(fields[3]) is None:
                raise ValueError(f"{sample_key}: invalid fragment row: {line!r}")
            observed += 1
            if observed >= 10_000:
                return observed
    if observed == 0:
        raise ValueError(f"{sample_key}: fragment file contains no rows")
    return observed


def metrics_barcode_bridge(
    cells: pd.DataFrame, metrics_path: Path, sample_key: str
) -> tuple[pd.DataFrame, dict[str, Any]]:
    usecols = [
        "barcode",
        "atac_barcode",
        "is_cell",
        "atac_fragments",
        "atac_peak_region_cutsites",
    ]
    metrics = pd.read_csv(metrics_path, usecols=usecols)
    metrics["barcode_sequence"] = (
        metrics["barcode"].astype(str).str.split("-", n=1).str[0]
    )
    called = metrics.loc[metrics["is_cell"] == 1].copy()
    if called["barcode_sequence"].duplicated().any():
        raise ValueError(f"{sample_key}: called metrics barcodes are not injective")
    selected = cells.loc[cells["sample_key"] == sample_key].copy()
    selected = selected.drop(
        columns=[
            "fragment_barcode",
            "atac_library_barcode",
            "metric_atac_fragments",
            "metric_atac_peak_region_cutsites",
        ],
        errors="ignore",
    )
    selected = selected.merge(
        called[
            [
                "barcode_sequence",
                "barcode",
                "atac_barcode",
                "atac_fragments",
                "atac_peak_region_cutsites",
            ]
        ],
        on="barcode_sequence",
        how="left",
        validate="one_to_one",
    )
    if selected["barcode"].isna().any():
        missing = selected.loc[selected["barcode"].isna(), "cell_id"].head().tolist()
        raise ValueError(f"{sample_key}: processed cells missing from called metrics: {missing}")
    if selected["barcode"].duplicated().any():
        raise ValueError(f"{sample_key}: common fragment barcode bridge is not injective")
    selected = selected.rename(
        columns={
            "barcode": "fragment_barcode",
            "atac_barcode": "atac_library_barcode",
            "atac_fragments": "metric_atac_fragments",
            "atac_peak_region_cutsites": "metric_atac_peak_region_cutsites",
        }
    )
    summary = {
        "metrics_rows": len(metrics),
        "cellranger_called_cells": len(called),
        "processed_annotated_cells": len(selected),
        "processed_cells_are_called_subset": True,
        "exact_one_to_one_common_fragment_barcode_bridge": True,
        "atac_library_barcode_retained_for_provenance": True,
    }
    return selected, summary


def normalize_fragment_suite(
    config: Mapping[str, Any], sample_keys: Optional[set[str]] = None
) -> None:
    root = processed_root(config)
    labels_path = root / "labels/source_broad7_v1/cells.tsv.gz"
    if not labels_path.is_file():
        raise FileNotFoundError(f"Run metadata preprocessing first: {labels_path}")
    cells = pd.read_csv(labels_path, sep="\t", low_memory=False)
    bridged = []
    records = []
    manifest = fragment_manifest(config)
    if sample_keys is not None:
        unknown = sample_keys - set(manifest)
        if unknown:
            raise ValueError(f"Unknown fragment sample keys: {sorted(unknown)}")
    selected_manifest = {
        key: value
        for key, value in manifest.items()
        if sample_keys is None or key in sample_keys
    }
    for sample_key, record in selected_manifest.items():
        gsm = str(record["gsm"])
        source_fragment = source_sample_file(config, gsm, "atac_fragments.tsv.gz")
        source_index = source_sample_file(config, gsm, "atac_fragments.tsv.gz.tbi")
        source_metrics = source_sample_file(config, gsm, "per_barcode_metrics.csv")
        destination = root / "normalized_fragments/GRCh38" / sample_key
        fragment = destination / "fragments.tsv.gz"
        index = Path(f"{fragment}.tbi")
        metrics = destination / "per_barcode_metrics.csv"
        retain_hard_link(source_fragment, fragment)
        retain_hard_link(source_index, index)
        retain_hard_link(source_metrics, metrics)
        with pysam.TabixFile(str(fragment), index=str(index)) as tabix:
            contigs = list(tabix.contigs)
            prefix_rows = validate_fragment_prefix(tabix, sample_key)
        selected, bridge = metrics_barcode_bridge(cells, metrics, sample_key)
        bridged.append(selected)
        records.append(
            {
                "sample_key": sample_key,
                "gsm": gsm,
                "fragment": str(fragment.relative_to(ROOT)),
                "fragment_bytes": fragment.stat().st_size,
                "index": str(index.relative_to(ROOT)),
                "index_bytes": index.stat().st_size,
                "metrics": str(metrics.relative_to(ROOT)),
                "metrics_bytes": metrics.stat().st_size,
                "tabix_contigs": len(contigs),
                "fragment_schema_prefix_rows": prefix_rows,
                **bridge,
            }
        )
        print(f"normalized {sample_key}: {bridge['processed_annotated_cells']} cells", flush=True)
    incoming = pd.concat(bridged, ignore_index=True).set_index("cell_id")
    reconciled = cells.set_index("cell_id")
    bridge_columns = [
        "fragment_barcode",
        "atac_library_barcode",
        "metric_atac_fragments",
        "metric_atac_peak_region_cutsites",
    ]
    for column in bridge_columns:
        if column not in reconciled:
            reconciled[column] = pd.NA
        reconciled.loc[incoming.index, column] = incoming[column]
    reconciled = reconciled.reset_index().sort_values("cell_index")
    if reconciled["cell_id"].tolist() != cells.sort_values("cell_index")["cell_id"].tolist():
        raise ValueError("Barcode reconciliation changed the H5AD cell axis")
    complete = set(selected_manifest) == set(manifest)
    if complete and reconciled["fragment_barcode"].isna().any():
        raise ValueError("Complete fragment normalization left cells without barcodes")
    atomic_tsv_gzip(reconciled, labels_path)
    write_lodo_folds(reconciled, root / "splits/broad7_lodo_v1")
    summary = {
        "schema_version": 1,
        "accession": config["accession"],
        "samples": len(records),
        "sample_keys": list(selected_manifest),
        "complete": complete,
        "source_bytes": sum(
            int(item["fragment_bytes"] + item["index_bytes"] + item["metrics_bytes"])
            for item in selected_manifest.values()
        ),
        "retention": "hard links preserve exact provider bytes without duplicate blocks",
        "processed_cells": len(reconciled),
        "processed_cells_with_exact_metrics_bridge": int(
            reconciled["fragment_barcode"].notna().sum()
        ),
        "all_processed_cells_have_exact_metrics_bridge": bool(
            reconciled["fragment_barcode"].notna().all()
        ),
        "files": records,
    }
    name = "fragment_suite.yaml" if complete else "pilot_fragment_suite.yaml"
    atomic_yaml(root / "source_audit" / name, summary)


def csr_selected_values(
    h5ad_path: Path, rows: list[int], columns: list[int]
) -> np.ndarray:
    output = np.zeros((len(rows), len(columns)), dtype=np.int64)
    column_lookup = {column: index for index, column in enumerate(columns)}
    with h5py.File(h5ad_path, "r") as handle:
        group = handle["layers/counts"]
        indptr = group["indptr"]
        indices = group["indices"]
        data = group["data"]
        for output_row, source_row in enumerate(rows):
            start, end = int(indptr[source_row]), int(indptr[source_row + 1])
            row_indices = indices[start:end]
            row_data = data[start:end]
            for source_column, value in zip(row_indices, row_data):
                destination = column_lookup.get(int(source_column))
                if destination is not None:
                    if value < 0 or not float(value).is_integer():
                        raise ValueError("Raw count layer contains non-integer values")
                    output[output_row, destination] = int(value)
    return output


def count_cut_sites(
    fragment_path: Path,
    index_path: Path,
    cells: pd.DataFrame,
    peaks: pd.DataFrame,
    right_offset: int,
) -> np.ndarray:
    result = np.zeros((len(cells), len(peaks)), dtype=np.int64)
    cell_lookup = {
        str(barcode): index for index, barcode in enumerate(cells["fragment_barcode"])
    }
    with pysam.TabixFile(str(fragment_path), index=str(index_path)) as tabix:
        for peak_index, peak in enumerate(peaks.itertuples(index=False)):
            query_start = max(0, int(peak.start) - 1)
            query_end = int(peak.end) + 1
            for line in tabix.fetch(str(peak.chromosome), query_start, query_end):
                fields = line.split("\t")
                cell_index = cell_lookup.get(fields[3])
                if cell_index is None:
                    continue
                left, right = int(fields[1]), int(fields[2]) + right_offset
                if int(peak.start) <= left < int(peak.end):
                    result[cell_index, peak_index] += 1
                if int(peak.start) <= right < int(peak.end):
                    result[cell_index, peak_index] += 1
    return result


def audit_pilot_matrix(config: Mapping[str, Any]) -> None:
    root = processed_root(config)
    cells = pd.read_csv(root / "labels/source_broad7_v1/cells.tsv.gz", sep="\t", low_memory=False)
    features = pd.read_csv(root / "feature_axes/source_axis_v1/features.tsv.gz", sep="\t")
    pilot = next(
        record
        for record in config["openproblems_fragments"]["files"]
        if record["gsm"] == config["pilot_gsm"]
    )
    sample_key = str(pilot["sample_key"])
    selected_cells = cells.loc[cells["sample_key"] == sample_key].head(64).copy()
    candidates = features.loc[
        (features["feature_type"] == "ATAC") & (features["chromosome"] == "chr1")
    ].head(4096)
    candidate_counts = csr_selected_values(
        source_h5ad(config),
        selected_cells["cell_index"].astype(int).tolist(),
        candidates["feature_index"].astype(int).tolist(),
    )
    nonzero = np.flatnonzero(candidate_counts.sum(axis=0) > 0)[:256]
    if len(nonzero) < 128:
        raise ValueError("Pilot matrix audit could not find 128 nonzero chr1 peaks")
    peaks = candidates.iloc[nonzero].reset_index(drop=True)
    expected = candidate_counts[:, nonzero]
    directory = root / "normalized_fragments/GRCh38" / sample_key
    fragment = directory / "fragments.tsv.gz"
    index = Path(f"{fragment}.tbi")
    candidates_audit = []
    for right_offset in (0, -1):
        observed = count_cut_sites(
            fragment, index, selected_cells, peaks, right_offset=right_offset
        )
        difference = observed - expected
        candidates_audit.append(
            {
                "right_cut_offset": right_offset,
                "depth_normalized_total_cuts": int(expected.sum()),
                "full_fragment_total_cuts": int(observed.sum()),
                "depth_normalized_to_full_ratio": float(expected.sum() / observed.sum()),
                "mismatched_entries": int(np.count_nonzero(difference)),
                "entries_below_depth_normalized_count": int(np.count_nonzero(difference < 0)),
                "entries_above_depth_normalized_count": int(np.count_nonzero(difference > 0)),
                "absolute_error": int(np.abs(difference).sum()),
                "maximum_absolute_error": int(np.abs(difference).max(initial=0)),
                "entrywise_contains_depth_normalized_counts": bool(np.all(difference >= 0)),
                "exact": bool(np.array_equal(observed, expected)),
            }
        )
    compatible = [
        item
        for item in candidates_audit
        if item["entrywise_contains_depth_normalized_counts"]
        and item["full_fragment_total_cuts"] > item["depth_normalized_total_cuts"]
    ]
    minimum_error = min(item["absolute_error"] for item in compatible) if compatible else None
    selected = [item for item in compatible if item["absolute_error"] == minimum_error]
    if len(selected) != 1 or selected[0]["right_cut_offset"] != 0:
        raise ValueError(f"Pilot matrix semantics are not uniquely supported: {candidates_audit}")
    audit = {
        "schema_version": 1,
        "accession": config["accession"],
        "pilot_gsm": pilot["gsm"],
        "sample_key": sample_key,
        "cells": len(selected_cells),
        "peaks": len(peaks),
        "matrix_entries": int(expected.size),
        "depth_normalized_total_cuts": int(expected.sum()),
        "full_fragment_total_cuts": selected[0]["full_fragment_total_cuts"],
        "depth_normalized_to_full_ratio": selected[0]["depth_normalized_to_full_ratio"],
        "matrix_relationship": "H5AD counts are a Cell Ranger ARC aggr depth-normalized subsample of the full per-sample fragments",
        "fragment_weight": "one per deduplicated row; readSupport ignored",
        "left_cut_offset": 0,
        "candidates": candidates_audit,
        "selected_right_cut_offset": selected[0]["right_cut_offset"],
        "selection_basis": "entrywise containment plus the unique minimum-error endpoint-boundary convention",
        "entrywise_subsample_reconstruction": True,
        "exact_reconstruction": False,
        "cell_ids": selected_cells["cell_id"].tolist(),
        "feature_ids": peaks["feature_id"].tolist(),
    }
    atomic_yaml(root / "source_audit/pilot_matrix_reconstruction.yaml", audit)
    print(
        f"validated depth-normalized containment for {audit['matrix_entries']} pilot matrix entries",
        flush=True,
    )


def preprocess_metadata(config: Mapping[str, Any]) -> None:
    path = ensure_source_h5ad(config)
    adata = ad.read_h5ad(path, backed="r")
    try:
        cells = cell_table(adata, config)
        features = feature_table(adata)
        details = {
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "obs_columns": [str(column) for column in adata.obs.columns],
            "var_columns": [str(column) for column in adata.var.columns],
            "layers": sorted(str(key) for key in adata.layers.keys()),
            "obsm": sorted(str(key) for key in adata.obsm.keys()),
            "varm": sorted(str(key) for key in adata.varm.keys()),
            "uns": sorted(str(key) for key in adata.uns.keys()),
        }
    finally:
        adata.file.close()
    root = processed_root(config)
    labels_path = root / "labels/source_broad7_v1/cells.tsv.gz"
    features_path = root / "feature_axes/source_axis_v1/features.tsv.gz"
    atomic_tsv_gzip(cells, labels_path)
    atomic_tsv_gzip(features, features_path)
    write_lodo_folds(cells, root / "splits/broad7_lodo_v1")
    modality_counts = features["feature_type"].value_counts().sort_index()
    summary = {
        "schema_version": 1,
        "accession": config["accession"],
        "source_h5ad": str(path.relative_to(ROOT)),
        "source_h5ad_bytes": path.stat().st_size,
        "source_h5ad_sha256": hash_file(path),
        "source_h5ad_gzip_sha256": hash_file(raw_h5ad_gzip(config)),
        **details,
        "cells": len(cells),
        "features": len(features),
        "sites": sorted(int(value) for value in cells["site"].unique()),
        "donors": sorted(int(value) for value in cells["donor"].unique()),
        "samples": cells["sample_key"].nunique(),
        "cell_types": cells["harmonized_cell_type"].nunique(),
        "cells_by_sample": {str(k): int(v) for k, v in cells["sample_key"].value_counts().sort_index().items()},
        "cells_by_donor": {str(k): int(v) for k, v in cells["donor"].value_counts().sort_index().items()},
        "cells_by_type": {str(k): int(v) for k, v in cells["harmonized_cell_type"].value_counts().sort_index().items()},
        "features_by_type": {str(k): int(v) for k, v in modality_counts.items()},
        "cell_axis_sha256": sequence_sha256(cells["cell_id"]),
        "feature_axis_sha256": sequence_sha256(features["feature_id"]),
        "outputs": {
            "labels": str(labels_path.relative_to(ROOT)),
            "features": str(features_path.relative_to(ROOT)),
            "lodo_root": str((root / "splits/broad7_lodo_v1").relative_to(ROOT)),
        },
    }
    atomic_yaml(root / "source_audit/processed_h5ad.yaml", summary)
    print(f"processed {len(cells)} cells and {len(features)} features", flush=True)


def fragment_from_alignment(record: pysam.AlignedSegment) -> Optional[tuple[str, int, int, str, int]]:
    if (
        not record.is_paired
        or not record.is_proper_pair
        or record.is_unmapped
        or record.mate_is_unmapped
        or record.is_secondary
        or record.is_supplementary
        or record.is_qcfail
        or record.is_duplicate
        or record.mapping_quality < 30
        or record.next_reference_id != record.reference_id
        or record.template_length <= 0
        or not record.has_tag("CB")
    ):
        return None
    start = int(record.reference_start) + 4
    end = int(record.reference_start) + int(record.template_length) - 5
    if start < 0 or end <= start:
        return None
    return str(record.reference_name), start, end, str(record.get_tag("CB")), 1


def pilot_bam_path(config: Mapping[str, Any]) -> Path:
    sample = next(
        sample for sample in config["atac_samples"] if sample["gsm"] == config["pilot_gsm"]
    )
    return (
        project_path(str(config["raw_directory"]))
        / "samples"
        / str(sample["gsm"])
        / "source_files"
        / str(sample["bam_filename"])
    )


def reconstruct_pilot(config: Mapping[str, Any], overwrite: bool = False) -> None:
    root = processed_root(config)
    labels_path = root / "labels/source_broad7_v1/cells.tsv.gz"
    if not labels_path.is_file():
        raise FileNotFoundError(f"Run metadata preprocessing first: {labels_path}")
    cells = pd.read_csv(labels_path, sep="\t", dtype=str)
    pilot = next(
        sample for sample in config["atac_samples"] if sample["gsm"] == config["pilot_gsm"]
    )
    selected = cells.loc[cells["gsm"] == pilot["gsm"]]
    if "fragment_barcode" not in selected:
        raise ValueError("Run fragment normalization to create the metrics barcode bridge")
    expected = set(selected["fragment_barcode"])
    if not expected:
        raise ValueError("Processed H5AD contains no cells for the pilot sample")
    sample_key = f"s{pilot['site']}d{pilot['donor']}"
    output = root / "reconstructed_fragments" / sample_key / "fragments.tsv.gz"
    index = Path(f"{output}.tbi")
    if output.exists() and index.exists() and not overwrite:
        print(f"reusing {output}", flush=True)
        return
    bam = pilot_bam_path(config)
    if not bam.is_file():
        raise FileNotFoundError(f"Acquire the pilot BAM first: {bam}")
    if bam.stat().st_size != int(pilot["bam_bytes"]):
        raise ValueError("Pilot BAM byte size does not match the frozen manifest")
    if hash_file(bam, "md5") != str(pilot["bam_md5"]):
        raise ValueError("Pilot BAM MD5 does not match NCBI")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.unlink(missing_ok=True)
    counts = Counter()
    seen = set()
    try:
        with pysam.AlignmentFile(bam, "rb") as alignments, pysam.BGZFile(
            str(temporary), "wb"
        ) as fragments:
            if alignments.header.to_dict().get("HD", {}).get("SO") != "coordinate":
                raise ValueError("Pilot BAM is not coordinate sorted")
            for alignment in alignments.fetch(until_eof=True):
                counts["alignment_records"] += 1
                fragment = fragment_from_alignment(alignment)
                if fragment is None:
                    continue
                counts["passing_fragments_all_barcodes"] += 1
                chromosome, start, end, barcode, support = fragment
                if barcode not in expected:
                    continue
                seen.add(barcode)
                counts["retained_fragments"] += 1
                fragments.write(
                    f"{chromosome}\t{start}\t{end}\t{barcode}\t{support}\n".encode()
                )
        missing = expected - seen
        if missing:
            raise ValueError(
                f"Pilot BAM is missing {len(missing)} processed-cell barcodes; "
                f"examples={sorted(missing)[:5]}"
            )
        os.replace(temporary, output)
        index.unlink(missing_ok=True)
        pysam.tabix_index(str(output), preset="bed", force=True)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    with pysam.TabixFile(str(output)) as tabix:
        contigs = list(tabix.contigs)
    audit = {
        "schema_version": 1,
        "accession": config["accession"],
        "pilot_gsm": pilot["gsm"],
        "pilot_run": pilot["run"],
        "sample_key": sample_key,
        "bam": str(bam.relative_to(ROOT)),
        "bam_bytes": bam.stat().st_size,
        "bam_md5": pilot["bam_md5"],
        "processed_cells_expected": len(expected),
        "processed_cells_seen": len(seen),
        "exact_barcode_reconciliation": seen == expected,
        "filters": {
            "proper_primary_pair": True,
            "minimum_mapq": 30,
            "exclude_qcfail_secondary_supplementary_duplicate": True,
            "corrected_barcode_tag": "CB",
            "one_record_per_pair": "positive template_length",
            "tn5_left_shift": 4,
            "tn5_right_shift": -5,
            "read_support": 1,
        },
        "counts": dict(counts),
        "fragments": str(output.relative_to(ROOT)),
        "fragments_bytes": output.stat().st_size,
        "fragments_sha256": hash_file(output),
        "tabix_index": str(index.relative_to(ROOT)),
        "tabix_sha256": hash_file(index),
        "tabix_contigs": contigs,
    }
    atomic_yaml(root / "source_audit/pilot_fragments.yaml", audit)
    print(
        f"retained {counts['retained_fragments']} fragments for {len(seen)} pilot cells",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=("metadata", "fragments", "matrix-audit", "bam-fallback", "all"),
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--sample-key",
        action="append",
        default=[],
        help="limit the fragments stage to declared sample keys",
    )
    args = parser.parse_args()
    config = load_config(args.config)
    if args.stage in {"metadata", "all"}:
        preprocess_metadata(config)
    if args.stage in {"fragments", "all"}:
        normalize_fragment_suite(config, sample_keys=set(args.sample_key) or None)
    if args.stage in {"matrix-audit", "all"}:
        audit_pilot_matrix(config)
    if args.stage == "bam-fallback":
        reconstruct_pilot(config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
