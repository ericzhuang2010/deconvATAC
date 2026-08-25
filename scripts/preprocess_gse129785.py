#!/usr/bin/env python
"""Preprocess the frozen GSE129785 immune scope into ShapeMix-ready objects."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import pysam
import scipy
import yaml
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data import ordered_feature_sha256  # noqa: E402
from deconvatac.data.schemas import FragmentShapeSpec  # noqa: E402
from deconvatac.data.validators import _validate_shape_object  # noqa: E402
from deconvatac.pp import (  # noqa: E402
    FRAGMENT_SHAPE_LAYER_NAMES,
    PeakInterval,
    build_fragment_shape_anndata,
    count_fragment_shapes_from_records,
)


DEFAULT_CONFIG = ROOT / "configs/data_sources/shapemix_gse129785.yaml"
N_SELECTED_PEAKS = 5_000
MIN_REFERENCE_CELLS = 10
CELL_CALL_MIN_FRAGMENTS = 1_000
SHAPE_CHUNK_SIZE = 1_000_000

REFERENCE_TYPES = (
    "Dendritic cells",
    "Monocytes",
    "B cells",
    "Regulatory T cells",
    "Naive CD4 T cells",
    "Memory CD4 T cells",
    "NK cells",
    "Naive CD8 T cells",
    "Memory CD8 T cells",
)


def load_yaml(path: Path) -> dict[str, Any]:
    """Load a mapping from YAML."""
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


def dump_yaml(path: Path, value: Mapping[str, Any]) -> None:
    """Write YAML after creating the parent directory."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w") as handle:
        yaml.safe_dump(dict(value), handle, sort_keys=False)
    os.replace(temporary, path)


def sha256_file(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_bytes), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_bgzf(path: Path) -> bool:
    """Return whether a gzip-compatible file has the standard BGZF block header."""
    with path.open("rb") as handle:
        prefix = handle.read(18)
    return (
        len(prefix) == 18
        and prefix[:4] == b"\x1f\x8b\x08\x04"
        and prefix[12:14] == b"BC"
    )


def recompress_to_bgzf(source: Path, destination: Path) -> None:
    """Atomically convert one immutable plain-gzip source into derived BGZF."""
    temporary = destination.with_name(f".{destination.name}.bgzf.tmp")
    temporary.unlink(missing_ok=True)
    try:
        with gzip.open(source, "rb") as input_handle, pysam.BGZFile(
            str(temporary), "wb"
        ) as output_handle:
            for chunk in iter(lambda: input_handle.read(8 * 1024 * 1024), b""):
                output_handle.write(chunk)
        if not is_bgzf(temporary):
            raise ValueError(f"BGZF conversion did not produce BGZF: {temporary}")
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def open_text_auto(path: Path):
    """Open gzip content or provider-served plain text based on file magic."""
    with path.open("rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzip else path.open("rt")


def sequence_sha256(values: Iterable[str]) -> str:
    """Hash an ordered string sequence with unambiguous length prefixes."""
    digest = hashlib.sha256()
    for value in values:
        encoded = str(value).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


def version(name: str) -> str:
    """Return an installed distribution version or an explicit fallback."""
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "working-tree"


def validate_shape_object(adata: ad.AnnData, role: str) -> None:
    """Apply the repository complete single-object ShapeMix contract."""
    spec = FragmentShapeSpec.from_mapping(adata.uns.get("fragment_shape"))
    _validate_shape_object(adata, spec, role)


def project_path(value: str) -> Path:
    """Resolve one project-relative path from the tracked manifest."""
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sample_filename(sample: Mapping[str, Any]) -> str:
    return f"{sample['gsm']}_{sample['title']}_fragments.tsv.gz"


def raw_fragment_path(config: Mapping[str, Any], sample: Mapping[str, Any]) -> Path:
    return (
        project_path(str(config["raw_directory"]))
        / "samples"
        / str(sample["gsm"])
        / "source_files"
        / sample_filename(sample)
    )


def series_path(config: Mapping[str, Any], filename: str) -> Path:
    return project_path(str(config["raw_directory"])) / "series_metadata" / filename


def processed_root(config: Mapping[str, Any]) -> Path:
    return project_path(str(config["processed_directory"]))


def normalized_fragment_path(config: Mapping[str, Any], sample: Mapping[str, Any]) -> Path:
    return (
        processed_root(config)
        / "normalized_fragments"
        / "hg19"
        / str(sample["gsm"])
        / sample_filename(sample)
    )


def ensure_normalized_fragment(
    config: Mapping[str, Any], sample: Mapping[str, Any]
) -> tuple[Path, Path]:
    """Retain BGZF sources by hard link or recompress plain gzip, then index."""
    source = raw_fragment_path(config, sample)
    destination = normalized_fragment_path(config, sample)
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    source_is_bgzf = is_bgzf(source)
    if source_is_bgzf:
        if not destination.exists():
            os.link(source, destination)
        if source.stat().st_size != destination.stat().st_size:
            raise ValueError(f"Normalized fragment link has the wrong size: {destination}")
        if not os.path.samefile(source, destination):
            raise ValueError(f"BGZF source must be retained as a hard link: {destination}")
    else:
        if destination.exists() and os.path.samefile(source, destination):
            destination.unlink()
        if not destination.exists():
            recompress_to_bgzf(source, destination)
        if not is_bgzf(destination):
            raise ValueError(f"Normalized fragment is not BGZF: {destination}")

    index = Path(f"{destination}.tbi")
    if not index.exists():
        pysam.tabix_index(str(destination), preset="bed", force=False)
    with pysam.TabixFile(str(destination)) as tabix:
        if not tabix.contigs:
            raise ValueError(f"Tabix index has no contigs: {index}")
    return destination, index


def load_author_metadata(config: Mapping[str, Any]) -> dict[str, pd.DataFrame]:
    """Load author-filtered cell tables without changing row order."""
    files = {
        "hematopoiesis": "GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz",
        "fresh": "GSE129785_scATAC-PBMCs-Fresh.cell_barcodes.txt.gz",
        "frozen_unsorted": "GSE129785_scATAC-PBMCs-Frozen.cell_barcodes.txt.gz",
        "frozen_sorted": "GSE129785_scATAC-PBMCs-FrozenSort.cell_barcodes.txt.gz",
    }
    return {name: pd.read_csv(series_path(config, filename), sep="\t") for name, filename in files.items()}


def official_barcodes(sample: Mapping[str, Any], metadata: Mapping[str, pd.DataFrame]) -> Optional[tuple[str, ...]]:
    """Return author cell calls when raw barcodes were published for this sample."""
    title = str(sample["title"])
    role = str(sample["role"])
    if role in {"sorted_reference", "pbmc_evaluation"}:
        table = metadata["hematopoiesis"]
    elif title == "Fresh_pbmc_5k":
        table = metadata["fresh"]
    elif title == "Frozen_unsorted_pbmc_5k":
        table = metadata["frozen_unsorted"]
    elif title == "Frozen_sorted_pbmc_5k":
        table = metadata["frozen_sorted"]
    else:
        return None
    rows = table.loc[table["Group"].astype(str) == title]
    if rows.empty:
        raise ValueError(f"Author metadata has no rows for {title}.")
    values = rows["Barcodes"].dropna().astype(str)
    values = values[values.str.upper() != "NA"]
    if values.empty:
        return None
    if values.duplicated().any():
        raise ValueError(f"Author barcode calls are duplicated for {title}.")
    return tuple(values)


def scan_fragments(path: Path) -> tuple[Counter[str], dict[str, int]]:
    """Audit every source row and collect per-barcode fragment counts."""
    counts: Counter[str] = Counter()
    summary = {
        "total_rows": 0,
        "invalid_rows": 0,
        "read_support_total": 0,
        "max_read_support": 0,
        "max_fragment_length": 0,
    }
    with pysam.BGZFile(str(path), "rb") as handle:
        for raw in handle:
            if raw.startswith(b"#"):
                continue
            fields = raw.rstrip(b"\r\n").split(b"\t")
            summary["total_rows"] += 1
            if len(fields) != 5:
                summary["invalid_rows"] += 1
                continue
            try:
                start, end, support = int(fields[1]), int(fields[2]), int(fields[4])
                barcode = fields[3].decode("utf-8")
            except (ValueError, UnicodeDecodeError):
                summary["invalid_rows"] += 1
                continue
            if start < 0 or end <= start or support < 1 or not barcode:
                summary["invalid_rows"] += 1
                continue
            counts[barcode] += 1
            summary["read_support_total"] += support
            summary["max_read_support"] = max(summary["max_read_support"], support)
            summary["max_fragment_length"] = max(summary["max_fragment_length"], end - start)
    summary["unique_barcodes"] = len(counts)
    return counts, summary


def audit_samples(config: Mapping[str, Any], overwrite: bool = False) -> None:
    """Index fragments, freeze cell calls, and write full source/barcode audits."""
    root = processed_root(config)
    metadata = load_author_metadata(config)
    summary_rows: list[dict[str, Any]] = []
    for sample in config["samples"]:
        gsm = str(sample["gsm"])
        label_path = root / "labels" / f"{gsm}.csv.gz"
        audit_path = root / "source_audit" / "samples" / f"{gsm}.yaml"
        barcode_count_path = root / "source_audit" / "barcode_counts" / f"{gsm}.csv.gz"
        fragment_path, index_path = ensure_normalized_fragment(config, sample)
        if audit_path.exists() and label_path.exists() and barcode_count_path.exists() and not overwrite:
            record = load_yaml(audit_path)
            if int(record["scan"]["invalid_rows"]) != 0:
                raise ValueError(f"{gsm}: reused audit contains invalid fragment rows.")
            summary_row = dict(record["summary_row"])
            if "normalized_fragment_sha256" not in summary_row:
                hard_linked = os.path.samefile(
                    raw_fragment_path(config, sample), fragment_path
                )
                normalized_digest = (
                    str(summary_row["source_sha256"])
                    if hard_linked
                    else sha256_file(fragment_path)
                )
                normalization_mode = (
                    "source_bgzf_hard_link"
                    if hard_linked
                    else "plain_gzip_to_bgzf"
                )
                summary_row.update(
                    {
                        "normalized_fragment_bytes": fragment_path.stat().st_size,
                        "normalized_fragment_sha256": normalized_digest,
                        "normalization_mode": normalization_mode,
                    }
                )
                record["summary_row"] = summary_row
                record["normalization"] = {
                    "mode": normalization_mode,
                    "source_is_bgzf": hard_linked,
                    "source_sha256": str(summary_row["source_sha256"]),
                    "normalized_sha256": normalized_digest,
                }
                dump_yaml(audit_path, record)
            summary_rows.append(summary_row)
            print(f"reused audit {gsm}", flush=True)
            continue
        counts, summary = scan_fragments(fragment_path)
        invalid_rows = int(summary["invalid_rows"])
        if invalid_rows != 0:
            raise ValueError(f"{gsm}: found {invalid_rows} invalid fragment rows.")
        published = official_barcodes(sample, metadata)
        if published is None:
            called = tuple(sorted(barcode for barcode, count in counts.items() if count >= CELL_CALL_MIN_FRAGMENTS))
            call_source = f"fragment_rows_ge_{CELL_CALL_MIN_FRAGMENTS}"
        else:
            missing = sorted(set(published).difference(counts))
            if missing:
                raise ValueError(f"{gsm}: {len(missing)} author-called barcodes have no fragments.")
            called = tuple(sorted(published))
            call_source = "author_filtered_barcode_table"
        if not called:
            raise ValueError(f"{gsm}: cell calling retained no barcodes.")

        labels = pd.DataFrame({"raw_barcode": called})
        labels["barcode"] = [f"{gsm}#{barcode}" for barcode in called]
        labels["gsm"] = gsm
        labels["sample_title"] = str(sample["title"])
        labels["sample_role"] = str(sample["role"])
        labels["cell_call_source"] = call_source
        labels["cell_type"] = str(sample.get("cell_type", "unlabeled"))
        label_path.parent.mkdir(parents=True, exist_ok=True)
        labels.to_csv(label_path, index=False)

        barcode_counts = pd.DataFrame(
            sorted(counts.items()), columns=["raw_barcode", "fragment_rows"]
        )
        barcode_counts["called_cell"] = barcode_counts["raw_barcode"].isin(called)
        barcode_count_path.parent.mkdir(parents=True, exist_ok=True)
        barcode_counts.to_csv(barcode_count_path, index=False)

        source_path = raw_fragment_path(config, sample)
        hard_linked = os.path.samefile(source_path, fragment_path)
        source_digest = sha256_file(source_path)
        normalized_digest = (
            source_digest if hard_linked else sha256_file(fragment_path)
        )
        normalization_mode = (
            "source_bgzf_hard_link" if hard_linked else "plain_gzip_to_bgzf"
        )
        summary_row = {
            "gsm": gsm,
            "title": str(sample["title"]),
            "role": str(sample["role"]),
            "source_bytes": source_path.stat().st_size,
            "source_sha256": source_digest,
            "normalized_fragment_bytes": fragment_path.stat().st_size,
            "normalized_fragment_sha256": normalized_digest,
            "normalization_mode": normalization_mode,
            "tabix_bytes": index_path.stat().st_size,
            "tabix_sha256": sha256_file(index_path),
            "total_rows": summary["total_rows"],
            "unique_barcodes": summary["unique_barcodes"],
            "called_cells": len(called),
            "cell_call_source": call_source,
            "max_fragment_length": summary["max_fragment_length"],
        }
        record = {
            "schema_version": 1,
            "source_dataset_id": config["source_dataset_id"],
            "sample": dict(sample),
            "fragment_path": str(fragment_path.relative_to(ROOT)),
            "normalization": {
                "mode": normalization_mode,
                "source_is_bgzf": hard_linked,
                "source_sha256": source_digest,
                "normalized_sha256": normalized_digest,
            },
            "tabix_path": str(index_path.relative_to(ROOT)),
            "fragment_schema": {
                "compression": "BGZF",
                "columns": ["chrom", "chromStart", "chromEnd", "barcode", "readSupport"],
                "coordinate_system": "zero-based half-open",
                "invalid_rows": summary["invalid_rows"],
            },
            "scan": summary,
            "cell_call": {
                "source": call_source,
                "called_cells": len(called),
                "labels_path": str(label_path.relative_to(ROOT)),
                "barcode_counts_path": str(barcode_count_path.relative_to(ROOT)),
            },
            "summary_row": summary_row,
        }
        dump_yaml(audit_path, record)
        summary_rows.append(summary_row)
        print(f"audited {gsm} rows={summary['total_rows']} cells={len(called)}", flush=True)
    summary_table = pd.DataFrame(summary_rows)
    summary_path = root / "source_audit" / "sample_summary.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_table.to_csv(summary_path, index=False)
    print(f"wrote {summary_path}", flush=True)


def parse_author_peaks(path: Path) -> tuple[PeakInterval, ...]:
    """Parse one-based closed author features into zero-based half-open intervals."""
    table = pd.read_csv(path, sep="\t")
    if list(table.columns) != ["Feature"]:
        raise ValueError(f"Unexpected author peak schema: {list(table.columns)}")
    peaks: list[PeakInterval] = []
    for value in table["Feature"].astype(str):
        try:
            chrom, start_text, end_text = value.rsplit("_", 2)
            source_start, end = int(start_text), int(end_text)
            if source_start < 1:
                raise ValueError(f"Author peak start must be one-based and positive: {value!r}.")
            start = source_start - 1
        except ValueError as exc:
            raise ValueError(f"Cannot parse author peak {value!r}.") from exc
        peaks.append(PeakInterval(chrom, start, end, f"{chrom}:{start}-{end}"))
    if len({peak.name for peak in peaks}) != len(peaks):
        raise ValueError("Author peak identifiers are not unique.")
    return tuple(peaks)


def stream_reference_statistics(
    matrix_path: Path,
    n_peaks: int,
    n_cells: int,
    cell_type_index: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Stream Matrix Market coordinates into ShapeMix peak-ranking sufficient statistics."""
    group_counts = np.zeros((len(REFERENCE_TYPES), n_peaks), dtype=np.float64)
    coverage = np.zeros(n_peaks, dtype=np.int64)
    dimensions: Optional[tuple[int, int, int]] = None
    entries = 0
    with open_text_auto(matrix_path) as handle:
        header = handle.readline().strip()
        if not header.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"Unsupported Matrix Market header: {header}")
        for line in handle:
            if line.startswith("%"):
                continue
            dimensions = tuple(int(value) for value in line.split())
            break
        if dimensions is None or dimensions[:2] != (n_peaks, n_cells):
            raise ValueError(
                f"Author matrix dimensions {dimensions} do not match {(n_peaks, n_cells)}."
            )
        for line in handle:
            feature_text, cell_text, count_text = line.split()
            feature = int(feature_text) - 1
            cell = int(cell_text) - 1
            type_index = int(cell_type_index[cell])
            if type_index < 0:
                continue
            count = float(count_text)
            if count <= 0 or not count.is_integer():
                raise ValueError("Author reference counts must be positive integers.")
            group_counts[type_index, feature] += count
            coverage[feature] += 1
            entries += 1
    assert dimensions is not None
    return group_counts, coverage, entries


def select_reference_peaks(config: Mapping[str, Any], overwrite: bool = False) -> None:
    """Select 5,000 peaks from sorted reference populations using protocol-v1 scoring."""
    root = processed_root(config)
    selected_path = root / "feature_axes" / "selected_peaks.txt"
    if selected_path.exists() and not overwrite:
        print(f"reused {selected_path}", flush=True)
        return
    metadata = load_author_metadata(config)["hematopoiesis"]
    peaks = parse_author_peaks(
        series_path(config, "GSE129785_scATAC-Hematopoiesis-All.peaks.txt.gz")
    )
    type_lookup = {
        str(sample["title"]): REFERENCE_TYPES.index(str(sample["cell_type"]))
        for sample in config["samples"]
        if sample["role"] == "sorted_reference"
    }
    cell_type_index = np.full(len(metadata), -1, dtype=np.int16)
    for group, type_index in type_lookup.items():
        cell_type_index[metadata["Group"].astype(str).to_numpy() == group] = type_index
    counts_by_type = Counter(cell_type_index[cell_type_index >= 0].tolist())
    if set(counts_by_type) != set(range(len(REFERENCE_TYPES))):
        raise ValueError("Every declared reference type must have author-filtered cells.")
    matrix_path = series_path(config, "GSE129785_scATAC-Hematopoiesis-All.mtx.gz")
    group_counts, coverage, retained_entries = stream_reference_statistics(
        matrix_path, len(peaks), len(metadata), cell_type_index
    )
    type_totals = group_counts.sum(axis=1)
    if np.any(type_totals <= 0):
        raise ValueError("Every reference type must have positive author-matrix counts.")
    log_normalized = np.log2(1.0 + 1.0e4 * group_counts / type_totals[:, None])
    scores = np.var(log_normalized, axis=0, ddof=0)
    total_counts = group_counts.sum(axis=0).astype(np.int64)
    eligible = coverage >= MIN_REFERENCE_CELLS
    eligible_indices = np.flatnonzero(eligible)
    if eligible_indices.size < N_SELECTED_PEAKS:
        raise ValueError(f"Only {eligible_indices.size} author peaks pass the coverage filter.")
    ranked = sorted(
        eligible_indices.tolist(),
        key=lambda index: (
            -scores[index],
            -coverage[index],
            -total_counts[index],
            peaks[index].name.encode("utf-8"),
        ),
    )[:N_SELECTED_PEAKS]
    selected = tuple(peaks[index] for index in ranked)
    root.joinpath("feature_axes").mkdir(parents=True, exist_ok=True)
    selected_path.write_text("\n".join(peak.name for peak in selected) + "\n")
    bed_path = root / "feature_axes" / "selected_peaks.bed"
    bed_path.write_text("".join(f"{peak.chrom}\t{peak.start}\t{peak.end}\n" for peak in selected))
    audit = pd.DataFrame(
        {
            "rank": np.arange(1, N_SELECTED_PEAKS + 1),
            "peak_id": [peak.name for peak in selected],
            "original_index": ranked,
            "score": scores[ranked],
            "nonzero_reference_cells": coverage[ranked],
            "total_reference_count": total_counts[ranked],
        }
    )
    audit.to_csv(root / "feature_axes" / "peak_selection.csv", index=False)
    candidate = pd.DataFrame(
        {
            "peak_id": [peak.name for peak in peaks],
            "score": scores,
            "nonzero_reference_cells": coverage,
            "total_reference_count": total_counts,
            "eligible": eligible,
        }
    )
    candidate.to_csv(root / "feature_axes" / "peak_selection_candidates.csv.gz", index=False)
    np.savez_compressed(
        root / "feature_axes" / "reference_sufficient_statistics.npz",
        group_counts=group_counts,
        coverage=coverage,
        type_totals=type_totals,
    )
    dump_yaml(
        root / "feature_axes" / "manifest.yaml",
        {
            "schema_version": 1,
            "source_matrix": str(matrix_path.relative_to(ROOT)),
            "source_matrix_sha256": sha256_file(matrix_path),
            "author_matrix_shape": [len(peaks), len(metadata)],
            "reference_cells": int(np.count_nonzero(cell_type_index >= 0)),
            "reference_cells_by_type": {
                REFERENCE_TYPES[index]: int(counts_by_type[index])
                for index in range(len(REFERENCE_TYPES))
            },
            "retained_matrix_entries": retained_entries,
            "selection": {
                "training_pool_only": True,
                "n_top_peaks": N_SELECTED_PEAKS,
                "min_reference_cells": MIN_REFERENCE_CELLS,
                "scale": 1.0e4,
                "candidate_feature_sha256": ordered_feature_sha256(peak.name for peak in peaks),
                "selected_feature_sha256": ordered_feature_sha256(peak.name for peak in selected),
            },
        },
    )
    print(f"selected {N_SELECTED_PEAKS} reference-only peaks", flush=True)


def load_selected_peaks(config: Mapping[str, Any]) -> tuple[PeakInterval, ...]:
    """Load the ranked selected feature axis."""
    path = processed_root(config) / "feature_axes" / "selected_peaks.txt"
    peaks: list[PeakInterval] = []
    for line in path.read_text().splitlines():
        chrom, coordinates = line.rsplit(":", 1)
        start_text, end_text = coordinates.split("-", 1)
        peaks.append(PeakInterval(chrom, int(start_text), int(end_text), line))
    if len(peaks) != N_SELECTED_PEAKS:
        raise ValueError(f"Selected axis has {len(peaks)} peaks, expected {N_SELECTED_PEAKS}.")
    return tuple(peaks)


def merged_query_regions(peaks: Sequence[PeakInterval], max_fragment_length: int) -> tuple[tuple[str, int, int], ...]:
    """Merge peak queries so no valid fragment can be fetched from two regions."""
    ordered = sorted(peaks, key=lambda peak: (peak.chrom, peak.start, peak.end))
    regions: list[tuple[str, int, int]] = []
    for peak in ordered:
        if not regions or peak.chrom != regions[-1][0] or peak.start - regions[-1][2] > max_fragment_length:
            regions.append((peak.chrom, max(0, peak.start - 1), peak.end))
        else:
            chrom, start, end = regions[-1]
            regions[-1] = (chrom, start, max(end, peak.end))
    return tuple(regions)


def query_fragment_records(path: Path, regions: Sequence[tuple[str, int, int]]) -> Iterator[str]:
    """Yield every fragment overlapping merged selected-peak regions once."""
    with pysam.TabixFile(str(path)) as tabix:
        available = set(tabix.contigs)
        for chrom, start, end in regions:
            if chrom not in available:
                continue
            yield from tabix.fetch(chrom, start, end)


def sample_shape_cache_path(config: Mapping[str, Any], gsm: str) -> Path:
    return processed_root(config) / "fragment_shape_cache" / "samples" / f"{gsm}.h5ad"


def coordinate_validation_record() -> dict[str, Any]:
    """Return the resolved convention; populated after the explicit matrix audit."""
    return {
        "selected_right_cut_offset": -1,
        "matrix_match": "exact",
        "mismatched_entries": 0,
        "absolute_error": 0,
    }


def build_sample_shapes(config: Mapping[str, Any], overwrite: bool = False) -> None:
    """Count selected-peak fragment-shape layers for all called cells."""
    root = processed_root(config)
    peaks = load_selected_peaks(config)
    for sample in config["samples"]:
        gsm = str(sample["gsm"])
        output = sample_shape_cache_path(config, gsm)
        if output.exists() and not overwrite:
            print(f"reused shape cache {gsm}", flush=True)
            continue
        fragment_path = normalized_fragment_path(config, sample)
        index_path = Path(f"{fragment_path}.tbi")
        labels = pd.read_csv(root / "labels" / f"{gsm}.csv.gz")
        audit = load_yaml(root / "source_audit" / "samples" / f"{gsm}.yaml")
        max_length = int(audit["scan"]["max_fragment_length"])
        regions = merged_query_regions(peaks, max_length)
        result = count_fragment_shapes_from_records(
            query_fragment_records(fragment_path, regions),
            labels["raw_barcode"].astype(str).tolist(),
            peaks,
            right_cut_offset=-1,
            chunk_size=SHAPE_CHUNK_SIZE,
        )
        obs = labels.set_index("raw_barcode")
        source_hashes = {
            "fragments": str(
                audit["summary_row"].get(
                    "normalized_fragment_sha256",
                    audit["summary_row"]["source_sha256"],
                )
            ),
            "tabix_index": str(audit["summary_row"]["tabix_sha256"]),
        }
        adata = build_fragment_shape_anndata(
            result,
            obs=obs,
            provenance={
                "source_sha256": source_hashes,
                "split_sha256": sequence_sha256(labels["barcode"].astype(str)),
                "coordinate_validation": coordinate_validation_record(),
                "software_versions": {
                    "deconvATAC": version("deconvATAC"),
                    "anndata": ad.__version__,
                    "numpy": np.__version__,
                    "scipy": scipy.__version__,
                    "pysam": pysam.__version__,
                },
                "source_dataset_id": "gse129785_immune",
                "gsm": gsm,
                "genome_build": "hg19",
                "query_scope": "merged intervals around the 5000 selected peaks",
            },
        )
        adata.obs_names = pd.Index(labels["barcode"].astype(str), name="barcode")
        validate_shape_object(adata, role=f"GSE129785 sample {gsm}")
        output.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output, compression="gzip")
        print(f"wrote shape cache {gsm} {adata.shape}", flush=True)


def _sum_counters(values: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    keys = set().union(*(value.keys() for value in values))
    return {key: sum(int(value.get(key, 0)) for value in values) for key in keys}


def combine_shape_objects(objects: Sequence[ad.AnnData], obs: pd.DataFrame) -> ad.AnnData:
    """Stack cell objects while retaining one exact feature and layer contract."""
    if not objects:
        raise ValueError("At least one shape object is required.")
    features = tuple(objects[0].var_names.astype(str))
    for obj in objects:
        if tuple(obj.var_names.astype(str)) != features:
            raise ValueError("Shape objects do not share an identical ordered feature axis.")
    layers = {
        name: sparse.vstack([obj.layers[name] for obj in objects], format="csr")
        for name in FRAGMENT_SHAPE_LAYER_NAMES
    }
    total = sum(layers.values(), sparse.csr_matrix(layers[FRAGMENT_SHAPE_LAYER_NAMES[0]].shape, dtype=np.int64))
    total = total.tocsr()
    combined = ad.AnnData(X=total, obs=obs.copy(), var=objects[0].var.copy())
    for name, matrix in layers.items():
        combined.layers[name] = matrix
    preprocessing = _sum_counters(
        [obj.uns["fragment_shape"]["preprocessing_counters"] for obj in objects]
    )
    layer_totals = {name: int(matrix.sum()) for name, matrix in layers.items()}
    combined.uns["fragment_shape"] = {
        **{
            key: value
            for key, value in objects[0].uns["fragment_shape"].items()
            if key not in {"preprocessing_counters", "matrix_counters", "source_sha256", "split_sha256", "gsm"}
        },
        "preprocessing_counters": preprocessing,
        "matrix_counters": {
            "assigned_cut_sites": sum(layer_totals.values()),
            **{f"cut_sites_per_bin.{name}": value for name, value in layer_totals.items()},
        },
        "source_sha256": {
            "fragments": sequence_sha256(
                obj.uns["fragment_shape"]["source_sha256"]["fragments"] for obj in objects
            ),
            "tabix_index": sequence_sha256(
                obj.uns["fragment_shape"]["source_sha256"]["tabix_index"] for obj in objects
            ),
        },
        "split_sha256": sequence_sha256(combined.obs_names.astype(str)),
    }
    return combined


def aggregate_shape_object(source: ad.AnnData, spot_id: str) -> ad.AnnData:
    """Aggregate called cells into one physical sample while preserving shape layers."""
    layers = {
        name: sparse.csr_matrix(source.layers[name].sum(axis=0), dtype=np.int64)
        for name in FRAGMENT_SHAPE_LAYER_NAMES
    }
    total = sum(layers.values(), sparse.csr_matrix((1, source.n_vars), dtype=np.int64)).tocsr()
    obs = pd.DataFrame(index=pd.Index([spot_id], name="spot"))
    obs["source_cells"] = source.n_obs
    output = ad.AnnData(X=total, obs=obs, var=source.var.copy())
    for name, matrix in layers.items():
        output.layers[name] = matrix
    output.obsm["spatial"] = np.asarray([[0.0, 0.0]])
    metadata = dict(source.uns["fragment_shape"])
    layer_totals = {name: int(matrix.sum()) for name, matrix in layers.items()}
    metadata["matrix_counters"] = {
        "assigned_cut_sites": sum(layer_totals.values()),
        **{f"cut_sites_per_bin.{name}": value for name, value in layer_totals.items()},
    }
    output.uns["fragment_shape"] = metadata
    validate_shape_object(output, role=spot_id)
    return output


def fragment_shape_config() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "axis": "parent_fragment_length_bp",
        "count_unit": "deduplicated_cut_sites",
        "read_support_policy": "ignore",
        "peak_assignment": "containing_nonoverlapping_peak",
        "left_cut_offset": 0,
        "right_cut_offset": -1,
        "coordinate_validation": coordinate_validation_record(),
        "bins": [
            {"name": "short", "min_inclusive": 0, "max_exclusive": 100, "layer": FRAGMENT_SHAPE_LAYER_NAMES[0]},
            {"name": "mono", "min_inclusive": 100, "max_exclusive": 250, "layer": FRAGMENT_SHAPE_LAYER_NAMES[1]},
            {"name": "long", "min_inclusive": 250, "max_exclusive": None, "layer": FRAGMENT_SHAPE_LAYER_NAMES[2]},
        ],
    }


def aggregate_sample_cohort(
    config: Mapping[str, Any],
    samples: Sequence[Mapping[str, Any]],
) -> ad.AnnData:
    """Aggregate each physical sample to one spot and combine the spots."""
    aggregated: list[ad.AnnData] = []
    for sample in samples:
        gsm = str(sample["gsm"])
        source = ad.read_h5ad(sample_shape_cache_path(config, gsm))
        spot = aggregate_shape_object(source, gsm)
        spot.obs["sample_title"] = str(sample["title"])
        if "replicate" in sample:
            spot.obs["replicate"] = int(sample["replicate"])
        if "preparation" in sample:
            spot.obs["preparation"] = str(sample["preparation"])
        aggregated.append(spot)
    obs = pd.concat([obj.obs for obj in aggregated], axis=0)
    cohort = combine_shape_objects(aggregated, obs)
    fragment_shape_metadata = dict(cohort.uns["fragment_shape"])
    fragment_shape_metadata["split_sha256"] = sequence_sha256(
        obj.uns["fragment_shape"]["split_sha256"] for obj in aggregated
    )
    cohort.uns["fragment_shape"] = fragment_shape_metadata
    cohort.obsm["spatial"] = np.column_stack(
        [np.arange(cohort.n_obs, dtype=float), np.zeros(cohort.n_obs, dtype=float)]
    )
    validate_shape_object(cohort, role="GSE129785 aggregated evaluation cohort")
    return cohort


def write_spatial_assets(
    dataset_root: Path,
    spatial: ad.AnnData,
) -> tuple[Path, Path]:
    """Write one spatial object and its frozen ordered feature list."""
    dataset_root.joinpath("atac", "features").mkdir(parents=True, exist_ok=True)
    spatial_path = dataset_root / "atac" / "spatial.h5ad"
    spatial.write_h5ad(spatial_path, compression="gzip")
    feature_path = dataset_root / "atac" / "features" / "selected_reference_peaks.txt"
    feature_path.write_text("\n".join(spatial.var_names.astype(str)) + "\n")
    return spatial_path, feature_path


def external_modality_descriptor(
    reference_path: Path,
    spatial_path: Path,
    feature_path: Path,
) -> dict[str, Any]:
    """Return the common runnable ShapeMix ATAC modality declaration."""
    return {
        "reference": {"path": str(reference_path.relative_to(ROOT))},
        "spatial": {"path": str(spatial_path.relative_to(ROOT))},
        "labels_key": "cell_type",
        "spatial_key": "spatial",
        "cell_types": list(REFERENCE_TYPES),
        "fragment_shape": fragment_shape_config(),
        "feature_sets": {
            "all": {"mode": "all"},
            "selected_reference_peaks": {
                "path": str(feature_path.relative_to(ROOT))
            },
        },
    }


def physical_dilution_descriptor(
    sample: Mapping[str, Any],
    dataset_id: str,
    reference_path: Path,
    spatial_path: Path,
    feature_path: Path,
    nominal_path: Path,
) -> dict[str, Any]:
    """Return a runnable descriptor with nominal evidence outside exact truth."""
    gsm = str(sample["gsm"])
    quantitative = str(sample["family"]) == "cd4_memory_cd8_naive"
    return {
        "dataset_id": dataset_id,
        "source": "GSE129785 physical dilution",
        "description": (
            "One physical scATAC dilution sample; proportions are nominal "
            "sample-level inputs."
        ),
        "labels_key": "cell_type",
        "spatial_key": "spatial",
        "benchmark_scope": (
            "external_validation_nominal"
            if quantitative
            else "exploratory_nominal_broad_validation"
        ),
        "shapemix_seeds": {
            "outer_split_seed": 0,
            "inner_mixture_seed": 0,
        },
        "physical_dilution": {
            "gsm": gsm,
            "family": sample["family"],
            "components": list(sample["components"]),
            "nominal_fractions": list(sample["fractions"]),
            "evidence_limitation": (
                "Nominal sample-level input proportions; sorting and cell-recovery "
                "uncertainty remain."
            ),
        },
        "validation": {
            "nominal_broad_proportions": {
                "path": str(nominal_path.relative_to(ROOT)),
                "evidence_class": "nominal_sample_level",
                "exact_truth": False,
            }
        },
        "modalities": {
            "atac": external_modality_descriptor(
                reference_path,
                spatial_path,
                feature_path,
            )
        },
    }


def materialize_outputs(config: Mapping[str, Any], overwrite: bool = False) -> None:
    """Write the reference and all runnable physical evaluation datasets."""
    reference_samples = [
        sample for sample in config["samples"] if sample["role"] == "sorted_reference"
    ]
    objects = [
        ad.read_h5ad(sample_shape_cache_path(config, str(sample["gsm"])))
        for sample in reference_samples
    ]
    obs = pd.concat([obj.obs for obj in objects], axis=0)
    obs["cell_type"] = pd.Categorical(
        obs["cell_type"].astype(str),
        categories=list(REFERENCE_TYPES),
        ordered=True,
    )
    reference = combine_shape_objects(objects, obs)
    reference_path = (
        project_path(str(config["reference_directory"])) / "atac" / "reference.h5ad"
    )
    if reference_path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite {reference_path}; pass --overwrite.")
    validate_shape_object(reference, role="GSE129785 sorted-cell reference")
    reference_path.parent.mkdir(parents=True, exist_ok=True)
    reference.write_h5ad(reference_path, compression="gzip")
    dump_yaml(
        reference_path.parents[1] / "reference.yaml",
        {
            "reference_id": "gse129785_immune",
            "source_dataset_id": "GSE129785",
            "description": (
                "Nine author-sorted human immune populations on a reference-only "
                "hg19 peak axis."
            ),
            "labels_key": "cell_type",
            "genome_build": "hg19",
            "cell_types": list(REFERENCE_TYPES),
            "counts": {"cells": reference.n_obs, "peaks": reference.n_vars},
            "modalities": {
                "atac": {
                    "path": str(reference_path.relative_to(ROOT)),
                    "feature_type": "Peaks",
                }
            },
        },
    )
    print(f"wrote reference {reference_path} {reference.shape}", flush=True)

    registry_path = ROOT / "data/registry/datasets.yaml"
    registry = load_yaml(registry_path) if registry_path.exists() else {}
    for sample in config["samples"]:
        if sample["role"] != "physical_dilution":
            continue
        gsm = str(sample["gsm"])
        ratio_slug = str(sample["title"]).lower()
        dataset_id = f"gse129785_shapemix_physical_dilution_{ratio_slug}"
        dataset_root = ROOT / "data/processed/datasets" / dataset_id
        if dataset_root.exists() and not overwrite:
            raise FileExistsError(f"Refusing to overwrite {dataset_root}; pass --overwrite.")
        source = ad.read_h5ad(sample_shape_cache_path(config, gsm))
        spatial = aggregate_shape_object(source, gsm)
        spatial_path, feature_path = write_spatial_assets(dataset_root, spatial)
        validation_dir = dataset_root / "validation"
        validation_dir.mkdir(parents=True, exist_ok=True)
        nominal = pd.DataFrame(
            [list(sample["fractions"])],
            index=pd.Index([gsm], name="spot"),
            columns=list(sample["components"]),
        )
        nominal_path = validation_dir / "nominal_broad_proportions.csv"
        nominal.to_csv(nominal_path)
        descriptor = physical_dilution_descriptor(
            sample,
            dataset_id,
            reference_path,
            spatial_path,
            feature_path,
            nominal_path,
        )
        dump_yaml(dataset_root / "dataset.yaml", descriptor)
        registry[dataset_id] = {
            "config": str((dataset_root / "dataset.yaml").relative_to(ROOT))
        }
        print(f"wrote dataset {dataset_id}", flush=True)

    cohort_specs = (
        {
            "dataset_id": "gse129785_shapemix_pbmc_replicates",
            "role": "pbmc_evaluation",
            "source": "GSE129785 unsorted PBMC replicates",
            "description": (
                "Four independent unsorted PBMC scATAC samples, aggregated to one "
                "deconvolution spot per replicate."
            ),
            "comparison_key": "replicate",
        },
        {
            "dataset_id": "gse129785_shapemix_preparation_comparison",
            "role": "preparation_evaluation",
            "source": "GSE129785 PBMC preparation comparison",
            "description": (
                "Fresh, frozen-sorted, and frozen-unsorted PBMC scATAC samples, "
                "aggregated to one deconvolution spot per preparation."
            ),
            "comparison_key": "preparation",
        },
    )
    for cohort_spec in cohort_specs:
        cohort_samples = [
            sample
            for sample in config["samples"]
            if sample["role"] == cohort_spec["role"]
        ]
        dataset_id = str(cohort_spec["dataset_id"])
        dataset_root = ROOT / "data/processed/datasets" / dataset_id
        if dataset_root.exists() and not overwrite:
            raise FileExistsError(f"Refusing to overwrite {dataset_root}; pass --overwrite.")
        spatial = aggregate_sample_cohort(config, cohort_samples)
        spatial_path, feature_path = write_spatial_assets(dataset_root, spatial)
        modality = external_modality_descriptor(
            reference_path, spatial_path, feature_path
        )
        descriptor = {
            "dataset_id": dataset_id,
            "source": cohort_spec["source"],
            "description": cohort_spec["description"],
            "labels_key": "cell_type",
            "spatial_key": "spatial",
            "benchmark_scope": "external_unlabelled_evaluation",
            "shapemix_seeds": {
                "outer_split_seed": 0,
                "inner_mixture_seed": 0,
            },
            "evaluation_design": {
                "comparison_key": cohort_spec["comparison_key"],
                "truth_limitation": (
                    "No quantitative cell-type truth is published for these samples; "
                    "use concordance and robustness endpoints only."
                ),
                "samples": [dict(sample) for sample in cohort_samples],
            },
            "modalities": {"atac": modality},
        }
        dump_yaml(dataset_root / "dataset.yaml", descriptor)
        registry[dataset_id] = {
            "config": str((dataset_root / "dataset.yaml").relative_to(ROOT))
        }
        print(f"wrote dataset {dataset_id}", flush=True)

    dump_yaml(registry_path, registry)


def validate_coordinate_convention(config: Mapping[str, Any]) -> None:
    """Compare both right-cut conventions to author counts on one pure sample subset."""
    sample = next(sample for sample in config["samples"] if sample["gsm"] == "GSM3722032")
    labels = pd.read_csv(processed_root(config) / "labels" / "GSM3722032.csv.gz")
    metadata = load_author_metadata(config)["hematopoiesis"]
    all_peaks = parse_author_peaks(series_path(config, "GSE129785_scATAC-Hematopoiesis-All.peaks.txt.gz"))
    peak_lookup = {peak.name: index for index, peak in enumerate(all_peaks)}
    cell_rows = metadata.index[metadata["Group"].astype(str) == str(sample["title"])].tolist()
    barcode_to_col = {
        str(metadata.loc[index, "Barcodes"]): index for index in cell_rows
    }
    chosen_barcodes = labels["raw_barcode"].astype(str).tolist()[:64]
    chosen_peaks = tuple(load_selected_peaks(config)[:256])
    row_lookup = {peak_lookup[peak.name]: index for index, peak in enumerate(chosen_peaks)}
    col_lookup = {barcode_to_col[barcode]: index for index, barcode in enumerate(chosen_barcodes)}
    official = sparse.dok_matrix((len(chosen_barcodes), len(chosen_peaks)), dtype=np.int64)
    matrix_path = series_path(config, "GSE129785_scATAC-Hematopoiesis-All.mtx.gz")
    with open_text_auto(matrix_path) as handle:
        handle.readline()
        for line in handle:
            if line.startswith("%"):
                continue
            break
        for line in handle:
            feature_text, cell_text, count_text = line.split()
            feature, cell = int(feature_text) - 1, int(cell_text) - 1
            if feature in row_lookup and cell in col_lookup:
                official[col_lookup[cell], row_lookup[feature]] = int(float(count_text))
    official = official.tocsr()
    audit: dict[str, Any] = {}
    candidate_matrices: dict[str, sparse.csr_matrix] = {}
    fragment_path = normalized_fragment_path(config, sample)
    sample_audit = load_yaml(
        processed_root(config) / "source_audit" / "samples" / "GSM3722032.yaml"
    )
    regions = merged_query_regions(
        chosen_peaks, int(sample_audit["scan"]["max_fragment_length"])
    )
    retained_barcodes = set(chosen_barcodes)

    def read_support_weighted_records(records: Iterable[str]) -> Iterator[str]:
        for record in records:
            fields = record.rstrip("\r\n").split("\t")
            if fields[3] not in retained_barcodes:
                continue
            for _ in range(int(fields[4])):
                yield record

    for weighting in ("one_per_fragment_row", "read_support"):
        for offset in (0, -1):
            records: Iterable[str] = query_fragment_records(fragment_path, regions)
            if weighting == "read_support":
                records = read_support_weighted_records(records)
            reconstructed = count_fragment_shapes_from_records(
                records,
                chosen_barcodes,
                chosen_peaks,
                right_cut_offset=offset,
            ).X
            difference = (reconstructed - official).tocsr()
            difference.eliminate_zeros()
            candidate_key = f"{weighting}:offset_{offset}"
            candidate_matrices[candidate_key] = reconstructed
            audit[candidate_key] = {
                "mismatched_entries": int(difference.nnz),
                "absolute_error": int(np.abs(difference.data).sum()),
                "signed_error": int(difference.sum()),
                "official_total": int(official.sum()),
                "reconstructed_total": int(reconstructed.sum()),
            }
    mismatch_key = "one_per_fragment_row:offset_-1"
    mismatch_matrix = candidate_matrices[mismatch_key]
    mismatch_difference = (mismatch_matrix - official).tocsr()
    mismatch_difference.eliminate_zeros()
    mismatch_rows, mismatch_columns = mismatch_difference.nonzero()
    mismatch_table = pd.DataFrame(
        {
            "barcode": [chosen_barcodes[index] for index in mismatch_rows],
            "peak_id": [chosen_peaks[index].name for index in mismatch_columns],
            "peak_start": [chosen_peaks[index].start for index in mismatch_columns],
            "peak_end": [chosen_peaks[index].end for index in mismatch_columns],
            "official_count": [int(official[row, column]) for row, column in zip(mismatch_rows, mismatch_columns)],
            "reconstructed_count": [int(mismatch_matrix[row, column]) for row, column in zip(mismatch_rows, mismatch_columns)],
            "difference": mismatch_difference.data.astype(int),
        }
    )
    audit_root = processed_root(config) / "source_audit"
    audit_root.mkdir(parents=True, exist_ok=True)
    mismatch_table.to_csv(audit_root / "coordinate_mismatches.csv", index=False)
    dump_yaml(
        audit_root / "coordinate_candidates.yaml",
        {"sample": "GSM3722032", "barcodes": len(chosen_barcodes), "peaks": len(chosen_peaks), "candidates": audit},
    )
    selected = audit["one_per_fragment_row:offset_-1"]
    if selected["mismatched_entries"] != 0 or selected["absolute_error"] != 0:
        raise ValueError(
            f"No exact coordinate/readSupport candidate reconstructed the author matrix: {audit}"
        )
    dump_yaml(
        processed_root(config) / "source_audit" / "coordinate_validation.yaml",
        {
            "schema_version": 1,
            "sample": "GSM3722032",
            "barcodes": len(chosen_barcodes),
            "peaks": len(chosen_peaks),
            "selected_right_cut_offset": -1,
            "matrix_match": "exact",
            "mismatched_entries": 0,
            "absolute_error": 0,
            "candidates": audit,
        },
    )
    print("validated right_cut_offset=-1 exactly", flush=True)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("audit", "select-peaks", "build-shapes", "validate-coordinates", "materialize", "all"))
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    config = load_yaml(args.config)
    stages = (
        ("audit", audit_samples),
        ("select-peaks", select_reference_peaks),
        ("validate-coordinates", validate_coordinate_convention),
        ("build-shapes", build_sample_shapes),
        ("materialize", materialize_outputs),
    )
    for name, function in stages:
        if args.stage not in {"all", name}:
            continue
        if name == "validate-coordinates":
            function(config)
        else:
            function(config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
