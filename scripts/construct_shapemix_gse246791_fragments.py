#!/usr/bin/env python3
"""Construct and validate GSE246791 fragments with author SnapATAC2 semantics."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import tempfile
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, TextIO

import h5py
import numpy as np
from scipy.stats import spearmanr
import yaml


ROOT = Path(__file__).resolve().parents[1]
CHUNK_BYTES = 8 * 1024 * 1024
BARCODE_REGEX = r"^(\w+):.+"
CHUNK_SIZE = 5_000_000


def project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def repository_path(path: Path) -> str:
    return str(path.absolute().relative_to(ROOT.absolute()))


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_h5ad_fragment_counts(path: Path) -> tuple[list[str], np.ndarray, dict[str, int]]:
    with h5py.File(path, "r") as handle:
        barcodes = [
            value.decode("ascii") if isinstance(value, bytes) else str(value)
            for value in handle["obs/index"][:]
        ]
        counts = np.asarray(handle["obs/n_fragment"][:], dtype=np.int64)
        reference_group = handle["uns/reference_sequences"]
        names = [
            value.decode("ascii") if isinstance(value, bytes) else str(value)
            for value in reference_group["reference_seq_name"][:]
        ]
        lengths = np.asarray(
            reference_group["reference_seq_length"][:], dtype=np.int64
        )
    if len(barcodes) != len(counts) or len(barcodes) != len(set(barcodes)):
        raise ValueError(f"Invalid H5AD barcode/count contract: {path}")
    if len(names) != len(lengths) or len(names) != len(set(names)):
        raise ValueError(f"Invalid H5AD reference-sequence contract: {path}")
    return barcodes, counts, dict(zip(names, lengths.tolist(), strict=True))


def open_fragment_text(path: Path) -> TextIO:
    """Open compiled-backend output by content, independent of its work suffix."""
    with path.open("rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzip else path.open("rt")


def filter_fragments(
    source: Path,
    destination: Path,
    *,
    barcodes: list[str],
    chrom_sizes: Mapping[str, int],
) -> dict[str, Any]:
    selected = set(barcodes)
    per_barcode: Counter[str] = Counter()
    raw_rows = 0
    kept_rows = 0
    unknown_barcode_rows = 0
    excluded_contig_rows = 0
    duplicate_support_total = 0
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.stem}.partial.tsv.gz")
    if destination.exists():
        raise FileExistsError(destination)
    if temporary.exists():
        raise FileExistsError(f"Unresolved fragment partial: {temporary}")
    try:
        with open_fragment_text(source) as input_handle, gzip.open(
            temporary, "xt", compresslevel=6
        ) as output_handle:
            for line_number, line in enumerate(input_handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                raw_rows += 1
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 5:
                    raise ValueError(
                        f"Expected five fragment fields at row {line_number}; found {len(fields)}"
                    )
                chrom, start_text, end_text, barcode, support_text = fields
                try:
                    start = int(start_text)
                    end = int(end_text)
                    support = int(support_text)
                except ValueError as exc:
                    raise ValueError(f"Non-integer fragment field at row {line_number}") from exc
                if support < 1:
                    raise ValueError(f"Invalid duplicate support at row {line_number}")
                if barcode not in selected:
                    unknown_barcode_rows += 1
                    continue
                chrom_length = chrom_sizes.get(chrom)
                if chrom_length is None:
                    excluded_contig_rows += 1
                    continue
                if start < 0 or end <= start or end > chrom_length:
                    raise ValueError(f"Out-of-bounds fragment at row {line_number}")
                output_handle.write(line if line.endswith("\n") else line + "\n")
                kept_rows += 1
                per_barcode[barcode] += 1
                duplicate_support_total += support
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    observed = np.asarray([per_barcode[barcode] for barcode in barcodes], dtype=np.int64)
    return {
        "raw_rows": raw_rows,
        "kept_rows": kept_rows,
        "unknown_barcode_rows": unknown_barcode_rows,
        "excluded_contig_rows": excluded_contig_rows,
        "duplicate_support_total": duplicate_support_total,
        "barcodes_expected": len(barcodes),
        "barcodes_observed": int(np.count_nonzero(observed)),
        "observed_counts": observed,
    }


def fragment_count_concordance(
    expected: np.ndarray,
    observed: np.ndarray,
    *,
    thresholds: Mapping[str, Any],
) -> dict[str, Any]:
    expected = np.asarray(expected, dtype=np.float64)
    observed = np.asarray(observed, dtype=np.float64)
    if expected.shape != observed.shape or expected.ndim != 1:
        raise ValueError("Expected and observed fragment-count vectors must align")
    if np.any(expected <= 0) or np.any(observed < 0):
        raise ValueError("Fragment counts must have positive expectations and nonnegative observations")
    absolute_relative_error = np.abs(observed - expected) / expected
    correlation = float(spearmanr(expected, observed).statistic)
    median_error = float(np.median(absolute_relative_error))
    p95_error = float(np.quantile(absolute_relative_error, 0.95))
    missing = int(np.count_nonzero(observed == 0))
    gates = {
        "all_h5ad_barcodes_present": (
            missing == 0
            if bool(thresholds["require_all_h5ad_barcodes_present"])
            else True
        ),
        "median_absolute_relative_error": (
            median_error <= float(thresholds["median_absolute_relative_error_max"])
        ),
        "p95_absolute_relative_error": (
            p95_error <= float(thresholds["p95_absolute_relative_error_max"])
        ),
        "spearman_r": correlation >= float(thresholds["spearman_r_min"]),
    }
    return {
        "cells": int(expected.size),
        "missing_barcodes": missing,
        "expected_total": int(expected.sum()),
        "observed_total": int(observed.sum()),
        "median_absolute_relative_error": median_error,
        "p95_absolute_relative_error": p95_error,
        "spearman_r": correlation,
        "thresholds": dict(thresholds),
        "gates": gates,
        "passed": all(gates.values()),
    }


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "configs/data_sources/shapemix_gse246791_fragment_reads.yaml",
    )
    parser.add_argument("--gsm", default="GSM7877011")
    parser.add_argument("--bam", type=Path)
    parser.add_argument("--h5ad", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cleanup-bam", action="store_true")
    return parser.parse_args()

def cleanup_alignment_bam(bam: Path, gsm: str, *, enabled: bool) -> None:
    if not enabled:
        return
    expected = (
        ROOT
        / "data/work/preprocessing/gse246791_mouse_brain_reference"
        / gsm
        / "alignment/name_sorted.bam"
    ).absolute()
    if bam.absolute() != expected:
        raise ValueError("--cleanup-bam is restricted to the canonical generated work BAM")
    if bam.is_symlink():
        raise ValueError(f"Refusing to remove a symlinked work BAM: {bam}")
    if bam.exists():
        if not bam.is_file():
            raise ValueError(f"Canonical work BAM is not a regular file: {bam}")
        removed_bytes = bam.stat().st_size
        bam.unlink()
        print(
            f"alignment_work_cleanup gsm={gsm} bytes={removed_bytes} status=completed",
            flush=True,
        )



def main() -> None:
    args = parse_args()
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run fragment construction through run_shapemix_low_impact.sh")
    if os.environ.get("RAYON_NUM_THREADS") != "1":
        raise RuntimeError("RAYON_NUM_THREADS=1 is required for SnapATAC2 containment")
    config = load_yaml(args.config)
    bam = args.bam or (
        ROOT
        / "data/work/preprocessing/gse246791_mouse_brain_reference"
        / args.gsm
        / "alignment/name_sorted.bam"
    )
    if args.h5ad is None:
        matches = sorted(
            (
                ROOT
                / "data/processed/shapemix/gse246791_mouse_brain_reference/source_audit/source_objects"
            ).glob(f"{args.gsm}_*.h5ad")
        )
        if len(matches) != 1:
            raise ValueError(f"Expected one selected H5AD for {args.gsm}; found {len(matches)}")
        h5ad = matches[0]
    else:
        h5ad = args.h5ad
    family_root = project_path(str(config["processed_directory"]))
    output = args.output or (
        family_root / "normalized_fragments" / args.gsm / "fragments.tsv.gz"
    )
    audit_output = output.parent / "manifest.yaml"
    if output.is_file() and audit_output.is_file():
        cleanup_alignment_bam(bam, args.gsm, enabled=args.cleanup_bam)
        print(f"fragments gsm={args.gsm} status=reused", flush=True)
        return
    if output.exists() or audit_output.exists():
        raise FileExistsError(f"Partial immutable fragment state for {args.gsm}")
    if not bam.is_file():
        raise FileNotFoundError(bam)
    work_root = (
        ROOT
        / "data/work/preprocessing/gse246791_mouse_brain_reference"
        / args.gsm
        / "fragments"
    )
    work_root.mkdir(parents=True, exist_ok=True)
    raw_fragments = work_root / "snapatac2.fragments.tsv.gz"
    if raw_fragments.exists():
        raise FileExistsError(f"Unresolved raw fragment intermediate: {raw_fragments}")
    candidate = work_root / "whitelisted.fragments.tsv.gz"
    if candidate.exists():
        raise FileExistsError(f"Unresolved whitelisted fragment candidate: {candidate}")
    import snapatac2 as snap

    started = time.monotonic()
    statistics = snap.pp.make_fragment_file(
        bam_file=bam,
        output_file=raw_fragments,
        is_paired=True,
        barcode_regex=BARCODE_REGEX,
        shift_left=4,
        shift_right=-5,
        min_mapq=30,
        chunk_size=CHUNK_SIZE,
    )
    barcodes, expected_counts, chrom_sizes = load_h5ad_fragment_counts(h5ad)
    filtering = filter_fragments(
        raw_fragments,
        candidate,
        barcodes=barcodes,
        chrom_sizes=chrom_sizes,
    )
    observed_counts = filtering.pop("observed_counts")
    thresholds = config["preprocessing_policy"]["representative_concordance_gate"]
    concordance = fragment_count_concordance(
        expected_counts,
        observed_counts,
        thresholds=thresholds,
    )
    record = {
        "schema_version": 1,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "gsm": args.gsm,
        "inputs": {
            "bam": repository_path(bam),
            "bam_sha256": file_digest(bam),
            "h5ad": repository_path(h5ad),
            "h5ad_sha256": file_digest(h5ad),
        },
        "method": {
            "snapatac2_version": snap.__version__,
            "is_paired": True,
            "barcode_regex": BARCODE_REGEX,
            "shift_left": 4,
            "shift_right": -5,
            "min_mapq": 30,
            "chunk_size": CHUNK_SIZE,
            "statistics": str(statistics),
        },
        "filtering": filtering,
        "concordance": concordance,
        "output": {
            "path": repository_path(output),
            "bytes": candidate.stat().st_size,
            "sha256": file_digest(candidate),
        },
        "elapsed_seconds": time.monotonic() - started,
        "resource_policy": {
            "rayon_threads": int(os.environ["RAYON_NUM_THREADS"]),
            "guard": "scripts/run_shapemix_low_impact.sh",
        },
    }
    if not concordance["passed"]:
        atomic_yaml(work_root / "concordance_failure.yaml", record)
        raise ValueError(
            f"Representative fragment-count concordance gate failed: {concordance['gates']}"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    os.replace(candidate, output)
    atomic_yaml(audit_output, record)
    raw_fragments.unlink()
    cleanup_alignment_bam(bam, args.gsm, enabled=args.cleanup_bam)
    print(
        f"fragments gsm={args.gsm} rows={filtering['kept_rows']} "
        f"spearman={concordance['spearman_r']:.6f} "
        f"seconds={record['elapsed_seconds']:.1f} status=completed",
        flush=True,
    )


if __name__ == "__main__":
    main()
