#!/usr/bin/env python3
"""Materialize the audited GSE216371 E13.5 ShapeMix reference."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import platform
import shutil
import struct
import subprocess
import tarfile
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Mapping

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
import yaml

from deconvatac.data import FragmentShapeSpec, ordered_feature_sha256
from deconvatac.data.validators import (
    validate_fragment_shape_feature_axis,
    validate_fragment_shape_spec,
)
from deconvatac.pp.fragment_shapes import (
    DEFAULT_FRAGMENT_LENGTH_BINS,
    FragmentShapeQC,
    FragmentShapeResult,
    PeakInterval,
    build_fragment_shape_anndata,
)
from scripts.prepare_shapemix_gse216371 import CCRE_VERSION, LABEL_VERSION


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs/data_sources/shapemix_gse216371_reference.yaml"
LOCK_PATH = ROOT / "configs/data_sources/shapemix_gse216371_reference_lock.yaml"
CONFIG = yaml.safe_load(CONFIG_PATH.read_text())
FAMILY_ROOT = ROOT / str(CONFIG["processed_directory"])
AUTHOR_LABEL_ROOT = FAMILY_ROOT / "labels" / LABEL_VERSION
CCRE_ROOT = FAMILY_ROOT / "feature_axes" / CCRE_VERSION
BROAD_LABEL_VERSION = "major_types_v1"
BROAD_LABEL_ROOT = FAMILY_ROOT / "labels" / BROAD_LABEL_VERSION
FRAGMENT_CACHE_ROOT = FAMILY_ROOT / "normalized_fragments" / BROAD_LABEL_VERSION
AXIS_ROOT = FAMILY_ROOT / "feature_axes" / BROAD_LABEL_VERSION
AUDIT_PATH = FAMILY_ROOT / "manifests/embryo_fragment_coordinate_audit.yaml"
REFERENCE_ROOT = (
    ROOT / "data/processed/references" / str(CONFIG["standardized_reference_id"])
)
WORK_ROOT = ROOT / "data/work/preprocessing/gse216371_reference"
STREAM_SOURCE = ROOT / "scripts/shapemix_gse216371_stream.cpp"
CHROM_SIZES = ROOT / "data/raw/sources/ucsc/mm10_initial/mm10.chrom.sizes"
N_TOP_PEAKS = 5_000
MIN_REFERENCE_CELLS_PER_PEAK = 10
FEATURE_CHUNK = 100_000
EVENT_CHUNK = 2_000_000
EVENT_DTYPE = np.dtype(
    [("cell", "<u4"), ("feature", "<u4"), ("layer", "u1")], align=False
)


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


def repository_path(path: Path) -> str:
    return path.absolute().relative_to(ROOT.absolute()).as_posix()


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_yaml(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}."
    )
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
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }


def ontology() -> tuple[tuple[str, ...], dict[str, str], int]:
    value = CONFIG.get("broad_ontology")
    if not isinstance(value, Mapping):
        raise ValueError(
            "GSE216371 broad_ontology must be frozen after the author-label audit "
            "and before fragment feature selection"
        )
    cell_types = tuple(str(item) for item in value.get("ordered_cell_types", ()))
    mapping_value = value.get("author_main_cluster_mapping")
    if (
        not cell_types
        or len(cell_types) != len(set(cell_types))
        or not isinstance(mapping_value, Mapping)
    ):
        raise ValueError("Invalid frozen GSE216371 broad ontology")
    mapping = {str(key): str(item) for key, item in mapping_value.items()}
    if not mapping or set(mapping.values()) != set(cell_types):
        raise ValueError("GSE216371 ontology mapping does not cover its output classes")
    minimum = int(value.get("minimum_cells_per_type", 100))
    if minimum < 1:
        raise ValueError("minimum_cells_per_type must be positive")
    return cell_types, mapping, minimum


def build_broad_labels() -> pd.DataFrame:
    output_path = BROAD_LABEL_ROOT / "cells.tsv.gz"
    manifest_path = BROAD_LABEL_ROOT / "manifest.yaml"
    if output_path.is_file() and manifest_path.is_file():
        return pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    if BROAD_LABEL_ROOT.exists():
        raise FileExistsError(f"Partial immutable broad-label directory: {BROAD_LABEL_ROOT}")
    author_path = AUTHOR_LABEL_ROOT / "cells.tsv.gz"
    author_manifest = AUTHOR_LABEL_ROOT / "manifest.yaml"
    if not author_path.is_file() or not author_manifest.is_file():
        raise FileNotFoundError("Run prepare_shapemix_gse216371.py first")
    cell_types, mapping, minimum = ontology()
    labels = pd.read_csv(author_path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "cell_id",
        "round4_barcode",
        "fragments",
        "author_main_cluster",
    }
    missing = sorted(required.difference(labels.columns))
    if missing:
        raise ValueError(f"Author E13.5 labels lack required columns: {missing}")
    if labels.empty or labels["cell_id"].duplicated().any():
        raise ValueError("Author E13.5 cell IDs are empty or duplicated")
    observed = set(labels["author_main_cluster"])
    missing_mapping = sorted(observed.difference(mapping))
    unused_mapping = sorted(set(mapping).difference(observed))
    if missing_mapping or unused_mapping:
        raise ValueError(
            f"Frozen ontology does not exactly match audited author clusters; "
            f"missing={missing_mapping} unused={unused_mapping}"
        )
    if (labels["round4_barcode"].str.len() == 0).any():
        raise ValueError("An E13.5 label has no Round4 barcode")
    expected = pd.to_numeric(labels["fragments"], errors="raise")
    if (
        not np.isfinite(expected).all()
        or (expected <= 0).any()
        or not np.equal(expected, np.floor(expected)).all()
    ):
        raise ValueError("Workbook fragment totals must be positive integers")
    labels["fragments"] = expected.astype(np.int64).astype(str)
    labels["cell_type"] = labels["author_main_cluster"].map(mapping)
    type_index = {cell_type: index for index, cell_type in enumerate(cell_types)}
    labels["cell_type_index"] = labels["cell_type"].map(type_index).astype(str)
    support = Counter(labels["cell_type"])
    failures = {
        cell_type: support[cell_type]
        for cell_type in cell_types
        if support[cell_type] < minimum
    }
    if failures:
        raise ValueError(f"Embryo broad-label support gate failed: {failures}")

    BROAD_LABEL_ROOT.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{BROAD_LABEL_VERSION}.", dir=BROAD_LABEL_ROOT.parent
        )
    )
    try:
        labels.to_csv(
            temporary / "cells.tsv.gz", sep="\t", index=False, compression="gzip"
        )
        with (temporary / "support.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("cell_type", "cells"))
            for cell_type in cell_types:
                writer.writerow((cell_type, support[cell_type]))
        atomic_yaml(
            temporary / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "completed_at": datetime.now(timezone.utc).isoformat(),
                "label_version": BROAD_LABEL_VERSION,
                "cell_types": list(cell_types),
                "minimum_cells_per_type": minimum,
                "retained_cells": len(labels),
                "support": {
                    cell_type: support[cell_type] for cell_type in cell_types
                },
                "mapping": mapping,
                "outcome_data_used": False,
                "inputs": {
                    "author_labels": repository_path(author_path),
                    "author_labels_sha256": file_digest(author_path),
                    "author_manifest": repository_path(author_manifest),
                },
                "outputs": {"cells": "cells.tsv.gz", "support": "support.tsv"},
            },
        )
        temporary.rename(BROAD_LABEL_ROOT)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(
        f"gse216371 broad_labels cells={len(labels)} types={len(cell_types)} "
        "status=completed",
        flush=True,
    )
    return labels


def locked_record(role: str) -> Mapping[str, Any]:
    records = [
        value
        for value in read_yaml(LOCK_PATH).get("files", ())
        if isinstance(value, Mapping) and value.get("role") == role
    ]
    if len(records) != 1:
        raise ValueError(f"Expected one locked source with role {role}")
    return records[0]


def expected_fragment_members() -> tuple[str, ...]:
    record = locked_record("geo_family_metadata")
    path = ROOT / str(record["path"])
    if (
        not path.is_file()
        or path.stat().st_size != int(record["bytes"])
        or file_digest(path) != str(record["sha256"])
    ):
        raise ValueError("Locked GSE216371 GEO metadata changed or is absent")
    members: list[str] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            prefix = "!Sample_supplementary_file_1 = "
            if not line.startswith(prefix):
                continue
            name = line[len(prefix) :].strip().rsplit("/", 1)[-1]
            if not name.endswith(".bed.gz"):
                raise ValueError(f"Unexpected GSE216371 supplementary member: {name}")
            members.append(name)
    if len(members) != 68 or len(members) != len(set(members)):
        raise ValueError(f"Expected 68 unique GEO fragment members; found {len(members)}")
    return tuple(members)


def archive_path() -> tuple[Path, Mapping[str, Any]]:
    record = locked_record("complete_author_fragment_bed_suite")
    path = ROOT / str(record["path"])
    if not path.is_file() or path.stat().st_size != int(record["bytes"]):
        raise ValueError("Locked GSE216371 fragment archive changed or is absent")
    if str(record["sha256"]) != "0db7906e74767a650942b0a8b357f408cd6d47f85ca6a914bfb63fc113658a5e":
        raise ValueError("Frozen GSE216371 archive digest declaration changed")
    return path, record


def compile_streamer() -> tuple[Path, dict[str, Any]]:
    compiler = shutil.which("g++")
    if compiler is None:
        raise FileNotFoundError("g++ is required for the GSE216371 fragment stream")
    source_sha256 = file_digest(STREAM_SOURCE)
    binary_root = WORK_ROOT / "bin"
    binary = binary_root / "shapemix_gse216371_stream"
    manifest_path = binary_root / "manifest.yaml"
    if binary.is_file() and manifest_path.is_file():
        manifest = read_yaml(manifest_path)
        if (
            manifest.get("source_sha256") == source_sha256
            and manifest.get("binary_sha256") == file_digest(binary)
        ):
            return binary, manifest
    binary_root.mkdir(parents=True, exist_ok=True)
    temporary = binary_root / ".shapemix_gse216371_stream.tmp"
    temporary.unlink(missing_ok=True)
    command = [
        compiler,
        "-O3",
        "-DNDEBUG",
        "-std=c++20",
        "-Wall",
        "-Wextra",
        "-Werror",
        str(STREAM_SOURCE),
        "-lz",
        "-o",
        str(temporary),
    ]
    subprocess.run(command, check=True, cwd=ROOT)
    os.replace(temporary, binary)
    compiler_version = subprocess.run(
        [compiler, "--version"],
        check=True,
        text=True,
        capture_output=True,
    ).stdout.splitlines()[0]
    manifest = {
        "schema_version": 1,
        "status": "complete",
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "source": repository_path(STREAM_SOURCE),
        "source_sha256": source_sha256,
        "binary": repository_path(binary),
        "binary_sha256": file_digest(binary),
        "compiler": compiler_version,
        "command": [
            value if value not in {str(STREAM_SOURCE), str(temporary)} else Path(value).name
            for value in command
        ],
        "threads": 1,
    }
    atomic_yaml(manifest_path, manifest)
    return binary, manifest


def parse_summary(path: Path) -> dict[str, Any]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any(set(row) != {"key", "value"} for row in rows):
        raise ValueError(f"Invalid streamer summary: {path}")
    result: dict[str, Any] = {}
    for row in rows:
        key, value = str(row["key"]), str(row["value"])
        if key in result:
            raise ValueError(f"Duplicate streamer summary key: {key}")
        result[key] = value if key == "mode" else int(value)
    return result


def write_member_hashes(path: Path, values: Mapping[str, str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("member", "compressed_sha256"))
        for member, digest in values.items():
            writer.writerow((member, digest))


def read_member_hashes(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("member", "compressed_sha256"):
            raise ValueError(f"Invalid member-hash table: {path}")
        result = {str(row["member"]): str(row["compressed_sha256"]) for row in reader}
    if len(result) != 68 or any(len(value) != 64 for value in result.values()):
        raise ValueError("GSE216371 member hashes are incomplete")
    return result


def stream_archive_to_helper(
    command: list[str],
    *,
    source: Path,
    expected_members: tuple[str, ...],
) -> dict[str, str]:
    process = subprocess.Popen(command, stdin=subprocess.PIPE, cwd=ROOT)
    if process.stdin is None:
        raise RuntimeError("Failed to open streamer stdin")
    member_hashes: dict[str, str] = {}
    error: BaseException | None = None
    try:
        with tarfile.open(source, "r:") as archive:
            members = archive.getmembers()
            names = tuple(member.name for member in members)
            if (
                len(members) != 68
                or set(names) != set(expected_members)
                or len(names) != len(set(names))
            ):
                raise ValueError("Tar member set does not exactly match GEO declarations")
            for member in members:
                pure = PurePosixPath(member.name)
                if (
                    not member.isfile()
                    or pure.is_absolute()
                    or len(pure.parts) != 1
                    or ".." in pure.parts
                ):
                    raise ValueError(f"Unsafe GSE216371 tar member: {member.name!r}")
                extracted = archive.extractfile(member)
                if extracted is None:
                    raise ValueError(f"Cannot stream tar member: {member.name}")
                marker = gzip.compress(
                    f"#shapemix_member\t{member.name}\n".encode(), mtime=0
                )
                process.stdin.write(marker)
                digest = hashlib.sha256()
                for chunk in iter(lambda: extracted.read(8 * 1024 * 1024), b""):
                    digest.update(chunk)
                    process.stdin.write(chunk)
                member_hashes[member.name] = digest.hexdigest()
        process.stdin.close()
        return_code = process.wait()
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, command)
    except BaseException as exc:
        error = exc
        try:
            process.stdin.close()
        except (BrokenPipeError, OSError):
            pass
        process.wait()
    if error is not None:
        raise error
    return member_hashes


def validate_fragment_concordance(
    labels: pd.DataFrame, totals_path: Path
) -> dict[str, Any]:
    totals = pd.read_csv(totals_path, sep="\t", dtype={"cell_id": str})
    if (
        tuple(totals.columns) != ("cell_id", "bed_rows", "read_support_sum")
        or totals["cell_id"].duplicated().any()
        or set(totals["cell_id"]) != set(labels["cell_id"])
    ):
        raise ValueError("Streamer cell-total axis does not match frozen E13.5 labels")
    expected = labels.set_index("cell_id")["fragments"].astype(np.int64)
    aligned = totals.set_index("cell_id").loc[expected.index]
    if (aligned[["bed_rows", "read_support_sum"]] < 0).any().any():
        raise ValueError("Streamer emitted negative fragment totals")
    row_match = bool(np.array_equal(aligned["bed_rows"].to_numpy(), expected.to_numpy()))
    support_match = bool(
        np.array_equal(aligned["read_support_sum"].to_numpy(), expected.to_numpy())
    )
    if int(row_match) + int(support_match) != 1:
        raise ValueError(
            "Exactly one deposited fragment-total convention must match every "
            f"retained E13.5 cell; bed_rows={row_match} read_support_sum={support_match}"
        )
    if (aligned["bed_rows"] <= 0).any():
        raise ValueError("A retained E13.5 cell has no deposited fragment row")
    return {
        "passed": True,
        "matching_convention": "bed_rows" if row_match else "read_support_sum",
        "bed_rows_match_every_cell": row_match,
        "read_support_sum_matches_every_cell": support_match,
        "cells_compared": len(expected),
        "total_bed_rows": int(aligned["bed_rows"].sum()),
        "total_read_support": int(aligned["read_support_sum"].sum()),
        "total_workbook_fragments": int(expected.sum()),
    }


def write_coordinate_audit(
    cache_manifest: Mapping[str, Any], concordance: Mapping[str, Any]
) -> None:
    atomic_yaml(
        AUDIT_PATH,
        {
            "schema_version": 1,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "passed": True,
            "genome_build": "mm10",
            "coordinate_system": "zero_based_half_open",
            "left_cut_offset": 0,
            "right_cut_offset": 0,
            "fragment_length": "end_minus_start",
            "read_support_policy": "ignore_for_model_counts",
            "semantic_match": "exact",
            "fragment_total_match": "exact",
            "concordance": dict(concordance),
            "source_archive": cache_manifest["source_archive"],
            "source_archive_sha256": cache_manifest["source_archive_sha256"],
            "source_members": cache_manifest["source_members"],
            "normalized_fragment_manifest": repository_path(
                FRAGMENT_CACHE_ROOT / "manifest.yaml"
            ),
            "outcome_data_used": False,
        },
    )


def build_fragment_statistics() -> Path:
    manifest_path = FRAGMENT_CACHE_ROOT / "manifest.yaml"
    fragments_path = FRAGMENT_CACHE_ROOT / "fragments.tsv.gz"
    if manifest_path.is_file() and fragments_path.is_file():
        manifest = read_yaml(manifest_path)
        if not manifest.get("concordance", {}).get("passed", False):
            raise ValueError("Stored embryo fragment concordance did not pass")
        write_coordinate_audit(manifest, manifest["concordance"])
        print("gse216371 fragment_statistics status=reused", flush=True)
        return fragments_path
    if FRAGMENT_CACHE_ROOT.exists():
        raise FileExistsError(
            f"Partial immutable embryo fragment cache: {FRAGMENT_CACHE_ROOT}"
        )
    labels = build_broad_labels()
    ccre_path = CCRE_ROOT / "candidate_ccres.tsv.gz"
    ccre_manifest = CCRE_ROOT / "manifest.yaml"
    if not ccre_path.is_file() or not ccre_manifest.is_file():
        raise FileNotFoundError("Run the GSE216371 cCRE audit first")
    if not CHROM_SIZES.is_file():
        raise FileNotFoundError(
            "Acquire the checksum-pinned UCSC mm10 chromosome sizes first"
        )
    binary, binary_manifest = compile_streamer()
    source, source_record = archive_path()
    expected_members = expected_fragment_members()

    FRAGMENT_CACHE_ROOT.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{BROAD_LABEL_VERSION}.", dir=FRAGMENT_CACHE_ROOT.parent
        )
    )
    try:
        summary_path = temporary / "stream_summary.tsv"
        statistics_path = temporary / "feature_statistics.bin"
        totals_path = temporary / "cell_totals.tsv"
        command = [
            str(binary),
            "--mode",
            "statistics",
            "--labels",
            str(BROAD_LABEL_ROOT / "cells.tsv.gz"),
            "--peaks",
            str(ccre_path),
            "--chrom-sizes",
            str(CHROM_SIZES),
            "--output-fragments",
            str(temporary / "fragments.tsv.gz"),
            "--output-statistics",
            str(statistics_path),
            "--output-cell-totals",
            str(totals_path),
            "--output-summary",
            str(summary_path),
        ]
        member_hashes = stream_archive_to_helper(
            command, source=source, expected_members=expected_members
        )
        if tuple(member_hashes) != expected_members:
            raise ValueError("Streamer did not preserve the GEO tar-member order")
        write_member_hashes(temporary / "source_member_sha256.tsv", member_hashes)
        summary = parse_summary(summary_path)
        if (
            summary.get("mode") != "statistics"
            or summary.get("cells") != len(labels)
            or summary.get("features")
            != int(read_yaml(ccre_manifest)["candidate_intervals"])
            or summary.get("valid_rows") != summary.get("total_rows")
        ):
            raise ValueError(f"Invalid embryo statistics summary: {summary}")
        concordance = validate_fragment_concordance(labels, totals_path)
        manifest = {
            "schema_version": 1,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "cache_version": BROAD_LABEL_VERSION,
            "source_archive": repository_path(source),
            "source_archive_bytes": int(source_record["bytes"]),
            "source_archive_sha256": str(source_record["sha256"]),
            "source_members": len(member_hashes),
            "source_member_sha256_table": "source_member_sha256.tsv",
            "labels": repository_path(BROAD_LABEL_ROOT / "cells.tsv.gz"),
            "candidate_ccres": repository_path(ccre_path),
            "candidate_ccres_sha256": file_digest(ccre_path),
            "chrom_sizes": repository_path(CHROM_SIZES),
            "chrom_sizes_sha256": file_digest(CHROM_SIZES),
            "streamer": binary_manifest,
            "command": [
                (
                    Path(value).name
                    if value == str(binary) or value.startswith(f"{temporary}{os.sep}")
                    else repository_path(Path(value))
                    if value.startswith(f"{ROOT}{os.sep}")
                    else value
                )
                for value in command
            ],
            "threads": 1,
            "counters": summary,
            "concordance": concordance,
            "coordinate_semantics": {
                "genome_build": "mm10",
                "source": "zero_based_half_open_parent_fragments",
                "left_cut_offset": 0,
                "right_cut_offset": 0,
                "fragment_length": "end_minus_start",
                "read_support_policy": "ignore",
            },
            "outputs": {
                "fragments": "fragments.tsv.gz",
                "fragments_bytes": (temporary / "fragments.tsv.gz").stat().st_size,
                "fragments_sha256": file_digest(temporary / "fragments.tsv.gz"),
                "feature_statistics": "feature_statistics.bin",
                "cell_totals": "cell_totals.tsv",
                "stream_summary": "stream_summary.tsv",
            },
        }
        atomic_yaml(temporary / "manifest.yaml", manifest)
        temporary.rename(FRAGMENT_CACHE_ROOT)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    write_coordinate_audit(manifest, concordance)
    print(
        f"gse216371 fragment_statistics rows={summary['total_rows']} "
        f"retained={summary['retained_fragments']} status=completed",
        flush=True,
    )
    return fragments_path


def load_statistics(path: Path) -> tuple[np.memmap, np.memmap]:
    with path.open("rb") as handle:
        if handle.read(8) != b"SM216C01":
            raise ValueError("Invalid GSE216371 feature-statistics magic")
        type_count, feature_count = struct.unpack("<QQ", handle.read(16))
    cell_types, _, _ = ontology()
    candidate_count = int(
        read_yaml(CCRE_ROOT / "manifest.yaml")["candidate_intervals"]
    )
    if type_count != len(cell_types) or feature_count != candidate_count:
        raise ValueError("Feature-statistics axes do not match frozen inputs")
    counts = np.memmap(
        path,
        mode="r",
        dtype="<u8",
        offset=24,
        shape=(type_count, feature_count),
    )
    coverage = np.memmap(
        path,
        mode="r",
        dtype="u1",
        offset=24 + type_count * feature_count * 8,
        shape=(feature_count,),
    )
    expected_bytes = 24 + type_count * feature_count * 8 + feature_count
    if path.stat().st_size != expected_bytes:
        raise ValueError("Feature-statistics byte count does not match its header")
    return counts, coverage


def load_candidate_rows(indices: Iterable[int]) -> dict[int, dict[str, str]]:
    requested = {int(value) for value in indices}
    result: dict[int, dict[str, str]] = {}
    with gzip.open(CCRE_ROOT / "candidate_ccres.tsv.gz", "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            index = int(row["source_index"])
            if index in requested:
                result[index] = dict(row)
    missing = sorted(requested.difference(result))
    if missing:
        raise ValueError(f"Candidate cCRE rows are missing: {missing[:10]}")
    return result


def rank_features(
    score: np.ndarray,
    coverage: np.ndarray,
    total: np.ndarray,
    n_top: int = N_TOP_PEAKS,
) -> list[int]:
    score = np.asarray(score, dtype=np.float64)
    coverage = np.asarray(coverage, dtype=np.uint8)
    total = np.asarray(total, dtype=np.uint64)
    if score.ndim != 1 or score.shape != coverage.shape or score.shape != total.shape:
        raise ValueError("Feature ranker arrays must be aligned one-dimensional vectors")
    eligible = np.flatnonzero(
        np.isfinite(score) & (coverage >= MIN_REFERENCE_CELLS_PER_PEAK)
    )
    if len(eligible) < n_top:
        raise ValueError(f"Only {len(eligible)} embryo cCREs pass the support gate")
    if np.any(total > np.iinfo(np.int64).max):
        raise ValueError("A feature total exceeds the deterministic signed rank range")
    signed_total = total.astype(np.int64)
    numeric_order = np.lexsort(
        (eligible, -signed_total[eligible], -score[eligible])
    )
    boundary = int(eligible[numeric_order[n_top - 1]])
    strict = eligible[
        (score[eligible] > score[boundary])
        | (
            (score[eligible] == score[boundary])
            & (total[eligible] > total[boundary])
        )
    ]
    tied = eligible[
        (score[eligible] == score[boundary])
        & (total[eligible] == total[boundary])
    ]
    remaining = n_top - len(strict)
    if remaining < 1 or remaining > len(tied):
        raise RuntimeError(
            "The deterministic embryo feature-rank boundary is inconsistent"
        )
    tied_rows = load_candidate_rows(tied)
    chosen = sorted(tied.tolist(), key=lambda index: tied_rows[index]["peak_id"])[
        :remaining
    ]
    selected = np.concatenate(
        (strict.astype(np.int64), np.asarray(chosen, dtype=np.int64))
    )
    selected_rows = load_candidate_rows(selected)
    if len(selected) != n_top or len(set(selected.tolist())) != n_top:
        raise RuntimeError("The embryo feature ranker did not select exactly n_top")
    return sorted(
        selected.tolist(),
        key=lambda index: (
            -score[index],
            -int(total[index]),
            selected_rows[index]["peak_id"],
        ),
    )


def build_feature_axis() -> pd.DataFrame:
    selected_path = AXIS_ROOT / "selected_500bp_intervals.tsv.gz"
    manifest_path = AXIS_ROOT / "manifest.yaml"
    if selected_path.is_file() and manifest_path.is_file():
        return pd.read_csv(selected_path, sep="\t")
    if AXIS_ROOT.exists():
        raise FileExistsError(f"Partial immutable embryo feature axis: {AXIS_ROOT}")
    build_fragment_statistics()
    counts, coverage = load_statistics(
        FRAGMENT_CACHE_ROOT / "feature_statistics.bin"
    )
    type_totals = np.asarray(counts.sum(axis=1, dtype=np.float64))
    if np.any(type_totals <= 0) or not np.isfinite(type_totals).all():
        raise ValueError("An embryo broad class has no candidate-cCRE cut counts")
    score = np.empty(counts.shape[1], dtype=np.float64)
    for start in range(0, counts.shape[1], FEATURE_CHUNK):
        stop = min(start + FEATURE_CHUNK, counts.shape[1])
        normalized = np.log2(
            1.0
            + 1.0e4
            * counts[:, start:stop].astype(np.float64)
            / type_totals[:, None]
        )
        score[start:stop] = np.var(normalized, axis=0, ddof=0)
    total = np.asarray(counts.sum(axis=0, dtype=np.uint64)).ravel()
    selected_indices = rank_features(score, coverage, total)
    source_rows = load_candidate_rows(selected_indices)
    records = []
    for rank, index in enumerate(selected_indices, start=1):
        row = source_rows[index]
        records.append(
            {
                "rank": rank,
                "source_feature_index": index,
                "peak_id": row["peak_id"],
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "score": float(score[index]),
                "nonzero_reference_cells_min_cap": int(coverage[index]),
                "total_reference_counts": int(total[index]),
                "ccre_type": row["ccre_type"],
                "nearest_tss": row["nearest_tss"],
            }
        )
    selected = pd.DataFrame.from_records(records)
    cell_types, _, _ = ontology()
    AXIS_ROOT.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{BROAD_LABEL_VERSION}.", dir=AXIS_ROOT.parent)
    )
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
                "completed_at": datetime.now(timezone.utc).isoformat(),
                "axis_version": BROAD_LABEL_VERSION,
                "candidate_intervals": counts.shape[1],
                "eligible_intervals": int(
                    np.count_nonzero(coverage >= MIN_REFERENCE_CELLS_PER_PEAK)
                ),
                "selected_intervals": N_TOP_PEAKS,
                "feature_sha256": ordered_feature_sha256(
                    selected["peak_id"].astype(str)
                ),
                "label_ontology": list(cell_types),
                "type_candidate_cut_totals": {
                    cell_type: int(value)
                    for cell_type, value in zip(
                        cell_types, type_totals.astype(np.uint64), strict=True
                    )
                },
                "selector": {
                    "minimum_nonzero_reference_cells": MIN_REFERENCE_CELLS_PER_PEAK,
                    "coverage_storage": "exact_until_threshold_then_capped",
                    "score": "population_variance_log2_1_plus_scaled_type_rate",
                    "scale": 1.0e4,
                    "tie_breaks": [
                        "score_desc",
                        "total_count_desc",
                        "peak_id_asc",
                    ],
                    "outcome_data_used": False,
                },
                "inputs": {
                    "labels": repository_path(BROAD_LABEL_ROOT / "cells.tsv.gz"),
                    "candidate_ccres": repository_path(
                        CCRE_ROOT / "candidate_ccres.tsv.gz"
                    ),
                    "feature_statistics": repository_path(
                        FRAGMENT_CACHE_ROOT / "feature_statistics.bin"
                    ),
                },
                "software_versions": software_versions(),
            },
        )
        temporary.rename(AXIS_ROOT)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print("gse216371 feature_axis peaks=5000 status=completed", flush=True)
    return selected


def push_sparse_run(
    levels: list[sparse.csr_matrix | None], run: sparse.csr_matrix
) -> None:
    level = 0
    while True:
        if level == len(levels):
            levels.append(run)
            return
        existing = levels[level]
        if existing is None:
            levels[level] = run
            return
        levels[level] = None
        run = (existing + run).tocsr()
        run.sum_duplicates()
        run.eliminate_zeros()
        run.sort_indices()
        level += 1


def event_layers(path: Path, shape: tuple[int, int]) -> dict[str, sparse.csr_matrix]:
    if path.stat().st_size % EVENT_DTYPE.itemsize:
        raise ValueError("Packed embryo shape-event file has a partial record")
    events = np.memmap(path, mode="r", dtype=EVENT_DTYPE)
    levels: list[list[sparse.csr_matrix | None]] = [[], [], []]
    for start in range(0, len(events), EVENT_CHUNK):
        chunk = events[start : start + EVENT_CHUNK]
        if (
            len(chunk)
            and (
                int(chunk["cell"].max()) >= shape[0]
                or int(chunk["feature"].max()) >= shape[1]
                or int(chunk["layer"].max()) >= 3
            )
        ):
            raise ValueError("Packed embryo shape event is outside its declared axes")
        for layer in range(3):
            selected = chunk["layer"] == layer
            if not np.any(selected):
                continue
            run = sparse.coo_matrix(
                (
                    np.ones(int(np.count_nonzero(selected)), dtype=np.int64),
                    (chunk["cell"][selected], chunk["feature"][selected]),
                ),
                shape=shape,
                dtype=np.int64,
            ).tocsr()
            run.sum_duplicates()
            run.sort_indices()
            push_sparse_run(levels[layer], run)
    result = {}
    for bin_, layer_levels in zip(DEFAULT_FRAGMENT_LENGTH_BINS, levels, strict=True):
        matrix = sparse.csr_matrix(shape, dtype=np.int64)
        for run in layer_levels:
            if run is not None:
                matrix = (matrix + run).tocsr()
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        matrix.sort_indices()
        result[bin_.layer] = matrix
    return result


def reference_shape_result(
    labels: pd.DataFrame,
    selected: pd.DataFrame,
    layers: Mapping[str, sparse.csr_matrix],
    summary: Mapping[str, Any],
) -> FragmentShapeResult:
    qc = FragmentShapeQC(
        total_rows=int(summary["total_rows"]),
        header_rows=int(summary["header_rows"]),
        invalid_schema_rows=0,
        invalid_coordinate_rows=0,
        unknown_barcodes=int(summary["unknown_barcodes"]),
        filtered_contigs=0,
        valid_rows=int(summary["valid_rows"]),
        retained_fragments=int(summary["retained_fragments"]),
        fragments_with_assigned_cut_sites=int(
            summary["fragments_with_assigned_cut_sites"]
        ),
        cut_sites_outside_peaks=int(summary["cut_sites_outside_peaks"]),
        assigned_cut_sites=int(summary["assigned_cut_sites"]),
        read_support_total=int(summary["read_support_total"]),
        cut_sites_per_bin={
            bin_.layer: int(summary[f"cut_sites_per_bin.{bin_.layer}"])
            for bin_ in DEFAULT_FRAGMENT_LENGTH_BINS
        },
    )
    peaks = tuple(
        PeakInterval(str(chrom), int(start), int(end), str(peak_id))
        for chrom, start, end, peak_id in selected[
            ["chrom", "start", "end", "peak_id"]
        ].itertuples(index=False, name=None)
    )
    return FragmentShapeResult(
        barcodes=tuple(labels["cell_id"].astype(str)),
        peaks=peaks,
        bins=DEFAULT_FRAGMENT_LENGTH_BINS,
        layers=dict(layers),
        qc=qc,
        right_cut_offset=0,
    )


def build_reference() -> Path:
    reference_path = REFERENCE_ROOT / "atac/reference.h5ad"
    manifest_path = REFERENCE_ROOT / "reference.yaml"
    if reference_path.is_file() and manifest_path.is_file():
        print("gse216371 reference status=reused", flush=True)
        return reference_path
    if REFERENCE_ROOT.exists():
        raise FileExistsError(f"Partial immutable embryo reference: {REFERENCE_ROOT}")
    labels = build_broad_labels()
    build_fragment_statistics()
    selected = build_feature_axis()
    binary, binary_manifest = compile_streamer()
    event_root = WORK_ROOT / "shape_events" / BROAD_LABEL_VERSION
    event_root.mkdir(parents=True, exist_ok=True)
    event_path = event_root / "events.bin"
    summary_path = event_root / "summary.tsv"
    command = [
        str(binary),
        "--mode",
        "shape",
        "--labels",
        str(BROAD_LABEL_ROOT / "cells.tsv.gz"),
        "--peaks",
        str(AXIS_ROOT / "selected_500bp_intervals.tsv.gz"),
        "--chrom-sizes",
        str(CHROM_SIZES),
        "--output-events",
        str(event_path),
        "--output-summary",
        str(summary_path),
    ]
    with (FRAGMENT_CACHE_ROOT / "fragments.tsv.gz").open("rb") as input_handle:
        subprocess.run(command, stdin=input_handle, check=True, cwd=ROOT)
    summary = parse_summary(summary_path)
    if (
        summary.get("mode") != "shape"
        or summary.get("cells") != len(labels)
        or summary.get("features") != N_TOP_PEAKS
        or summary.get("unknown_barcodes") != 0
        or summary.get("valid_rows") != summary.get("retained_fragments")
        or event_path.stat().st_size // EVENT_DTYPE.itemsize
        != summary.get("event_records")
    ):
        raise ValueError(f"Invalid embryo shape-event summary: {summary}")
    layers = event_layers(event_path, (len(labels), len(selected)))
    result = reference_shape_result(labels, selected, layers, summary)
    cache_manifest = read_yaml(FRAGMENT_CACHE_ROOT / "manifest.yaml")
    source_hashes = read_member_hashes(
        FRAGMENT_CACHE_ROOT / str(cache_manifest["source_member_sha256_table"])
    )
    obs = labels.set_index(labels["cell_id"].astype(str)).copy()
    obs.index.name = "cell_id"
    var = selected.set_index("peak_id")[
        ["chrom", "start", "end", "ccre_type", "nearest_tss"]
    ].copy()
    reference = build_fragment_shape_anndata(
        result,
        obs=obs,
        var=var,
        provenance={
            "split_sha256": file_digest(BROAD_LABEL_ROOT / "cells.tsv.gz"),
            "source_sha256": source_hashes,
            "coordinate_validation": {
                "selected_right_cut_offset": 0,
                "matrix_match": "not_available",
                "validation_method": "deposited_bed_half_open_parent_fragments",
                "semantic_match": "exact",
                "fragment_total_match": "exact",
                "audit": repository_path(AUDIT_PATH),
            },
            "software_versions": software_versions(),
        },
    )
    validate_fragment_shape_spec(
        FragmentShapeSpec.from_mapping(reference.uns["fragment_shape"])
    )
    validate_fragment_shape_feature_axis(reference, "GSE216371 reference")
    cell_types, _, _ = ontology()
    if set(reference.obs["cell_type"].astype(str)) != set(cell_types):
        raise ValueError("Assembled embryo reference label universe changed")

    REFERENCE_ROOT.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{CONFIG['standardized_reference_id']}.",
            dir=REFERENCE_ROOT.parent,
        )
    )
    try:
        (temporary / "atac").mkdir()
        reference.write_h5ad(temporary / "atac/reference.h5ad", compression="gzip")
        atomic_yaml(
            temporary / "reference.yaml",
            {
                "schema_version": 1,
                "reference_id": CONFIG["standardized_reference_id"],
                "source_dataset_id": CONFIG["source_dataset_id"],
                "description": (
                    "Whole-embryo E13.5 SPATAC reference on a reference-only "
                    "5,000-author-cCRE mm10 axis."
                ),
                "labels_key": "cell_type",
                "genome_build": "mm10",
                "cell_types": list(cell_types),
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
                    "source_lock": repository_path(LOCK_PATH),
                    "label_manifest": repository_path(
                        BROAD_LABEL_ROOT / "manifest.yaml"
                    ),
                    "feature_manifest": repository_path(AXIS_ROOT / "manifest.yaml"),
                    "coordinate_audit": repository_path(AUDIT_PATH),
                    "fragment_manifest": repository_path(
                        FRAGMENT_CACHE_ROOT / "manifest.yaml"
                    ),
                    "streamer_manifest": binary_manifest,
                    "shape_event_summary": summary,
                    "shape_event_sha256": file_digest(event_path),
                },
            },
        )
        temporary.rename(REFERENCE_ROOT)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(
        f"gse216371 reference cells={reference.n_obs} peaks={reference.n_vars} "
        "status=completed",
        flush=True,
    )
    return reference_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("labels", "fragment-statistics", "feature-axis", "reference", "all"),
        default="all",
    )
    return parser.parse_args()


def main() -> None:
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run this materializer through run_shapemix_low_impact.sh")
    args = parse_args()
    if args.stage in {"labels", "all"}:
        build_broad_labels()
    if args.stage in {"fragment-statistics", "all"}:
        build_fragment_statistics()
    if args.stage in {"feature-axis", "all"}:
        build_feature_axis()
    if args.stage in {"reference", "all"}:
        build_reference()


if __name__ == "__main__":
    main()
