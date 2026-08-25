#!/usr/bin/env python
"""Preprocess GSE205055 or GSE263333 spatial deposits for ShapeMix."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import tarfile
import tempfile
import urllib.parse
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Mapping, Optional

import anndata as ad
import numpy as np
import pandas as pd
import pysam
import yaml
from scipy import io as scipy_io
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
CHUNK_BYTES = 8 * 1024 * 1024
GSM_PATTERN = re.compile(r"^(GSM\d+)_")
POSITION_COLUMNS = (
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pixel_row_fullres",
    "pixel_col_fullres",
)


def load_yaml(path: Path) -> dict[str, Any]:
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


def project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_member_name(name: str) -> bool:
    path = PurePosixPath(name)
    return bool(name) and not path.is_absolute() and ".." not in path.parts


def accession_from_filename(path: Path) -> str:
    match = GSM_PATTERN.match(path.name)
    if match is None:
        raise ValueError(f"Supplementary filename lacks a GSM prefix: {path.name}")
    return match.group(1)


def all_declared(config: Mapping[str, Any], key: str) -> set[str]:
    return {
        str(gsm)
        for group in config["sample_groups"]
        for gsm in group.get(key, [])
    }


def classify_payload(path: Path, config: Mapping[str, Any]) -> str:
    name = path.name
    lower = name.lower()
    gsm = accession_from_filename(path)
    if lower.endswith("_spatial.tar.gz") or "_spatial_" in lower:
        return "spatial"
    if lower.endswith(".fragments.tsv.gz") or lower.endswith("_fragments.tsv.gz"):
        if "h3k27ac" in lower:
            return "h3k27ac_fragments"
        if "h3k27me3" in lower:
            return "h3k27me3_fragments"
        if "h3k4me3" in lower:
            return "h3k4me3_fragments"
        if "_atac." in lower or gsm in all_declared(config, "atac"):
            return "atac_fragments"
        return "other_epigenome_fragments"
    if lower.endswith(".tsv.gz") and "matrix" in lower:
        for mark in ("h3k27ac", "h3k27me3", "h3k4me3"):
            if mark in lower:
                return f"{mark}_dense_matrix"
        if gsm in all_declared(config, "rna"):
            return "rna_dense_matrix"
    if lower.endswith("_rna.tar.gz"):
        return "rna_matrix_bundle"
    if lower.endswith("_protein.tar.gz"):
        return "protein_matrix_bundle"
    return "unclassified"


def source_archive(config: Mapping[str, Any]) -> Path:
    return (
        project_path(str(config["raw_directory"]))
        / "processed_downloads"
        / str(config["archive"]["filename"])
    )


def extracted_root(config: Mapping[str, Any]) -> Path:
    return project_path(str(config["processed_directory"])) / "extracted_payload"


def iter_payloads(config: Mapping[str, Any]) -> tuple[Path, ...]:
    root = extracted_root(config)
    return tuple(sorted(path for path in root.glob("GSM*/*") if path.is_file()))


def expected_supplementary_basenames(config: Mapping[str, Any]) -> set[str]:
    metadata_root = project_path(str(config["raw_directory"])) / "series_metadata"
    names: set[str] = set()
    for record in config["series_metadata"]:
        accession = str(record["accession"])
        path = metadata_root / f"{accession}_family.soft.gz"
        if not path.is_file():
            raise FileNotFoundError(path)
        with gzip.open(path, "rt", errors="replace") as handle:
            for line in handle:
                if not line.startswith("!Sample_supplementary_file"):
                    continue
                _, separator, value = line.rstrip("\r\n").partition(" = ")
                if not separator or not value or value.upper() == "NONE":
                    continue
                url_path = urllib.parse.urlparse(value).path
                name = urllib.parse.unquote(PurePosixPath(url_path).name)
                if not name:
                    raise ValueError(f"Malformed supplementary-file URL in {path}: {value!r}")
                names.add(name)
    if not names:
        raise ValueError(f"No supplementary files declared by metadata for {config['accession']}")
    return names


def extract_source_archive(config: Mapping[str, Any], overwrite: bool = False) -> None:
    source = source_archive(config)
    if not source.is_file():
        raise FileNotFoundError(source)
    target_root = extracted_root(config)
    target_root.mkdir(parents=True, exist_ok=True)
    records = []
    with tarfile.open(source, "r:") as archive:
        for member in archive:
            if member.isdir():
                continue
            if not safe_member_name(member.name) or not member.isfile():
                raise ValueError(f"Unsafe or unsupported archive member: {member.name!r}")
            name = PurePosixPath(member.name).name
            match = GSM_PATTERN.match(name)
            if match is None:
                raise ValueError(f"Archive member lacks a GSM prefix: {member.name!r}")
            destination = target_root / match.group(1) / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists() and not overwrite:
                if destination.stat().st_size != member.size:
                    raise ValueError(f"Existing extracted member has wrong size: {destination}")
            else:
                stream = archive.extractfile(member)
                if stream is None:
                    raise ValueError(f"Could not read archive member: {member.name!r}")
                temporary = destination.with_name(f".{destination.name}.tmp")
                temporary.unlink(missing_ok=True)
                try:
                    with temporary.open("wb") as output:
                        shutil.copyfileobj(stream, output, length=CHUNK_BYTES)
                    if temporary.stat().st_size != member.size:
                        raise IOError(f"Extracted member has wrong size: {member.name!r}")
                    os.replace(temporary, destination)
                except BaseException:
                    temporary.unlink(missing_ok=True)
                    raise
            records.append(
                {
                    "gsm": match.group(1),
                    "filename": name,
                    "tar_member": member.name,
                    "bytes": member.size,
                    "path": str(destination.relative_to(ROOT)),
                    "assay_role": classify_payload(destination, config),
                }
            )
    if not records:
        raise ValueError(f"No files were extracted from {source}")
    expected = expected_supplementary_basenames(config)
    observed = {str(record["filename"]) for record in records}
    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        raise ValueError(
            f"Archive does not match GEO family metadata; missing={missing}, extra={extra}"
        )
    frame = pd.DataFrame.from_records(records).sort_values(["gsm", "filename"])
    audit_root = project_path(str(config["processed_directory"])) / "source_audit"
    output = audit_root / "payload_inventory.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, sep="\t", index=False)
    atomic_yaml(
        audit_root / "extraction.yaml",
        {
            "schema_version": 1,
            "accession": config["accession"],
            "source_archive": str(source.relative_to(ROOT)),
            "source_archive_bytes": source.stat().st_size,
            "expected_supplementary_files": len(expected),
            "observed_archive_files": len(observed),
            "exact_metadata_match": True,
            "filenames": sorted(observed),
        },
    )
    print(f"extracted {len(records)} payloads to {target_root}", flush=True)


def validate_nested_payload(path: Path) -> tuple[str, Optional[int]]:
    if path.name.endswith(".tar.gz"):
        members = 0
        with tarfile.open(path, "r:gz") as archive:
            for member in archive:
                if not safe_member_name(member.name):
                    raise ValueError(f"Unsafe nested archive member: {member.name!r}")
                if member.issym() or member.islnk() or member.isdev():
                    raise ValueError(f"Unsupported nested archive member: {member.name!r}")
                members += 1
        if members == 0:
            raise ValueError(f"Nested archive is empty: {path}")
        return "passed_nested_tar_gzip_stream", members
    if path.name.endswith(".gz"):
        with gzip.open(path, "rb") as handle:
            for _ in iter(lambda: handle.read(CHUNK_BYTES), b""):
                pass
        return "passed_full_gzip_stream", None
    raise ValueError(f"Unexpected uncompressed supplementary payload: {path}")


def audit_payloads(config: Mapping[str, Any]) -> None:
    records = []
    for path in iter_payloads(config):
        role = classify_payload(path, config)
        if role == "unclassified":
            raise ValueError(f"Unclassified supplementary payload: {path}")
        integrity, members = validate_nested_payload(path)
        records.append(
            {
                "gsm": accession_from_filename(path),
                "filename": path.name,
                "role": role,
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "integrity": integrity,
                "nested_members": members,
            }
        )
        print(f"audited {path.name}", flush=True)
    output = project_path(str(config["processed_directory"])) / "source_audit/payloads.yaml"
    atomic_yaml(
        output,
        {
            "schema_version": 1,
            "accession": config["accession"],
            "payloads": len(records),
            "files": records,
        },
    )


def fragment_assay(role: str) -> str:
    if role == "atac_fragments":
        return "atac"
    return role.removesuffix("_fragments")


def normalize_one_fragment(path: Path, role: str, config: Mapping[str, Any], overwrite: bool) -> dict[str, Any]:
    gsm = accession_from_filename(path)
    assay = fragment_assay(role)
    family_root = project_path(str(config["processed_directory"]))
    if assay == "atac":
        output_root = family_root / "normalized_atac_fragments" / gsm
    else:
        output_root = family_root / "validation_modalities/epigenome" / gsm / assay
    output_root.mkdir(parents=True, exist_ok=True)
    destination = output_root / path.name
    temporary = destination.with_name(f".{destination.name}.bgzf.tmp")
    if destination.exists() and not overwrite:
        raise FileExistsError(
            f"Normalized fragment exists; rerun without normalize or pass --overwrite: {destination}"
        )
    temporary.unlink(missing_ok=True)
    rows = 0
    invalid_rows = 0
    read_support_total = 0
    source_column_counts: Counter[int] = Counter()
    barcodes: set[str] = set()
    length_bins = {"lt_100": 0, "100_249": 0, "ge_250": 0}
    max_length = 0
    last_chrom: Optional[str] = None
    last_start = -1
    closed_chroms: set[str] = set()
    coordinate_sorted = True
    adjacent_duplicates = 0
    previous: Optional[tuple[str, int, int, str]] = None
    try:
        with gzip.open(path, "rt") as source, pysam.BGZFile(str(temporary), "wb") as output:
            for line_number, line in enumerate(source, start=1):
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\r\n").split("\t")
                source_column_counts[len(fields)] += 1
                if len(fields) not in {4, 5}:
                    invalid_rows += 1
                    raise ValueError(f"{path}: row {line_number} has {len(fields)} columns")
                chrom, start_text, end_text, barcode = fields[:4]
                support_text = fields[4] if len(fields) == 5 else "1"
                try:
                    start, end, support = int(start_text), int(end_text), int(support_text)
                except ValueError as exc:
                    invalid_rows += 1
                    raise ValueError(f"{path}: non-integer fragment field at row {line_number}") from exc
                if not chrom or start < 0 or end <= start or not barcode or support < 1:
                    invalid_rows += 1
                    raise ValueError(f"{path}: invalid fragment at row {line_number}")
                if last_chrom is None:
                    last_chrom = chrom
                elif chrom != last_chrom:
                    closed_chroms.add(last_chrom)
                    if chrom in closed_chroms:
                        coordinate_sorted = False
                    last_chrom = chrom
                    last_start = -1
                if start < last_start:
                    coordinate_sorted = False
                last_start = start
                key = (chrom, start, end, barcode)
                if key == previous:
                    adjacent_duplicates += 1
                previous = key
                length = end - start
                if length < 100:
                    length_bins["lt_100"] += 1
                elif length < 250:
                    length_bins["100_249"] += 1
                else:
                    length_bins["ge_250"] += 1
                rows += 1
                read_support_total += support
                max_length = max(max_length, length)
                barcodes.add(barcode)
                output.write(f"{chrom}\t{start}\t{end}\t{barcode}\t{support}\n".encode())
        if rows == 0:
            raise ValueError(f"No fragment rows in {path}")
        if not coordinate_sorted:
            raise ValueError(f"Source fragment file is not coordinate sorted: {path}")
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    index = Path(f"{destination}.tbi")
    if index.exists() and overwrite:
        index.unlink()
    pysam.tabix_index(str(destination), preset="bed", force=False)
    barcode_path = (
        family_root / "source_audit/barcodes" / f"{gsm}__{path.name}.txt"
    )
    barcode_path.parent.mkdir(parents=True, exist_ok=True)
    barcode_path.write_text("\n".join(sorted(barcodes)) + "\n")
    return {
        "gsm": gsm,
        "source": str(path.relative_to(ROOT)),
        "source_sha256": sha256_file(path),
        "role": role,
        "assay": assay,
        "normalized": str(destination.relative_to(ROOT)),
        "normalized_bytes": destination.stat().st_size,
        "normalized_sha256": sha256_file(destination),
        "tabix": str(index.relative_to(ROOT)),
        "tabix_bytes": index.stat().st_size,
        "tabix_sha256": sha256_file(index),
        "rows": rows,
        "invalid_rows": invalid_rows,
        "source_column_counts": dict(sorted(source_column_counts.items())),
        "implicit_support_one": set(source_column_counts) == {4},
        "unique_barcodes": len(barcodes),
        "read_support_total": read_support_total,
        "max_fragment_length": max_length,
        "length_bins": length_bins,
        "coordinate_sorted": coordinate_sorted,
        "adjacent_duplicate_rows": adjacent_duplicates,
        "barcode_path": str(barcode_path.relative_to(ROOT)),
        "barcode_sha256": sha256_file(barcode_path),
    }


def normalize_fragments(config: Mapping[str, Any], overwrite: bool = False) -> None:
    records = []
    for path in iter_payloads(config):
        role = classify_payload(path, config)
        if not role.endswith("_fragments"):
            continue
        record = normalize_one_fragment(path, role, config, overwrite)
        records.append(record)
        print(f"normalized {path.name}: {record['rows']} rows", flush=True)
    if not records:
        raise ValueError("No fragment payloads were found")
    atomic_yaml(
        project_path(str(config["processed_directory"])) / "source_audit/fragments.yaml",
        {
            "schema_version": 1,
            "accession": config["accession"],
            "fragment_files": len(records),
            "atac_fragment_files": sum(record["assay"] == "atac" for record in records),
            "files": records,
        },
    )


def extract_nested_archive(source: Path, target: Path, overwrite: bool) -> tuple[Path, ...]:
    target.mkdir(parents=True, exist_ok=True)
    outputs = []
    with tarfile.open(source, "r:gz") as archive:
        for member in archive:
            if member.isdir():
                continue
            if not safe_member_name(member.name) or not member.isfile():
                raise ValueError(f"Unsafe or unsupported nested member: {member.name!r}")
            relative = Path(*PurePosixPath(member.name).parts)
            destination = target / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists() and not overwrite:
                if destination.stat().st_size != member.size:
                    raise ValueError(f"Existing nested file has wrong size: {destination}")
            else:
                stream = archive.extractfile(member)
                if stream is None:
                    raise ValueError(f"Could not read nested member: {member.name!r}")
                temporary = destination.with_name(f".{destination.name}.tmp")
                with temporary.open("wb") as output:
                    shutil.copyfileobj(stream, output, length=CHUNK_BYTES)
                os.replace(temporary, destination)
            outputs.append(destination)
    return tuple(outputs)


def parse_positions(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, header=None, names=POSITION_COLUMNS, dtype={"barcode": str})
    if frame.shape[1] != 6 or frame.empty:
        raise ValueError(f"Invalid tissue-position table: {path}")
    if frame["barcode"].isna().any() or frame["barcode"].duplicated().any():
        raise ValueError(f"Missing or duplicated spatial barcodes: {path}")
    numeric = POSITION_COLUMNS[1:]
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    return frame


def preprocess_spatial_assets(config: Mapping[str, Any], overwrite: bool = False) -> None:
    family_root = project_path(str(config["processed_directory"]))
    records = []
    for path in iter_payloads(config):
        if classify_payload(path, config) != "spatial":
            continue
        gsm = accession_from_filename(path)
        target = family_root / "spatial_coordinates" / gsm / "deposited"
        files = extract_nested_archive(path, target, overwrite)
        positions = [candidate for candidate in files if candidate.name == "tissue_positions_list.csv"]
        scales = [candidate for candidate in files if candidate.name == "scalefactors_json.json"]
        if len(positions) != 1 or len(scales) != 1:
            raise ValueError(f"Spatial archive lacks one positions/scales pair: {path}")
        frame = parse_positions(positions[0])
        with scales[0].open() as handle:
            scale_factors = json.load(handle)
        canonical = family_root / "spatial_coordinates" / gsm / "coordinates.csv"
        frame.to_csv(canonical, index=False)
        records.append(
            {
                "gsm": gsm,
                "source": str(path.relative_to(ROOT)),
                "pixels": len(frame),
                "in_tissue_pixels": int((frame["in_tissue"].astype(float) > 0).sum()),
                "canonical_coordinates": str(canonical.relative_to(ROOT)),
                "canonical_coordinates_sha256": sha256_file(canonical),
                "scale_factors": scale_factors,
                "deposited_files": [str(item.relative_to(ROOT)) for item in files],
            }
        )
        print(f"spatial {gsm}: {len(frame)} pixels", flush=True)
    atomic_yaml(
        family_root / "source_audit/spatial.yaml",
        {
            "schema_version": 1,
            "accession": config["accession"],
            "spatial_archives": len(records),
            "files": records,
        },
    )


def unique_names(values: Iterable[str]) -> list[str]:
    counts: Counter[str] = Counter()
    output = []
    for value in values:
        index = counts[value]
        output.append(value if index == 0 else f"{value}-{index}")
        counts[value] += 1
    return output


def read_dense_gene_by_pixel(path: Path) -> ad.AnnData:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if not header or any(not barcode for barcode in header):
            raise ValueError(f"Dense RNA header has missing barcodes: {path}")
        rows: list[int] = []
        columns: list[int] = []
        values: list[int] = []
        genes: list[str] = []
        for row_index, line in enumerate(handle):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(header) + 1:
                raise ValueError(f"Dense RNA row {row_index + 2} has the wrong width: {path}")
            genes.append(fields[0])
            for column_index, text in enumerate(fields[1:]):
                value = int(text)
                if value < 0:
                    raise ValueError(f"Negative RNA count in {path}")
                if value:
                    rows.append(column_index)
                    columns.append(row_index)
                    values.append(value)
    matrix = sparse.coo_matrix(
        (np.asarray(values, dtype=np.int64), (rows, columns)),
        shape=(len(header), len(genes)),
        dtype=np.int64,
    ).tocsr()
    obs = pd.DataFrame(index=pd.Index(header, name="pixel"))
    var = pd.DataFrame(
        {"source_name": genes},
        index=pd.Index(unique_names(genes), name="feature"),
    )
    return ad.AnnData(X=matrix, obs=obs, var=var)


def read_10x_bundle(root: Path) -> ad.AnnData:
    matrices = list(root.rglob("matrix.mtx.gz")) + list(root.rglob("matrix.mtx"))
    features = list(root.rglob("features.tsv.gz")) + list(root.rglob("features.tsv"))
    barcodes = list(root.rglob("barcodes.tsv.gz")) + list(root.rglob("barcodes.tsv"))
    if len(matrices) != 1 or len(features) != 1 or len(barcodes) != 1:
        raise ValueError(f"Expected one 10x matrix/features/barcodes triplet under {root}")
    matrix = scipy_io.mmread(matrices[0]).tocsr()
    opener = gzip.open if features[0].name.endswith(".gz") else open
    with opener(features[0], "rt") as handle:
        feature_rows = [line.rstrip("\r\n").split("\t") for line in handle if line.strip()]
    opener = gzip.open if barcodes[0].name.endswith(".gz") else open
    with opener(barcodes[0], "rt") as handle:
        barcode_values = [line.rstrip("\r\n").split("\t")[0] for line in handle if line.strip()]
    if matrix.shape != (len(feature_rows), len(barcode_values)):
        raise ValueError(f"10x bundle axes do not match matrix dimensions under {root}")
    names = [row[1] if len(row) > 1 else row[0] for row in feature_rows]
    ids = [row[0] for row in feature_rows]
    obs = pd.DataFrame(index=pd.Index(barcode_values, name="pixel"))
    var = pd.DataFrame(
        {"feature_id": ids, "source_name": names},
        index=pd.Index(unique_names(names), name="feature"),
    )
    return ad.AnnData(X=matrix.T.tocsr().astype(np.int64), obs=obs, var=var)


def preprocess_matrices(config: Mapping[str, Any], overwrite: bool = False) -> None:
    family_root = project_path(str(config["processed_directory"]))
    records = []
    for path in iter_payloads(config):
        role = classify_payload(path, config)
        if not (role in {"rna_dense_matrix", "rna_matrix_bundle", "protein_matrix_bundle"}
                or role.endswith("_dense_matrix")):
            continue
        gsm = accession_from_filename(path)
        if role == "protein_matrix_bundle":
            modality = "protein"
            output_root = family_root / "validation_modalities/protein" / gsm
            output_name = "protein.h5ad"
        elif role in {"rna_dense_matrix", "rna_matrix_bundle"}:
            modality = "rna"
            output_root = family_root / "validation_modalities/rna" / gsm
            output_name = "rna.h5ad"
        else:
            modality = role.removesuffix("_dense_matrix")
            output_root = family_root / "validation_modalities/epigenome" / gsm / modality
            output_name = "matrix.h5ad"
        output_root.mkdir(parents=True, exist_ok=True)
        if role.endswith("_dense_matrix"):
            matrix = read_dense_gene_by_pixel(path)
            extracted_files: tuple[Path, ...] = ()
        else:
            bundle_root = output_root / "deposited"
            extracted_files = extract_nested_archive(path, bundle_root, overwrite)
            matrix = read_10x_bundle(bundle_root)
        matrix.uns["source"] = {
            "accession": config["accession"],
            "gsm": gsm,
            "source_path": str(path.relative_to(ROOT)),
            "source_sha256": sha256_file(path),
            "modality": modality,
            "truth_status": "orthogonal_validation_not_composition_truth",
        }
        output = output_root / output_name
        if output.exists() and not overwrite:
            raise FileExistsError(output)
        matrix.write_h5ad(output, compression="gzip")
        records.append(
            {
                "gsm": gsm,
                "role": role,
                "modality": modality,
                "source": str(path.relative_to(ROOT)),
                "output": str(output.relative_to(ROOT)),
                "output_bytes": output.stat().st_size,
                "output_sha256": sha256_file(output),
                "pixels": matrix.n_obs,
                "features": matrix.n_vars,
                "nonzero": int(matrix.X.nnz),
                "total_counts": int(matrix.X.sum()),
                "deposited_files": [str(item.relative_to(ROOT)) for item in extracted_files],
            }
        )
        print(f"matrix {gsm} {modality}: {matrix.shape}", flush=True)
    atomic_yaml(
        family_root / "source_audit/matrices.yaml",
        {
            "schema_version": 1,
            "accession": config["accession"],
            "matrix_payloads": len(records),
            "files": records,
        },
    )


def read_lines(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def canonical_fragment_barcodes(values: set[str], config: Mapping[str, Any]) -> set[str]:
    policy = config.get("barcode_policy", {})
    suffix = policy.get("fragment_terminal_suffix_to_strip")
    if not suffix:
        return values
    suffix = str(suffix)
    if policy.get("require_suffix_on_every_fragment_barcode", False):
        unexpected = sorted(value for value in values if not value.endswith(suffix))
        if unexpected:
            raise ValueError(
                f"Fragment barcodes do not all have required suffix {suffix!r}: "
                f"{unexpected[:5]}"
            )
    canonical = {
        value[: -len(suffix)] if value.endswith(suffix) else value for value in values
    }
    if len(canonical) != len(values):
        raise ValueError("Fragment barcode normalization is not injective")
    return canonical


def barcode_sources_for_gsm(
    family_root: Path, gsm: str, config: Mapping[str, Any]
) -> dict[str, set[str]]:
    sources: dict[str, set[str]] = {}
    for path in sorted((family_root / "source_audit/barcodes").glob(f"{gsm}__*.txt")):
        sources[f"fragments:{path.name}"] = canonical_fragment_barcodes(
            read_lines(path), config
        )
    coordinate_path = family_root / "spatial_coordinates" / gsm / "coordinates.csv"
    if coordinate_path.exists():
        frame = pd.read_csv(coordinate_path, dtype={"barcode": str})
        sources["coordinates"] = set(frame["barcode"])
    for modality in ("rna", "protein"):
        matrix_path = family_root / "validation_modalities" / modality / gsm / f"{modality}.h5ad"
        if matrix_path.exists():
            matrix = ad.read_h5ad(matrix_path, backed="r")
            sources[f"matrix:{modality}"] = set(map(str, matrix.obs_names))
            matrix.file.close()
    epigenome_root = family_root / "validation_modalities/epigenome" / gsm
    for matrix_path in sorted(epigenome_root.glob("*/matrix.h5ad")):
        matrix = ad.read_h5ad(matrix_path, backed="r")
        sources[f"matrix:epigenome:{matrix_path.parent.name}"] = set(map(str, matrix.obs_names))
        matrix.file.close()
    return sources


def align_modalities(config: Mapping[str, Any]) -> None:
    family_root = project_path(str(config["processed_directory"]))
    output_root = family_root / "cross_modality_alignment"
    output_root.mkdir(parents=True, exist_ok=True)
    summaries = []
    for group in config["sample_groups"]:
        gsms = []
        for key in (
            "atac",
            "rna",
            "protein",
            "validation_epigenome",
            "validation_epigenome_matrix",
        ):
            gsms.extend(map(str, group.get(key, [])))
        gsms = list(dict.fromkeys(gsms))
        sources: dict[str, set[str]] = {}
        for gsm in gsms:
            for label, values in barcode_sources_for_gsm(family_root, gsm, config).items():
                sources[f"{gsm}:{label}"] = values
        intersections = []
        labels = sorted(sources)
        for left_index, left in enumerate(labels):
            for right in labels[left_index + 1 :]:
                overlap = sources[left] & sources[right]
                union = sources[left] | sources[right]
                intersections.append(
                    {
                        "left": left,
                        "right": right,
                        "overlap": len(overlap),
                        "left_only": len(sources[left] - sources[right]),
                        "right_only": len(sources[right] - sources[left]),
                        "jaccard": float(len(overlap) / len(union)) if union else 1.0,
                    }
                )
        record = {
            "schema_version": 1,
            "group": group["id"],
            "declared_samples": gsms,
            "barcode_sources": {label: len(values) for label, values in sources.items()},
            "barcode_policy": config.get("barcode_policy", {"canonical_form": "identity"}),
            "pairwise": intersections,
            "exact_global_match": bool(sources) and len({frozenset(values) for values in sources.values()}) == 1,
            "truth_status": "orthogonal_validation_not_exact_composition_truth",
        }
        atomic_yaml(output_root / f"{group['id']}.yaml", record)
        summaries.append(record)
    atomic_yaml(
        family_root / "source_audit/alignment.yaml",
        {
            "schema_version": 1,
            "accession": config["accession"],
            "groups": summaries,
        },
    )


def finalize(config: Mapping[str, Any]) -> None:
    family_root = project_path(str(config["processed_directory"]))
    fragments = load_yaml(family_root / "source_audit/fragments.yaml")
    spatial = load_yaml(family_root / "source_audit/spatial.yaml")
    matrices = load_yaml(family_root / "source_audit/matrices.yaml")
    alignment = load_yaml(family_root / "source_audit/alignment.yaml")
    atac_files = [record for record in fragments["files"] if record["assay"] == "atac"]
    gates = {
        "source_archive_acquired": source_archive(config).is_file(),
        "all_fragment_rows_valid": all(record["invalid_rows"] == 0 for record in fragments["files"]),
        "all_fragment_files_coordinate_sorted": all(record["coordinate_sorted"] for record in fragments["files"]),
        "atac_files_separated_from_validation_epigenome": bool(atac_files),
        "spatial_assets_preprocessed": spatial["spatial_archives"] > 0,
        "orthogonal_matrices_preprocessed": matrices["matrix_payloads"] > 0,
        "reference_feature_axis_ready": False,
        "shape_layers_materialized": False,
        "runnable_dataset_registered": False,
    }
    record = {
        "schema_version": 1,
        "source_dataset_id": config["source_dataset_id"],
        "accession": config["accession"],
        "status": "source_preprocessing_complete_reference_gate_pending",
        "payload_files": len(iter_payloads(config)),
        "fragment_files": fragments["fragment_files"],
        "atac_fragment_files": len(atac_files),
        "spatial_archives": spatial["spatial_archives"],
        "matrix_payloads": matrices["matrix_payloads"],
        "sample_groups": len(alignment["groups"]),
        "validation_gates": gates,
        "reference_candidates": config.get("reference_candidates", {}),
        "blocking_gate": (
            "A compatible labeled fragment-level single-cell ATAC reference must be "
            "frozen before reference-only peak selection and ShapeMix layer construction."
        ),
        "truth_policy": "real spatial qualitative validation; no exact cell proportions",
    }
    atomic_yaml(family_root / "manifests/preprocessing.yaml", record)
    print(f"wrote preprocessing manifest for {config['accession']}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=("extract", "audit", "normalize", "spatial", "matrices", "align", "finalize", "all"),
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_yaml(args.config)
    if config.get("accession") not in {"GSE205055", "GSE263333"}:
        raise ValueError(f"Unsupported spatial manifest: {args.config}")
    stages = (
        ("extract", lambda: extract_source_archive(config, args.overwrite)),
        ("audit", lambda: audit_payloads(config)),
        ("normalize", lambda: normalize_fragments(config, args.overwrite)),
        ("spatial", lambda: preprocess_spatial_assets(config, args.overwrite)),
        ("matrices", lambda: preprocess_matrices(config, args.overwrite)),
        ("align", lambda: align_modalities(config)),
        ("finalize", lambda: finalize(config)),
    )
    for name, function in stages:
        if args.stage in {"all", name}:
            function()


if __name__ == "__main__":
    main()
