#!/usr/bin/env python3
"""Prepare and index the checksum-pinned UCSC mm10 reference for ShapeMix."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import yaml


ROOT = Path(__file__).resolve().parents[1]
CHUNK_BYTES = 8 * 1024 * 1024
BWA_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
PRIMARY_MM10 = {
    "chr1": 195471971,
    "chr2": 182113224,
    "chr3": 160039680,
    "chr4": 156508116,
    "chr5": 151834684,
    "chr6": 149736546,
    "chr7": 145441459,
    "chr8": 129401213,
    "chr9": 124595110,
    "chr10": 130694993,
    "chr11": 122082543,
    "chr12": 120129022,
    "chr13": 120421639,
    "chr14": 124902244,
    "chr15": 104043685,
    "chr16": 98207768,
    "chr17": 94987271,
    "chr18": 90702639,
    "chr19": 61431566,
    "chrX": 171031299,
    "chrY": 91744698,
}


def repository_path(path: Path) -> str:
    return str(path.absolute().relative_to(ROOT.absolute()))


def project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


def file_digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_chrom_sizes(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not fields[0]:
                raise ValueError(f"Malformed chromosome-size row {line_number}: {path}")
            if fields[0] in result:
                raise ValueError(f"Duplicate chromosome {fields[0]!r}: {path}")
            result[fields[0]] = int(fields[1])
    if not result:
        raise ValueError(f"Chromosome-size file is empty: {path}")
    return result


def parse_fasta_lengths(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    current: str | None = None
    length = 0
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip(b"\r\n")
            if line.startswith(b">"):
                if current is not None:
                    result[current] = length
                identifier = line[1:].split(None, 1)[0]
                if not identifier:
                    raise ValueError(f"Empty FASTA identifier at line {line_number}")
                current = identifier.decode("ascii")
                if current in result:
                    raise ValueError(f"Duplicate FASTA sequence {current!r}")
                length = 0
            else:
                if current is None:
                    raise ValueError("FASTA sequence appears before its first header")
                if not line or re.fullmatch(rb"[A-Za-z*.-]+", line) is None:
                    raise ValueError(f"Malformed FASTA sequence at line {line_number}")
                length += len(line)
    if current is not None:
        result[current] = length
    if not result:
        raise ValueError(f"FASTA contains no sequences: {path}")
    return result


def validate_reference_contract(
    fasta_lengths: Mapping[str, int], chrom_sizes: Mapping[str, int]
) -> None:
    if dict(fasta_lengths) != dict(chrom_sizes):
        missing = sorted(set(chrom_sizes).difference(fasta_lengths))
        extra = sorted(set(fasta_lengths).difference(chrom_sizes))
        different = sorted(
            name
            for name in set(fasta_lengths).intersection(chrom_sizes)
            if fasta_lengths[name] != chrom_sizes[name]
        )
        raise ValueError(
            "FASTA/chrom-size mismatch: "
            f"missing={missing} extra={extra} different={different}"
        )
    observed_primary = {
        name: fasta_lengths[name] for name in PRIMARY_MM10 if name in fasta_lengths
    }
    if observed_primary != PRIMARY_MM10:
        raise ValueError(
            "Reference does not reproduce the deposited GSE246791 primary-contig contract"
        )


def source_record(lock: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    matches = [
        record
        for record in lock.get("files", [])
        if isinstance(record, Mapping) and record.get("name") == name
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one lock record for {name}; found {len(matches)}")
    return matches[0]


def decompress_fasta(source: Path, destination: Path) -> None:
    if destination.exists():
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".partial")
    if temporary.exists():
        raise FileExistsError(f"Unresolved partial decompression: {temporary}")
    try:
        with gzip.open(source, "rb") as input_handle, temporary.open("xb") as output_handle:
            for chunk in iter(lambda: input_handle.read(CHUNK_BYTES), b""):
                output_handle.write(chunk)
            output_handle.flush()
            os.fsync(output_handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def bwa_version(bwa: Path) -> str:
    completed = subprocess.run(
        [str(bwa)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    match = re.search(r"Version:\s*(\S+)", completed.stdout)
    if match is None:
        raise ValueError(f"Could not determine BWA version from {bwa}")
    return match.group(1)


def build_bwa_index(fasta: Path, bwa: Path, work_root: Path) -> dict[str, Any]:
    outputs = [Path(f"{fasta}{suffix}") for suffix in BWA_SUFFIXES]
    existing = [path.is_file() and path.stat().st_size > 0 for path in outputs]
    if any(existing) and not all(existing):
        raise FileExistsError("Incomplete final BWA index exists; refusing to overwrite it")
    if not all(existing):
        work_root.mkdir(parents=True, exist_ok=True)
        scratch = Path(tempfile.mkdtemp(prefix="bwa-index-", dir=work_root))
        prefix = scratch / "mm10.fa"
        subprocess.run(
            [str(bwa), "index", "-p", str(prefix), str(fasta)],
            check=True,
        )
        staged = [Path(f"{prefix}{suffix}") for suffix in BWA_SUFFIXES]
        if any(not path.is_file() or path.stat().st_size == 0 for path in staged):
            raise ValueError(f"BWA did not produce a complete nonempty index under {scratch}")
        for staged_path, output_path in zip(staged, outputs, strict=True):
            os.replace(staged_path, output_path)
    return {
        "version": bwa_version(bwa),
        "executable": str(bwa.absolute()),
        "command": [str(bwa), "index", "-p", "<scratch>/mm10.fa", repository_path(fasta)],
        "files": [
            {
                "path": repository_path(path),
                "bytes": path.stat().st_size,
                "sha256": file_digest(path, "sha256"),
            }
            for path in outputs
        ],
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
        default=ROOT / "configs/data_sources/shapemix_ucsc_mm10_initial.yaml",
    )
    parser.add_argument(
        "--lock",
        type=Path,
        default=ROOT / "configs/data_sources/shapemix_ucsc_mm10_initial_lock.yaml",
    )
    parser.add_argument("--bwa", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run this preprocessing stage through run_shapemix_low_impact.sh")
    config = load_yaml(args.config)
    lock = load_yaml(args.lock)
    raw_root = project_path(str(config["raw_directory"]))
    processed_root = project_path(str(config["processed_directory"]))
    fasta_source = raw_root / "mm10.fa.gz"
    sizes_source = raw_root / "mm10.chrom.sizes"
    for path in (fasta_source, sizes_source):
        if not path.is_file():
            raise FileNotFoundError(path)
        record = source_record(lock, path.name)
        if path.stat().st_size != int(record["bytes"]):
            raise ValueError(f"Locked source byte count changed: {path}")
        if file_digest(path, "sha256") != str(record["sha256"]):
            raise ValueError(f"Locked source SHA-256 changed: {path}")
    fasta = processed_root / "mm10.fa"
    decompress_fasta(fasta_source, fasta)
    lengths = parse_fasta_lengths(fasta)
    sizes = parse_chrom_sizes(sizes_source)
    validate_reference_contract(lengths, sizes)
    index = None
    if args.bwa is not None:
        bwa = args.bwa.absolute()
        if not bwa.is_file():
            raise FileNotFoundError(bwa)
        index = build_bwa_index(
            fasta,
            bwa,
            ROOT / "data/work/preprocessing/mm10_ucsc_initial/bwa_index",
        )
    record = {
        "schema_version": 1,
        "reference_id": "mm10_ucsc_initial",
        "assembly": "GRCm38 initial release",
        "assembly_accession": "GCA_000001635.2",
        "genome_build": "mm10",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_config": repository_path(args.config),
        "source_lock": repository_path(args.lock),
        "fasta": {
            "path": repository_path(fasta),
            "bytes": fasta.stat().st_size,
            "sha256": file_digest(fasta, "sha256"),
            "sequences": len(lengths),
        },
        "chrom_sizes": {
            "path": repository_path(sizes_source),
            "bytes": sizes_source.stat().st_size,
            "sha256": file_digest(sizes_source, "sha256"),
        },
        "validation": {
            "fasta_matches_chrom_sizes": True,
            "gse246791_primary_contigs_match": True,
            "primary_contigs": len(PRIMARY_MM10),
        },
        "bwa_index": index,
    }
    atomic_yaml(processed_root / "reference.yaml", record)
    print(
        f"prepared {repository_path(fasta)} sequences={len(lengths)} "
        f"bwa_index={'ready' if index is not None else 'not_requested'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
