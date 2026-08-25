#!/usr/bin/env python
"""Acquire and validate the gated GSE194122 ShapeMix source scope."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import os
import subprocess
import tempfile
import urllib.request
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Optional

import pysam
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/data_sources/shapemix_gse194122.yaml"
CHUNK_BYTES = 8 * 1024 * 1024


@dataclass(frozen=True)
class Resource:
    name: str
    role: str
    url: str
    destination: Path
    staging_path: Path
    expected_bytes: int
    payload: str
    expected_md5: Optional[str] = None
    expected_etag: Optional[str] = None
    gsm: Optional[str] = None
    run: Optional[str] = None
    sample_key: Optional[str] = None


@dataclass(frozen=True)
class DownloadResult:
    name: str
    role: str
    url: str
    path: str
    bytes: int
    sha256: str
    md5: str
    integrity: str
    etag: Optional[str]
    gsm: Optional[str]
    run: Optional[str]
    sample_key: Optional[str]


def load_config(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict) or value.get("accession") != "GSE194122":
        raise ValueError(f"Not a GSE194122 source manifest: {path}")
    return value


def project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def pilot_sample(config: Mapping[str, Any]) -> Mapping[str, Any]:
    matches = [
        sample
        for sample in config["atac_samples"]
        if str(sample["gsm"]) == str(config["pilot_gsm"])
    ]
    if len(matches) != 1:
        raise ValueError("pilot_gsm must select exactly one ATAC sample")
    return matches[0]


def resources_from_config(
    config: Mapping[str, Any],
    include_pilot_bam: bool,
    include_fragments: bool = False,
    fragment_sample_keys: Optional[set[str]] = None,
) -> tuple[Resource, ...]:
    raw = project_path(str(config["raw_directory"]))
    work = project_path(str(config["work_directory"]))
    records = (
        (config["series_file"], raw / "series_metadata", work / "series_metadata", "gzip"),
        (config["processed_object"], raw / "processed_downloads", work, "gzip"),
    )
    resources = [
        Resource(
            name=str(record["filename"]),
            role=str(record["role"]),
            url=str(record["url"]),
            destination=destination / str(record["filename"]),
            staging_path=staging / f"{record['filename']}.part",
            expected_bytes=int(record["expected_bytes"]),
            payload=payload,
        )
        for record, destination, staging, payload in records
    ]
    if include_pilot_bam:
        sample = pilot_sample(config)
        name = str(sample["bam_filename"])
        resources.append(
            Resource(
                name=name,
                role="pilot_author_cellranger_arc_atac_bam",
                url=str(sample["bam_url"]),
                destination=(
                    raw
                    / "samples"
                    / str(sample["gsm"])
                    / "source_files"
                    / name
                ),
                staging_path=work / "samples" / str(sample["gsm"]) / f"{name}.part",
                expected_bytes=int(sample["bam_bytes"]),
                expected_md5=str(sample["bam_md5"]),
                payload="bam",
                gsm=str(sample["gsm"]),
                run=str(sample["run"]),
            )
        )
    if include_fragments:
        suite = config["openproblems_fragments"]
        base_url = str(suite["base_url"])
        roles = {
            "atac_fragments.tsv.gz": "author_cellranger_arc_fragments",
            "atac_fragments.tsv.gz.tbi": "author_cellranger_arc_fragment_index",
            "per_barcode_metrics.csv": "author_cellranger_arc_barcode_metrics",
        }
        payloads = {
            "atac_fragments.tsv.gz": "gzip",
            "atac_fragments.tsv.gz.tbi": "tabix",
            "per_barcode_metrics.csv": "csv",
        }
        for record in suite["files"]:
            sample_key = str(record["sample_key"])
            if fragment_sample_keys is not None and sample_key not in fragment_sample_keys:
                continue
            gsm = str(record["gsm"])
            for name, byte_key, etag_key in (
                ("atac_fragments.tsv.gz", "fragment_bytes", "fragment_etag"),
                ("atac_fragments.tsv.gz.tbi", "index_bytes", "index_etag"),
                ("per_barcode_metrics.csv", "metrics_bytes", "metrics_etag"),
            ):
                resources.append(
                    Resource(
                        name=f"{sample_key}/{name}",
                        role=roles[name],
                        url=f"{base_url}/{sample_key}/{name}",
                        destination=raw / "samples" / gsm / "source_files" / name,
                        staging_path=work / "samples" / gsm / f"{name}.part",
                        expected_bytes=int(record[byte_key]),
                        expected_etag=str(record[etag_key]),
                        payload=payloads[name],
                        gsm=gsm,
                        sample_key=sample_key,
                    )
                )
    return tuple(resources)


def hash_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def remote_identity(url: str, timeout: int) -> tuple[int, Optional[str]]:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=timeout) as response:
        value = response.headers.get("Content-Length")
        etag = response.headers.get("ETag")
    if value is None or int(value) <= 0:
        raise ValueError(f"Source has no positive Content-Length: {url}")
    return int(value), etag.strip('"') if etag else None


def download(resource: Resource, timeout: int, aria2c: Optional[Path] = None) -> Path:
    official, etag = remote_identity(resource.url, timeout)
    if official != resource.expected_bytes:
        raise ValueError(
            f"Official byte size changed for {resource.name}: "
            f"{official} != {resource.expected_bytes}"
        )
    if resource.expected_etag is not None and etag != resource.expected_etag:
        raise ValueError(
            f"Official ETag changed for {resource.name}: {etag} != {resource.expected_etag}"
        )
    resource.destination.parent.mkdir(parents=True, exist_ok=True)
    if resource.destination.exists():
        if resource.destination.stat().st_size != official:
            raise IOError(f"Existing immutable source has wrong size: {resource.destination}")
        return resource.destination
    resource.staging_path.parent.mkdir(parents=True, exist_ok=True)
    staged = resource.staging_path.stat().st_size if resource.staging_path.exists() else 0
    if staged > official:
        raise IOError(f"Oversized staging file will not be overwritten: {resource.staging_path}")
    if staged < official:
        if aria2c is None:
            command = [
                "wget",
                "-q",
                "--continue",
                f"--timeout={timeout}",
                "--tries=20",
                "--output-document",
                str(resource.staging_path),
                resource.url,
            ]
        else:
            command = [
                str(aria2c),
                "--continue=true",
                "--allow-overwrite=true",
                "--auto-file-renaming=false",
                "--file-allocation=none",
                "--max-connection-per-server=4",
                "--split=4",
                "--min-split-size=20M",
                "--max-tries=20",
                "--retry-wait=3",
                f"--timeout={timeout}",
                "--summary-interval=0",
                "--console-log-level=warn",
                "--show-console-readout=false",
                "--dir",
                str(resource.staging_path.parent),
                "--out",
                resource.staging_path.name,
                resource.url,
            ]
        subprocess.run(command, check=True)
    if resource.staging_path.stat().st_size != official:
        raise IOError(f"Incomplete transfer: {resource.staging_path}")
    return resource.staging_path


def validate_bam(path: Path) -> str:
    with pysam.AlignmentFile(path, "rb") as handle:
        header = handle.header.to_dict()
        if header.get("HD", {}).get("SO") != "coordinate":
            raise ValueError("Pilot BAM is not coordinate-sorted")
        references = set(handle.references)
        if not {"chr1", "chrX", "chrM"}.issubset(references):
            raise ValueError("Pilot BAM does not use the expected GRCh38 chr-style axis")
        observed = 0
        has_cb = False
        for record in handle.fetch(until_eof=True):
            observed += 1
            has_cb = has_cb or record.has_tag("CB")
            if observed >= 10_000:
                break
        if observed == 0 or not has_cb:
            raise ValueError("Pilot BAM prefix has no records with corrected CB tags")
    return "passed_header_reference_and_10000_record_cb_audit"


def validate_and_promote(
    resource: Resource, timeout: int, aria2c: Optional[Path] = None
) -> DownloadResult:
    source = download(resource, timeout, aria2c=aria2c)
    if resource.payload == "gzip":
        with gzip.open(source, "rb") as handle:
            for _ in iter(lambda: handle.read(CHUNK_BYTES), b""):
                pass
        integrity = "passed_full_gzip_stream"
    elif resource.payload == "bam":
        integrity = validate_bam(source)
    elif resource.payload == "tabix":
        with source.open("rb") as handle:
            if handle.read(2) != b"\x1f\x8b":
                raise ValueError(f"Tabix index is not gzip/BGZF: {source}")
        integrity = "passed_tabix_gzip_signature"
    elif resource.payload == "csv":
        with source.open("rt", errors="strict") as handle:
            header = handle.readline().rstrip("\r\n").split(",")
        if not header or "barcode" not in {column.lower() for column in header}:
            raise ValueError(f"Per-barcode metrics has no barcode column: {source}")
        integrity = "passed_csv_header"
    else:
        raise ValueError(f"Unsupported payload: {resource.payload}")
    md5 = hash_file(source, "md5")
    if resource.expected_md5 is not None and md5 != resource.expected_md5:
        raise ValueError(f"MD5 mismatch for {resource.name}: {md5}")
    sha256 = hash_file(source, "sha256")
    if source == resource.staging_path:
        os.replace(source, resource.destination)
    _, etag = remote_identity(resource.url, timeout)
    return DownloadResult(
        name=resource.name,
        role=resource.role,
        url=resource.url,
        path=str(resource.destination.relative_to(ROOT)),
        bytes=resource.expected_bytes,
        sha256=sha256,
        md5=md5,
        integrity=integrity,
        etag=etag,
        gsm=resource.gsm,
        run=resource.run,
        sample_key=resource.sample_key,
    )


def atomic_yaml(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(descriptor, "w") as handle:
            yaml.safe_dump(dict(value), handle, sort_keys=False)
        os.replace(name, path)
    except BaseException:
        Path(name).unlink(missing_ok=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--include-pilot-bam", action="store_true")
    parser.add_argument("--include-fragments", action="store_true")
    parser.add_argument(
        "--fragment-sample",
        action="append",
        default=[],
        help="limit fragment-suite acquisition to a sample key; repeat as needed",
    )
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument(
        "--aria2c",
        type=Path,
        help="optional explicit aria2c executable for four-range resumable transfers",
    )
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument(
        "--lock-output",
        type=Path,
        default=ROOT / "configs/data_sources/shapemix_gse194122_lock.yaml",
    )
    args = parser.parse_args()
    config = load_config(args.config)
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    if args.aria2c is not None and not args.aria2c.is_file():
        raise FileNotFoundError(args.aria2c)
    selected_fragment_samples = set(args.fragment_sample) or None
    if selected_fragment_samples is not None:
        declared = {
            str(record["sample_key"])
            for record in config["openproblems_fragments"]["files"]
        }
        unknown = selected_fragment_samples - declared
        if unknown:
            raise ValueError(f"Unknown fragment sample keys: {sorted(unknown)}")
    resources = resources_from_config(
        config,
        args.include_pilot_bam,
        include_fragments=args.include_fragments,
        fragment_sample_keys=selected_fragment_samples,
    )
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        pending = {
            executor.submit(
                validate_and_promote, resource, args.timeout, args.aria2c
            ): resource
            for resource in resources
        }
        for future in concurrent.futures.as_completed(pending):
            result = future.result()
            results.append(result)
            print(f"validated {result.name} {result.bytes} {result.sha256}", flush=True)
    order = {resource.name: index for index, resource in enumerate(resources)}
    results.sort(key=lambda item: order[item.name])
    inventory = config["atac_samples"]
    lock = {
        "schema_version": 1,
        "source_dataset_id": config["source_dataset_id"],
        "accession": config["accession"],
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "config": str(args.config.resolve().relative_to(ROOT.resolve())),
        "pilot_bam_included": args.include_pilot_bam,
        "fragment_suite_included": args.include_fragments,
        "fragment_samples": sorted(selected_fragment_samples) if selected_fragment_samples else "all",
        "validated_bytes": sum(item.bytes for item in results),
        "files": [asdict(item) for item in results],
        "atac_remote_inventory": {
            "samples": len(inventory),
            "total_bam_bytes": sum(int(item["bam_bytes"]) for item in inventory),
            "policy": "NCBI BAMs are a frozen fallback and are not acquired when the author fragment suite validates",
        },
    }
    atomic_yaml(args.lock_output, lock)
    audit = project_path(str(config["processed_directory"])) / "source_audit/downloads.yaml"
    atomic_yaml(audit, lock)
    print(f"wrote {args.lock_output}", flush=True)


if __name__ == "__main__":
    main()
