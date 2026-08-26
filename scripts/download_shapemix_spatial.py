#!/usr/bin/env python
"""Download and validate a frozen ShapeMix spatial GEO source family."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import os
import subprocess
import tarfile
import tempfile
import urllib.request
import zipfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Mapping, Optional

import yaml


ROOT = Path(__file__).resolve().parents[1]
CHUNK_BYTES = 8 * 1024 * 1024


@dataclass(frozen=True)
class Resource:
    name: str
    role: str
    accession: str
    url: str
    destination: Path
    staging_path: Path
    expected_bytes: Optional[int] = None
    expected_etag: Optional[str] = None
    expected_md5: Optional[str] = None
    payload: str = "gzip"


@dataclass(frozen=True)
class DownloadResult:
    name: str
    role: str
    accession: str
    url: str
    path: str
    bytes: int
    remote_etag: Optional[str]
    provider_md5: Optional[str]
    observed_md5: Optional[str]
    sha256: str
    integrity: str
    tar_members: Optional[int]


def load_config(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    accession = str(value.get("accession", ""))
    if accession not in {
        "GSE205055",
        "GSE263333",
        "GSE216371",
        "GSE246791",
        "GSE244618",
        "UCSC_mm10_initial",
    }:
        raise ValueError(f"Unsupported ShapeMix source manifest: {path}")
    return value


def project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def geo_series_stem(accession: str) -> str:
    if not accession.startswith("GSE") or len(accession) < 7:
        raise ValueError(f"Invalid GEO Series accession: {accession}")
    return f"{accession[:-3]}nnn"


def metadata_url(accession: str) -> str:
    stem = geo_series_stem(accession)
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{accession}/soft/"
        f"{accession}_family.soft.gz"
    )


def resources_from_config(config: Mapping[str, Any]) -> tuple[Resource, ...]:
    raw_root = project_path(str(config["raw_directory"]))
    work_root = project_path(str(config["work_directory"]))
    if "resources" in config:
        resources = []
        for record in config["resources"]:
            relative = Path(str(record["destination"]))
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"Unsafe resource destination: {relative}")
            resources.append(
                Resource(
                    name=str(record["name"]),
                    role=str(record["role"]),
                    accession=str(record.get("gsm", config["accession"])),
                    url=str(record["url"]),
                    destination=raw_root / relative,
                    staging_path=work_root / relative.parent / f"{relative.name}.part",
                    expected_bytes=int(record["expected_bytes"]),
                    expected_etag=(
                        str(record["expected_etag"])
                        if record.get("expected_etag") is not None
                        else None
                    ),
                    expected_md5=(
                        str(record["expected_md5"]).lower()
                        if record.get("expected_md5") is not None
                        else None
                    ),
                    payload=str(record.get("payload", "gzip")),
                )
            )
        names = [resource.name for resource in resources]
        destinations = [resource.destination for resource in resources]
        if len(names) != len(set(names)) or len(destinations) != len(set(destinations)):
            raise ValueError("Resolved source names and destinations must be unique")
        return tuple(resources)
    resources: list[Resource] = []
    for record in config["series_metadata"]:
        accession = str(record["accession"])
        name = f"{accession}_family.soft.gz"
        resources.append(
            Resource(
                name=name,
                role=str(record["role"]),
                accession=accession,
                url=metadata_url(accession),
                destination=raw_root / "series_metadata" / name,
                staging_path=work_root / "series_metadata" / f"{name}.part",
            )
        )
    archive = config["archive"]
    name = str(archive["filename"])
    resources.append(
        Resource(
            name=name,
            role=str(archive["role"]),
            accession=str(config["accession"]),
            url=str(archive["url"]),
            destination=raw_root / "processed_downloads" / name,
            staging_path=work_root / f"{name}.part",
            expected_bytes=int(archive["expected_bytes"]),
            payload="tar",
        )
    )
    names = [resource.name for resource in resources]
    if len(names) != len(set(names)):
        raise ValueError("Resolved source filenames are not unique")
    return tuple(resources)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_etag(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    normalized = value.strip()
    if normalized.startswith("W/"):
        normalized = normalized[2:].strip()
    return normalized.strip('"')


def remote_metadata(url: str, timeout: int) -> tuple[int, Optional[str]]:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=timeout) as response:
        value = response.headers.get("Content-Length")
        etag = normalize_etag(response.headers.get("ETag"))
    if value is None or int(value) <= 0:
        raise ValueError(f"Source did not publish a positive Content-Length: {url}")
    return int(value), etag


def safe_tar_member_name(name: str) -> bool:
    path = PurePosixPath(name)
    return bool(name) and not path.is_absolute() and ".." not in path.parts


def validate_tar(path: Path) -> int:
    count = 0
    with tarfile.open(path, "r:") as archive:
        for member in archive:
            if not safe_tar_member_name(member.name):
                raise ValueError(f"Unsafe tar member: {member.name!r}")
            if member.issym() or member.islnk() or member.isdev():
                raise ValueError(f"Unsupported tar member type: {member.name!r}")
            count += 1
    if count == 0:
        raise ValueError(f"Archive contains no members: {path}")
    return count


def validate_gzip(path: Path) -> None:
    with gzip.open(path, "rb") as handle:
        for _ in iter(lambda: handle.read(CHUNK_BYTES), b""):
            pass


def validate_gzip_schema(path: Path, name: str) -> None:
    lower = name.lower()
    if lower.endswith(".h5ad.gz"):
        with gzip.open(path, "rb") as handle:
            if handle.read(8) != b"\x89HDF\r\n\x1a\n":
                raise ValueError(f"Gzip payload is not an HDF5/H5AD object: {path}")
        return
    minimum_columns = None
    if lower.endswith(".bedpe.gz"):
        minimum_columns = 6
    elif lower.endswith(".bed.gz"):
        minimum_columns = 3
    elif lower.endswith(".tsv.gz"):
        minimum_columns = 2
    if minimum_columns is not None:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
            first = next(
                (
                    line
                    for line in handle
                    if line.strip() and not line.startswith("#")
                ),
                None,
            )
        if (
            first is None
            or len(first.rstrip("\r\n").split("\t")) < minimum_columns
        ):
            raise ValueError(f"Compressed table has no valid first row: {path}")
    elif lower.endswith(".soft.gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            found = any(
                line.startswith("^SERIES")
                for _, line in zip(range(1000), handle)
            )
        if not found:
            raise ValueError(f"GEO SOFT payload lacks a SERIES record: {path}")


def validate_xlsx(path: Path) -> None:
    with zipfile.ZipFile(path) as workbook:
        bad_member = workbook.testzip()
        if bad_member is not None:
            raise ValueError(f"Corrupt XLSX member {bad_member!r} in {path}")
        if "[Content_Types].xml" not in workbook.namelist():
            raise ValueError(f"XLSX lacks [Content_Types].xml: {path}")


def download_once(resource: Resource, expected_bytes: int, timeout: int) -> None:
    resource.staging_path.parent.mkdir(parents=True, exist_ok=True)
    observed = resource.staging_path.stat().st_size if resource.staging_path.exists() else 0
    if observed > expected_bytes:
        raise IOError(
            f"Staging file is larger than the official source and will not be overwritten: "
            f"{resource.staging_path}"
        )
    if observed < expected_bytes:
        subprocess.run(
            [
                "wget",
                "-q",
                "--continue",
                f"--timeout={timeout}",
                "--tries=20",
                "--output-document",
                str(resource.staging_path),
                resource.url,
            ],
            check=True,
        )
    final_size = resource.staging_path.stat().st_size
    if final_size != expected_bytes:
        raise IOError(
            f"Incomplete transfer for {resource.name}: {final_size} != {expected_bytes}"
        )


def acquire(resource: Resource, timeout: int) -> DownloadResult:
    official_bytes, remote_etag = remote_metadata(resource.url, timeout)
    if resource.expected_bytes is not None and official_bytes != resource.expected_bytes:
        raise ValueError(
            f"Official size changed for {resource.name}: "
            f"{official_bytes} != frozen {resource.expected_bytes}"
        )
    if resource.expected_etag is not None:
        if remote_etag != normalize_etag(resource.expected_etag):
            raise ValueError(
                f"Official ETag changed for {resource.name}: "
                f"{remote_etag!r} != frozen {normalize_etag(resource.expected_etag)!r}"
            )
    resource.destination.parent.mkdir(parents=True, exist_ok=True)
    source = resource.destination
    if source.exists() and source.stat().st_size != official_bytes:
        raise IOError(f"Existing immutable source has an unexpected size: {source}")
    if not source.exists():
        download_once(resource, official_bytes, timeout)
        source = resource.staging_path

    tar_members: Optional[int] = None
    if resource.payload == "tar":
        tar_members = validate_tar(source)
        integrity = "passed_tar_stream_and_safe_member_audit"
    elif resource.payload == "gzip":
        validate_gzip(source)
        validate_gzip_schema(source, resource.name)
        integrity = "passed_full_gzip_stream_and_schema"
    elif resource.payload == "xlsx":
        validate_xlsx(source)
        integrity = "passed_xlsx_zip_crc_and_manifest"
    elif resource.payload in {"text", "csv", "tsv"}:
        integrity = "passed_exact_byte_transfer"
    else:
        raise ValueError(
            f"Unsupported payload type for {resource.name}: {resource.payload}"
        )
    digest = sha256_file(source)
    observed_md5 = md5_file(source) if resource.expected_md5 is not None else None
    if observed_md5 != resource.expected_md5:
        raise ValueError(
            f"Provider MD5 changed for {resource.name}: "
            f"{observed_md5!r} != frozen {resource.expected_md5!r}"
        )
    if source == resource.staging_path:
        os.replace(source, resource.destination)
    return DownloadResult(
        name=resource.name,
        role=resource.role,
        accession=resource.accession,
        url=resource.url,
        path=str(resource.destination.relative_to(ROOT)),
        bytes=official_bytes,
        remote_etag=remote_etag,
        provider_md5=resource.expected_md5,
        observed_md5=observed_md5,
        sha256=digest,
        integrity=integrity,
        tar_members=tar_members,
    )


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
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--lock-output", type=Path)
    parser.add_argument("--audit-output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    config = load_config(args.config)
    worker_limit = min(2, int(config.get("download_workers_max", 2)))
    if args.workers > worker_limit:
        raise ValueError(f"--workers must be <= {worker_limit} for this campaign")
    resources = resources_from_config(config)
    results: list[DownloadResult] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        pending = {
            executor.submit(acquire, resource, args.timeout): resource for resource in resources
        }
        for future in concurrent.futures.as_completed(pending):
            result = future.result()
            results.append(result)
            print(f"validated {result.name} {result.bytes} {result.sha256}", flush=True)
    order = {resource.name: index for index, resource in enumerate(resources)}
    results.sort(key=lambda result: order[result.name])
    config_path = args.config.resolve().relative_to(ROOT.resolve())
    record = {
        "schema_version": 1,
        "source_dataset_id": config["source_dataset_id"],
        "accession": config["accession"],
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "config": str(config_path),
        "resources": len(results),
        "total_bytes": sum(result.bytes for result in results),
        "files": [asdict(result) for result in results],
    }
    default_lock = args.config.with_name(f"{args.config.stem}_lock.yaml")
    lock_output = args.lock_output or default_lock
    audit_filename = str(config.get("source_audit_filename", "downloads.yaml"))
    if Path(audit_filename).name != audit_filename:
        raise ValueError(f"Unsafe source-audit filename: {audit_filename}")
    default_audit = (
        project_path(str(config["processed_directory"]))
        / "source_audit"
        / audit_filename
    )
    audit_output = args.audit_output or default_audit
    atomic_yaml(lock_output, record)
    atomic_yaml(audit_output, record)
    print(f"wrote {lock_output}", flush=True)
    print(f"wrote {audit_output}", flush=True)


if __name__ == "__main__":
    main()

