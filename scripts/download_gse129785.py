#!/usr/bin/env python
"""Download and validate the frozen ShapeMix GSE129785 acquisition scope."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import os
import subprocess
import tempfile
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/data_sources/shapemix_gse129785.yaml"
CHUNK_BYTES = 8 * 1024 * 1024


@dataclass(frozen=True)
class Resource:
    """One immutable file in the frozen acquisition scope."""

    name: str
    role: str
    url: str
    destination: Path
    staging_path: Path
    gsm: Optional[str] = None


@dataclass(frozen=True)
class DownloadResult:
    """Validated source identity recorded after atomic promotion."""

    name: str
    role: str
    url: str
    path: str
    gsm: Optional[str]
    bytes: int
    sha256: str
    gzip_integrity: str
    fragment_schema: str


def load_config(path: Path) -> dict[str, Any]:
    """Load the tracked acquisition manifest."""
    with path.open() as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict) or config.get("accession") != "GSE129785":
        raise ValueError(f"Not a GSE129785 source manifest: {path}")
    return config


def project_path(value: str) -> Path:
    """Resolve a project-relative manifest path without hiding its identity."""
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def resources_from_config(config: dict[str, Any], scope: str) -> list[Resource]:
    """Resolve series and sample files into deterministic local destinations."""
    raw_root = project_path(config["raw_directory"])
    work_root = project_path(config["work_directory"])
    resources: list[Resource] = []
    if scope in {"all", "series"}:
        for record in config["series_files"]:
            name = str(record["filename"])
            resources.append(
                Resource(
                    name=name,
                    role=str(record["role"]),
                    url=str(record["url"]),
                    destination=raw_root / "series_metadata" / name,
                    staging_path=work_root / "series_metadata" / f"{name}.part",
                )
            )
    if scope in {"all", "samples"}:
        template = str(config["download_policy"]["fragment_url_template"])
        for sample in config["samples"]:
            gsm = str(sample["gsm"])
            title = str(sample["title"])
            name = f"{gsm}_{title}_fragments.tsv.gz"
            resources.append(
                Resource(
                    name=name,
                    role=str(sample["role"]),
                    url=template.format(gsm=gsm, title=title),
                    destination=raw_root / "samples" / gsm / "source_files" / name,
                    staging_path=work_root / "samples" / gsm / f"{name}.part",
                    gsm=gsm,
                )
            )
    names = [resource.name for resource in resources]
    if len(names) != len(set(names)):
        raise ValueError("Resolved source filenames are not unique.")
    return resources


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gzip_integrity(path: Path) -> None:
    """Read every decompressed byte so CRC and EOF failures cannot be missed."""
    with gzip.open(path, "rb") as handle:
        for _ in iter(lambda: handle.read(CHUNK_BYTES), b""):
            pass

def payload_integrity(path: Path, resource: Resource) -> str:
    """Validate gzip sources and GEO Matrix Market files served as plain text."""
    with path.open("rb") as handle:
        prefix = handle.read(32)
    if prefix.startswith(b"\x1f\x8b"):
        gzip_integrity(path)
        return "passed_full_gzip_stream"
    if resource.gsm is None and "count_matrix" in resource.role and prefix.startswith(b"%%MatrixMarket"):
        return "passed_plain_matrix_market_header_and_size"
    raise gzip.BadGzipFile(f"Unexpected payload signature for {resource.name}: {prefix[:8]!r}")



def validate_fragment_schema(path: Path, sample_rows: int = 10_000) -> None:
    """Validate a representative prefix against the five-column data contract."""
    observed = 0
    with gzip.open(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 5 or any(value == "" for value in fields):
                raise ValueError(f"{path}: fragment row {line_number} is not five-column TSV.")
            try:
                start, end, support = int(fields[1]), int(fields[2]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"{path}: fragment row {line_number} has non-integer fields.") from exc
            if start < 0 or end <= start or support < 1:
                raise ValueError(f"{path}: fragment row {line_number} has invalid coordinates/support.")
            observed += 1
            if observed >= sample_rows:
                break
    if observed == 0:
        raise ValueError(f"{path}: no fragment data rows were found.")


def remote_size(url: str, timeout: int) -> int:
    """Return the required official Content-Length."""
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=timeout) as response:
        value = response.headers.get("Content-Length")
    if value is None or int(value) <= 0:
        raise ValueError(f"Source did not publish a positive Content-Length: {url}")
    return int(value)


def _download_once(resource: Resource, expected_bytes: int, timeout: int) -> None:
    """Resume one source with wget and require the exact official byte size."""
    resource.staging_path.parent.mkdir(parents=True, exist_ok=True)
    existing = resource.staging_path.stat().st_size if resource.staging_path.exists() else 0
    if existing > expected_bytes:
        resource.staging_path.unlink()
        existing = 0
    if existing == expected_bytes:
        return
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
    observed = resource.staging_path.stat().st_size
    if observed != expected_bytes:
        raise IOError(
            f"Incomplete transfer for {resource.name}: observed {observed}, expected {expected_bytes}."
        )


def validate_and_promote(resource: Resource, expected_bytes: int, timeout: int) -> DownloadResult:
    """Acquire one file, validate it, and atomically publish immutable raw bytes."""
    resource.destination.parent.mkdir(parents=True, exist_ok=True)
    source = resource.destination
    if source.exists() and source.stat().st_size != expected_bytes:
        raise IOError(
            f"Existing raw file has unexpected size and will not be overwritten: {source}"
        )
    if not source.exists():
        _download_once(resource, expected_bytes, timeout)
        source = resource.staging_path
    integrity = payload_integrity(source, resource)
    schema = "not_applicable"
    if resource.gsm is not None:
        validate_fragment_schema(source)
        schema = "passed_10000_prefix_rows"
    digest = sha256_file(source)
    if source == resource.staging_path:
        os.replace(source, resource.destination)
        source = resource.destination
    return DownloadResult(
        name=resource.name,
        role=resource.role,
        url=resource.url,
        path=str(resource.destination.relative_to(ROOT)),
        gsm=resource.gsm,
        bytes=expected_bytes,
        sha256=digest,
        gzip_integrity=integrity,
        fragment_schema=schema,
    )


def download_resource(resource: Resource, timeout: int) -> DownloadResult:
    """Resolve source size, resume interruptions, and restart corrupt staging payloads."""
    expected = remote_size(resource.url, timeout)
    last_error: Optional[BaseException] = None
    for _ in range(3):
        try:
            return validate_and_promote(resource, expected, timeout)
        except (EOFError, gzip.BadGzipFile) as exc:
            last_error = exc
            if resource.destination.exists():
                raise
            resource.staging_path.unlink(missing_ok=True)
        except (OSError, subprocess.CalledProcessError, urllib.error.URLError) as exc:
            last_error = exc
    assert last_error is not None
    raise last_error


def atomic_yaml(path: Path, value: dict[str, Any]) -> None:
    """Write a YAML record using an adjacent temporary file and atomic rename."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(descriptor, "w") as handle:
            yaml.safe_dump(value, handle, sort_keys=False)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--scope", choices=("all", "series", "samples"), default="all")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=ROOT / "data/processed/shapemix/gse129785_immune/source_audit/downloads.yaml",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be positive.")
    config = load_config(args.config)
    resources = resources_from_config(config, args.scope)
    results: list[DownloadResult] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_resource = {
            executor.submit(download_resource, resource, args.timeout): resource
            for resource in resources
        }
        for future in concurrent.futures.as_completed(future_to_resource):
            resource = future_to_resource[future]
            result = future.result()
            results.append(result)
            print(f"validated {result.name} {result.bytes} {result.sha256}", flush=True)

    order = {resource.name: index for index, resource in enumerate(resources)}
    results.sort(key=lambda result: order[result.name])
    total_bytes = sum(result.bytes for result in results)
    audit = {
        "schema_version": 1,
        "source_dataset_id": config["source_dataset_id"],
        "accession": config["accession"],
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "config": str(args.config.resolve().relative_to(ROOT.resolve())),
        "scope": args.scope,
        "resources": len(results),
        "total_bytes": total_bytes,
        "files": [asdict(result) for result in results],
    }
    atomic_yaml(args.audit_output, audit)
    print(f"wrote {args.audit_output}", flush=True)
    print(f"validated_bytes {total_bytes}", flush=True)


if __name__ == "__main__":
    main()
