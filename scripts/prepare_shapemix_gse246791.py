#!/usr/bin/env python
"""Prepare audited source objects for the GSE246791 ShapeMix reference."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import platform
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import h5py
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs/data_sources/shapemix_gse246791_reference.yaml"
LOCK_PATH = ROOT / "configs/data_sources/shapemix_gse246791_reference_lock.yaml"
HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"
SOURCE_ROLE = "selected_sample_h5ad_with_raw_fragments"


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a mapping.")
    return value


CONFIG = _read_yaml(CONFIG_PATH)
LOCK = _read_yaml(LOCK_PATH)
RAW_ROOT = ROOT / str(CONFIG["raw_directory"])
FAMILY_ROOT = ROOT / str(CONFIG["processed_directory"])
SOURCE_OBJECT_ROOT = FAMILY_ROOT / "source_audit" / "source_objects"
MANIFEST_PATH = SOURCE_OBJECT_ROOT / "manifest.yaml"
SCHEMA_ROOT = SOURCE_OBJECT_ROOT / "schema"


@dataclass(frozen=True)
class SelectedSource:
    gsm: str
    sample: str
    major_region: str
    name: str
    source_path: Path
    expected_source_bytes: int
    source_sha256: str
    output_path: Path


def _repository_path(path: Path) -> str:
    return path.absolute().relative_to(ROOT.absolute()).as_posix()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_yaml(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
    )
    try:
        with os.fdopen(descriptor, "w") as handle:
            yaml.safe_dump(dict(value), handle, sort_keys=False)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def _git_revision() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _lock_by_name(lock: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    files = lock.get("files")
    if not isinstance(files, list):
        raise ValueError("The GSE246791 source lock has no file records.")
    result: dict[str, Mapping[str, Any]] = {}
    for record in files:
        if not isinstance(record, Mapping):
            raise ValueError("Every source-lock file record must be a mapping.")
        name = str(record["name"])
        if name in result:
            raise ValueError(f"Duplicate source-lock record: {name}")
        result[name] = record
    if len(result) != int(lock["resources"]):
        raise ValueError("The GSE246791 source lock is incomplete.")
    return result


def selected_sources(
    config: Mapping[str, Any] = CONFIG,
    lock: Mapping[str, Any] = LOCK,
) -> list[SelectedSource]:
    locked = _lock_by_name(lock)
    records: list[SelectedSource] = []
    for resource in config.get("resources", []):
        if resource.get("role") != SOURCE_ROLE:
            continue
        name = str(resource["name"])
        if not name.endswith(".h5ad.gz"):
            raise ValueError(f"Selected source is not an H5AD gzip: {name}")
        lock_record = locked.get(name)
        if lock_record is None:
            raise ValueError(f"Selected source is absent from the lock: {name}")
        if lock_record.get("role") != SOURCE_ROLE:
            raise ValueError(f"Selected source role disagrees with lock: {name}")
        configured_bytes = int(resource["expected_bytes"])
        if configured_bytes != int(lock_record["bytes"]):
            raise ValueError(f"Selected source byte count disagrees with lock: {name}")
        records.append(
            SelectedSource(
                gsm=str(resource["gsm"]),
                sample=str(resource["sample"]),
                major_region=str(resource["major_region"]),
                name=name,
                source_path=ROOT / str(lock_record["path"]),
                expected_source_bytes=configured_bytes,
                source_sha256=str(lock_record["sha256"]),
                output_path=SOURCE_OBJECT_ROOT / name.removesuffix(".gz"),
            )
        )
    expected_samples = int(config["selection"]["expected_samples"])
    expected_regions = [str(value) for value in config["selection"]["major_regions"]]
    if len(records) != expected_samples:
        raise ValueError(
            f"Expected {expected_samples} selected H5ADs, found {len(records)}."
        )
    if [record.major_region for record in records] != expected_regions:
        raise ValueError("Selected H5AD order does not match the frozen region order.")
    if len({record.gsm for record in records}) != len(records):
        raise ValueError("Selected GSE246791 GSM accessions must be unique.")
    return records


def _hdf5_validation(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        signature = handle.read(len(HDF5_SIGNATURE))
    if signature != HDF5_SIGNATURE:
        raise ValueError(f"Not an HDF5 file: {path}")
    with h5py.File(path, "r") as handle:
        root_keys = sorted(handle.keys())
        missing = sorted({"X", "obs", "var"}.difference(root_keys))
        if missing:
            raise ValueError(f"H5AD is missing root objects {missing}: {path}")
        root_attributes = {
            str(key): _attribute_value(value)
            for key, value in sorted(handle.attrs.items())
        }
    return {
        "hdf5_signature": "passed",
        "h5py_open": "passed",
        "required_h5ad_root_objects": "passed",
        "root_keys": root_keys,
        "root_attributes": root_attributes,
    }


def _attribute_value(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, list):
        return [_attribute_value(item) for item in value]
    return repr(value)


def _record_for_output(source: SelectedSource) -> dict[str, Any]:
    source_stat = source.source_path.stat()
    if source_stat.st_size != source.expected_source_bytes:
        raise ValueError(f"Source byte count changed: {source.source_path}")
    observed_source_sha256 = _sha256(source.source_path)
    if observed_source_sha256 != source.source_sha256:
        raise ValueError(f"Source SHA-256 changed: {source.source_path}")
    validation = _hdf5_validation(source.output_path)
    return {
        "gsm": source.gsm,
        "sample": source.sample,
        "major_region": source.major_region,
        "input": _repository_path(source.source_path),
        "input_bytes": source_stat.st_size,
        "input_sha256": source.source_sha256,
        "output": _repository_path(source.output_path),
        "output_bytes": source.output_path.stat().st_size,
        "output_sha256": _sha256(source.output_path),
        "validation": validation,
    }


def decompress_source(source: SelectedSource) -> dict[str, Any]:
    if source.output_path.is_file():
        print(f"gse246791 decompress gsm={source.gsm} status=validating_existing")
        return _record_for_output(source)
    source_stat = source.source_path.stat()
    if source_stat.st_size != source.expected_source_bytes:
        raise ValueError(f"Source byte count changed: {source.source_path}")
    observed_source_sha256 = _sha256(source.source_path)
    if observed_source_sha256 != source.source_sha256:
        raise ValueError(f"Source SHA-256 changed: {source.source_path}")
    temporary = Path(f"{source.output_path}.partial")
    if temporary.exists():
        raise FileExistsError(f"Unresolved partial source object: {temporary}")
    source.output_path.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    try:
        with gzip.open(source.source_path, "rb") as input_handle, temporary.open(
            "xb"
        ) as output_handle:
            for chunk in iter(lambda: input_handle.read(8 * 1024 * 1024), b""):
                output_handle.write(chunk)
                digest.update(chunk)
            output_handle.flush()
            os.fsync(output_handle.fileno())
        _hdf5_validation(temporary)
        os.replace(temporary, source.output_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    record = {
        "gsm": source.gsm,
        "sample": source.sample,
        "major_region": source.major_region,
        "input": _repository_path(source.source_path),
        "input_bytes": source_stat.st_size,
        "input_sha256": source.source_sha256,
        "output": _repository_path(source.output_path),
        "output_bytes": source.output_path.stat().st_size,
        "output_sha256": digest.hexdigest(),
        "validation": _hdf5_validation(source.output_path),
    }
    print(
        f"gse246791 decompress gsm={source.gsm} "
        f"bytes={record['output_bytes']} status=completed",
        flush=True,
    )
    return record


def _manifest(records: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    records = sorted(records, key=lambda value: str(value["major_region"]))
    expected = len(selected_sources())
    return {
        "schema_version": 1,
        "stage": "decompressed_author_h5ad_source_objects",
        "status": "complete" if len(records) == expected else "in_progress",
        "source_dataset_id": CONFIG["source_dataset_id"],
        "accession": CONFIG["accession"],
        "genome_build": CONFIG["genome_build"],
        "expected_objects": expected,
        "completed_objects": len(records),
        "parameters": {
            "operation": "gzip_decompression_without_content_transformation",
            "chunk_bytes": 8 * 1024 * 1024,
            "random_seed": None,
        },
        "inputs": {
            "config": _repository_path(CONFIG_PATH),
            "config_sha256": _sha256(CONFIG_PATH),
            "source_lock": _repository_path(LOCK_PATH),
            "source_lock_sha256": _sha256(LOCK_PATH),
        },
        "code": {
            "script": _repository_path(Path(__file__)),
            "script_sha256": _sha256(Path(__file__)),
            "git_revision": _git_revision(),
        },
        "software_versions": {
            "python": platform.python_version(),
            "h5py": h5py.__version__,
            "pyyaml": yaml.__version__,
        },
        "objects": list(records),
    }


def _existing_manifest_records() -> dict[str, Mapping[str, Any]]:
    if not MANIFEST_PATH.is_file():
        return {}
    manifest = _read_yaml(MANIFEST_PATH)
    values = manifest.get("objects", [])
    if not isinstance(values, list):
        raise ValueError("The processed source-object manifest is malformed.")
    return {str(value["gsm"]): value for value in values}


def decompress_sources(sources: Iterable[SelectedSource]) -> None:
    completed = _existing_manifest_records()
    for source in sources:
        completed[source.gsm] = decompress_source(source)
        _write_yaml(MANIFEST_PATH, _manifest(completed.values()))


def schema_snapshot(source: SelectedSource) -> dict[str, Any]:
    if not source.output_path.is_file():
        raise FileNotFoundError(f"Decompressed source object is absent: {source.output_path}")
    objects: list[dict[str, Any]] = []
    with h5py.File(source.output_path, "r") as handle:
        root_attributes = {
            str(key): _attribute_value(value)
            for key, value in sorted(handle.attrs.items())
        }

        def visit(name: str, value: h5py.Group | h5py.Dataset) -> None:
            record: dict[str, Any] = {
                "path": name,
                "kind": "dataset" if isinstance(value, h5py.Dataset) else "group",
            }
            if isinstance(value, h5py.Dataset):
                record["shape"] = list(value.shape)
                record["dtype"] = str(value.dtype)
                record["chunks"] = list(value.chunks) if value.chunks else None
                record["compression"] = value.compression
            if value.attrs:
                record["attributes"] = {
                    str(key): _attribute_value(attribute)
                    for key, attribute in sorted(value.attrs.items())
                }
            objects.append(record)

        handle.visititems(visit)
    return {
        "schema_version": 1,
        "gsm": source.gsm,
        "sample": source.sample,
        "major_region": source.major_region,
        "source_object": _repository_path(source.output_path),
        "source_object_sha256": _sha256(source.output_path),
        "root_attributes": root_attributes,
        "objects": objects,
    }


def inspect_sources(sources: Iterable[SelectedSource]) -> None:
    for source in sources:
        output = SCHEMA_ROOT / f"{source.gsm}.yaml"
        _write_yaml(output, schema_snapshot(source))
        print(
            f"gse246791 inspect gsm={source.gsm} output={_repository_path(output)} "
            "status=completed",
            flush=True,
        )


def _select_gsms(values: list[str] | None) -> list[SelectedSource]:
    sources = selected_sources()
    if not values:
        return sources
    requested = set(values)
    selected = [source for source in sources if source.gsm in requested]
    missing = sorted(requested.difference(source.gsm for source in selected))
    if missing:
        raise ValueError(f"Unknown selected GSM accession(s): {', '.join(missing)}")
    return selected


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("decompress", "inspect", "all"),
        default="all",
    )
    parser.add_argument(
        "--gsm",
        action="append",
        help="Restrict the stage to one or more selected GSM accessions.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    sources = _select_gsms(args.gsm)
    if args.stage in {"decompress", "all"}:
        decompress_sources(sources)
    if args.stage in {"inspect", "all"}:
        inspect_sources(sources)


if __name__ == "__main__":
    main()
