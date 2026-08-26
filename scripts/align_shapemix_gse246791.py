#!/usr/bin/env python3
"""Stream barcode-restored GSE246791 read pairs into one-thread BWA-MEM."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import queue
import re
import subprocess
import tempfile
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import pysam
import yaml

from scripts.preprocess_gse246791_fragment_reads import iter_normalized_pairs


ROOT = Path(__file__).resolve().parents[1]
CHUNK_BYTES = 8 * 1024 * 1024
BARCODE_QNAME_RE = re.compile(r"^[ACGT]{32}:")
QUEUE_RECORDS = 512


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


def file_digest(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sample_pair(config: Mapping[str, Any], gsm: str) -> dict[str, Any]:
    records = [
        record
        for record in config.get("resources", [])
        if isinstance(record, Mapping) and str(record.get("gsm")) == gsm
    ]
    if len(records) != 2:
        raise ValueError(f"Expected exactly two FASTQ resources for {gsm}; found {len(records)}")
    by_read = {int(record["read"]): record for record in records}
    if set(by_read) != {1, 2}:
        raise ValueError(f"Expected read 1 and read 2 resources for {gsm}")
    first, second = by_read[1], by_read[2]
    shared_fields = ("srr", "srx", "sample", "major_region")
    for field in shared_fields:
        if str(first[field]) != str(second[field]):
            raise ValueError(f"Mate metadata mismatch for {gsm}: {field}")
    raw_root = project_path(str(config["raw_directory"]))
    return {
        "gsm": gsm,
        "srr": str(first["srr"]),
        "srx": str(first["srx"]),
        "sample": str(first["sample"]),
        "major_region": str(first["major_region"]),
        "read1": raw_root / str(first["destination"]),
        "read2": raw_root / str(second["destination"]),
        "expected_bytes1": int(first["expected_bytes"]),
        "expected_bytes2": int(second["expected_bytes"]),
        "expected_md5_1": str(first["expected_md5"]),
        "expected_md5_2": str(second["expected_md5"]),
    }


def validate_inputs(sample: Mapping[str, Any], reference: Path, bwa: Path) -> None:
    for read_number in (1, 2):
        path = Path(sample[f"read{read_number}"])
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_size != int(sample[f"expected_bytes{read_number}"]):
            raise ValueError(f"FASTQ byte count changed: {path}")
    if not reference.is_file():
        raise FileNotFoundError(reference)
    if not bwa.is_file() or not os.access(bwa, os.X_OK):
        raise FileNotFoundError(f"BWA is missing or not executable: {bwa}")
    missing_index = [
        suffix
        for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")
        if not Path(f"{reference}{suffix}").is_file()
    ]
    if missing_index:
        raise FileNotFoundError(f"Reference lacks BWA index files: {missing_index}")


def _writer(
    file_descriptor: int,
    records: "queue.Queue[bytes | None]",
    errors: "queue.Queue[BaseException]",
) -> None:
    try:
        with os.fdopen(file_descriptor, "wb", buffering=1024 * 1024) as handle:
            while True:
                value = records.get()
                if value is None:
                    break
                handle.write(value)
    except BaseException as exc:
        errors.put(exc)


def _safe_put(
    records: "queue.Queue[bytes | None]",
    value: bytes | None,
    *,
    process: subprocess.Popen[bytes],
    writer: threading.Thread,
) -> None:
    while True:
        if not writer.is_alive():
            raise RuntimeError("FASTQ writer stopped before the input stream completed")
        if process.poll() is not None:
            raise RuntimeError(f"BWA stopped before the input stream completed: {process.returncode}")
        try:
            records.put(value, timeout=1.0)
            return
        except queue.Full:
            continue


def stream_bwa_sam(
    *,
    read1: Path,
    read2: Path,
    srr: str,
    reference: Path,
    bwa: Path,
    sam_output: Path,
    stderr_output: Path,
) -> dict[str, Any]:
    """Validate pairs and stream normalized records without materializing FASTQs."""
    read_fd1, write_fd1 = os.pipe()
    read_fd2, write_fd2 = os.pipe()
    queue1: "queue.Queue[bytes | None]" = queue.Queue(maxsize=QUEUE_RECORDS)
    queue2: "queue.Queue[bytes | None]" = queue.Queue(maxsize=QUEUE_RECORDS)
    errors: "queue.Queue[BaseException]" = queue.Queue()
    command = [
        str(bwa),
        "mem",
        "-t",
        "1",
        str(reference),
        f"/dev/fd/{read_fd1}",
        f"/dev/fd/{read_fd2}",
    ]
    pair_count = 0
    first_ordinal: int | None = None
    last_ordinal: int | None = None
    unique_h5_barcodes: set[str] = set()
    sam_output.parent.mkdir(parents=True, exist_ok=True)
    stderr_output.parent.mkdir(parents=True, exist_ok=True)
    if sam_output.exists():
        raise FileExistsError(f"Unresolved or existing SAM output: {sam_output}")
    with sam_output.open("xb") as sam_handle, stderr_output.open("xb") as stderr_handle:
        process = subprocess.Popen(
            command,
            stdout=sam_handle,
            stderr=stderr_handle,
            pass_fds=(read_fd1, read_fd2),
        )
        os.close(read_fd1)
        os.close(read_fd2)
        writer1 = threading.Thread(
            target=_writer, args=(write_fd1, queue1, errors), name="normalized-read1"
        )
        writer2 = threading.Thread(
            target=_writer, args=(write_fd2, queue2, errors), name="normalized-read2"
        )
        writer1.start()
        writer2.start()
        try:
            with gzip.open(read1, "rb") as input1, gzip.open(read2, "rb") as input2:
                for record1, record2, header in iter_normalized_pairs(
                    input1, input2, expected_srr=srr, limit=None
                ):
                    pair_count += 1
                    if first_ordinal is None:
                        first_ordinal = header.ordinal
                    last_ordinal = header.ordinal
                    if len(unique_h5_barcodes) < 1_000_000:
                        unique_h5_barcodes.add(header.barcode)
                    _safe_put(
                        queue1,
                        b"".join(record1),
                        process=process,
                        writer=writer1,
                    )
                    _safe_put(
                        queue2,
                        b"".join(record2),
                        process=process,
                        writer=writer2,
                    )
            _safe_put(queue1, None, process=process, writer=writer1)
            _safe_put(queue2, None, process=process, writer=writer2)
            writer1.join()
            writer2.join()
            return_code = process.wait()
            if not errors.empty():
                raise errors.get()
            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, command)
        except BaseException:
            process.terminate()
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            for records in (queue1, queue2):
                try:
                    records.put_nowait(None)
                except queue.Full:
                    pass
            writer1.join(timeout=30)
            writer2.join(timeout=30)
            raise
    if pair_count == 0:
        raise ValueError("No read pairs were streamed to BWA")
    return {
        "command": command,
        "threads": 1,
        "read_pairs": pair_count,
        "first_ordinal": first_ordinal,
        "last_ordinal": last_ordinal,
        "unique_barcodes_capped": len(unique_h5_barcodes),
        "unique_barcode_count_is_lower_bound": len(unique_h5_barcodes) == 1_000_000,
    }


def sort_and_validate_bam(sam: Path, bam: Path) -> dict[str, Any]:
    if bam.exists():
        raise FileExistsError(f"Final BAM already exists: {bam}")
    temporary = bam.with_suffix(bam.suffix + ".partial")
    if temporary.exists():
        raise FileExistsError(f"Unresolved BAM partial: {temporary}")
    bam.parent.mkdir(parents=True, exist_ok=True)
    pysam.sort("-n", "-o", str(temporary), str(sam))
    alignments = 0
    checked = 0
    with pysam.AlignmentFile(temporary, "rb") as handle:
        header = handle.header.to_dict()
        for record in handle.fetch(until_eof=True):
            alignments += 1
            if checked < 10_000:
                checked += 1
                if BARCODE_QNAME_RE.match(record.query_name or "") is None:
                    raise ValueError(f"BAM QNAME does not retain a 32-base barcode: {record.query_name!r}")
    if alignments == 0:
        raise ValueError("BWA produced an empty BAM")
    os.replace(temporary, bam)
    return {
        "alignments": alignments,
        "qnames_checked": checked,
        "barcode_first_qname_audit": "passed",
        "header_sort_order": header.get("HD", {}).get("SO"),
        "bytes": bam.stat().st_size,
        "sha256": file_digest(bam),
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
    parser.add_argument(
        "--reference",
        type=Path,
        default=ROOT / "data/processed/references/mm10_ucsc_initial/mm10.fa",
    )
    parser.add_argument(
        "--bwa",
        type=Path,
        default=ROOT / ".venv-shapemix-fragments/bin/bwa",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--audit-output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run alignment through scripts/run_shapemix_low_impact.sh")
    config = load_yaml(args.config)
    sample = sample_pair(config, args.gsm)
    reference = args.reference.absolute()
    bwa = args.bwa.absolute()
    validate_inputs(sample, reference, bwa)
    work_root = (
        ROOT
        / "data/work/preprocessing/gse246791_mouse_brain_reference"
        / args.gsm
        / "alignment"
    )
    output = args.output or work_root / "name_sorted.bam"
    audit_output = args.audit_output or (
        project_path(str(config["processed_directory"]))
        / "source_audit"
        / "alignments"
        / f"{args.gsm}.yaml"
    )
    if output.is_file() and audit_output.is_file():
        print(f"alignment gsm={args.gsm} status=reused", flush=True)
        return
    if output.exists() or audit_output.exists():
        raise FileExistsError(f"Partial immutable alignment state for {args.gsm}")
    sam = work_root / "alignment.sam.partial"
    stderr = work_root / "bwa.stderr.log"
    started = time.monotonic()
    stream = stream_bwa_sam(
        read1=Path(sample["read1"]),
        read2=Path(sample["read2"]),
        srr=str(sample["srr"]),
        reference=reference,
        bwa=bwa,
        sam_output=sam,
        stderr_output=stderr,
    )
    bam_validation = sort_and_validate_bam(sam, output)
    sam.unlink()
    record = {
        "schema_version": 1,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "gsm": args.gsm,
        "srr": sample["srr"],
        "sample": sample["sample"],
        "major_region": sample["major_region"],
        "inputs": {
            "read1": repository_path(Path(sample["read1"])),
            "read1_bytes": Path(sample["read1"]).stat().st_size,
            "read1_sha256": file_digest(Path(sample["read1"])),
            "read2": repository_path(Path(sample["read2"])),
            "read2_bytes": Path(sample["read2"]).stat().st_size,
            "read2_sha256": file_digest(Path(sample["read2"])),
            "reference": repository_path(reference),
            "reference_sha256": file_digest(reference),
        },
        "header_normalization": {
            "rule": "restore ENA second header field as barcode-first native QNAME",
            "materialized_normalized_fastqs": False,
            "validation": "every pair checked before streaming",
        },
        "alignment": stream,
        "bam": {
            "path": repository_path(output),
            **bam_validation,
        },
        "stderr_log": repository_path(stderr),
        "elapsed_seconds": time.monotonic() - started,
        "resource_policy": {
            "bwa_threads": 1,
            "pysam_sort_extra_threads": 0,
            "guard": "scripts/run_shapemix_low_impact.sh",
        },
    }
    atomic_yaml(audit_output, record)
    print(
        f"alignment gsm={args.gsm} pairs={stream['read_pairs']} "
        f"seconds={record['elapsed_seconds']:.1f} status=completed",
        flush=True,
    )


if __name__ == "__main__":
    main()
