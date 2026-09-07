#!/usr/bin/env python3
"""Audit and normalize ENA FASTQ headers for GSE246791 fragment recovery.

ENA preserves the authors' original read name as the second whitespace-delimited
FASTQ header field while prepending an ``SRR.accession`` identifier.  BWA retains
only the first field as QNAME.  The original barcode-first name must therefore be
restored before alignment so SnapATAC2's production barcode regex can recover the
32-base cell barcode.

Raw FASTQs are never modified.  The reusable functions here are also used by the
streaming alignment stage, which writes normalized records directly to BWA.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import BinaryIO, Iterator


ROOT = Path(__file__).resolve().parents[1]
ACCESSION_RE = re.compile(rb"^@(SRR[0-9]+)\.([1-9][0-9]*)$")
BARCODE_RE = re.compile(rb"^[ACGT]{32}$")


@dataclass(frozen=True)
class NormalizedHeader:
    accession: str
    ordinal: int
    barcode: str
    native_qname: bytes


@dataclass(frozen=True)
class PairAudit:
    read_pairs: int
    first_ordinal: int | None
    last_ordinal: int | None
    barcode_length: int
    read1_length_min: int | None
    read1_length_max: int | None
    read2_length_min: int | None
    read2_length_max: int | None


def parse_rewritten_header(
    header: bytes, *, expected_srr: str, expected_mate: int
) -> NormalizedHeader:
    """Return the original barcode-first QNAME after strict ENA validation."""
    if expected_mate not in (1, 2):
        raise ValueError("expected_mate must be 1 or 2")
    line = header.rstrip(b"\r\n")
    fields = line.split()
    if len(fields) != 2:
        raise ValueError(
            "FASTQ header must have the ENA accession field and one native QNAME"
        )
    accession_match = ACCESSION_RE.fullmatch(fields[0])
    if accession_match is None:
        raise ValueError(f"Malformed ENA FASTQ identifier: {fields[0]!r}")
    observed_srr = accession_match.group(1).decode("ascii")
    if observed_srr != expected_srr:
        raise ValueError(
            f"Unexpected run accession {observed_srr}; expected {expected_srr}"
        )
    native_qname = fields[1]
    mate_suffix = f"/{expected_mate}".encode("ascii")
    if not native_qname.endswith(mate_suffix):
        raise ValueError(
            f"Native QNAME lacks expected mate suffix {mate_suffix!r}: "
            f"{native_qname!r}"
        )
    barcode = native_qname.split(b":", 1)[0]
    if BARCODE_RE.fullmatch(barcode) is None:
        raise ValueError(f"Native QNAME lacks an exact 32-base A/C/G/T barcode: {barcode!r}")
    return NormalizedHeader(
        accession=observed_srr,
        ordinal=int(accession_match.group(2)),
        barcode=barcode.decode("ascii"),
        native_qname=native_qname,
    )


def iter_fastq_records(handle: BinaryIO, *, label: str) -> Iterator[tuple[bytes, ...]]:
    record_number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        record_number += 1
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not sequence or not plus or not quality:
            raise ValueError(f"Truncated FASTQ record {record_number} in {label}")
        sequence_value = sequence.rstrip(b"\r\n")
        quality_value = quality.rstrip(b"\r\n")
        if not plus.startswith(b"+"):
            raise ValueError(f"Invalid plus line in FASTQ record {record_number} of {label}")
        if len(sequence_value) != len(quality_value):
            raise ValueError(
                f"Sequence/quality length mismatch in record {record_number} of {label}"
            )
        yield header, sequence, plus, quality


def iter_normalized_pairs(
    read1: BinaryIO,
    read2: BinaryIO,
    *,
    expected_srr: str,
    limit: int | None = None,
) -> Iterator[tuple[tuple[bytes, ...], tuple[bytes, ...], NormalizedHeader]]:
    """Yield validated paired records with barcode-first normalized headers."""
    iterator1 = iter_fastq_records(read1, label="read 1")
    iterator2 = iter_fastq_records(read2, label="read 2")
    previous_ordinal: int | None = None
    pair_number = 0
    while limit is None or pair_number < limit:
        record1 = next(iterator1, None)
        record2 = next(iterator2, None)
        if record1 is None and record2 is None:
            return
        if record1 is None or record2 is None:
            raise ValueError("Paired FASTQs contain different numbers of records")
        pair_number += 1
        header1 = parse_rewritten_header(
            record1[0], expected_srr=expected_srr, expected_mate=1
        )
        header2 = parse_rewritten_header(
            record2[0], expected_srr=expected_srr, expected_mate=2
        )
        if header1.ordinal != header2.ordinal:
            raise ValueError(
                f"Mate ordinal mismatch at pair {pair_number}: "
                f"{header1.ordinal} != {header2.ordinal}"
            )
        if header1.native_qname[:-2] != header2.native_qname[:-2]:
            raise ValueError(f"Native mate QNAME mismatch at pair {pair_number}")
        if header1.barcode != header2.barcode:
            raise ValueError(f"Mate barcode mismatch at pair {pair_number}")
        if previous_ordinal is not None and header1.ordinal != previous_ordinal + 1:
            raise ValueError(
                f"Nonconsecutive ENA ordinal at pair {pair_number}: "
                f"{header1.ordinal} after {previous_ordinal}"
            )
        previous_ordinal = header1.ordinal
        normalized1 = (b"@" + header1.native_qname + b"\n", *record1[1:])
        normalized2 = (b"@" + header2.native_qname + b"\n", *record2[1:])
        yield normalized1, normalized2, header1


def audit_pair(
    read1_path: Path,
    read2_path: Path,
    *,
    expected_srr: str,
    limit: int | None,
) -> PairAudit:
    count = 0
    first_ordinal: int | None = None
    last_ordinal: int | None = None
    read1_min: int | None = None
    read1_max: int | None = None
    read2_min: int | None = None
    read2_max: int | None = None
    with gzip.open(read1_path, "rb") as read1, gzip.open(read2_path, "rb") as read2:
        for record1, record2, header in iter_normalized_pairs(
            read1, read2, expected_srr=expected_srr, limit=limit
        ):
            count += 1
            first_ordinal = header.ordinal if first_ordinal is None else first_ordinal
            last_ordinal = header.ordinal
            length1 = len(record1[1].rstrip(b"\r\n"))
            length2 = len(record2[1].rstrip(b"\r\n"))
            read1_min = length1 if read1_min is None else min(read1_min, length1)
            read1_max = length1 if read1_max is None else max(read1_max, length1)
            read2_min = length2 if read2_min is None else min(read2_min, length2)
            read2_max = length2 if read2_max is None else max(read2_max, length2)
    return PairAudit(
        read_pairs=count,
        first_ordinal=first_ordinal,
        last_ordinal=last_ordinal,
        barcode_length=32,
        read1_length_min=read1_min,
        read1_length_max=read1_max,
        read2_length_min=read2_min,
        read2_length_max=read2_max,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--srr", required=True)
    parser.add_argument("--limit", type=int, default=100_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.limit < 1:
        parser.error("--limit must be at least 1")
    return args


def main() -> None:
    args = parse_args()
    audit = audit_pair(
        args.read1.resolve(),
        args.read2.resolve(),
        expected_srr=args.srr,
        limit=args.limit,
    )
    payload = {
        "schema_version": 1,
        "header_interpretation": (
            "restore the second ENA header field as barcode-first native QNAME"
        ),
        "expected_srr": args.srr,
        "limit": args.limit,
        **asdict(audit),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(serialized, end="")
    else:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_suffix(output.suffix + ".partial")
        temporary.write_text(serialized)
        temporary.replace(output)


if __name__ == "__main__":
    main()
