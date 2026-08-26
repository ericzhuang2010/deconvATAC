import gzip
import io
from pathlib import Path

import pytest

from scripts.preprocess_gse246791_fragment_reads import (
    audit_pair,
    iter_normalized_pairs,
    parse_rewritten_header,
)


BARCODE = "ACGTACGTACGTACGTACGTACGTACGTACGT"
NATIVE = f"{BARCODE}:7001113:869:H53K3BCX2:1:1106:1185:1857"


def record(ordinal: int, mate: int, *, native: str = NATIVE) -> bytes:
    return (
        f"@SRR26585986.{ordinal} {native}/{mate}\n"
        "ACGT\n"
        "+\n"
        "IIII\n"
    ).encode()


def test_parse_rewritten_header_restores_barcode_first_native_qname():
    parsed = parse_rewritten_header(
        f"@SRR26585986.17 {NATIVE}/1\n".encode(),
        expected_srr="SRR26585986",
        expected_mate=1,
    )
    assert parsed.ordinal == 17
    assert parsed.barcode == BARCODE
    assert parsed.native_qname == f"{NATIVE}/1".encode()


@pytest.mark.parametrize(
    ("header", "message"),
    [
        (f"@SRR26585986.1\n", "one native QNAME"),
        (f"@SRR1.1 {NATIVE}/1\n", "Unexpected run accession"),
        (f"@SRR26585986.1 {NATIVE}/2\n", "mate suffix"),
        (f"@SRR26585986.1 SHORT:instrument/1\n", "32-base"),
    ],
)
def test_parse_rewritten_header_fails_closed(header: str, message: str):
    with pytest.raises(ValueError, match=message):
        parse_rewritten_header(
            header.encode(), expected_srr="SRR26585986", expected_mate=1
        )


def test_iter_normalized_pairs_validates_mates_and_restores_headers():
    read1 = io.BytesIO(record(1, 1) + record(2, 1))
    read2 = io.BytesIO(record(1, 2) + record(2, 2))
    pairs = list(
        iter_normalized_pairs(read1, read2, expected_srr="SRR26585986")
    )
    assert len(pairs) == 2
    assert pairs[0][0][0] == f"@{NATIVE}/1\n".encode()
    assert pairs[0][1][0] == f"@{NATIVE}/2\n".encode()


def test_iter_normalized_pairs_rejects_ordinal_and_count_mismatch():
    with pytest.raises(ValueError, match="ordinal mismatch"):
        list(
            iter_normalized_pairs(
                io.BytesIO(record(1, 1)),
                io.BytesIO(record(2, 2)),
                expected_srr="SRR26585986",
            )
        )
    with pytest.raises(ValueError, match="different numbers"):
        list(
            iter_normalized_pairs(
                io.BytesIO(record(1, 1) + record(2, 1)),
                io.BytesIO(record(1, 2)),
                expected_srr="SRR26585986",
            )
        )


def test_audit_pair_reads_gzip_without_creating_normalized_fastqs(tmp_path: Path):
    read1 = tmp_path / "r1.fastq.gz"
    read2 = tmp_path / "r2.fastq.gz"
    with gzip.open(read1, "wb") as handle:
        handle.write(record(8, 1) + record(9, 1))
    with gzip.open(read2, "wb") as handle:
        handle.write(record(8, 2) + record(9, 2))
    audit = audit_pair(
        read1, read2, expected_srr="SRR26585986", limit=100
    )
    assert audit.read_pairs == 2
    assert audit.first_ordinal == 8
    assert audit.last_ordinal == 9
    assert audit.read1_length_min == audit.read1_length_max == 4
    assert audit.read2_length_min == audit.read2_length_max == 4
