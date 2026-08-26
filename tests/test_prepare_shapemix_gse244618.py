from pathlib import Path

import pytest

from scripts.prepare_shapemix_gse244618 import (
    CELL_TYPES,
    _sample_records,
    broad_label,
    parse_bedpe_record,
)


def test_gse244618_broad_ontology_is_frozen_in_output_order():
    assert CELL_TYPES == (
        "Inhibitory neurons",
        "Excitatory neurons",
        "Astrocytes",
        "Oligodendrocytes",
        "OPCs",
        "Microglia",
    )
    assert broad_label("GABA", "PVALB") == "Inhibitory neurons"
    assert broad_label("GLUT", "SUB") == "Excitatory neurons"
    assert broad_label("NonN", "ACBGM") == "Astrocytes"
    assert broad_label("NonN", "ASCNT") == "Astrocytes"
    assert broad_label("NonN", "ASCT") == "Astrocytes"
    assert broad_label("NonN", "OGC") == "Oligodendrocytes"
    assert broad_label("NonN", "OPC") == "OPCs"
    assert broad_label("NonN", "MGC") == "Microglia"
    assert broad_label("NonN", "EC") is None
    assert broad_label("NonN", "PER") is None
    assert broad_label("NonN", "SMC") is None


def test_bedpe_parser_uses_strand_aware_five_prime_cut_sites():
    forward_reverse = parse_bedpe_record(
        "chr1\t10\t20\tchr1\t50\t60\tBARCODE:read\t60\t+\t-\n"
    )
    reverse_forward = parse_bedpe_record(
        "chr1\t50\t60\tchr1\t10\t20\tBARCODE:read\t60\t-\t+\n"
    )

    assert forward_reverse == reverse_forward
    assert forward_reverse.chrom == "chr1"
    assert forward_reverse.start == 10
    assert forward_reverse.end == 60
    assert forward_reverse.length == 50
    assert forward_reverse.barcode == "BARCODE"
    assert forward_reverse.read_support == 1


@pytest.mark.parametrize(
    "line, message",
    [
        (
            "chr1\t10\t20\tchr2\t50\t60\tBC:read\t60\t+\t-\n",
            "same chromosome",
        ),
        (
            "chr1\t10\t20\tchr1\t50\t60\tBC:read\t60\t+\t+\n",
            "opposite strands",
        ),
        (
            "chr1\t20\t10\tchr1\t50\t60\tBC:read\t60\t+\t-\n",
            "invalid",
        ),
        (
            "chr1\t10\t20\tchr1\t0\t10\tBC:read\t60\t+\t-\n",
            "positive fragment",
        ),
    ],
)
def test_bedpe_parser_fails_closed_on_invalid_parent_fragments(line, message):
    with pytest.raises(ValueError, match=message):
        parse_bedpe_record(line)


def test_frozen_gse244618_subset_is_nine_unique_existing_bedpe_sources():
    records = _sample_records()

    assert len(records) == 9
    assert len({record["sample"] for record in records}) == 9
    assert {record["donor"] for record in records} == {1, 2, 4}
    assert {record["region"] for record in records} == {"HiT", "HiB", "Sub"}
    assert all(isinstance(record["path"], Path) and record["path"].is_file() for record in records)
