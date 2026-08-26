from pathlib import Path

import pytest

from scripts.prepare_shapemix_mm10_reference import (
    PRIMARY_MM10,
    parse_chrom_sizes,
    parse_fasta_lengths,
    validate_reference_contract,
)


def test_parse_fasta_lengths_and_chrom_sizes(tmp_path: Path):
    fasta = tmp_path / "tiny.fa"
    fasta.write_bytes(b">chr1 description\nACGT\nNN\n>chr2\nAAA\n")
    sizes = tmp_path / "tiny.chrom.sizes"
    sizes.write_text("chr1\t6\nchr2\t3\n")
    assert parse_fasta_lengths(fasta) == {"chr1": 6, "chr2": 3}
    assert parse_chrom_sizes(sizes) == {"chr1": 6, "chr2": 3}


def test_validate_reference_contract_accepts_exact_mm10_primary_plus_scaffolds():
    lengths = {**PRIMARY_MM10, "chrM": 16299, "chrUn_GL456359": 22974}
    validate_reference_contract(lengths, lengths)


def test_validate_reference_contract_rejects_sequence_size_disagreement():
    fasta = {**PRIMARY_MM10}
    sizes = {**PRIMARY_MM10, "chrM": 16299}
    with pytest.raises(ValueError, match="FASTA/chrom-size mismatch"):
        validate_reference_contract(fasta, sizes)


def test_validate_reference_contract_rejects_primary_length_change():
    changed = {**PRIMARY_MM10, "chr1": PRIMARY_MM10["chr1"] - 1}
    with pytest.raises(ValueError, match="primary-contig"):
        validate_reference_contract(changed, changed)


def test_parsers_reject_duplicate_or_malformed_records(tmp_path: Path):
    duplicate_fasta = tmp_path / "duplicate.fa"
    duplicate_fasta.write_text(">chr1\nA\n>chr1\nC\n")
    with pytest.raises(ValueError, match="Duplicate FASTA"):
        parse_fasta_lengths(duplicate_fasta)

    malformed_sizes = tmp_path / "bad.sizes"
    malformed_sizes.write_text("chr1\n")
    with pytest.raises(ValueError, match="Malformed chromosome-size"):
        parse_chrom_sizes(malformed_sizes)


def test_primary_contract_exactly_matches_deposited_h5ad_lengths():
    assert PRIMARY_MM10 == {
        "chr1": 195471971, "chr2": 182113224, "chr3": 160039680,
        "chr4": 156508116, "chr5": 151834684, "chr6": 149736546,
        "chr7": 145441459, "chr8": 129401213, "chr9": 124595110,
        "chr10": 130694993, "chr11": 122082543, "chr12": 120129022,
        "chr13": 120421639, "chr14": 124902244, "chr15": 104043685,
        "chr16": 98207768, "chr17": 94987271, "chr18": 90702639,
        "chr19": 61431566, "chrX": 171031299, "chrY": 91744698,
    }
