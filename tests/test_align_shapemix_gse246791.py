from pathlib import Path

import pysam
import pytest

from scripts.align_shapemix_gse246791 import sample_pair, sort_and_validate_bam
from scripts.download_shapemix_spatial import load_config


ROOT = Path(__file__).resolve().parents[1]
BARCODE = "ACGTACGTACGTACGTACGTACGTACGTACGT"


def test_representative_sample_pair_resolves_immutable_raw_paths():
    config = load_config(
        ROOT / "configs/data_sources/shapemix_gse246791_fragment_reads.yaml"
    )
    sample = sample_pair(config, "GSM7877011")
    assert sample["srr"] == "SRR26585986"
    assert sample["major_region"] == "STR"
    assert str(sample["read1"]).endswith(
        "data/raw/sources/ncbi_sra/GSE246791/samples/GSM7877011/raw_reads/SRR26585986_1.fastq.gz"
    )
    assert str(sample["read2"]).endswith("SRR26585986_2.fastq.gz")


def write_sam(path: Path, query_name: str) -> None:
    path.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "@SQ\tSN:chr1\tLN:1000\n"
        f"{query_name}\t99\tchr1\t11\t60\t4M\t=\t31\t24\tACGT\tIIII\n"
        f"{query_name}\t147\tchr1\t31\t60\t4M\t=\t11\t-24\tTGCA\tIIII\n"
    )


def test_sort_and_validate_bam_accepts_barcode_first_qnames(tmp_path: Path):
    sam = tmp_path / "input.sam"
    bam = tmp_path / "output.bam"
    write_sam(sam, f"{BARCODE}:instrument:1")
    result = sort_and_validate_bam(sam, bam)
    assert result["alignments"] == 2
    assert result["barcode_first_qname_audit"] == "passed"
    with pysam.AlignmentFile(bam, "rb") as handle:
        assert next(handle.fetch(until_eof=True)).query_name.startswith(BARCODE)


def test_sort_and_validate_bam_rejects_barcode_lost_qnames(tmp_path: Path):
    sam = tmp_path / "input.sam"
    bam = tmp_path / "output.bam"
    write_sam(sam, "SRR26585986.1")
    with pytest.raises(ValueError, match="does not retain"):
        sort_and_validate_bam(sam, bam)


def test_sample_pair_rejects_unknown_gsm():
    config = load_config(
        ROOT / "configs/data_sources/shapemix_gse246791_fragment_reads.yaml"
    )
    with pytest.raises(ValueError, match="exactly two"):
        sample_pair(config, "GSM_DOES_NOT_EXIST")
