import gzip
import os
from pathlib import Path

import pysam
import pytest

from scripts.align_shapemix_gse246791 import (
    alignment_commands,
    stream_bwa_name_sorted_bam,
    sample_pair,
    sort_and_validate_bam,
    validate_and_promote_bam,
    validate_cpu_affinity,
)
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


def test_cpu_affinity_is_fail_closed_at_two_cpus(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: {6, 7})
    assert validate_cpu_affinity() == [6, 7]
    monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: {5, 6, 7})
    with pytest.raises(RuntimeError, match="cap it to 2"):
        validate_cpu_affinity()


def test_alignment_commands_stream_one_thread_bwa_into_one_thread_name_sort(
    tmp_path: Path,
):
    bwa_command, sort_command = alignment_commands(
        reference=tmp_path / "mm10.fa",
        bwa=tmp_path / "bwa",
        sorter_python=tmp_path / "python",
        sorter_script=tmp_path / "sorter.py",
        read_fd1=11,
        read_fd2=12,
        bam_output=tmp_path / "name_sorted.bam.partial",
    )
    assert bwa_command[1:4] == ["mem", "-t", "1"]
    assert bwa_command[-2:] == ["/dev/fd/11", "/dev/fd/12"]
    assert sort_command == [
        str(tmp_path / "python"),
        str(tmp_path / "sorter.py"),
        "--output",
        str(tmp_path / "name_sorted.bam.partial"),
    ]
    assert "alignment.sam" not in " ".join(sort_command)


def test_streaming_alignment_pipe_materializes_only_name_sorted_bam(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    sorter_python = ROOT / ".venv-shapemix-fragments/bin/python"
    sorter_script = ROOT / "scripts/sort_shapemix_bam_stream.py"
    if not sorter_python.is_file():
        pytest.skip("Pinned pysam environment is not installed")
    monkeypatch.setenv("DECONVATAC_RESOURCE_GUARD", "1")
    read1 = tmp_path / "read1.fastq.gz"
    read2 = tmp_path / "read2.fastq.gz"
    native = f"{BARCODE}:instrument:lane:tile:x:y"
    with gzip.open(read1, "wb") as handle:
        handle.write(
            (
                f"@SRR1.1 {native}/1\n"
                "ACGT\n"
                "+\n"
                "IIII\n"
            ).encode()
        )
    with gzip.open(read2, "wb") as handle:
        handle.write(
            (
                f"@SRR1.1 {native}/2\n"
                "TGCA\n"
                "+\n"
                "IIII\n"
            ).encode()
        )
    fake_bwa = tmp_path / "fake_bwa"
    fake_bwa.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "for path in sys.argv[-2:]:\n"
        "    pathlib.Path(path).read_bytes()\n"
        f"query = '{BARCODE}:instrument:lane:tile:x:y'\n"
        "sys.stdout.write('@HD\\tVN:1.6\\tSO:unsorted\\n'\n"
        "                 '@SQ\\tSN:chr1\\tLN:1000\\n'\n"
        "                 f'{query}\\t99\\tchr1\\t11\\t60\\t4M\\t=\\t31\\t24\\tACGT\\tIIII\\n'\n"
        "                 f'{query}\\t147\\tchr1\\t31\\t60\\t4M\\t=\\t11\\t-24\\tTGCA\\tIIII\\n')\n"
    )
    fake_bwa.chmod(0o755)
    partial = tmp_path / "name_sorted.bam.partial"
    final = tmp_path / "name_sorted.bam"
    result = stream_bwa_name_sorted_bam(
        read1=read1,
        read2=read2,
        srr="SRR1",
        reference=tmp_path / "unused.fa",
        bwa=fake_bwa,
        sorter_python=sorter_python,
        sorter_script=sorter_script,
        bam_output=partial,
        bwa_stderr_output=tmp_path / "bwa.stderr.log",
        sort_stderr_output=tmp_path / "pysam.stderr.log",
    )
    validation = validate_and_promote_bam(partial, final)
    assert result["read_pairs"] == 1
    assert result["materialized_sam"] is False
    assert validation["alignments"] == 2
    assert final.is_file()
    assert not (tmp_path / "alignment.sam").exists()


def test_sample_pair_rejects_unknown_gsm():
    config = load_config(
        ROOT / "configs/data_sources/shapemix_gse246791_fragment_reads.yaml"
    )
    with pytest.raises(ValueError, match="exactly two"):
        sample_pair(config, "GSM_DOES_NOT_EXIST")
