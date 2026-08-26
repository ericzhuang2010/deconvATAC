import gzip
from pathlib import Path

import h5py
import numpy as np
import pytest

from scripts.construct_shapemix_gse246791_fragments import (
    filter_fragments,
    fragment_count_concordance,
    load_h5ad_fragment_counts,
)


BARCODE1 = "ACGTACGTACGTACGTACGTACGTACGTACGT"
BARCODE2 = "TGCATGCATGCATGCATGCATGCATGCATGCA"
THRESHOLDS = {
    "require_all_h5ad_barcodes_present": True,
    "median_absolute_relative_error_max": 0.02,
    "p95_absolute_relative_error_max": 0.10,
    "spearman_r_min": 0.98,
}


def write_h5ad_contract(path: Path) -> None:
    with h5py.File(path, "w") as handle:
        obs = handle.create_group("obs")
        obs.create_dataset("index", data=np.asarray([BARCODE1, BARCODE2], dtype="S32"))
        obs.create_dataset("n_fragment", data=np.asarray([2, 1], dtype=np.uint64))
        uns = handle.create_group("uns")
        reference = uns.create_group("reference_sequences")
        reference.create_dataset("reference_seq_name", data=np.asarray(["chr1"], dtype="S4"))
        reference.create_dataset("reference_seq_length", data=np.asarray([100], dtype=np.uint64))


def test_load_h5ad_fragment_counts_reads_exact_contract(tmp_path: Path):
    path = tmp_path / "source.h5ad"
    write_h5ad_contract(path)
    barcodes, counts, chrom_sizes = load_h5ad_fragment_counts(path)
    assert barcodes == [BARCODE1, BARCODE2]
    np.testing.assert_array_equal(counts, [2, 1])
    assert chrom_sizes == {"chr1": 100}


def test_filter_fragments_keeps_whitelist_and_h5ad_contigs(tmp_path: Path):
    source = tmp_path / "raw.tsv.gz"
    output = tmp_path / "final.tsv.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(f"chr1\t1\t10\t{BARCODE1}\t1\n")
        handle.write(f"chr1\t20\t30\t{BARCODE1}\t2\n")
        handle.write(f"chr1\t40\t50\t{BARCODE2}\t1\n")
        handle.write(f"chrM\t1\t10\t{BARCODE1}\t1\n")
        handle.write("chr1\t1\t10\tUNKNOWN\t1\n")
    result = filter_fragments(
        source,
        output,
        barcodes=[BARCODE1, BARCODE2],
        chrom_sizes={"chr1": 100},
    )
    np.testing.assert_array_equal(result.pop("observed_counts"), [2, 1])
    assert result["raw_rows"] == 5
    assert result["kept_rows"] == 3
    assert result["excluded_contig_rows"] == 1
    assert result["unknown_barcode_rows"] == 1
    with gzip.open(output, "rt") as handle:
        assert len(handle.readlines()) == 3


def test_filter_fragments_accepts_plain_compiled_backend_output(tmp_path: Path):
    source = tmp_path / "backend-output.tsv.gz"
    output = tmp_path / "standardized.tsv.gz"
    source.write_text(f"chr1\t1\t10\t{BARCODE1}\t1\n")
    result = filter_fragments(
        source,
        output,
        barcodes=[BARCODE1],
        chrom_sizes={"chr1": 100},
    )
    np.testing.assert_array_equal(result.pop("observed_counts"), [1])
    assert result["kept_rows"] == 1
    with gzip.open(output, "rt") as handle:
        assert handle.read() == f"chr1\t1\t10\t{BARCODE1}\t1\n"


def test_fragment_count_concordance_passes_exact_and_fails_missing():
    exact = fragment_count_concordance(
        np.asarray([100, 200, 300]),
        np.asarray([100, 200, 300]),
        thresholds=THRESHOLDS,
    )
    assert exact["passed"]
    missing = fragment_count_concordance(
        np.asarray([100, 200, 300]),
        np.asarray([100, 0, 300]),
        thresholds=THRESHOLDS,
    )
    assert not missing["passed"]
    assert not missing["gates"]["all_h5ad_barcodes_present"]


@pytest.mark.parametrize(
    "line",
    [
        f"chr1\t1\t10\t{BARCODE1}\n",
        f"chr1\t10\t1\t{BARCODE1}\t1\n",
        f"chr1\t1\t101\t{BARCODE1}\t1\n",
        f"chr1\t1\t10\t{BARCODE1}\t0\n",
    ],
)
def test_filter_fragments_rejects_malformed_selected_rows(tmp_path: Path, line: str):
    source = tmp_path / "raw.tsv.gz"
    output = tmp_path / "final.tsv.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(line)
    with pytest.raises(ValueError):
        filter_fragments(
            source,
            output,
            barcodes=[BARCODE1],
            chrom_sizes={"chr1": 100},
        )
