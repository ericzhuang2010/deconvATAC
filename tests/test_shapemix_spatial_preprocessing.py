import gzip
from pathlib import Path

import numpy as np

from scripts.preprocess_shapemix_spatial import (
    canonical_fragment_barcodes,
    classify_payload,
    expected_supplementary_basenames,
    parse_positions,
    read_dense_gene_by_pixel,
    safe_member_name,
)


def _config():
    return {
        "sample_groups": [
            {
                "atac": ["GSM1"],
                "rna": ["GSM2"],
                "protein": [],
                "validation_epigenome": ["GSM1"],
                "validation_epigenome_matrix": ["GSM3"],
            }
        ]
    }


def test_payload_classifier_never_mistakes_histone_for_atac():
    config = _config()
    assert classify_payload(Path("GSM1_section_ATAC.fragments.tsv.gz"), config) == "atac_fragments"
    assert classify_payload(Path("GSM1_section_H3K27me3.fragments.tsv.gz"), config) == "h3k27me3_fragments"
    assert classify_payload(Path("GSM2_section_RNA.tar.gz"), config) == "rna_matrix_bundle"
    assert classify_payload(Path("GSM2_section_spatial_RNA.tar.gz"), config) == "spatial"
    assert classify_payload(Path("GSM3_section_H3K27ac_matrix.tsv.gz"), config) == "h3k27ac_dense_matrix"


def test_safe_member_name_rejects_archive_escapes():
    assert safe_member_name("spatial/tissue_positions_list.csv")
    assert not safe_member_name("../escape")
    assert not safe_member_name("/absolute")


def test_fragment_barcode_suffix_mapping_is_explicit_and_injective():
    config = {
        "barcode_policy": {
            "fragment_terminal_suffix_to_strip": "-1",
            "require_suffix_on_every_fragment_barcode": True,
        }
    }

    assert canonical_fragment_barcodes({"pixel-a-1", "pixel-b-1"}, config) == {
        "pixel-a",
        "pixel-b",
    }


def test_supplementary_manifest_parser_deduplicates_subseries_metadata(tmp_path: Path):
    metadata = tmp_path / "series_metadata"
    metadata.mkdir()
    url_a = "ftp://ftp.ncbi.nlm.nih.gov/path/GSM1_file.tsv.gz"
    url_b = "ftp://ftp.ncbi.nlm.nih.gov/path/GSM2%20file.tar.gz"
    for accession, urls in (("GSE1", (url_a,)), ("GSE2", (url_a, url_b))):
        with gzip.open(metadata / f"{accession}_family.soft.gz", "wt") as handle:
            for index, url in enumerate(urls, start=1):
                handle.write(f"!Sample_supplementary_file_{index} = {url}\n")
    config = {
        "accession": "GSE1",
        "raw_directory": str(tmp_path),
        "series_metadata": [{"accession": "GSE1"}, {"accession": "GSE2"}],
    }

    assert expected_supplementary_basenames(config) == {
        "GSM1_file.tsv.gz",
        "GSM2 file.tar.gz",
    }


def test_positions_parser_accepts_integer_and_float_deposits(tmp_path: Path):
    path = tmp_path / "tissue_positions_list.csv"
    path.write_text("bc1,1,0,1,10,20\nbc2,1.0,2,3,30.0,40.0\n")
    frame = parse_positions(path)

    assert frame["barcode"].tolist() == ["bc1", "bc2"]
    assert frame["pixel_col_fullres"].tolist() == [20.0, 40.0]


def test_dense_gene_by_pixel_parser_builds_sparse_transposed_matrix(tmp_path: Path):
    path = tmp_path / "GSM2_matrix.tsv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("pixel_a\tpixel_b\n")
        handle.write("Gene1\t0\t2\n")
        handle.write("Gene1\t3\t0\n")
    matrix = read_dense_gene_by_pixel(path)

    assert matrix.shape == (2, 2)
    assert matrix.obs_names.tolist() == ["pixel_a", "pixel_b"]
    assert matrix.var_names.tolist() == ["Gene1", "Gene1-1"]
    np.testing.assert_array_equal(matrix.X.toarray(), [[0, 3], [2, 0]])
