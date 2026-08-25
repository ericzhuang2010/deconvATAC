import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pysam

from scripts.preprocess_gse194122 import (
    feature_table,
    fragment_from_alignment,
    count_cut_sites,
    csr_selected_values,
    metrics_barcode_bridge,
    barcode_sequence,
    parse_barcode,
    parse_sample_key,
)


def test_processed_cell_identifier_parsing_is_explicit():
    assert parse_sample_key("AAACCCAAGGTCCAGA-1-s3d10") == ("s3d10", 3, 10)
    assert parse_sample_key("site3_donor10") is None
    assert parse_barcode("AAACCCAAGGTCCAGA-1-s3d10") == "AAACCCAAGGTCCAGA-1"
    assert parse_barcode("AAACCCAAGGTCCAGA-9-s3d6") == "AAACCCAAGGTCCAGA-9"
    assert barcode_sequence("AAACCCAAGGTCCAGA-9-s3d6") == "AAACCCAAGGTCCAGA"
    assert parse_barcode("AAACCCAAGGTCCAGA-s2d4") == "AAACCCAAGGTCCAGA"
    assert parse_barcode("not-a-cell") is None


def test_feature_table_requires_parseable_peak_coordinates():
    matrix = ad.AnnData(
        X=np.zeros((1, 2)),
        var=pd.DataFrame(
            {"feature_types": ["GEX", "ATAC"]},
            index=["ENSG000001", "chr1-100-200"],
        ),
    )

    result = feature_table(matrix)

    assert result["chromosome"].tolist() == ["", "chr1"]
    assert result["start"].tolist() == [pd.NA, 100]
    assert result["end"].tolist() == [pd.NA, 200]


def _alignment(*, duplicate=False, mapq=60, template_length=150):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    )
    record = pysam.AlignedSegment(header)
    record.query_name = "read1"
    record.query_sequence = "A" * 50
    record.flag = 0x1 | 0x2 | (0x400 if duplicate else 0)
    record.reference_id = 0
    record.reference_start = 100
    record.mapping_quality = mapq
    record.cigar = ((0, 50),)
    record.next_reference_id = 0
    record.next_reference_start = 200
    record.template_length = template_length
    record.set_tag("CB", "AAACCCAAGGTCCAGA-1")
    return record


def test_arc_pair_becomes_one_tn5_shifted_deduplicated_fragment():
    assert fragment_from_alignment(_alignment()) == (
        "chr1",
        104,
        245,
        "AAACCCAAGGTCCAGA-1",
        1,
    )
    assert fragment_from_alignment(_alignment(duplicate=True)) is None
    assert fragment_from_alignment(_alignment(mapq=20)) is None
    assert fragment_from_alignment(_alignment(template_length=-150)) is None


def test_metrics_bridge_maps_processed_sequence_to_common_fragment_barcode(tmp_path):
    cells = pd.DataFrame(
        {
            "cell_index": [0, 1],
            "cell_id": ["AAAAAAAAAAAAAAAA-9-s3d6", "CCCCCCCCCCCCCCCC-s3d6"],
            "sample_key": ["s3d6", "s3d6"],
            "barcode_sequence": ["AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC"],
        }
    )
    metrics = tmp_path / "metrics.csv"
    pd.DataFrame(
        {
            "barcode": ["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1", "GGGGGGGGGGGGGGGG-1"],
            "atac_barcode": ["TTTTTTTTTTTTTTTT-1", "ACACACACACACACAC-1", "CACACACACACACACA-1"],
            "is_cell": [1, 1, 0],
            "atac_fragments": [100, 200, 0],
            "atac_peak_region_cutsites": [50, 75, 0],
        }
    ).to_csv(metrics, index=False)

    bridged, summary = metrics_barcode_bridge(cells, metrics, "s3d6")

    assert bridged["fragment_barcode"].tolist() == [
        "AAAAAAAAAAAAAAAA-1",
        "CCCCCCCCCCCCCCCC-1",
    ]
    assert bridged["atac_library_barcode"].tolist() == [
        "TTTTTTTTTTTTTTTT-1",
        "ACACACACACACACAC-1",
    ]
    assert summary["exact_one_to_one_common_fragment_barcode_bridge"] is True


def test_sparse_raw_count_reader_selects_only_requested_rows_and_columns(tmp_path):
    path = tmp_path / "source.h5ad"
    with h5py.File(path, "w") as handle:
        group = handle.create_group("layers/counts")
        group.create_dataset("indptr", data=np.array([0, 2, 3], dtype=np.int32))
        group.create_dataset("indices", data=np.array([0, 3, 2], dtype=np.int32))
        group.create_dataset("data", data=np.array([5, 7, 11], dtype=np.float32))

    observed = csr_selected_values(path, rows=[1, 0], columns=[3, 2])

    np.testing.assert_array_equal(observed, [[0, 11], [7, 0]])


def test_right_endpoint_zero_is_distinguishable_from_minus_one(tmp_path):
    fragments = tmp_path / "fragments.tsv.gz"
    with pysam.BGZFile(str(fragments), "wb") as handle:
        handle.write(b"chr1\t100\t200\tAAAAAAAAAAAAAAAA-1\t1\n")
    pysam.tabix_index(str(fragments), preset="bed")
    index = tmp_path / "fragments.tsv.gz.tbi"
    cells = pd.DataFrame({"fragment_barcode": ["AAAAAAAAAAAAAAAA-1"]})
    peaks = pd.DataFrame(
        {
            "chromosome": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [150, 250],
        }
    )

    np.testing.assert_array_equal(
        count_cut_sites(fragments, index, cells, peaks, right_offset=0), [[1, 1]]
    )
    np.testing.assert_array_equal(
        count_cut_sites(fragments, index, cells, peaks, right_offset=-1), [[1, 0]]
    )
