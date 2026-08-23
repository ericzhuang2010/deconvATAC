import sys
import types

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from deconvatac.data import ordered_feature_sha256
from deconvatac.pp import (
    DEFAULT_FRAGMENT_LENGTH_BINS,
    FRAGMENT_SHAPE_LAYER_NAMES,
    FragmentLengthBin,
    FragmentParseError,
    FragmentRecord,
    PeakInterval,
    build_fragment_shape_anndata,
    build_peak_index,
    count_fragment_shapes,
    count_fragment_shapes_from_records,
    fragment_length_bin,
    parse_fragment_length_bins,
    parse_fragment_line,
    read_peaks_bed,
)
from deconvatac.pp.fragment_shapes import _ChunkedSparseAccumulator


def _matrix_values(result):
    return {name: result.layers[name].toarray() for name in FRAGMENT_SHAPE_LAYER_NAMES}


def test_fragment_length_boundaries_and_typed_bin_parser():
    assert [fragment_length_bin(length) for length in (99, 100, 249, 250)] == [0, 1, 1, 2]

    counted = count_fragment_shapes_from_records(
        [f"chr1\t10\t{10 + length}\tcell\t1" for length in (99, 100, 249, 250)],
        ["cell"],
        [PeakInterval("chr1", 0, 1000, "peak")],
        right_cut_offset=-1,
    )
    assert [int(counted.layers[name].sum()) for name in FRAGMENT_SHAPE_LAYER_NAMES] == [2, 4, 2]

    mappings = [
        {"name": bin_.name, "min_inclusive": bin_.min_inclusive, "max_exclusive": bin_.max_exclusive, "layer": bin_.layer}
        for bin_ in DEFAULT_FRAGMENT_LENGTH_BINS
    ]
    assert parse_fragment_length_bins(mappings) == DEFAULT_FRAGMENT_LENGTH_BINS
    assert isinstance(parse_fragment_length_bins(mappings)[0], FragmentLengthBin)

    with pytest.raises(ValueError, match="contiguous"):
        parse_fragment_length_bins(
            [
                {"name": "a", "min_inclusive": 0, "max_exclusive": 99, "layer": "a"},
                {"name": "b", "min_inclusive": 100, "max_exclusive": None, "layer": "b"},
            ]
        )


def test_parse_exact_five_column_schema_and_headers():
    record = parse_fragment_line("chr1\t5\t105\tAA-1\t17\n")
    assert record is not None
    assert (record.chrom, record.start, record.end, record.barcode, record.read_support) == (
        "chr1",
        5,
        105,
        "AA-1",
        17,
    )
    assert record.length == 100
    assert parse_fragment_line("# header") is None

    with pytest.raises(FragmentParseError) as missing:
        parse_fragment_line("chr1\t5\t105\tAA-1")
    assert missing.value.category == "schema"
    with pytest.raises(FragmentParseError) as coordinate:
        parse_fragment_line("chr1\tfive\t105\tAA-1\t1")
    assert coordinate.value.category == "coordinates"
    with pytest.raises(FragmentParseError) as support:
        parse_fragment_line("chr1\t5\t105\tAA-1\tmany")
    assert support.value.category == "schema"
    with pytest.raises(FragmentParseError, match="positive integer"):
        parse_fragment_line("chr1\t5\t105\tAA-1\t0")


def test_qc_counts_malformed_unknown_filtered_and_outside_rows():
    lines = [
        "# embedded header\n",
        "chr1\t10\t60\tknown\t3\n",  # two assigned cuts
        "chr1\t80\t180\tknown\t4\n",  # left assigned, right outside
        "chr1\t10\t60\tunknown\t5\n",
        "chrUn\t10\t60\tknown\t6\n",
        "chrUn\t10\t60\tunknown\t7\n",  # counted in both exclusion categories
        "chr1\tbad\t60\tknown\t8\n",
        "chr1\t60\t60\tknown\t9\n",
        "chr1\t10\t60\tknown\n",
        "chr1\t10\t60\tknown\tbad\n",
    ]
    result = count_fragment_shapes_from_records(
        lines,
        ["known"],
        [PeakInterval("chr1", 0, 100, "peak")],
        right_cut_offset=-1,
        chunk_size=2,
    )

    qc = result.qc
    assert qc.header_rows == 1
    assert qc.total_rows == 9
    assert qc.invalid_coordinate_rows == 2
    assert qc.invalid_schema_rows == 2
    assert qc.valid_rows == 5
    assert qc.read_support_total == 25
    assert qc.unknown_barcodes == 2
    assert qc.filtered_contigs == 2
    assert qc.retained_fragments == 2
    assert qc.fragments_with_assigned_cut_sites == 2
    assert qc.assigned_cut_sites == 3
    assert qc.cut_sites_outside_peaks == 1
    assert sum(qc.cut_sites_per_bin.values()) == qc.assigned_cut_sites
    assert int(result.X.sum()) == 3


def test_programmatic_fragment_records_receive_runtime_validation():
    result = count_fragment_shapes_from_records(
        [
            FragmentRecord("chr1", 10, 60, "known", 1),
            FragmentRecord("chr1", 10.5, 60, "known", 1),
            FragmentRecord("chr1", 10, 60, "known", 0),
        ],
        ["known"],
        [PeakInterval("chr1", 0, 100, "peak")],
        right_cut_offset=0,
    )

    assert result.qc.total_rows == 3
    assert result.qc.valid_rows == 1
    assert result.qc.invalid_coordinate_rows == 1
    assert result.qc.invalid_schema_rows == 1
    assert int(result.X.sum()) == 2


def test_read_support_is_qc_only_and_never_weights_counts():
    peaks = [PeakInterval("chr1", 0, 1000, "peak")]
    low = count_fragment_shapes_from_records(
        ["chr1\t10\t60\tcell\t1"], ["cell"], peaks, right_cut_offset=-1
    )
    high = count_fragment_shapes_from_records(
        ["chr1\t10\t60\tcell\t99"], ["cell"], peaks, right_cut_offset=-1
    )

    assert np.array_equal(low.X.toarray(), high.X.toarray())
    assert int(low.X.sum()) == int(high.X.sum()) == 2
    assert low.qc.read_support_total == 1
    assert high.qc.read_support_total == 99


def test_right_cut_coordinate_conventions_are_explicit():
    peaks = [
        PeakInterval("chr1", 0, 10, "left_peak"),
        PeakInterval("chr1", 10, 20, "right_peak"),
    ]
    line = ["chr1\t5\t10\tcell\t1"]
    end_minus_one = count_fragment_shapes_from_records(
        line, ["cell"], peaks, right_cut_offset=-1
    )
    end_coordinate = count_fragment_shapes_from_records(
        line, ["cell"], peaks, right_cut_offset=0
    )

    assert end_minus_one.X.toarray().tolist() == [[2, 0]]
    assert end_coordinate.X.toarray().tolist() == [[1, 1]]
    assert end_minus_one.fragment_shape_metadata()["right_cut_offset"] == -1
    assert end_coordinate.fragment_shape_metadata()["right_cut_offset"] == 0
    with pytest.raises(ValueError, match="0 or -1"):
        count_fragment_shapes_from_records(line, ["cell"], peaks, right_cut_offset=1)
    with pytest.raises(TypeError, match="integer"):
        count_fragment_shapes_from_records(line, ["cell"], peaks, right_cut_offset=None)


def test_chunk_parity_feature_order_and_determinism():
    # Deliberately rank p2 before p1; lookup is genomic, output stays ranked.
    peaks = [
        PeakInterval("chr1", 100, 200, "p2"),
        PeakInterval("chr1", 0, 100, "p1"),
    ]
    lines = [
        "chr1\t5\t105\tb\t1",
        "chr1\t7\t106\ta\t2",
        "chr1\t10\t110\ta\t3",
        "chr1\t20\t269\tb\t4",
        "chr1\t30\t280\tb\t5",
    ]
    tiny = count_fragment_shapes_from_records(
        iter(lines), ["b", "a"], peaks, right_cut_offset=-1, chunk_size=1
    )
    large = count_fragment_shapes_from_records(
        iter(lines), ["b", "a"], peaks, right_cut_offset=-1, chunk_size=100
    )
    repeated = count_fragment_shapes_from_records(
        iter(lines), ["b", "a"], peaks, right_cut_offset=-1, chunk_size=2
    )

    assert [peak.name for peak in tiny.peaks] == ["p2", "p1"]
    for name in FRAGMENT_SHAPE_LAYER_NAMES:
        assert sparse.isspmatrix_csr(tiny.layers[name])
        assert np.array_equal(tiny.layers[name].toarray(), large.layers[name].toarray())
        assert np.array_equal(tiny.layers[name].toarray(), repeated.layers[name].toarray())
        assert np.array_equal(tiny.layers[name].indptr, repeated.layers[name].indptr)
        assert np.array_equal(tiny.layers[name].indices, repeated.layers[name].indices)
        assert np.array_equal(tiny.layers[name].data, repeated.layers[name].data)
    assert tiny.qc.to_dict() == large.qc.to_dict() == repeated.qc.to_dict()


def test_log_structured_accumulator_many_tiny_flushes():
    shape = (31, 37)
    n_layers = 3
    events = [((index * 7) % shape[0], (index * 11) % shape[1], index % n_layers) for index in range(2049)]
    tiny = _ChunkedSparseAccumulator(shape, n_layers=n_layers, chunk_size=1)
    single_flush = _ChunkedSparseAccumulator(shape, n_layers=n_layers, chunk_size=len(events))

    event_counts = [0] * n_layers
    for row, column, layer in events:
        event_counts[layer] += 1
        tiny.add(row, column, layer)
        single_flush.add(row, column, layer)

    expected_levels = tuple(
        tuple(level for level in range(count.bit_length()) if count & (1 << level))
        for count in event_counts
    )
    assert tiny.occupied_levels == expected_levels
    assert all(len(levels) <= count.bit_length() for levels, count in zip(tiny.occupied_levels, event_counts))
    assert tiny.merge_count == sum(count - len(levels) for count, levels in zip(event_counts, expected_levels))
    assert single_flush.occupied_levels == ((0,), (0,), (0,))

    tiny_matrices = tiny.finish()
    single_flush_matrices = single_flush.finish()
    assert tiny.occupied_levels == ((), (), ())
    assert tiny.merge_count == sum(count - 1 for count in event_counts)
    for observed, expected in zip(tiny_matrices, single_flush_matrices):
        assert sparse.isspmatrix_csr(observed)
        assert np.array_equal(observed.indptr, expected.indptr)
        assert np.array_equal(observed.indices, expected.indices)
        assert np.array_equal(observed.data, expected.data)


def test_peak_validation_and_bed_order(tmp_path):
    bed = tmp_path / "ranked.bed"
    bed.write_text("chr1\t100\t200\tignored_name\nchr1\t0\t100\n")
    peaks = read_peaks_bed(bed)
    assert [peak.name for peak in peaks] == ["chr1:100-200", "chr1:0-100"]
    assert build_peak_index(peaks).assign("chr1", 50) == 1

    overlap = tmp_path / "overlap.bed"
    overlap.write_text("chr1\t0\t100\nchr1\t99\t150\n")
    with pytest.raises(ValueError, match="overlap"):
        read_peaks_bed(overlap)
    with pytest.raises(ValueError, match="unique"):
        build_peak_index(["chr1:0-100", "chr1:0-100"])


def test_anndata_layers_conservation_metadata_and_h5ad_round_trip(tmp_path):
    result = count_fragment_shapes_from_records(
        [
            "chr1\t10\t109\tcell_b\t2",
            "chr1\t200\t300\tcell_a\t3",
            "chr1\t400\t650\tcell_a\t4",
        ],
        ["cell_b", "cell_a"],
        [PeakInterval("chr1", 0, 1000, "chr1:0-1000")],
        right_cut_offset=0,
        chunk_size=2,
    )
    obs = pd.DataFrame({"cell_type": ["B", "A"]}, index=["cell_b", "cell_a"])
    adata = build_fragment_shape_anndata(
        result,
        obs=obs,
        provenance={
            "source_sha256": {"fragments": "a" * 64, "tabix_index": "b" * 64},
            "split_sha256": "c" * 64,
            "coordinate_validation": {
                "selected_right_cut_offset": 0,
                "matrix_match": "exact",
                "mismatched_entries": 0,
                "absolute_error": 0,
            },
            "software_versions": {"deconvatac": "0.0.1", "pysam": "0.24.0"},
        },
    )

    assert sparse.isspmatrix_csr(adata.X)
    layer_sum = sparse.csr_matrix(adata.shape, dtype=np.int64)
    for name in FRAGMENT_SHAPE_LAYER_NAMES:
        assert sparse.isspmatrix_csr(adata.layers[name])
        layer_sum = (layer_sum + adata.layers[name]).tocsr()
    assert (adata.X != layer_sum).nnz == 0
    metadata = adata.uns["fragment_shape"]
    assert metadata["read_support_policy"] == "ignore"
    assert metadata["right_cut_offset"] == 0
    expected_feature_sha256 = ordered_feature_sha256(["chr1:0-1000"])
    assert metadata["feature_sha256"] == expected_feature_sha256
    assert metadata["matrix_counters"] == {
        "assigned_cut_sites": 6,
        "cut_sites_per_bin.fragment_length_lt_100": 2,
        "cut_sites_per_bin.fragment_length_100_249": 2,
        "cut_sites_per_bin.fragment_length_ge_250": 2,
    }
    assert metadata["preprocessing_counters"]["assigned_cut_sites"] == 6
    assert list(metadata["bins"]) == ["0", "1", "2"]
    assert [metadata["bins"][str(index)]["name"] for index in range(3)] == ["short", "mono", "long"]
    assert "max_exclusive" not in metadata["bins"]["2"]

    path = tmp_path / "shape.h5ad"
    adata.write_h5ad(path)
    restored = ad.read_h5ad(path)
    assert sparse.isspmatrix_csr(restored.X)
    assert restored.uns["fragment_shape"]["right_cut_offset"] == 0
    assert restored.uns["fragment_shape"]["feature_sha256"] == expected_feature_sha256
    assert set(restored.uns["fragment_shape"]["bins"]) == {"0", "1", "2"}
    from deconvatac.data import (
        DeconvolutionInput,
        FragmentShapeSpec,
        validate_deconvolution_input,
        validate_fragment_shape_spec,
    )

    restored_spec = FragmentShapeSpec.from_mapping(restored.uns["fragment_shape"])
    validate_fragment_shape_spec(restored_spec)
    restored_sum = sum(
        (restored.layers[name] for name in FRAGMENT_SHAPE_LAYER_NAMES),
        sparse.csr_matrix(restored.shape, dtype=np.int64),
    ).tocsr()
    assert (restored.X != restored_sum).nnz == 0

    spatial = restored.copy()
    spatial.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 0.0]])
    shape_input = DeconvolutionInput(
        dataset_id="shape_builder_round_trip",
        modality="atac",
        feature_set="all",
        reference=restored,
        spatial=spatial,
        labels_key="cell_type",
        truth=pd.DataFrame(
            [[0.0, 1.0], [1.0, 0.0]],
            index=spatial.obs_names,
            columns=["A", "B"],
        ),
        fragment_shape=FragmentShapeSpec.from_mapping(restored.uns["fragment_shape"]),
        cell_types=["A", "B"],
    )
    validate_deconvolution_input(shape_input)

    with pytest.raises(ValueError, match="cannot override.*feature_sha256"):
        result.fragment_shape_metadata({"feature_sha256": "0" * 64})


def test_tabix_wrapper_counts_header_and_imports_pysam_lazily(monkeypatch):
    class FakeTabixFile:
        header = ("# first", "# second")

        def __init__(self, path):
            self.path = path

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        def fetch(self, *args):
            assert not args
            return iter(["chr1\t0\t99\tcell\t8"])

    fake_pysam = types.ModuleType("pysam")
    fake_pysam.TabixFile = FakeTabixFile
    monkeypatch.setitem(sys.modules, "pysam", fake_pysam)

    result = count_fragment_shapes(
        "unused.tsv.gz",
        ["cell"],
        [PeakInterval("chr1", 0, 200, "peak")],
        right_cut_offset=-1,
    )
    assert result.qc.header_rows == 2
    assert result.qc.total_rows == 1
    assert result.qc.read_support_total == 8
    assert int(result.X.sum()) == 2
