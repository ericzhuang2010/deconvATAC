from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import scripts.materialize_shapemix_gse216371_reference as materializer


def test_embryo_ontology_is_fail_closed_until_explicitly_frozen(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(materializer, "CONFIG", {})
    with pytest.raises(ValueError, match="must be frozen"):
        materializer.ontology()


def test_embryo_ontology_preserves_declared_output_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        materializer,
        "CONFIG",
        {
            "broad_ontology": {
                "ordered_cell_types": ["Neural", "Mesenchymal"],
                "author_main_cluster_mapping": {
                    "Brain": "Neural",
                    "Mesenchyme": "Mesenchymal",
                },
                "minimum_cells_per_type": 20,
            }
        },
    )
    assert materializer.ontology() == (
        ("Neural", "Mesenchymal"),
        {"Brain": "Neural", "Mesenchyme": "Mesenchymal"},
        20,
    )


def test_embryo_ranker_uses_peak_id_for_final_tie(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = {
        0: {"peak_id": "chr1:0-500"},
        1: {"peak_id": "chr2:0-500"},
        2: {"peak_id": "chr1:500-1000"},
        3: {"peak_id": "chr3:0-500"},
    }
    monkeypatch.setattr(
        materializer,
        "load_candidate_rows",
        lambda indices: {int(index): rows[int(index)] for index in indices},
    )
    selected = materializer.rank_features(
        np.asarray([3.0, 2.0, 2.0, 1.0]),
        np.asarray([10, 10, 10, 10]),
        np.asarray([5, 5, 5, 5]),
        n_top=2,
    )
    assert selected == [0, 2]


def test_embryo_ranker_prefers_higher_unsigned_total_before_identifier(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = {
        0: {"peak_id": "chr1:0-500"},
        1: {"peak_id": "chr1:500-1000"},
        2: {"peak_id": "chr2:0-500"},
    }
    monkeypatch.setattr(
        materializer,
        "load_candidate_rows",
        lambda indices: {int(index): rows[int(index)] for index in indices},
    )
    selected = materializer.rank_features(
        np.asarray([3.0, 2.0, 2.0]),
        np.asarray([10, 10, 10]),
        np.asarray([1, 5, 10], dtype=np.uint64),
        n_top=2,
    )
    assert selected == [0, 2]


def test_canonical_mm10_contigs_rejects_config_drift(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        materializer,
        "CONFIG",
        {"preprocessing_policy": {"canonical_contigs": ["chr1", "chrX"]}},
    )
    with pytest.raises(ValueError, match="canonical_contigs"):
        materializer.canonical_mm10_contigs()


def test_fragment_total_concordance_requires_exactly_one_global_convention(
    tmp_path,
) -> None:
    labels = pd.DataFrame(
        {"cell_id": ["cellA", "cellB"], "fragments": ["1", "2"]}
    )
    totals = tmp_path / "totals.tsv"
    totals.write_text(
        "cell_id\tbed_rows\tread_support_sum\t"
        "excluded_canonical_invalid_coordinate_rows\t"
        "excluded_canonical_invalid_coordinate_read_support\t"
        "excluded_noncanonical_contig_rows\t"
        "excluded_noncanonical_contig_read_support\n"
        "cellA\t1\t3\t0\t0\t0\t0\n"
        "cellB\t2\t4\t0\t0\t0\t0\n"
    )
    observed = materializer.validate_fragment_concordance(labels, totals)
    assert observed["passed"] is True
    assert observed["matching_convention"] == "canonical_bed_rows"
    assert observed["cells_compared"] == 2

    totals.write_text(
        "cell_id\tbed_rows\tread_support_sum\t"
        "excluded_canonical_invalid_coordinate_rows\t"
        "excluded_canonical_invalid_coordinate_read_support\t"
        "excluded_noncanonical_contig_rows\t"
        "excluded_noncanonical_contig_read_support\n"
        "cellA\t1\t1\t0\t0\t0\t0\n"
        "cellB\t2\t2\t0\t0\t0\t0\n"
    )
    with pytest.raises(ValueError, match="Exactly one"):
        materializer.validate_fragment_concordance(labels, totals)


def test_fragment_total_concordance_includes_excluded_boundary_rows(
    tmp_path,
) -> None:
    labels = pd.DataFrame(
        {"cell_id": ["cellA", "cellB"], "fragments": ["2", "2"]}
    )
    totals = tmp_path / "totals.tsv"
    totals.write_text(
        "cell_id\tbed_rows\tread_support_sum\t"
        "excluded_canonical_invalid_coordinate_rows\t"
        "excluded_canonical_invalid_coordinate_read_support\t"
        "excluded_noncanonical_contig_rows\t"
        "excluded_noncanonical_contig_read_support\n"
        "cellA\t1\t5\t1\t2\t0\t0\n"
        "cellB\t2\t8\t0\t0\t0\t0\n"
    )
    observed = materializer.validate_fragment_concordance(labels, totals)
    assert observed["matching_convention"] == "canonical_bed_rows"
    assert observed["excluded_canonical_invalid_coordinate_rows"] == 1
    assert observed["total_canonical_source_bed_rows"] == 4
    assert observed["total_source_bed_rows"] == 4


def test_fragment_total_concordance_excludes_noncanonical_rows(
    tmp_path,
) -> None:
    labels = pd.DataFrame(
        {"cell_id": ["cellA", "cellB"], "fragments": ["1", "2"]}
    )
    totals = tmp_path / "totals.tsv"
    totals.write_text(
        "cell_id\tbed_rows\tread_support_sum\t"
        "excluded_canonical_invalid_coordinate_rows\t"
        "excluded_canonical_invalid_coordinate_read_support\t"
        "excluded_noncanonical_contig_rows\t"
        "excluded_noncanonical_contig_read_support\n"
        "cellA\t1\t4\t0\t0\t3\t9\n"
        "cellB\t2\t5\t0\t0\t0\t0\n"
    )
    observed = materializer.validate_fragment_concordance(labels, totals)
    assert observed["matching_convention"] == "canonical_bed_rows"
    assert observed["excluded_noncanonical_contig_rows"] == 3
    assert observed["total_canonical_source_bed_rows"] == 3
    assert observed["total_source_bed_rows"] == 6


def test_completed_failed_scan_is_preserved_and_incomplete_scan_is_removed(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    work_root = tmp_path / "work"
    monkeypatch.setattr(materializer, "WORK_ROOT", work_root)

    completed = tmp_path / ".major_types_v1.completed"
    completed.mkdir()
    (completed / "stream_summary.tsv").write_text("key\tvalue\nmode\tstatistics\n")
    preserved = materializer.preserve_failed_fragment_statistics(
        completed, ValueError("concordance gate failed")
    )
    assert preserved is not None
    assert preserved.parent == work_root / "failed_fragment_statistics"
    assert not completed.exists()
    assert (preserved / "failure.txt").read_text() == (
        "ValueError: concordance gate failed\n"
    )

    incomplete = tmp_path / ".major_types_v1.incomplete"
    incomplete.mkdir()
    assert (
        materializer.preserve_failed_fragment_statistics(
            incomplete, RuntimeError("stream interrupted")
        )
        is None
    )
    assert not incomplete.exists()


def test_event_layers_aggregate_duplicates_in_bounded_chunks(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    events = np.asarray(
        [
            (0, 0, 0),
            (0, 0, 0),
            (1, 1, 1),
            (1, 1, 2),
            (1, 1, 2),
            (1, 1, 2),
        ],
        dtype=materializer.EVENT_DTYPE,
    )
    path = tmp_path / "events.bin"
    events.tofile(path)
    monkeypatch.setattr(materializer, "EVENT_CHUNK", 2)
    layers = materializer.event_layers(path, (2, 2))
    np.testing.assert_array_equal(
        layers["fragment_length_lt_100"].toarray(), [[2, 0], [0, 0]]
    )
    np.testing.assert_array_equal(
        layers["fragment_length_100_249"].toarray(), [[0, 0], [0, 1]]
    )
    np.testing.assert_array_equal(
        layers["fragment_length_ge_250"].toarray(), [[0, 0], [0, 3]]
    )
