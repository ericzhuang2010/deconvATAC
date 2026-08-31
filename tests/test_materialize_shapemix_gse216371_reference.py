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


def test_fragment_total_concordance_requires_exactly_one_global_convention(
    tmp_path,
) -> None:
    labels = pd.DataFrame(
        {"cell_id": ["cellA", "cellB"], "fragments": ["1", "2"]}
    )
    totals = tmp_path / "totals.tsv"
    totals.write_text(
        "cell_id\tbed_rows\tread_support_sum\n"
        "cellA\t1\t3\n"
        "cellB\t2\t4\n"
    )
    observed = materializer.validate_fragment_concordance(labels, totals)
    assert observed["passed"] is True
    assert observed["matching_convention"] == "bed_rows"
    assert observed["cells_compared"] == 2

    totals.write_text(
        "cell_id\tbed_rows\tread_support_sum\n"
        "cellA\t1\t1\n"
        "cellB\t2\t2\n"
    )
    with pytest.raises(ValueError, match="Exactly one"):
        materializer.validate_fragment_concordance(labels, totals)


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
