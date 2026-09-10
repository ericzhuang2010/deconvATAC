from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scripts.summarize_shapemix_real_spatial import (
    append_marker_evidence_rows,
    marker_score,
    morans_i,
    neighbor_edges,
    read_prediction,
    safe_correlation,
    spatial_signal_totals,
)


def test_real_spatial_neighbor_graph_is_undirected_and_deterministic() -> None:
    coordinates = np.asarray([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    edges = neighbor_edges(coordinates, neighbors=1)

    assert edges.tolist() == [[0, 1], [1, 2], [2, 3]]
    assert morans_i(np.asarray([0.0, 1.0, 2.0, 3.0]), edges) > 0


def test_real_spatial_marker_score_uses_library_normalized_present_features() -> None:
    matrix = ad.AnnData(
        X=sparse.csr_matrix([[10, 0, 5], [0, 10, 5], [5, 5, 5]], dtype=np.int64),
        obs=pd.DataFrame(index=["a", "b", "c"]),
        var=pd.DataFrame(index=["GeneA", "GeneB", "Other"]),
    )
    score, present, missing = marker_score(
        matrix,
        ["genea", "GENEB", "absent"],
        scale_factor=10_000.0,
        minimum_features=2,
    )

    assert present == ["genea", "GENEB"]
    assert missing == ["absent"]
    assert score.index.tolist() == ["a", "b", "c"]
    assert np.isfinite(score.to_numpy()).all()


def test_real_spatial_constant_map_correlation_is_explicitly_undefined() -> None:
    value, constant = safe_correlation(
        np.asarray([1.0, 1.0, 1.0]),
        np.asarray([0.0, 1.0, 2.0]),
        "spearman",
    )
    assert np.isnan(value)
    assert constant is True


def write_prediction(path: Path, values: list[list[float]]) -> Path:
    results = path / "results"
    results.mkdir(parents=True)
    pd.DataFrame(
        values,
        index=["informative", "zero_signal"],
        columns=["type_a", "type_b"],
    ).to_csv(results / "proportions.csv")
    return path


def test_real_spatial_signal_totals_identify_exact_zero_input() -> None:
    spatial = ad.AnnData(
        X=sparse.csr_matrix([[1, 2], [0, 0]], dtype=np.int64),
        obs=pd.DataFrame(index=["informative", "zero_signal"]),
        var=pd.DataFrame(index=["peak_a", "peak_b"]),
    )

    totals = spatial_signal_totals(spatial)

    assert totals.to_dict() == {"informative": 3.0, "zero_signal": 0.0}


def test_real_spatial_prediction_allows_zero_only_for_zero_signal_input(
    tmp_path: Path,
) -> None:
    run_dir = write_prediction(tmp_path / "allowed", [[0.25, 0.75], [0.0, 0.0]])

    value = read_prediction(
        run_dir,
        ["type_a", "type_b"],
        zero_signal_spots=["zero_signal"],
    )

    assert value.loc["zero_signal"].sum() == 0.0


def test_real_spatial_prediction_rejects_unexplained_zero_row(tmp_path: Path) -> None:
    run_dir = write_prediction(tmp_path / "rejected", [[0.25, 0.75], [0.0, 0.0]])

    with pytest.raises(ValueError, match="Invalid prediction proportions"):
        read_prediction(run_dir, ["type_a", "type_b"])


def test_real_spatial_prediction_rejects_nonunit_nonzero_row(tmp_path: Path) -> None:
    run_dir = write_prediction(tmp_path / "rejected", [[0.25, 0.25], [0.0, 0.0]])

    with pytest.raises(ValueError, match="Invalid prediction proportions"):
        read_prediction(
            run_dir,
            ["type_a", "type_b"],
            zero_signal_spots=["zero_signal"],
        )


def test_real_spatial_unsupported_marker_panel_is_an_explicit_unavailable_row() -> None:
    matrix = ad.AnnData(
        X=sparse.csr_matrix([[1], [2], [3]], dtype=np.int64),
        obs=pd.DataFrame(index=["a", "b", "c"]),
        var=pd.DataFrame(index=["GeneA"]),
    )
    predictions = {
        "nnls": pd.DataFrame(
            [[1.0], [1.0], [1.0]],
            index=["a", "b", "c"],
            columns=["type_a"],
        )
    }
    rows: list[dict[str, object]] = []

    append_marker_evidence_rows(
        rows,
        predictions,
        matrix,
        ["GeneA", "GeneB"],
        dataset_id="dataset",
        cell_type="type_a",
        evidence_type="rna_marker_score",
        evidence_gsm="GSM1",
        assay="rna",
        evidence_class="independent_orthogonal_validation",
        scale_factor=10_000.0,
        minimum_features=2,
    )

    assert len(rows) == 1
    assert rows[0]["status"] == "unavailable_insufficient_present_marker_features"
    assert rows[0]["aligned_spots"] == 0
    assert np.isnan(rows[0]["spearman_r"])
    assert rows[0]["present_markers"] == "GeneA"
    assert rows[0]["missing_markers"] == "GeneB"
