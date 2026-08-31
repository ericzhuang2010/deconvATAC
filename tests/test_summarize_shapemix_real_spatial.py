from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scripts.summarize_shapemix_real_spatial import (
    marker_score,
    morans_i,
    neighbor_edges,
    safe_correlation,
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
