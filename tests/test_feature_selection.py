import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from deconvatac.data import ordered_feature_sha256
from deconvatac.pp import highly_accessible_peaks, highly_variable_peaks, select_reference_peaks


def _reference_adata(values, labels, peak_ids):
    return ad.AnnData(
        X=sparse.csr_matrix(np.asarray(values, dtype=np.int64)),
        obs=pd.DataFrame({"cell_type": labels}),
        var=pd.DataFrame(index=peak_ids),
    )


def test_highly_variable_peaks_adds_boolean_var_column():
    adata = ad.AnnData(
        X=sparse.csr_matrix(
            np.array(
                [
                    [10, 0, 1],
                    [8, 0, 1],
                    [0, 9, 1],
                    [0, 7, 1],
                ],
                dtype=float,
            )
        ),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}),
        var=pd.DataFrame(index=["peak1", "peak2", "peak3"]),
    )

    highly_variable_peaks(adata, cluster_key="cell_type", n_top_features=2)

    assert "highly_variable" in adata.var
    assert adata.var["highly_variable"].dtype == bool
    assert int(adata.var["highly_variable"].sum()) == 2


def test_highly_accessible_peaks_adds_boolean_var_column():
    adata = ad.AnnData(
        X=sparse.csr_matrix(
            np.array(
                [
                    [1, 0, 1],
                    [1, 0, 0],
                    [1, 1, 0],
                ],
                dtype=float,
            )
        ),
        var=pd.DataFrame(index=["peak1", "peak2", "peak3"]),
    )

    highly_accessible_peaks(adata, n_top_features=1)

    assert "highly_accessible" in adata.var
    assert adata.var["highly_accessible"].dtype == bool
    assert list(adata.var_names[adata.var["highly_accessible"]]) == ["peak1"]


def test_reference_peak_selector_uses_exact_protocol_score_and_audit_values():
    adata = _reference_adata(
        [
            [4, 0, 1],
            [2, 1, 0],
            [0, 3, 1],
            [1, 3, 0],
        ],
        ["A", "A", "B", "B"],
        ["peak1", "peak2", "peak3"],
    )

    result = select_reference_peaks(adata, "cell_type", n_top_peaks=3, min_reference_cells=1)

    group_counts = np.array([[6, 1, 1], [1, 6, 1]], dtype=float)
    expected_log = np.log2(1.0 + 1e4 * group_counts / group_counts.sum(axis=1, keepdims=True))
    expected_scores = np.var(expected_log, axis=0, ddof=0)
    np.testing.assert_allclose(result.candidate_scores, expected_scores, rtol=0, atol=0)
    np.testing.assert_array_equal(result.candidate_nonzero_reference_cells, [3, 3, 2])
    np.testing.assert_array_equal(result.candidate_total_reference_counts, [7, 7, 2])
    assert result.peak_ids == ("peak1", "peak2", "peak3")
    np.testing.assert_array_equal(result.indices, [0, 1, 2])
    assert result.candidate_feature_sha256 == ordered_feature_sha256(adata.var_names)
    assert result.selected_feature_sha256 == ordered_feature_sha256(result.peak_ids)

    audit = result.to_frame()
    assert audit.columns.tolist() == [
        "rank",
        "peak_id",
        "original_index",
        "score",
        "nonzero_reference_cells",
        "total_reference_count",
    ]
    assert audit["rank"].tolist() == [1, 2, 3]


def test_reference_peak_selector_enforces_minimum_cell_coverage():
    adata = _reference_adata(
        [[1, 0], [0, 1], [1, 0]],
        ["A", "A", "B"],
        ["eligible", "one_cell_only"],
    )

    result = select_reference_peaks(adata, "cell_type", n_top_peaks=1, min_reference_cells=2)
    assert result.peak_ids == ("eligible",)
    np.testing.assert_array_equal(result.eligible_mask, [True, False])

    with pytest.raises(ValueError, match=r"Only 1 peaks.*2 are required"):
        select_reference_peaks(adata, "cell_type", n_top_peaks=2, min_reference_cells=2)


def test_reference_peak_selector_applies_complete_stable_tie_key():
    # With one type every population-variance score is exactly zero.  This
    # isolates coverage, total-count, and finally bytewise peak-ID ordering.
    adata = _reference_adata(
        [
            [1, 1, 1, 3, 1],
            [0, 0, 0, 3, 1],
            [0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0],
        ],
        ["A", "A", "A", "A"],
        ["peak_z", "peak_b", "peak_a", "peak_total", "peak_coverage"],
    )

    result = select_reference_peaks(adata, "cell_type", n_top_peaks=5, min_reference_cells=1)

    assert result.peak_ids == ("peak_coverage", "peak_total", "peak_a", "peak_b", "peak_z")
    np.testing.assert_array_equal(result.indices, [4, 3, 2, 1, 0])
    np.testing.assert_array_equal(result.scores, np.zeros(5))
    assert result.to_frame()["peak_id"].tolist() == list(result.peak_ids)


def test_reference_peak_selector_is_invariant_to_reference_row_order():
    adata = _reference_adata(
        [[3, 0, 1], [1, 2, 0], [0, 4, 1], [2, 1, 0], [0, 2, 3]],
        ["A", "A", "B", "B", "B"],
        ["peak3", "peak1", "peak2"],
    )
    permutation = [4, 1, 3, 0, 2]

    original = select_reference_peaks(adata, "cell_type", n_top_peaks=3, min_reference_cells=1)
    reordered = select_reference_peaks(adata[permutation].copy(), "cell_type", n_top_peaks=3, min_reference_cells=1)

    assert reordered.peak_ids == original.peak_ids
    np.testing.assert_array_equal(reordered.indices, original.indices)
    np.testing.assert_array_equal(reordered.scores, original.scores)
    np.testing.assert_array_equal(reordered.nonzero_reference_cells, original.nonzero_reference_cells)
    np.testing.assert_array_equal(reordered.total_reference_counts, original.total_reference_counts)


def test_reference_peak_selector_is_independent_of_declared_type_order():
    adata = _reference_adata(
        [[3, 0, 1], [1, 2, 0], [0, 4, 1], [2, 1, 0]],
        ["A", "A", "B", "B"],
        ["peak3", "peak1", "peak2"],
    )

    forward = select_reference_peaks(
        adata, "cell_type", cell_types=["A", "B"], n_top_peaks=3, min_reference_cells=1
    )
    reverse = select_reference_peaks(
        adata, "cell_type", cell_types=["B", "A"], n_top_peaks=3, min_reference_cells=1
    )

    assert forward.cell_types == reverse.cell_types == ("A", "B")
    assert forward.peak_ids == reverse.peak_ids
    np.testing.assert_array_equal(forward.scores, reverse.scores)


def test_reference_peak_selector_fails_when_a_type_has_zero_total_count():
    adata = _reference_adata(
        [[2, 1], [0, 0]],
        ["A", "B"],
        ["peak1", "peak2"],
    )

    with pytest.raises(ValueError, match=r"positive total training count.*B"):
        select_reference_peaks(adata, "cell_type", n_top_peaks=1, min_reference_cells=1)
