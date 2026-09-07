from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from scripts.materialize_shapemix_gse246791_reference import (
    fragment_source_identifier,
    parse_interval,
    positive_reference_feature_mask,
    rank_selected_indices,
)


def test_gse246791_combined_source_identifier_names_the_fragment_file() -> None:
    identifier = fragment_source_identifier("GSM7876902")
    assert identifier == "GSM7876902_fragments.tsv.gz"
    assert "/" not in identifier


def test_gse246791_positive_reference_support_filter_is_order_preserving() -> None:
    layers = {
        "short": sparse.csr_matrix([[1, 0, 0], [0, 0, 2]]),
        "long": sparse.csr_matrix([[0, 0, 3], [0, 0, 0]]),
    }
    keep, totals = positive_reference_feature_mask(layers)
    np.testing.assert_array_equal(keep, [True, False, True])
    np.testing.assert_array_equal(totals, [1, 0, 5])


def test_gse246791_deposited_interval_parser_is_fail_closed() -> None:
    assert parse_interval("chr1:0-500") == ("chr1", 0, 500)
    assert parse_interval("chrY:91744500-91744698") == (
        "chrY",
        91_744_500,
        91_744_698,
    )
    with pytest.raises(ValueError, match="Invalid deposited interval"):
        parse_interval("chr1:500-500")


def test_gse246791_ranker_uses_peak_id_as_the_final_tie_break() -> None:
    score = np.asarray([3.0, 2.0, 2.0, 1.0])
    coverage = np.asarray([10, 10, 10, 10])
    total = np.asarray([5, 5, 5, 5])
    identifiers = {0: "chr1:0-500", 1: "chr2:0-500", 2: "chr1:500-1000", 3: "chr3:0-500"}

    selected = rank_selected_indices(
        score,
        coverage,
        total,
        n_top=2,
        id_loader=lambda indices: {
            int(index): identifiers[int(index)] for index in indices
        },
    )

    assert selected == [0, 2]
