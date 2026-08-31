from __future__ import annotations

import numpy as np
import pytest

from scripts.materialize_shapemix_gse246791_reference import (
    parse_interval,
    rank_selected_indices,
)


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
