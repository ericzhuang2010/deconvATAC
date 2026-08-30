from __future__ import annotations

import numpy as np
import pytest

from scripts.prepare_shapemix_reference_marker_features import rank_marker_indices


def test_reference_marker_ranker_is_reference_only_and_identifier_deterministic() -> None:
    mean_type = np.asarray([4.0, 2.0, 2.0, 1.0])
    mean_rest = np.asarray([1.0, 1.0, 1.0, 2.0])
    coverage = np.asarray([10, 10, 10, 100])
    totals = np.asarray([40.0, 20.0, 20.0, 100.0])
    names = ["peak-a", "peak-c", "peak-b", "peak-d"]

    assert rank_marker_indices(
        mean_type,
        mean_rest,
        coverage,
        totals,
        names,
        n_markers=2,
    ) == [0, 2]


def test_reference_marker_ranker_fails_when_specific_support_is_insufficient() -> None:
    with pytest.raises(ValueError, match="marker-support gate"):
        rank_marker_indices(
            np.asarray([2.0, 0.5]),
            np.asarray([1.0, 1.0]),
            np.asarray([9, 100]),
            np.asarray([2.0, 2.0]),
            ["a", "b"],
            n_markers=1,
        )
