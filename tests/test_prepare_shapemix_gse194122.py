from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.prepare_shapemix_gse194122 import (
    BROAD_LABELS,
    CELL_TYPES,
    CONDITIONS,
    DONORS,
    INNER_MIXTURE_SEEDS,
    N_TOP_PEAKS,
    _dataset_id,
    _fragment_records,
    _load_labels,
    _rank_fold,
)


def test_gse194122_broad_ontology_and_fold_support_are_frozen() -> None:
    labels = _load_labels()

    assert len(labels) == 69_249
    assert len(BROAD_LABELS) == 22
    assert tuple(dict.fromkeys(CELL_TYPES)) == CELL_TYPES
    assert set(labels["cell_type"]) == set(CELL_TYPES)
    assert set(labels["donor"]) == set(DONORS)


def test_gse194122_fold_dataset_ids_cover_the_frozen_40_units() -> None:
    identifiers = {
        _dataset_id(donor, condition, seed)
        for donor in DONORS
        for condition in CONDITIONS
        for seed in INNER_MIXTURE_SEEDS
    }

    assert len(identifiers) == 40
    assert all(identifier.startswith("gse194122_shapemix_broad7_lodo_") for identifier in identifiers)


def test_gse194122_ranker_emits_exact_deterministic_training_only_axis() -> None:
    n_candidates = 6_000
    candidates = pd.DataFrame(
        {
            "feature_index": np.arange(n_candidates),
            "feature_id": [f"chr1-{index * 10}-{index * 10 + 5}" for index in range(n_candidates)],
            "feature_type": "ATAC",
            "chromosome": "chr1",
            "start": np.arange(n_candidates) * 10,
            "end": np.arange(n_candidates) * 10 + 5,
        }
    )
    rng = np.random.default_rng(20260824)
    summed = rng.integers(
        1,
        100,
        size=(len(DONORS) * len(CELL_TYPES), n_candidates),
        dtype=np.int64,
    )
    coverage = np.ones_like(summed, dtype=np.int32)

    first = _rank_fold(1, candidates, summed, coverage)
    second = _rank_fold(1, candidates, summed, coverage)

    assert len(first) == N_TOP_PEAKS
    assert first["feature_id"].is_unique
    pd.testing.assert_frame_equal(first, second)


def test_gse194122_fragment_records_bind_all_source_hashes() -> None:
    records = _fragment_records()

    assert len(records) == 13
    assert len({record["sample_key"] for record in records}) == 13
    assert all(len(record["sha256"]) == 64 for record in records)
