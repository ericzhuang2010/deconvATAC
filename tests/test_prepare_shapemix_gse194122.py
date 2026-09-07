from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import scripts.prepare_shapemix_gse194122 as gse194122

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
    _repository_path,
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


def test_repository_path_preserves_lexical_path_through_symlink(
    tmp_path, monkeypatch
) -> None:
    repository = tmp_path / "repository"
    storage = tmp_path / "storage"
    repository.mkdir()
    storage.mkdir()
    (repository / "data").symlink_to(storage, target_is_directory=True)
    monkeypatch.setattr("scripts.prepare_shapemix_gse194122.ROOT", repository)

    assert _repository_path(repository / "data" / "artifact.h5ad") == (
        "data/artifact.h5ad"
    )


def test_complete_feature_axis_is_validated_and_reused(
    tmp_path, monkeypatch, capsys
) -> None:
    axis_root = tmp_path / "feature_axes"
    donors = (1, 2)
    n_top_peaks = 3
    union_ids = ["peak_a", "peak_b", "peak_c", "peak_d"]
    union = pd.DataFrame(
        {
            "feature_id": union_ids,
            "chromosome": ["chr1"] * len(union_ids),
            "start": np.arange(len(union_ids)) * 10,
            "end": np.arange(len(union_ids)) * 10 + 5,
        }
    )
    union_dir = axis_root / "union"
    union_dir.mkdir(parents=True)
    union.to_csv(
        union_dir / "peaks.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    (union_dir / "selected_peaks.txt").write_text("\n".join(union_ids) + "\n")

    fold_records = []
    for donor, selected_ids in ((1, union_ids[:3]), (2, union_ids[1:])):
        donor_dir = axis_root / f"donor_{donor}"
        donor_dir.mkdir()
        table_path = donor_dir / "peak_selection.tsv.gz"
        union[union["feature_id"].isin(selected_ids)].to_csv(
            table_path, sep="\t", index=False, compression="gzip"
        )
        (donor_dir / "selected_peaks.txt").write_text(
            "\n".join(selected_ids) + "\n"
        )
        fold_records.append(
            {
                "donor": donor,
                "selected_peaks": n_top_peaks,
                "selected_feature_sha256": (
                    gse194122.ordered_feature_sha256(selected_ids)
                ),
                "peak_selection_sha256": gse194122._sha256(table_path),
            }
        )
    gse194122._write_yaml(
        axis_root / "manifest.yaml",
        {
            "schema_version": 1,
            "status": "complete",
            "label_ontology": ["type_a", "type_b"],
            "donors": list(donors),
            "union_peaks": len(union_ids),
            "selector": {"n_top_peaks": n_top_peaks},
            "folds": fold_records,
        },
    )

    monkeypatch.setattr(gse194122, "AXIS_ROOT", axis_root)
    monkeypatch.setattr(gse194122, "DONORS", donors)
    monkeypatch.setattr(gse194122, "CELL_TYPES", ("type_a", "type_b"))
    monkeypatch.setattr(gse194122, "N_TOP_PEAKS", n_top_peaks)
    monkeypatch.setattr(
        gse194122,
        "_aggregate_donor_type_counts",
        lambda *_args, **_kwargs: pytest.fail(
            "complete axes must not be recomputed"
        ),
    )

    gse194122.build_feature_axes()

    assert capsys.readouterr().out == "feature_axes status=reused\n"


def test_incomplete_existing_feature_axis_fails_closed(
    tmp_path, monkeypatch
) -> None:
    axis_root = tmp_path / "feature_axes"
    axis_root.mkdir()
    monkeypatch.setattr(gse194122, "AXIS_ROOT", axis_root)

    with pytest.raises(FileExistsError, match="incomplete"):
        gse194122.build_feature_axes()


def test_existing_fold_split_with_current_schema_is_reusable(tmp_path) -> None:
    expected = pd.DataFrame(
        {
            "cell_id": ["cell-1", "cell-2"],
            "sample_key": ["s1d1", "s1d2"],
            "site": [1, 1],
            "donor": [1, 2],
            "author_cell_type": ["B1 B", "NK"],
            "cell_type": ["B/plasma", "NK/ILC"],
            "partition": ["heldout", "training"],
        }
    )
    path = tmp_path / "cells.tsv.gz"
    expected.to_csv(path, sep="\t", index=False, compression="gzip")

    gse194122._validate_existing_fold_split(path, expected)


def test_existing_fold_split_with_legacy_schema_fails_closed(tmp_path) -> None:
    expected = pd.DataFrame(
        {
            "cell_id": ["cell-1"],
            "sample_key": ["s1d1"],
            "site": [1],
            "donor": [1],
            "author_cell_type": ["B1 B"],
            "cell_type": ["B/plasma"],
            "partition": ["heldout"],
        }
    )
    path = tmp_path / "cells.tsv.gz"
    pd.DataFrame(
        {
            "cell_id": ["cell-1"],
            "role": ["heldout"],
            "sample_key": ["s1d1"],
            "site": [1],
            "donor": [1],
            "harmonized_cell_type": ["B1 B"],
        }
    ).to_csv(path, sep="\t", index=False, compression="gzip")

    with pytest.raises(ValueError, match="current schema/content"):
        gse194122._validate_existing_fold_split(path, expected)
