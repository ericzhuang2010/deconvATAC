from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import sparse
import yaml

from deconvatac.pp.fragment_shapes import (
    FragmentShapeQC,
    FragmentShapeResult,
    PeakInterval,
)
from scripts.prepare_shapemix_pbmc_sensitivity import (
    FIVE_BINS,
    OBSERVED_NK_FRACTION,
    TWO_BINS,
    _design_rows,
    _probabilities,
    derive_two_bin_result,
)


ROOT = Path(__file__).resolve().parents[1]


def test_frozen_design_matches_experiment_inventory() -> None:
    rows = _design_rows()
    experiment = yaml.safe_load(
        (ROOT / "configs/experiments/shapemix_pbmc_stress_v1.yaml").read_text()
    )
    assert [row.dataset_id for row in rows] == experiment["datasets"]
    assert len(rows) == len({row.dataset_id for row in rows}) == 40
    assert {
        factor: sum(row.factor == factor for row in rows)
        for factor in {row.factor for row in rows}
    } == {
        "anchor": 2,
        "depth": 6,
        "cells": 6,
        "rare_nk": 6,
        "features": 4,
        "subtype": 4,
        "reference_support": 8,
        "bins": 4,
    }


def test_controlled_rare_probability_preserves_total_and_target() -> None:
    for row in _design_rows():
        probabilities = _probabilities(row)
        assert np.isclose(sum(probabilities.values()), 1.0)
        if row.rare_nk_fraction is not None:
            assert probabilities["NK"] == row.rare_nk_fraction
        elif row.factor == "anchor":
            assert np.isclose(probabilities["NK"], OBSERVED_NK_FRACTION)


def test_two_bin_derivation_is_exact_sum_of_five_bins() -> None:
    layers = {
        bin_spec.layer: sparse.csr_matrix(np.asarray([[index + 1, index + 2]]))
        for index, bin_spec in enumerate(FIVE_BINS)
    }
    five = FragmentShapeResult(
        barcodes=("cell",),
        peaks=(
            PeakInterval("chr1", 0, 10, "chr1:0-10"),
            PeakInterval("chr1", 20, 30, "chr1:20-30"),
        ),
        bins=FIVE_BINS,
        layers=layers,
        qc=FragmentShapeQC(
            assigned_cut_sites=sum(int(layer.sum()) for layer in layers.values()),
            cut_sites_per_bin={
                name: int(layer.sum()) for name, layer in layers.items()
            },
        ),
        right_cut_offset=0,
    )
    two = derive_two_bin_result(five)
    assert two.bins == TWO_BINS
    np.testing.assert_array_equal(two.X.toarray(), five.X.toarray())
    np.testing.assert_array_equal(
        two.layers[TWO_BINS[0].layer].toarray(),
        (layers[FIVE_BINS[0].layer] + layers[FIVE_BINS[1].layer]).toarray(),
    )
    assert sum(two.qc.cut_sites_per_bin.values()) == int(five.X.sum())


def test_sensitivity_references_are_canonical_processed_objects() -> None:
    template = yaml.safe_load(
        (ROOT / "configs/datasets/shapemix_pbmc_stress_v1.yaml").read_text()
    )
    assert template["reference_root"].startswith("data/processed/references/")
    assert template["output_root"] == "data/processed/datasets"
    assert template["reusable_root"].startswith("data/processed/shapemix/")
