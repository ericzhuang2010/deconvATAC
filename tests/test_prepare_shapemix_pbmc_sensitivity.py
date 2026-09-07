from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
import yaml

from deconvatac.data import FragmentShapeSpec, validate_fragment_shape_spec
from deconvatac.pp.fragment_shapes import (
    FragmentShapeQC,
    FragmentShapeResult,
    PeakInterval,
    build_fragment_shape_anndata,
)
from scripts.prepare_shapemix_pbmc_sensitivity import (
    FIVE_BINS,
    OBSERVED_NK_FRACTION,
    TWO_BINS,
    _design_rows,
    _probabilities,
    derive_two_bin_result,
)
from scripts.regenerate_shapemix_pbmc_simulations import subset_shape_cells


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


def test_alternate_bins_use_schema_v2_and_can_be_subset() -> None:
    barcodes = ("reference_cell", "heldout_cell")
    peaks = (
        PeakInterval("chr1", 0, 10, "chr1:0-10"),
        PeakInterval("chr1", 20, 30, "chr1:20-30"),
    )
    layers = {
        bin_spec.layer: sparse.csr_matrix(
            np.asarray([[index + 1, index + 2], [index + 2, index + 1]])
        )
        for index, bin_spec in enumerate(FIVE_BINS)
    }
    result = FragmentShapeResult(
        barcodes=barcodes,
        peaks=peaks,
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
    obs = pd.DataFrame(
        {"split_pool": ["reference", "heldout"]},
        index=pd.Index(barcodes, name="barcode"),
    )
    var = pd.DataFrame(index=pd.Index([peak.name for peak in peaks], name="peak"))
    adata = build_fragment_shape_anndata(result, obs=obs, var=var)

    spec = FragmentShapeSpec.from_mapping(adata.uns["fragment_shape"])
    assert spec.schema_version == 2
    assert [
        (item.name, item.min_inclusive, item.max_exclusive, item.layer)
        for item in spec.bins
    ] == [
        (item.name, item.min_inclusive, item.max_exclusive, item.layer)
        for item in FIVE_BINS
    ]
    validate_fragment_shape_spec(spec)

    selected = subset_shape_cells(
        adata,
        observation_mask=adata.obs["split_pool"].astype(str).eq("reference"),
    )
    selected_spec = FragmentShapeSpec.from_mapping(selected.uns["fragment_shape"])
    assert selected_spec.schema_version == 2
    assert selected_spec.bins == spec.bins
    assert selected.n_obs == 1
    expected = sum(
        (selected.layers[layer] for layer in selected_spec.layer_names),
        sparse.csr_matrix(selected.shape, dtype=np.int64),
    )
    assert (selected.X != expected).nnz == 0


def test_sensitivity_references_are_canonical_processed_objects() -> None:
    template = yaml.safe_load(
        (ROOT / "configs/datasets/shapemix_pbmc_stress_v1.yaml").read_text()
    )
    assert template["reference_root"].startswith("data/processed/references/")
    assert template["output_root"] == "data/processed/datasets"
    assert template["reusable_root"].startswith("data/processed/shapemix/")
