from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
import yaml

from deconvatac.pp.fragment_shapes import FragmentShapeQC, FragmentShapeResult, PeakInterval
from scripts.prepare_shapemix_pbmc_sensitivity import FIVE_BINS
from scripts.prepare_shapemix_pbmc_sensitivity_v2 import (
    CAMPAIGN_ID,
    THREE_BINS,
    _derive_three_bin_result,
    _design_rows_v2,
    _reference_rows,
)
from deconvatac.pp.fragment_shapes import build_fragment_shape_anndata


ROOT = Path(__file__).resolve().parents[1]


def test_v2_design_matches_frozen_experiment_without_reusing_v1_ids() -> None:
    experiment = yaml.safe_load(
        (ROOT / "configs/experiments/shapemix_pbmc_stress_v2.yaml").read_text()
    )
    rows = _design_rows_v2()
    assert experiment["sensitivity_campaign"] == CAMPAIGN_ID
    assert [row.dataset_id for row in rows] == experiment["datasets"]
    assert len(rows) == len({row.dataset_id for row in rows}) == 40
    assert all("stress_v2" in row.dataset_id for row in rows)
    assert experiment["continue_on_error"] is False


def test_v2_three_bin_derivation_exactly_sums_the_five_bin_recount() -> None:
    layers = {
        item.layer: sparse.csr_matrix(np.asarray([[index + 1, index + 2]]))
        for index, item in enumerate(FIVE_BINS)
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
            assigned_cut_sites=sum(int(value.sum()) for value in layers.values()),
            cut_sites_per_bin={name: int(value.sum()) for name, value in layers.items()},
        ),
        right_cut_offset=0,
    )
    three = _derive_three_bin_result(five)
    assert three.bins == THREE_BINS
    np.testing.assert_array_equal(three.X.toarray(), five.X.toarray())
    np.testing.assert_array_equal(
        three.layers[THREE_BINS[0].layer].toarray(),
        (layers[FIVE_BINS[0].layer] + layers[FIVE_BINS[1].layer]).toarray(),
    )
    np.testing.assert_array_equal(
        three.layers[THREE_BINS[1].layer].toarray(),
        (layers[FIVE_BINS[2].layer] + layers[FIVE_BINS[3].layer]).toarray(),
    )


def test_reference_support_barcode_caps_are_deterministic_and_nested() -> None:
    barcodes = tuple(f"cell_{index}" for index in range(12))
    labels = ["A"] * 6 + ["B"] * 6
    result = FragmentShapeResult(
        barcodes=barcodes,
        peaks=(PeakInterval("chr1", 0, 10, "chr1:0-10"),),
        bins=FIVE_BINS,
        layers={
            item.layer: sparse.csr_matrix(np.ones((12, 1), dtype=np.int64))
            for item in FIVE_BINS
        },
        qc=FragmentShapeQC(
            assigned_cut_sites=12 * len(FIVE_BINS),
            cut_sites_per_bin={item.layer: 12 for item in FIVE_BINS},
        ),
        right_cut_offset=0,
    )
    obs = pd.DataFrame({"cell_type": labels}, index=pd.Index(barcodes))
    var = pd.DataFrame(index=pd.Index(["chr1:0-10"]))
    reference = build_fragment_shape_anndata(result, obs=obs, var=var)
    rows_two = _reference_rows(reference, ("A", "B"), 2)
    rows_three = _reference_rows(reference, ("A", "B"), 3)
    assert len(rows_two) == 4
    assert len(rows_three) == 6
    assert set(rows_two).issubset(rows_three)
    for label in ("A", "B"):
        expected = sorted(
            reference.obs_names[reference.obs["cell_type"] == label],
            key=lambda name: (hashlib.sha256(name.encode()).hexdigest(), name),
        )[:2]
        observed = reference.obs_names[rows_two][
            reference.obs.iloc[rows_two]["cell_type"].to_numpy() == label
        ].tolist()
        assert observed == expected


def test_v2_template_declares_additive_support_safe_amendment() -> None:
    template = yaml.safe_load(
        (ROOT / "configs/datasets/shapemix_pbmc_stress_v2.yaml").read_text()
    )
    assert template["template_id"] == CAMPAIGN_ID
    assert template["amendment"]["one_common_axis_across_factors"] is True
    assert template["amendment"]["prediction_values_or_accuracy_metrics_inspected"] is False
    assert template["amendment"]["no_v1_overwrite"] is True
    assert template["reference_root"].startswith("data/processed/references/")
