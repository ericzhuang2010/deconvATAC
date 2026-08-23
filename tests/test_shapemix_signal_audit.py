from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse

from deconvatac.data import ordered_feature_sha256
from scripts.audit_shapemix_signal import audit_shape_signal, write_signal_audit


LAYERS = (
    "fragment_length_lt_100",
    "fragment_length_100_249",
    "fragment_length_ge_250",
)


def _reference(outer_split_seed: int = 0) -> ad.AnnData:
    layers = (
        np.array([[4, 0], [3, 1], [0, 4], [1, 3]], dtype=np.int64),
        np.array([[1, 2], [1, 2], [3, 0], [2, 1]], dtype=np.int64),
        np.array([[0, 2], [1, 1], [2, 0], [2, 1]], dtype=np.int64),
    )
    matrices = [sparse.csr_matrix(values) for values in layers]
    x = sum(matrices, sparse.csr_matrix((4, 2), dtype=np.int64)).tocsr()
    result = ad.AnnData(
        X=x,
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B"]},
            index=["a2", "a1", "b2", "b1"],
        ),
        var=pd.DataFrame(index=["p1", "p2"]),
    )
    for name, matrix in zip(LAYERS, matrices):
        result.layers[name] = matrix
    result.uns["fragment_shape"] = {
        "schema_version": 1,
        "axis": "parent_fragment_length_bp",
        "count_unit": "deduplicated_cut_sites",
        "read_support_policy": "ignore",
        "peak_assignment": "containing_nonoverlapping_peak",
        "left_cut_offset": 0,
        "right_cut_offset": 0,
        "bins": [
            {"name": "short", "min_inclusive": 0, "max_exclusive": 100, "layer": LAYERS[0]},
            {"name": "mono", "min_inclusive": 100, "max_exclusive": 250, "layer": LAYERS[1]},
            {"name": "long", "min_inclusive": 250, "max_exclusive": None, "layer": LAYERS[2]},
        ],
        "feature_sha256": ordered_feature_sha256(result.var_names),
        "split_sha256": "a" * 64,
    }
    result.uns["shapemix_preparation"] = {
        "pool": "reference",
        "split_sha256": "a" * 64,
        "outer_split_seed": outer_split_seed,
    }
    result.obs["split_pool"] = "reference"
    return result


def test_signal_audit_reports_required_scopes_and_is_row_order_invariant() -> None:
    reference = _reference(outer_split_seed=7)
    peaks, summary = audit_shape_signal(
        reference,
        cell_types=["A", "B"],
        split_seed=7,
    )
    reordered_peaks, reordered_summary = audit_shape_signal(
        reference[[3, 1, 2, 0]].copy(),
        cell_types=["A", "B"],
        split_seed=7,
    )

    pd.testing.assert_frame_equal(peaks, reordered_peaks)
    assert summary["scope"] == "training_reference_only"
    assert summary["global_bin_counts"] == reordered_summary["global_bin_counts"]
    assert summary["split_half_reproducibility"] == reordered_summary["split_half_reproducibility"]
    assert summary["cells"] == 4
    assert summary["peaks"] == 2
    assert set(summary["technical_confounding"]) == {
        "positive_count_cells",
        *(f"spearman_log_depth_vs_fraction.{name}" for name in LAYERS),
    }
    assert list(peaks["peak_id"]) == ["p1", "p2"]
    assert (peaks["between_type_generalized_jsd_bits"] >= 0).all()


def test_signal_audit_rejects_heldout_or_mismatched_type_universe() -> None:
    heldout = _reference()
    heldout.uns["shapemix_preparation"]["pool"] = "heldout"
    heldout.obs["split_pool"] = "heldout"
    with pytest.raises(ValueError, match="reference pool"):
        audit_shape_signal(heldout, cell_types=["A", "B"], split_seed=0)

    with pytest.raises(ValueError, match="exactly match"):
        audit_shape_signal(_reference(), cell_types=["A", "C"], split_seed=0)


def test_signal_audit_writes_hashed_outputs_without_implicit_overwrite(tmp_path: Path) -> None:
    peaks, summary = audit_shape_signal(_reference(), cell_types=["A", "B"], split_seed=0)
    peak_path, summary_path = write_signal_audit(peaks, summary, tmp_path)

    restored = yaml.safe_load(summary_path.read_text())
    assert peak_path.exists()
    assert len(restored["files"]["signal_audit.csv"]["sha256"]) == 64
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        write_signal_audit(peaks, summary, tmp_path)
