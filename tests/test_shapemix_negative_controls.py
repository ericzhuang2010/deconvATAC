from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest
import yaml
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(SCRIPTS))

import run_shapemix_negative_controls as controls


CONFIG_PATH = ROOT / "configs" / "experiments" / "shapemix_negative_controls.yaml"


def test_checked_in_control_config_freezes_protocol_v1_contract() -> None:
    config = controls.load_negative_control_config(CONFIG_PATH)
    shape_aware, count_only = controls.load_ablation_configs(config)

    assert config.controls == controls.KNOWN_CONTROLS
    assert config.proportion_max_abs_tolerance == 1.0e-6
    assert config.one_bin_shape_log_likelihood_abs_tolerance == 0.0
    assert config.poisson_log_likelihood_abs_tolerance == 1.0e-10
    assert shape_aware.use_shape is True
    assert count_only.use_shape is False
    shape_aware.validate_protocol_v1()
    count_only.validate_protocol_v1()


def test_permutation_stream_is_repeatable_and_never_selects_identity() -> None:
    first, seed_tuple, draw_index = controls.deterministic_cell_type_permutation(3, 0)
    repeat, repeat_seed, repeat_draw = controls.deterministic_cell_type_permutation(3, 0)
    redrawn, redrawn_seed, redrawn_draw = controls.deterministic_cell_type_permutation(
        3, 2203
    )

    np.testing.assert_array_equal(first, [1, 0, 2])
    np.testing.assert_array_equal(repeat, first)
    assert seed_tuple == repeat_seed == (20260822, 0, 23)
    assert draw_index == repeat_draw == 1
    # The first draw for frozen outer seed 2203 is identity.  The declared
    # first-nonidentity policy deterministically selects the next stream draw.
    np.testing.assert_array_equal(redrawn, [1, 0, 2])
    assert redrawn_seed == (20260822, 2203, 23)
    assert redrawn_draw == 2
    assert not np.array_equal(redrawn, np.arange(3))

    with pytest.raises(ValueError, match="at least two"):
        controls.deterministic_cell_type_permutation(1, 0)


def test_poisson_factorization_golden_check_is_hashed_and_within_tolerance() -> None:
    result = controls.poisson_factorization_golden_check(1.0e-10)

    assert result["status"] == "pass"
    assert result["absolute_difference"] <= 1.0e-10
    assert len(result["counts_sha256"]) == 64
    assert len(result["rates_sha256"]) == 64


def test_primary_dataset_override_requires_positive_direction_attestation(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        controls,
        "load_deconvolution_input",
        lambda *args, **kwargs: pytest.fail("gate must fail before loading primary data"),
    )
    with pytest.raises(ValueError, match="positive preregistered primary direction"):
        controls.main(
            [
                "--config",
                str(CONFIG_PATH),
                "--dataset-id",
                "a_frozen_primary_dataset",
            ]
        )


class _TruthBomb:
    def __getattribute__(self, name: str):  # pragma: no cover - fires only on leakage
        raise AssertionError(f"Truth was accessed during fitting: {name}")


def _fake_fit(proportions: np.ndarray, *, shape_log_likelihood: float):
    restart = SimpleNamespace(steps=7, stopping_reason="objective_patience")
    return SimpleNamespace(
        proportions=np.asarray(proportions, dtype=np.float64),
        selected_restart=0,
        count_log_likelihood=-10.0,
        shape_log_likelihood=shape_log_likelihood,
        abundance_log_prior=-2.0,
        total_log_objective=-12.0 + shape_log_likelihood,
        restart_diagnostics=(restart,),
    )


def test_dataset_controls_preserve_count_inputs_and_do_not_read_truth(monkeypatch) -> None:
    shape_aware, count_only = controls.load_ablation_configs(
        controls.load_negative_control_config(CONFIG_PATH)
    )
    cell_types = ("B", "T", "NK")
    peak_names = ("peak_1", "peak_2")
    spot_names = ("spot_1", "spot_2")
    layer_names = ("short", "mono", "long")
    layers = {
        "short": sparse.csr_matrix([[2, 0], [1, 3]], dtype=np.int64),
        "mono": sparse.csr_matrix([[0, 2], [2, 0]], dtype=np.int64),
        "long": sparse.csr_matrix([[1, 1], [0, 2]], dtype=np.int64),
    }
    spatial = ad.AnnData(
        X=sum(layers.values(), sparse.csr_matrix((2, 2), dtype=np.int64)),
        layers=layers,
    )
    spatial.obs_names = list(spot_names)
    spatial.var_names = list(peak_names)
    bins = tuple(SimpleNamespace(name=name) for name in layer_names)
    data = SimpleNamespace(
        dataset_id="toy",
        modality="atac",
        feature_set="all",
        spatial=spatial,
        reference=object(),
        labels_key="cell_type",
        cell_types=list(cell_types),
        fragment_shape=SimpleNamespace(layer_names=layer_names, bins=bins),
        metadata={
            "dataset_config": {
                "simulation": {"outer_split_seed": 0, "inner_mixture_seed": 0}
            }
        },
        truth=_TruthBomb(),
    )
    A = np.asarray([[2.0, 0.5], [0.8, 2.1], [1.1, 1.3]])
    omega = np.asarray(
        [
            [[0.7, 0.2, 0.1], [0.2, 0.5, 0.3]],
            [[0.1, 0.6, 0.3], [0.5, 0.2, 0.3]],
            [[0.3, 0.2, 0.5], [0.1, 0.3, 0.6]],
        ]
    )
    u_peak = np.asarray([[0.4, 0.3, 0.3], [0.25, 0.35, 0.4]])
    signatures = SimpleNamespace(
        A=A,
        omega=omega,
        u_peak=u_peak,
        phi_ref=4.0,
        content_sha256="a" * 64,
        diagnostics=SimpleNamespace(feature_sha256="b" * 64),
        dispersion=SimpleNamespace(fold_membership_sha256="c" * 64),
    )
    monkeypatch.setattr(controls, "estimate_reference_signatures", lambda *a, **k: signatures)

    calls = []
    base = np.asarray([[0.5, 0.3, 0.2], [0.1, 0.6, 0.3]])
    permuted_result = np.asarray([[0.3, 0.5, 0.2], [0.2, 0.4, 0.4]])

    def fake_fit(shape_layers, *, A, omega, config, **kwargs):
        calls.append(
            {
                "layers": shape_layers,
                "A": np.asarray(A).copy(),
                "omega": np.asarray(omega).copy(),
                "config": config,
                "kwargs": kwargs,
            }
        )
        assert "truth" not in kwargs
        if len(shape_layers) == 1:
            return _fake_fit(base, shape_log_likelihood=0.0)
        homogeneous = np.array_equal(
            omega, np.broadcast_to(np.asarray(omega)[0:1], np.asarray(omega).shape)
        )
        if not config.use_shape or homogeneous:
            return _fake_fit(base, shape_log_likelihood=-4.0 if config.use_shape else 0.0)
        return _fake_fit(permuted_result, shape_log_likelihood=-5.0)

    evidence, predictions = controls.run_dataset_controls(
        data,
        shape_aware,
        count_only,
        controls.KNOWN_CONTROLS,
        fit_function=fake_fit,
    )

    assert evidence["truth_used_in_fitting"] is False
    assert evidence["all_acceptance_checks_passed"] is True
    assert evidence["controls"]["homogenized_omega"]["status"] == "pass"
    assert evidence["controls"]["one_bin"]["status"] == "pass"
    assert evidence["controls"]["poisson_factorization"]["status"] == "pass"
    permutation = evidence["controls"]["permuted_omega"]["permutation_indices"]
    assert permutation == [1, 0, 2]
    assert evidence["controls"]["permuted_omega"]["A_preserved"] is True
    assert set(predictions) == {
        "count_only",
        "homogenized_omega",
        "permuted_omega",
        "one_bin_count_only",
        "one_bin_shape",
    }
    assert all(np.array_equal(call["A"], A) for call in calls)
    collapsed_call_layers = [call["layers"] for call in calls if len(call["layers"]) == 1]
    expected_total = sum(layers.values(), sparse.csr_matrix((2, 2))).toarray()
    assert len(collapsed_call_layers) == 2
    for single_layer in collapsed_call_layers:
        np.testing.assert_array_equal(single_layer[0].toarray(), expected_total)


def test_compact_outputs_are_covered_by_sha256_manifest(tmp_path: Path) -> None:
    control_config = controls.load_negative_control_config(CONFIG_PATH)
    shape_aware, count_only = controls.load_ablation_configs(control_config)
    spatial = ad.AnnData(np.ones((2, 1)))
    spatial.obs_names = ["s1", "s2"]
    spatial.var_names = ["p1"]
    dataset_config_path = tmp_path / "dataset.yaml"
    dataset_config_path.write_text("dataset_id: toy\n")
    data = SimpleNamespace(
        dataset_id="toy",
        spatial=spatial,
        cell_types=["B", "T"],
        metadata={
            "dataset_config": {
                "dataset_id": "toy",
                "config_path": str(dataset_config_path),
            }
        },
    )
    evidence = {
        "schema_version": 1,
        "protocol_id": controls.PROTOCOL_ID,
        "all_acceptance_checks_passed": True,
    }
    predictions = {"count_only": np.asarray([[0.6, 0.4], [0.2, 0.8]])}

    manifest_path = controls.write_control_outputs(
        tmp_path / "outputs",
        data,
        evidence,
        predictions,
        control_config,
        shape_aware,
        count_only,
        overwrite=False,
    )

    manifest = yaml.safe_load(manifest_path.read_text())
    written_evidence = yaml.safe_load(
        (manifest_path.parent / "control_evidence.yaml").read_text()
    )
    executed_hashes = written_evidence["code_and_protocol"]["executed_code_sha256"]
    assert "src/deconvatac/shapemix/map.py" in executed_hashes
    assert "src/deconvatac/shapemix/likelihood.py" in executed_hashes
    assert set(manifest["outputs"]) == {
        "control_proportions.csv",
        "control_evidence.yaml",
    }
    for filename, metadata in manifest["outputs"].items():
        payload = (manifest_path.parent / filename).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == metadata["sha256"]
        assert len(payload) == metadata["bytes"]
    with pytest.raises(FileExistsError, match="already contains files"):
        controls.write_control_outputs(
            manifest_path.parent,
            data,
            evidence,
            predictions,
            control_config,
            shape_aware,
            count_only,
            overwrite=False,
        )
