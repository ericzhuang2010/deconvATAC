from __future__ import annotations

import hashlib
import re
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "data" / "registry" / "datasets.yaml"
PRIMARY_PATH = ROOT / "configs" / "experiments" / "shapemix_primary_ablation.yaml"
BASELINES_PATH = ROOT / "configs" / "experiments" / "shapemix_baselines.yaml"
STRESS_PATH = ROOT / "configs" / "experiments" / "shapemix_stress_tests.yaml"
SHAPE_CONFIG_PATH = ROOT / "configs" / "methods" / "shapemix.yaml"
COUNT_CONFIG_PATH = ROOT / "configs" / "methods" / "shapemix_count_only.yaml"

OUTER_SEEDS = [1103, 2203, 3301, 4409, 5501]
INNER_SEEDS = [101, 211]
CONDITIONS = ["observed_abundance", "equal_celltype"]
DATASET_PREFIX = "pbmc_granulocyte_sorted_10k_shapemix"
DATASET_PATTERN = re.compile(
    rf"^{DATASET_PREFIX}_(?P<condition>observed_abundance|equal_celltype)"
    r"_split_(?P<outer>\d+)_mix_(?P<inner>\d+)$"
)


def _read_yaml(path: Path) -> dict:
    with path.open() as handle:
        return yaml.safe_load(handle)


def _primary_dataset_ids() -> list[str]:
    return [
        f"{DATASET_PREFIX}_{condition}_split_{outer}_mix_{inner}"
        for outer in OUTER_SEEDS
        for inner in INNER_SEEDS
        for condition in CONDITIONS
    ]


def _resolved_project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_primary_config_is_the_complete_frozen_paired_ablation() -> None:
    experiment = _read_yaml(PRIMARY_PATH)

    assert experiment["benchmark_scope"] == "primary_one_donor_conditional_resampling"
    assert experiment["protocol"] == {
        "name": "shapemix_atac_benchmark",
        "version": 1,
        "frozen_date": "2026-08-22",
        "outer_split_seeds": OUTER_SEEDS,
        "inner_mixture_seeds": INNER_SEEDS,
        "conditions": CONDITIONS,
        "pairing_unit": ["outer_split_seed", "inner_mixture_seed", "condition"],
        "resampling_unit": "outer_split_seed",
    }
    assert experiment["datasets"] == _primary_dataset_ids()
    assert len(experiment["datasets"]) == 5 * 2 * 2
    assert experiment["modalities"] == ["atac"]
    assert experiment["feature_sets"] == {"atac": ["highly_variable"]}
    assert experiment["method_runs"] == [
        {
            "id": "shapemix_length",
            "method": "shapemix",
            "config": "configs/methods/shapemix.yaml",
        },
        {
            "id": "shapemix_count_only",
            "method": "shapemix",
            "config": "configs/methods/shapemix_count_only.yaml",
        },
    ]
    assert experiment["metrics"] == ["rmse_v1", "jsd_v2"]
    assert experiment["skip_missing"] is False
    assert experiment["continue_on_error"] is True
    assert experiment["overwrite"] is False
    assert experiment["run_id_template"].endswith("{method_run_id}")

    shape = _read_yaml(SHAPE_CONFIG_PATH)
    count = _read_yaml(COUNT_CONFIG_PATH)
    assert shape["method"] == count["method"] == "shapemix"
    assert set(shape["params"]) == set(count["params"])
    assert {
        key
        for key in shape["params"]
        if shape["params"][key] != count["params"][key]
    } == {"use_shape"}
    assert shape["params"]["use_shape"] is True
    assert count["params"]["use_shape"] is False


def test_every_primary_dataset_is_registered_complete_and_shape_ready() -> None:
    registry = _read_yaml(REGISTRY_PATH)
    reference_paths_by_outer: dict[int, set[str]] = {seed: set() for seed in OUTER_SEEDS}
    selected_peak_hashes_by_outer: dict[int, set[str]] = {
        seed: set() for seed in OUTER_SEEDS
    }

    for dataset_id in _primary_dataset_ids():
        assert dataset_id in registry
        config_path = _resolved_project_path(registry[dataset_id]["config"])
        assert config_path.is_file()
        dataset = _read_yaml(config_path)
        match = DATASET_PATTERN.fullmatch(dataset_id)
        assert match is not None
        outer = int(match.group("outer"))
        inner = int(match.group("inner"))
        condition = match.group("condition")

        assert dataset["dataset_id"] == dataset_id
        assert dataset["benchmark_scope"] == "primary"
        assert dataset["simulation"]["primary_dataset"] is True
        assert dataset["simulation"]["depth_retain_probability"] is None
        assert dataset["simulation"]["outer_split_seed"] == outer
        assert dataset["simulation"]["inner_mixture_seed"] == inner
        assert dataset["simulation"]["condition"] == condition

        atac = dataset["modalities"]["atac"]
        assert len(atac["truth"]["cell_types"]) == 16
        assert [bin_spec["layer"] for bin_spec in atac["fragment_shape"]["bins"]] == [
            "fragment_length_lt_100",
            "fragment_length_100_249",
            "fragment_length_ge_250",
        ]
        assert atac["feature_sets"]["all"] == {"mode": "all"}
        selected_path = _resolved_project_path(
            atac["feature_sets"]["highly_variable"]["path"]
        )
        reference_path = _resolved_project_path(atac["reference"]["path"])
        spatial_path = _resolved_project_path(atac["spatial"]["path"])
        truth_path = _resolved_project_path(atac["truth"]["path"])
        for path in (selected_path, reference_path, spatial_path, truth_path):
            assert path.is_file()
        assert len(selected_path.read_text().splitlines()) == 5000

        manifest_path = _resolved_project_path(dataset["simulation"]["manifest"])
        manifest = _read_yaml(manifest_path)
        assert manifest["status"] == "complete"
        assert manifest["benchmark_scope"] == "primary"
        assert manifest["simulation"]["num_spots"] == 1024
        assert manifest["simulation"]["num_features"] == 5000
        assert manifest["checks"] == {
            "x_equals_layer_sum": True,
            "truth_rows_sum_to_one": True,
            "reference_heldout_disjoint": True,
            "ordered_feature_hash_matches": True,
        }
        assert manifest["outputs"]["selected_peaks"]["sha256"] == _sha256(
            selected_path
        )

        reference_paths_by_outer[outer].add(str(reference_path))
        selected_peak_hashes_by_outer[outer].add(_sha256(selected_path))

    # Each nested pair and both conditions for one outer split use the exact
    # same training reference and reference-only selected peak order.
    assert all(len(paths) == 1 for paths in reference_paths_by_outer.values())
    assert all(len(hashes) == 1 for hashes in selected_peak_hashes_by_outer.values())


def test_baseline_config_runs_only_nnls_and_marks_optional_methods_as_gated() -> None:
    primary = _read_yaml(PRIMARY_PATH)
    baselines = _read_yaml(BASELINES_PATH)

    assert baselines["datasets"] == primary["datasets"]
    assert baselines["protocol"]["input_matrix"] == "X"
    assert baselines["protocol"]["input_contract"] == (
        "exact_layer_sum_on_the_primary_selected_peak_axis"
    )
    assert baselines["method_runs"] == [
        {
            "id": "nnls",
            "method": "nnls",
            "config": "configs/methods/nnls.yaml",
        }
    ]
    assert baselines["metrics"] == ["rmse_v1", "jsd_v2"]

    gates = baselines["gated_baselines"]["resource_or_dependency_gates"]
    assert {gate["id"] for gate in gates} == {
        "cell2location",
        "rctd",
        "spatialdwls",
    }
    assert {run["id"] for run in baselines["method_runs"]}.isdisjoint(
        gate["id"] for gate in gates
    )
    assert all(gate["status"].startswith("gated_") for gate in gates)
    assert all(gate["reason"] for gate in gates)
    assert all(gate["activation_requirements"] for gate in gates)
    assert "not runner jobs" in baselines["gated_baselines"]["policy"]

    for record in [*baselines["method_runs"], *gates]:
        method_config = _read_yaml(_resolved_project_path(record["config"]))
        assert method_config["method"] == record["method"]

    nnls_params = _read_yaml(ROOT / "configs" / "methods" / "nnls.yaml")["params"]
    assert nnls_params["layer_ref"] is None
    assert nnls_params["layer_spatial"] is None


def test_stress_config_is_an_explicit_non_executable_dataset_gate() -> None:
    stress = _read_yaml(STRESS_PATH)

    assert stress["benchmark_scope"] == "gated_secondary_sensitivity"
    assert stress["execution_gate"]["executable"] is False
    assert stress["execution_gate"]["status"] == (
        "gated_no_versioned_stress_datasets"
    )
    assert stress["execution_gate"]["activation_requirements"]
    assert stress["datasets"] == []
    assert stress["method_runs"] == [
        {
            "id": "shapemix_length",
            "method": "shapemix",
            "config": "configs/methods/shapemix.yaml",
        },
        {
            "id": "shapemix_count_only",
            "method": "shapemix",
            "config": "configs/methods/shapemix_count_only.yaml",
        },
    ]
    conditions = stress["stress_conditions"]
    assert {condition["id"] for condition in conditions} == {
        "depth_thinning",
        "cells_per_spot",
        "rare_cell_enrichment_depletion",
        "subtype_challenge",
        "feature_count",
        "fragment_length_cutoffs",
    }
    assert all(condition["status"] == "gated_dataset_not_built" for condition in conditions)
    assert all(condition["dataset_ids"] == [] for condition in conditions)
    assert all(condition["required_variation"] for condition in conditions)
    assert stress["metrics"] == ["rmse_v1", "jsd_v2"]
