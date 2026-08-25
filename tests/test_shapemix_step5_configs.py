from __future__ import annotations

from pathlib import Path

import yaml

from deconvatac.shapemix.config import (
    ShapeMixConfig,
    validate_nested_ablation_configs,
)


ROOT = Path(__file__).resolve().parents[1]
SHAPE_CONFIG = ROOT / "configs" / "methods" / "shapemix.yaml"
COUNT_CONFIG = ROOT / "configs" / "methods" / "shapemix_count_only.yaml"
CUDA_SHAPE_CONFIG = ROOT / "configs" / "methods" / "shapemix_cuda.yaml"
CUDA_COUNT_CONFIG = ROOT / "configs" / "methods" / "shapemix_count_only_cuda.yaml"
SMOKE_EXPERIMENT = ROOT / "configs" / "experiments" / "shapemix_smoke.yaml"
ALL_METHODS_EXPERIMENT = (
    ROOT / "configs" / "experiments" / "all_methods_all_atac_datasets.yaml"
)
SMOKE_DATASET = (
    "pbmc_granulocyte_sorted_10k_shapemix_"
    "equal_celltype_split_000_mix_000_smoke"
)


def _read_yaml(path: Path) -> dict:
    with path.open() as handle:
        return yaml.safe_load(handle)


def test_method_configs_are_full_frozen_nested_ablation_configs() -> None:
    shape_mapping = _read_yaml(SHAPE_CONFIG)
    count_mapping = _read_yaml(COUNT_CONFIG)
    shape = ShapeMixConfig.from_mapping(shape_mapping)
    count = ShapeMixConfig.from_mapping(count_mapping)

    validate_nested_ablation_configs(shape, count)
    shape.validate_protocol_v1()
    count.validate_protocol_v1()
    assert shape_mapping == {"method": "shapemix", "params": shape.to_dict()}
    assert count_mapping == {"method": "shapemix", "params": count.to_dict()}
    assert set(shape_mapping["params"]) == set(count_mapping["params"])
    differences = {
        key
        for key in shape_mapping["params"]
        if shape_mapping["params"][key] != count_mapping["params"][key]
    }
    assert differences == {"use_shape"}

    # This stronger textual check prevents comments, ordering, or omitted
    # defaults from making the two arms differ outside the one nested switch.
    normalized_shape = SHAPE_CONFIG.read_text().replace(
        "use_shape: true", "use_shape: <ablation>"
    )
    normalized_count = COUNT_CONFIG.read_text().replace(
        "use_shape: false", "use_shape: <ablation>"
    )
    assert normalized_shape == normalized_count
    assert not any(
        "fragment" in key or "layer" in key for key in shape_mapping["params"]
    )


def test_cuda_method_configs_preserve_the_nested_protocol_contract() -> None:
    shape_mapping = _read_yaml(CUDA_SHAPE_CONFIG)
    count_mapping = _read_yaml(CUDA_COUNT_CONFIG)
    shape = ShapeMixConfig.from_mapping(shape_mapping)
    count = ShapeMixConfig.from_mapping(count_mapping)

    validate_nested_ablation_configs(shape, count)
    shape.validate_protocol_v1()
    count.validate_protocol_v1()
    assert shape.is_protocol_v1
    assert count.is_protocol_v1
    assert shape.device == count.device == "cuda:0"
    assert {
        key
        for key in shape_mapping["params"]
        if shape_mapping["params"][key] != count_mapping["params"][key]
    } == {"use_shape"}
    assert CUDA_SHAPE_CONFIG.read_text().replace(
        "use_shape: true", "use_shape: <ablation>"
    ) == CUDA_COUNT_CONFIG.read_text().replace(
        "use_shape: false", "use_shape: <ablation>"
    )


def test_smoke_experiment_names_both_variants_and_is_fail_fast_development_only() -> None:
    experiment = _read_yaml(SMOKE_EXPERIMENT)

    assert experiment["benchmark_scope"] == "development_smoke"
    assert "development" in experiment["run_group"]
    assert "smoke" in experiment["run_group"]
    assert experiment["datasets"] == [SMOKE_DATASET]
    assert experiment["modalities"] == ["atac"]
    assert experiment["feature_sets"] == {"atac": ["highly_variable"]}
    assert "methods" not in experiment
    assert "method_configs" not in experiment
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
    assert experiment["metrics"] == ["rmse", "jsd"]
    assert experiment["output_root"].startswith("results/development")
    assert experiment["run_id_template"].endswith("{method_run_id}")
    assert experiment["skip_missing"] is False
    assert experiment["continue_on_error"] is False
    assert experiment["overwrite"] is False

    for method_run in experiment["method_runs"]:
        resolved = ROOT / method_run["config"]
        assert resolved.is_file()
        assert _read_yaml(resolved)["method"] == method_run["method"]


def test_shapemix_is_not_added_to_generic_all_atac_experiment() -> None:
    generic = _read_yaml(ALL_METHODS_EXPERIMENT)
    assert "shapemix" not in generic.get("methods", [])
    assert "shapemix" not in generic.get("method_configs", {})
    assert "shapemix" not in ALL_METHODS_EXPERIMENT.read_text().lower()
