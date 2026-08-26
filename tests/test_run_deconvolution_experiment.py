import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

from scripts.run_deconvolution import (
    collect_environment,
    method_config_sha256,
    resolve_experiment_jobs,
    resolve_method_runs,
    run_experiment,
    select_experiment_shard,
)


def _write_toy_dataset(tmp_path: Path, *, declare_cell_types: bool = True) -> Path:
    reference = ad.AnnData(
        X=np.array([[5, 0], [4, 1], [0, 5]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B"]}, index=["c1", "c2", "c3"]),
        var=pd.DataFrame({"highly_variable": [True, True]}, index=["f1", "f2"]),
    )
    spatial = ad.AnnData(
        X=np.array([[9, 1], [1, 8]], dtype=float),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["f1", "f2"]),
    )
    spatial.obsm["spatial"] = np.array([[0, 0], [1, 0]], dtype=float)

    truth = pd.DataFrame([[0.9, 0.1], [0.2, 0.8]], index=["s1", "s2"], columns=["A", "B"])

    reference_path = tmp_path / "reference.h5ad"
    spatial_path = tmp_path / "spatial.h5ad"
    truth_path = tmp_path / "truth.csv"
    reference.write_h5ad(reference_path)
    spatial.write_h5ad(spatial_path)
    truth.to_csv(truth_path)

    truth_spec = {"path": str(truth_path)}
    if declare_cell_types:
        truth_spec["cell_types"] = ["A", "B"]

    dataset_config = {
        "dataset_id": "toy",
        "labels_key": "cell_type",
        "spatial_key": "spatial",
        "modalities": {
            "atac": {
                "reference": {"path": str(reference_path)},
                "spatial": {"path": str(spatial_path)},
                "truth": truth_spec,
                "feature_sets": {
                    "highly_variable": {"var_column": "highly_variable"},
                    "all": {"mode": "all"},
                },
            }
        },
    }
    config_path = tmp_path / "toy.yaml"
    config_path.write_text(yaml.safe_dump(dataset_config))

    registry_path = tmp_path / "datasets.yaml"
    registry_path.write_text(yaml.safe_dump({"toy": {"config": str(config_path)}}))
    return registry_path


def _write_toy_dataset_without_truth(tmp_path: Path) -> Path:
    registry_path = _write_toy_dataset(tmp_path)
    config_path = tmp_path / "toy.yaml"
    dataset_config = yaml.safe_load(config_path.read_text())
    modality_config = dataset_config["modalities"]["atac"]
    modality_config.pop("truth")
    modality_config["cell_types"] = ["A", "B"]
    config_path.write_text(
        yaml.safe_dump(dataset_config, sort_keys=False)
    )
    return registry_path


def test_shapemix_environment_metadata_is_plain_yaml_safe_versions():
    environment = collect_environment("shapemix")
    assert {"torch", "pysam"}.issubset(environment["packages"])
    assert all(isinstance(version, str) for version in environment["packages"].values())
    assert yaml.safe_load(yaml.safe_dump(environment)) == environment


def test_run_experiment_writes_batch_comparison(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "experiment.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "toy_batch",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["highly_variable"],
                "methods": ["nnls"],
                "method_configs": {"nnls": {"method": "nnls", "params": {}}},
                "metrics": ["rmse", "jsd"],
                "skip_missing": True,
            }
        )
    )

    comparison_path = run_experiment(
        experiment_config_path=experiment_path,
        registry=registry_path,
        output_root_override=tmp_path / "results",
        overwrite=True,
    )

    assert comparison_path == tmp_path / "results" / "toy_batch" / "comparison.csv"
    comparison = pd.read_csv(comparison_path)
    assert set(comparison["metric"]) == {"rmse_v1", "jsd_v2"}
    assert set(comparison["metric_version"]) == {"v1", "v2"}
    assert set(comparison["evaluation_contract_version"]) == {"v2"}
    assert comparison["cell_type_universe"].map(
        lambda value: tuple(json.loads(value)) == ("A", "B")
    ).all()
    assert set(comparison["status"]) == {"success"}
    assert set(comparison["method"]) == {"nnls"}
    assert set(comparison["method_run_id"]) == {"nnls"}
    batch_dir = tmp_path / "results" / "toy_batch"
    runs = pd.read_csv(batch_dir / "runs.csv")
    assert runs.loc[0, "method"] == "nnls"
    assert runs.loc[0, "method_run_id"] == "nnls"
    run_dir = (
        tmp_path
        / "results"
        / "toy_batch"
        / "toy__atac__highly_variable__nnls"
    )
    assert (run_dir / "results" / "proportions.csv").exists()
    metadata = yaml.safe_load((run_dir / "run.yaml").read_text())
    assert metadata["method"] == "nnls"
    assert metadata["method_run_id"] == "nnls"


def test_prediction_only_experiment_writes_predictions_without_truth_metrics(tmp_path):
    registry_path = _write_toy_dataset_without_truth(tmp_path)
    experiment_path = tmp_path / "prediction_only.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "prediction_only",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "methods": ["nnls"],
                "method_configs": {"nnls": {"method": "nnls", "params": {}}},
                "evaluation_mode": "prediction_only",
                "metrics": [],
            }
        )
    )

    comparison_path = run_experiment(
        experiment_path,
        registry=registry_path,
        output_root_override=tmp_path / "results",
    )
    batch_dir = comparison_path.parent
    run_dir = batch_dir / "toy__atac__all__nnls"
    comparison = pd.read_csv(comparison_path)
    metadata = yaml.safe_load((run_dir / "run.yaml").read_text())
    inputs = yaml.safe_load((run_dir / "inputs.yaml").read_text())

    assert comparison.loc[0, "status"] == "success"
    assert comparison.loc[0, "evaluation_mode"] == "prediction_only"
    assert pd.isna(comparison.loc[0, "metric"])
    assert pd.isna(comparison.loc[0, "value"])
    assert not (run_dir / "results" / "truth.csv").exists()
    assert "truth" not in inputs
    assert metadata["proportion_evaluation"]["evidence"] == "prediction_only"
    batch_manifest = yaml.safe_load((batch_dir / "batch_manifest.yaml").read_text())
    assert batch_manifest["metrics"] == []


def test_named_method_runs_have_distinct_directories_and_full_provenance(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    file_config = {"method": "nnls", "params": {"seed": 7, "device": "cpu", "dtype": "float64"}}
    file_config_path = tmp_path / "nnls_file.yaml"
    file_config_path.write_text(yaml.safe_dump(file_config, sort_keys=False))
    inline_config = {"method": "nnls", "params": {"seed": 8, "device": "cpu", "dtype": "float32"}}
    experiment_path = tmp_path / "variants.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "variants",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["highly_variable"],
                "method_runs": [
                    {"id": "nnls_file", "method": "nnls", "config": str(file_config_path)},
                    {"id": "nnls_inline", "method": "nnls", "config": inline_config},
                ],
                "metrics": ["rmse"],
            },
            sort_keys=False,
        )
    )

    comparison_path = run_experiment(
        experiment_path,
        registry=registry_path,
        output_root_override=tmp_path / "results",
    )
    batch_dir = comparison_path.parent
    expected_ids = {"nnls_file", "nnls_inline"}
    assert {path.name for path in batch_dir.iterdir() if path.is_dir()} == {
        f"toy__atac__highly_variable__{variant}" for variant in expected_ids
    }

    runs = pd.read_csv(batch_dir / "runs.csv")
    comparison = pd.read_csv(comparison_path)
    assert set(runs["method"]) == {"nnls"}
    assert set(runs["method_run_id"]) == expected_ids
    assert set(comparison["method"]) == {"nnls"}
    assert set(comparison["method_run_id"]) == expected_ids
    assert set(runs["method_config_sha256"]) == {
        method_config_sha256(file_config),
        method_config_sha256(inline_config),
    }
    assert runs.loc[runs["method_run_id"] == "nnls_file", "method_config_source_sha256"].item() == (
        hashlib.sha256(file_config_path.read_bytes()).hexdigest()
    )
    assert pd.isna(
        runs.loc[runs["method_run_id"] == "nnls_inline", "method_config_source_path"].item()
    )

    for variant, expected_config in (("nnls_file", file_config), ("nnls_inline", inline_config)):
        run_dir = batch_dir / f"toy__atac__highly_variable__{variant}"
        metadata = yaml.safe_load((run_dir / "run.yaml").read_text())
        assert metadata["method"] == "nnls"
        assert metadata["method_run_id"] == variant
        assert metadata["method_config"] == expected_config
        assert metadata["method_config_sha256"] == method_config_sha256(expected_config)
        assert metadata["seed"] == expected_config["params"]["seed"]
        assert metadata["device"] == "cpu"
        assert metadata["dtype"] == expected_config["params"]["dtype"]
        assert metadata["determinism"] == {
            "seed": expected_config["params"]["seed"],
            "device": "cpu",
            "dtype": expected_config["params"]["dtype"],
            "source": "resolved_method_config",
        }
        assert metadata["environment"]["packages"] == metadata["software_versions"]
        assert "torch" not in metadata["software_versions"]
        assert "pysam" not in metadata["software_versions"]
        environment = (run_dir / "environment.txt").read_text()
        assert "torch=" not in environment
        assert "pysam=" not in environment
        assert json.loads(
            comparison.loc[comparison["method_run_id"] == variant, "method_config"].item()
        ) == expected_config


@pytest.mark.parametrize(
    "schema",
    [
        {},
        {"methods": ["nnls"], "method_runs": [{"id": "named", "method": "nnls"}]},
    ],
)
def test_mixed_or_missing_method_schema_fails_before_output_creation(tmp_path, schema):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "invalid_schema.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "must_not_exist",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                **schema,
            }
        )
    )
    output_root = tmp_path / "results"
    with pytest.raises(ValueError, match="exactly one of 'methods' or 'method_runs'"):
        run_experiment(experiment_path, registry_path, output_root_override=output_root)
    assert not (output_root / "must_not_exist").exists()


def test_duplicate_variant_ids_are_rejected():
    config = {
        "method_runs": [
            {"id": "same", "method": "nnls", "config": {"method": "nnls"}},
            {"id": "same", "method": "nnls", "config": {"method": "nnls"}},
        ]
    }
    with pytest.raises(ValueError, match="Duplicate method run id.*same"):
        resolve_method_runs(config)


def test_method_runs_rejects_legacy_method_configs_ambiguity():
    config = {
        "method_runs": [{"id": "named", "method": "nnls", "config": {"method": "nnls"}}],
        "method_configs": {"nnls": {"method": "nnls"}},
    }
    with pytest.raises(ValueError, match="method_configs is only valid"):
        resolve_method_runs(config)


def test_resolved_config_method_must_match_registry_method():
    config = {
        "method_runs": [
            {"id": "wrong", "method": "nnls", "config": {"method": "rctd", "params": {}}}
        ]
    }
    with pytest.raises(ValueError, match="does not match method 'nnls'"):
        resolve_method_runs(config)


def test_custom_template_run_id_collision_is_preflighted_before_output_creation(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "collision.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "collision",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "method_runs": [
                    {"id": "first", "method": "nnls", "config": {"method": "nnls"}},
                    {"id": "second", "method": "nnls", "config": {"method": "nnls"}},
                ],
                "run_id_template": "{dataset}__{method}",
            },
            sort_keys=False,
        )
    )
    output_root = tmp_path / "results"
    with pytest.raises(ValueError, match="Duplicate resolved run ID.*toy__nnls"):
        run_experiment(experiment_path, registry_path, output_root_override=output_root)
    assert not (output_root / "collision").exists()


@pytest.mark.parametrize("template", ["../{method_run_id}", "/absolute/{method_run_id}"])
def test_path_escaping_run_ids_are_rejected_before_output_creation(tmp_path, template):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "escaping.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "escaping",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "method_runs": [
                    {"id": "variant", "method": "nnls", "config": {"method": "nnls"}}
                ],
                "run_id_template": template,
            }
        )
    )
    output_root = tmp_path / "results"
    with pytest.raises(ValueError, match="single relative directory name"):
        run_experiment(experiment_path, registry_path, output_root_override=output_root)
    assert not (output_root / "escaping").exists()


def test_failure_rows_and_manifest_include_method_variant_and_config(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    failing_config = {"method": "nnls", "params": {"layer_ref": "not_a_layer", "seed": 3}}
    experiment_path = tmp_path / "failure.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "failure",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "method_runs": [
                    {"id": "nnls_expected_failure", "method": "nnls", "config": failing_config}
                ],
                "metrics": ["rmse"],
                "continue_on_error": True,
            },
            sort_keys=False,
        )
    )

    comparison_path = run_experiment(
        experiment_path,
        registry_path,
        output_root_override=tmp_path / "results",
    )
    batch_dir = comparison_path.parent
    for filename in ("runs.csv", "failures.csv", "comparison.csv"):
        frame = pd.read_csv(batch_dir / filename)
        assert frame.loc[0, "method"] == "nnls"
        assert frame.loc[0, "method_run_id"] == "nnls_expected_failure"
        assert frame.loc[0, "method_config_sha256"] == method_config_sha256(failing_config)
        assert json.loads(frame.loc[0, "method_config"]) == failing_config
        assert frame.loc[0, "status"] == "failed"


def test_metric_contract_failure_marks_completed_method_run_failed(tmp_path):
    registry_path = _write_toy_dataset(tmp_path, declare_cell_types=False)
    experiment_path = tmp_path / "missing_universe.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "missing_universe",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "methods": ["nnls"],
                "metrics": ["rmse_v1"],
                "continue_on_error": True,
            }
        )
    )

    comparison_path = run_experiment(
        experiment_path,
        registry_path,
        output_root_override=tmp_path / "results",
    )
    batch_dir = comparison_path.parent
    assert (batch_dir / "toy__atac__all__nnls" / "results" / "proportions.csv").exists()
    for filename in ("runs.csv", "failures.csv", "comparison.csv"):
        frame = pd.read_csv(batch_dir / filename)
        assert frame.loc[0, "status"] == "failed"
    assert "No declared cell-type universe" in pd.read_csv(batch_dir / "failures.csv").loc[0, "error"]


def test_metric_alias_collision_fails_before_output_creation(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "duplicate_metrics.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "duplicate_metrics",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "methods": ["nnls"],
                "metrics": ["rmse", "rmse_v1"],
            }
        )
    )

    output_root = tmp_path / "results"
    with pytest.raises(ValueError, match="duplicate canonical endpoint.*rmse_v1"):
        run_experiment(experiment_path, registry_path, output_root_override=output_root)
    assert not (output_root / "duplicate_metrics").exists()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_success_records_campaign_code_output_and_resource_provenance(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "provenance.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "provenance",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "methods": ["nnls"],
                "metrics": ["rmse_v1"],
            }
        )
    )

    comparison_path = run_experiment(
        experiment_path,
        registry_path,
        output_root_override=tmp_path / "results",
    )
    batch_dir = comparison_path.parent
    run_dir = batch_dir / "toy__atac__all__nnls"
    metadata = yaml.safe_load((run_dir / "run.yaml").read_text())
    provenance = metadata["execution_provenance"]

    assert metadata["status"] == "success"
    assert metadata["execution_action"] == "executed"
    assert provenance["schema_version"] == 1
    assert provenance["experiment_config"] == {
        "path": str(experiment_path),
        "source_sha256": _sha256(experiment_path),
        "resolved_sha256": provenance["experiment_config"]["resolved_sha256"],
        "resolved_hash_encoding": "canonical_json_sorted_keys_compact_utf8",
    }
    assert len(provenance["experiment_config"]["resolved_sha256"]) == 64
    assert provenance["benchmark_protocol"]["path"] == "docs/ShapeMix/benchmark_protocol.md"
    assert len(provenance["benchmark_protocol"]["source_sha256"]) == 64
    assert provenance["registry"] == {
        "path": str(registry_path),
        "source_sha256": _sha256(registry_path),
    }
    code = provenance["code"]
    canonical_files = json.dumps(
        code["files"],
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    assert code["manifest_sha256"] == hashlib.sha256(canonical_files).hexdigest()
    assert "scripts/run_deconvolution.py" in code["files"]
    assert "scripts/summarize_shapemix_benchmark.py" in code["files"]
    assert "src/deconvatac/methods/shapemix.py" in code["files"]
    assert code["git_commit"]
    assert isinstance(code["git_worktree_dirty"], bool)

    performance = metadata["performance"]
    assert performance["wall_runtime_seconds"] > 0
    assert performance["process_cpu_seconds"] >= 0
    assert performance["average_process_cpu_cores"] >= 0
    assert performance["average_process_cpu_percent_of_one_core"] >= 0
    assert performance["cpu_measurement"]["source"] == "psutil_process_tree_cpu_times"
    assert performance["resource_preflight"]["passed"] is True
    assert performance["peak_rss_bytes"] > 0
    assert performance["peak_rss_mb"] > 0
    assert performance["scope"] == "load_signature_fit_diagnostics_standard_write_and_evaluation"
    assert performance["rss_measurement"]["source"] == "psutil_sampled_process_tree"
    assert performance["rss_measurement"]["ru_maxrss"]["units"] in {
        "bytes",
        "kibibytes",
    }
    assert set(provenance["compute_environment"]["thread_environment"]) == {
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    }
    assert metadata["environment"]["compute"] == provenance["compute_environment"]

    manifest_record = metadata["output_manifest"]
    output_manifest_path = run_dir / manifest_record["path"]
    assert manifest_record["sha256"] == _sha256(output_manifest_path)
    output_manifest = yaml.safe_load(output_manifest_path.read_text())
    assert output_manifest["exclusions"] == ["run.yaml", "output_sha256.yaml"]
    assert set(output_manifest["files"]).issuperset(
        {
            "environment.txt",
            "inputs.yaml",
            "results/diagnostics.json",
            "results/proportions.csv",
            "results/truth.csv",
        }
    )
    for relative, digest in output_manifest["files"].items():
        assert digest == _sha256(run_dir / relative)

    batch_manifest = yaml.safe_load((batch_dir / "batch_manifest.yaml").read_text())
    assert batch_manifest["status"] == "completed"
    assert batch_manifest["execution_provenance"] == provenance
    for table_name in ("runs.csv", "comparison.csv"):
        table = pd.read_csv(batch_dir / table_name)
        assert table.loc[0, "experiment_config_source_sha256"] == _sha256(experiment_path)
        assert table.loc[0, "registry_source_sha256"] == _sha256(registry_path)
        assert table.loc[0, "code_manifest_sha256"] == code["manifest_sha256"]
        assert table.loc[0, "output_manifest_sha256"] == manifest_record["sha256"]
        assert table.loc[0, "wall_runtime_seconds"] > 0
        assert table.loc[0, "peak_rss_bytes"] > 0


def test_resume_reuses_only_fully_validated_success_without_rewriting_run(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "resume.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "resume",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "methods": ["nnls"],
                "metrics": ["rmse_v1"],
            }
        )
    )
    output_root = tmp_path / "results"
    comparison_path = run_experiment(experiment_path, registry_path, output_root)
    run_dir = comparison_path.parent / "toy__atac__all__nnls"
    run_yaml_before = (run_dir / "run.yaml").read_bytes()
    proportions_before = (run_dir / "results" / "proportions.csv").read_bytes()

    resumed_comparison = run_experiment(
        experiment_path,
        registry_path,
        output_root,
        resume=True,
    )

    assert resumed_comparison == comparison_path
    assert (run_dir / "run.yaml").read_bytes() == run_yaml_before
    assert (run_dir / "results" / "proportions.csv").read_bytes() == proportions_before
    runs = pd.read_csv(comparison_path.parent / "runs.csv")
    assert runs.loc[0, "status"] == "success"
    assert runs.loc[0, "execution_action"] == "resumed"


def test_resume_detects_tampered_output_before_batch_mutation(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "tamper.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "tamper",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "methods": ["nnls"],
                "metrics": ["rmse_v1"],
            }
        )
    )
    output_root = tmp_path / "results"
    comparison_path = run_experiment(experiment_path, registry_path, output_root)
    batch_dir = comparison_path.parent
    proportions_path = batch_dir / "toy__atac__all__nnls" / "results" / "proportions.csv"
    proportions_path.write_bytes(proportions_path.read_bytes() + b"\n")
    batch_manifest_before = (batch_dir / "batch_manifest.yaml").read_bytes()

    with pytest.raises(ValueError, match="Cannot resume run.*output SHA256 mismatch"):
        run_experiment(experiment_path, registry_path, output_root, resume=True)
    assert (batch_dir / "batch_manifest.yaml").read_bytes() == batch_manifest_before


def test_resource_gate_is_disabled_without_launcher_environment(monkeypatch):
    import scripts.run_deconvolution as runner

    monkeypatch.delenv(runner.RESOURCE_GUARD_ENV, raising=False)

    assert runner._resource_gate_state() == {"enabled": False, "passed": True}


def test_resource_gate_waits_and_rechecks_until_safe(monkeypatch):
    import scripts.run_deconvolution as runner

    states = iter(
        (
            {"enabled": True, "passed": False, "one_minute_load": 12.0},
            {"enabled": True, "passed": True, "one_minute_load": 1.0},
        )
    )
    sleeps = []
    monkeypatch.setattr(runner, "_resource_gate_state", lambda: next(states))
    monkeypatch.setattr(runner.time, "sleep", sleeps.append)

    selected = runner._wait_for_resource_gate()

    assert selected["passed"] is True
    assert sleeps == [runner.RESOURCE_RECHECK_SECONDS]


def test_failed_run_is_explicitly_finalized_and_never_reused(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "failed_resume.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "failed_resume",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "method_runs": [
                    {
                        "id": "bad",
                        "method": "nnls",
                        "config": {"method": "nnls", "params": {"layer_ref": "missing"}},
                    }
                ],
                "metrics": ["rmse_v1"],
                "continue_on_error": True,
            }
        )
    )
    output_root = tmp_path / "results"
    comparison_path = run_experiment(experiment_path, registry_path, output_root)
    run_dir = comparison_path.parent / "toy__atac__all__bad"
    metadata = yaml.safe_load((run_dir / "run.yaml").read_text())
    assert metadata["status"] == "failed"
    assert metadata["error"]
    assert metadata["performance"]["wall_runtime_seconds"] > 0
    assert (run_dir / "output_sha256.yaml").exists()
    failures = pd.read_csv(comparison_path.parent / "failures.csv")
    assert failures.loc[0, "experiment_config_source_sha256"] == _sha256(experiment_path)
    assert failures.loc[0, "output_manifest_sha256"] == metadata["output_manifest"]["sha256"]

    with pytest.raises(ValueError, match="Cannot resume run.*status is 'failed'"):
        run_experiment(experiment_path, registry_path, output_root, resume=True)


def test_failed_run_persists_and_hashes_method_diagnostics(tmp_path, monkeypatch):
    import scripts.run_deconvolution as runner

    class FailureDiagnostics:
        def to_dict(self):
            return {
                "success": False,
                "stopping_reason": "max_steps_without_convergence",
                "restarts": [{"restart_index": 0, "steps": 17}],
            }

    class FailingMethod:
        def __init__(self, **params):
            del params

        def run(self, data):
            del data
            error = RuntimeError("controlled optimizer failure")
            error.diagnostics = FailureDiagnostics()
            raise error

    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "failure_diagnostics.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "failure_diagnostics",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "methods": ["nnls"],
                "metrics": ["rmse_v1"],
                "continue_on_error": True,
            }
        )
    )
    monkeypatch.setattr(runner, "get_method", lambda _: FailingMethod)

    comparison_path = runner.run_experiment(
        experiment_path,
        registry_path,
        tmp_path / "results",
    )
    run_dir = comparison_path.parent / "toy__atac__all__nnls"
    diagnostics_path = run_dir / "failure_diagnostics.json"
    assert json.loads(diagnostics_path.read_text()) == FailureDiagnostics().to_dict()
    output_manifest = yaml.safe_load((run_dir / "output_sha256.yaml").read_text())
    assert "failure_diagnostics.json" in output_manifest["files"]


def test_resume_and_overwrite_are_mutually_exclusive_before_output(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "modes.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "modes",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "methods": ["nnls"],
            }
        )
    )
    output_root = tmp_path / "results"
    with pytest.raises(ValueError, match="resume and overwrite are mutually exclusive"):
        run_experiment(
            experiment_path,
            registry_path,
            output_root,
            overwrite=True,
            resume=True,
        )
    assert not (output_root / "modes").exists()


def test_resume_requires_stable_configured_run_group(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "no_group.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all"],
                "methods": ["nnls"],
            }
        )
    )
    output_root = tmp_path / "results"
    with pytest.raises(ValueError, match="resume requires experiment config field 'run_group'"):
        run_experiment(experiment_path, registry_path, output_root, resume=True)
    assert not output_root.exists()


def test_deterministic_shards_have_complete_union_no_overlap_and_keep_arms_together(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment = {
        "datasets": ["toy"],
        "modalities": ["atac"],
        "feature_sets": ["all", "highly_variable"],
        "method_runs": [
            {"id": "left", "method": "nnls", "config": {"method": "nnls"}},
            {"id": "right", "method": "nnls", "config": {"method": "nnls"}},
        ],
    }
    jobs = resolve_experiment_jobs(
        experiment,
        registry_path,
        "{dataset}__{modality}__{feature_set}__{method_run_id}",
    )
    all_ids = {job["run_id"] for job in jobs}
    shard_ids: list[set[str]] = []
    arm_shards: dict[tuple[str, str, str], set[int]] = {}
    for index in range(2):
        selected, metadata = select_experiment_shard(
            jobs,
            shard_index=index,
            shard_count=2,
        )
        assert metadata["run_group_suffix"] == f"__shard_0{index}_of_02"
        shard_ids.append({job["run_id"] for job in selected})
        for job in selected:
            unit = (job["dataset"], job["modality"], job["feature_set"])
            arm_shards.setdefault(unit, set()).add(index)
    assert shard_ids[0].isdisjoint(shard_ids[1])
    assert shard_ids[0] | shard_ids[1] == all_ids
    assert all(indices in ({0}, {1}) for indices in arm_shards.values())
    for unit in arm_shards:
        unit_jobs = [
            job for job in jobs if (job["dataset"], job["modality"], job["feature_set"]) == unit
        ]
        assert {job["method_run_id"] for job in unit_jobs} == {"left", "right"}


@pytest.mark.parametrize(
    ("index", "count", "message"),
    [
        (0, None, "provided together"),
        (None, 2, "provided together"),
        (-1, 2, "0 <= shard_index"),
        (2, 2, "0 <= shard_index"),
        (0, 0, "positive integer"),
    ],
)
def test_invalid_shard_arguments_fail_closed(index, count, message):
    jobs = [
        {
            "dataset": "d",
            "modality": "atac",
            "feature_set": "all",
            "run_id": "r",
        }
    ]
    with pytest.raises(ValueError, match=message):
        select_experiment_shard(jobs, shard_index=index, shard_count=count)


def test_sharded_experiment_uses_distinct_groups_and_shared_full_campaign_hash(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "sharded.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "sharded",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["all", "highly_variable"],
                "method_runs": [
                    {"id": "left", "method": "nnls", "config": {"method": "nnls"}},
                    {"id": "right", "method": "nnls", "config": {"method": "nnls"}},
                ],
                "metrics": ["rmse_v1"],
            }
        )
    )
    output_root = tmp_path / "results"
    comparisons = [
        run_experiment(
            experiment_path,
            registry_path,
            output_root,
            shard_index=index,
            shard_count=2,
        )
        for index in range(2)
    ]
    assert [path.parent.name for path in comparisons] == [
        "sharded__shard_00_of_02",
        "sharded__shard_01_of_02",
    ]
    manifests = [
        yaml.safe_load((comparison.parent / "batch_manifest.yaml").read_text())
        for comparison in comparisons
    ]
    assert [manifest["execution_provenance"]["shard"]["index"] for manifest in manifests] == [
        0,
        1,
    ]
    assert len(
        {
            manifest["execution_provenance"]["experiment_config"]["resolved_sha256"]
            for manifest in manifests
        }
    ) == 1
    assert len(
        {
            manifest["execution_provenance"]["code"]["manifest_sha256"]
            for manifest in manifests
        }
    ) == 1
    run_ids = [set(pd.read_csv(path.parent / "runs.csv")["run_id"]) for path in comparisons]
    assert run_ids[0].isdisjoint(run_ids[1])
    assert len(run_ids[0] | run_ids[1]) == 4
