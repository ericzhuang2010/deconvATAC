import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from deconvatac.metrics import evaluate_proportion_metric
from scripts.summarize_shapemix_benchmark import (
    FROZEN_CELL_TYPES,
    FROZEN_SHAPEMIX_PARAMS,
    PRIMARY_ARMS,
    PRIMARY_CONDITIONS,
    PRIMARY_INNER_MIXTURE_SEEDS,
    PRIMARY_METRICS,
    PRIMARY_OUTER_SPLIT_SEEDS,
    RARE_REFERENCE_TYPES,
    average_inner_effects,
    bootstrap_mean_interval,
    build_paired_effects,
    compute_rare_cell_metrics,
    exact_two_sided_sign_flip_pvalue,
    summarize_benchmark,
    summarize_outer_effects,
    validate_paired_run_contract,
    validate_metric_versions,
    validate_primary_design,
)


def _config(use_shape: bool) -> dict:
    return {
        "method": "shapemix",
        "params": {
            **FROZEN_SHAPEMIX_PARAMS,
            "use_shape": use_shape,
        },
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_hash(value: dict) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def _frozen_records() -> pd.DataFrame:
    rows = []
    for condition in PRIMARY_CONDITIONS:
        for outer in PRIMARY_OUTER_SPLIT_SEEDS:
            for inner in PRIMARY_INNER_MIXTURE_SEEDS:
                dataset_id = f"dataset_{condition}_{outer}_{inner}"
                for arm in PRIMARY_ARMS:
                    rows.append(
                        {
                            "source_run_group": "/tmp/group",
                            "run_id": f"{dataset_id}_{arm}",
                            "dataset_id": dataset_id,
                            "outer_split_seed": outer,
                            "inner_mixture_seed": inner,
                            "condition": condition,
                            "method": "shapemix",
                            "method_run_id": arm,
                            "method_config": _config(arm == "shapemix_length"),
                            "benchmark_scope": "primary",
                            "cell_type_universe": json.dumps(list(FROZEN_CELL_TYPES)),
                            "_cell_types": FROZEN_CELL_TYPES,
                            "status": "success",
                            "error": None,
                            "run_dir": f"/tmp/{dataset_id}_{arm}",
                        }
                    )
    return pd.DataFrame(rows)


def _metric_rows(records: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for record in records.to_dict("records"):
        if record["status"] != "success":
            continue
        base = 0.2 + record["outer_split_seed"] / 100_000 + record["inner_mixture_seed"] / 1_000_000
        for metric in PRIMARY_METRICS:
            rows.append(
                {
                    "run_id": record["run_id"],
                    "metric": metric,
                    "value": base - 0.01 if record["method_run_id"] == "shapemix_length" else base,
                }
            )
    return pd.DataFrame(rows)


def test_primary_grid_requires_exact_two_arm_pairing():
    records = _frozen_records()
    validate_primary_design(records)

    missing_arm = records.drop(records.index[0])
    with pytest.raises(ValueError, match="exactly both ShapeMix arms"):
        validate_primary_design(missing_arm)

    duplicate = pd.concat([records, records.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="exactly both ShapeMix arms"):
        validate_primary_design(duplicate)


def test_paired_contract_rejects_nonprotocol_config_even_when_arms_match():
    records = _frozen_records()
    records["status"] = "failed"  # signature files are irrelevant to this config test
    configs = records["method_config"].map(lambda value: json.loads(json.dumps(value)))
    first_pair = (
        (records["outer_split_seed"] == PRIMARY_OUTER_SPLIT_SEEDS[0])
        & (records["inner_mixture_seed"] == PRIMARY_INNER_MIXTURE_SEEDS[0])
        & (records["condition"] == PRIMARY_CONDITIONS[0])
    )
    for index in records.index[first_pair]:
        configs.iloc[index]["params"]["learning_rate"] = 0.04
    records["method_config"] = configs
    with pytest.raises(ValueError, match="exact frozen protocol-v1"):
        validate_paired_run_contract(records)


def test_signature_hash_must_match_across_conditions_and_inner_seeds(tmp_path):
    records = _frozen_records()
    records["status"] = "failed"
    selected = records[
        (records["outer_split_seed"] == PRIMARY_OUTER_SPLIT_SEEDS[0])
        & (records["method_run_id"] == "shapemix_length")
    ].iloc[[0, -1]]
    for index, (row_index, record) in enumerate(selected.iterrows()):
        run_dir = tmp_path / f"signature_{index}"
        native = run_dir / "results" / "raw_method_output"
        native.mkdir(parents=True)
        (native / "signature_summary.yaml").write_text(
            yaml.safe_dump({"signature": {"content_sha256": str(index + 1) * 64}})
        )
        records.loc[row_index, "run_dir"] = str(run_dir)
        records.loc[row_index, "status"] = "success"
    with pytest.raises(ValueError, match="different signatures across inner seeds or conditions"):
        validate_paired_run_contract(records)


def test_metric_versions_reject_legacy_aliases_and_stale_metadata():
    universe = json.dumps(list(FROZEN_CELL_TYPES), separators=(",", ":"))
    comparison = pd.DataFrame(
        [
            {
                "run_id": "run",
                "status": "success",
                "metric": "rmse_v1",
                "metric_name": "rmse",
                "metric_version": "v1",
                "cell_type_universe": universe,
                "value": 0.1,
            },
            {
                "run_id": "run",
                "status": "success",
                "metric": "jsd_v2",
                "metric_name": "jsd",
                "metric_version": "v2",
                "cell_type_universe": universe,
                "value": 0.2,
            },
        ]
    )
    validate_metric_versions(comparison, expected_cell_types_by_run={"run": FROZEN_CELL_TYPES})

    legacy = comparison.copy()
    legacy.loc[legacy["metric"] == "jsd_v2", "metric"] = "jsd"
    with pytest.raises(ValueError, match="Noncanonical"):
        validate_metric_versions(legacy)

    stale = comparison.copy()
    stale.loc[stale["metric"] == "jsd_v2", "metric_version"] = "v1"
    with pytest.raises(ValueError, match="version v2"):
        validate_metric_versions(stale)


def test_failed_arm_propagates_through_pair_outer_and_formal_summary():
    records = _frozen_records()
    failed = (
        (records["condition"] == "observed_abundance")
        & (records["outer_split_seed"] == PRIMARY_OUTER_SPLIT_SEEDS[0])
        & (records["inner_mixture_seed"] == PRIMARY_INNER_MIXTURE_SEEDS[0])
        & (records["method_run_id"] == "shapemix_length")
    )
    records.loc[failed, "status"] = "failed"
    records.loc[failed, "error"] = "optimizer did not converge"
    paired = build_paired_effects(records, _metric_rows(records))

    unavailable = paired[
        (paired["condition"] == "observed_abundance")
        & (paired["outer_split_seed"] == PRIMARY_OUTER_SPLIT_SEEDS[0])
        & (paired["inner_mixture_seed"] == PRIMARY_INNER_MIXTURE_SEEDS[0])
    ]
    assert set(unavailable["pair_status"]) == {"unavailable"}
    assert unavailable["unavailable_reason"].str.contains("optimizer did not converge").all()

    outer = average_inner_effects(paired)
    affected_outer = outer[
        (outer["condition"] == "observed_abundance")
        & (outer["outer_split_seed"] == PRIMARY_OUTER_SPLIT_SEEDS[0])
    ]
    assert set(affected_outer["outer_status"]) == {"unavailable"}
    assert affected_outer["delta_outer"].isna().all()

    summary = summarize_outer_effects(outer, bootstrap_replicates=100)
    affected_summary = summary[summary["condition"] == "observed_abundance"]
    assert set(affected_summary["analysis_status"]) == {"incomplete"}
    assert affected_summary["mean_outer_effect"].isna().all()
    assert not affected_summary["directional_support_rule_met"].any()


def test_nested_effects_and_resampling_are_deterministic():
    records = _frozen_records()
    paired = build_paired_effects(records, _metric_rows(records))
    assert np.allclose(paired["delta_length_minus_count_only"], -0.01)

    outer = average_inner_effects(paired)
    assert np.allclose(outer["delta_outer"], -0.01)
    assert set(outer["n_inner_available"]) == {2}

    first = bootstrap_mean_interval([-0.04, -0.03, -0.02, -0.01, 0.01], replicates=1000)
    second = bootstrap_mean_interval([-0.04, -0.03, -0.02, -0.01, 0.01], replicates=1000)
    assert first == second
    assert exact_two_sided_sign_flip_pvalue([-1, -2, -3, -4, -5]) == 0.0625

    summary = summarize_outer_effects(outer, bootstrap_replicates=100)
    assert set(summary["analysis_status"]) == {"complete"}
    assert (summary["n_outer_improved"] == 5).all()
    assert summary["directional_support_rule_met"].all()
    assert summary["reporting_scope"].str.contains("not biological replication").all()


def test_rare_cell_threshold_is_inclusive_and_undefined_types_are_retained(tmp_path):
    run_dir = tmp_path / "run"
    (run_dir / "results").mkdir(parents=True)
    spots = [f"spot_{index}" for index in range(4)]
    truth = pd.DataFrame(0.0, index=spots, columns=FROZEN_CELL_TYPES)
    truth["CD14 Mono"] = 1.0
    truth.loc[["spot_1", "spot_3"], "cDC"] = 0.1
    truth.loc[["spot_1", "spot_3"], "CD14 Mono"] = 0.9
    truth.loc["spot_2", "Treg"] = 0.1
    truth.loc["spot_2", "CD14 Mono"] = 0.9

    prediction = pd.DataFrame(0.0, index=spots, columns=FROZEN_CELL_TYPES)
    prediction["CD14 Mono"] = 1.0
    prediction.loc["spot_0", "cDC"] = 0.01
    prediction.loc["spot_0", "CD14 Mono"] = 0.99
    prediction.loc["spot_1", "cDC"] = 0.009
    prediction.loc["spot_1", "CD14 Mono"] = 0.991
    truth.to_csv(run_dir / "results" / "truth.csv")
    prediction.to_csv(run_dir / "results" / "proportions.csv")

    record = {
        "source_run_group": str(tmp_path),
        "run_id": "run",
        "dataset_id": "dataset",
        "outer_split_seed": 1103,
        "inner_mixture_seed": 101,
        "condition": "observed_abundance",
        "method_run_id": "shapemix_length",
        "status": "success",
        "run_dir": str(run_dir),
        "_cell_types": FROZEN_CELL_TYPES,
    }
    rare = compute_rare_cell_metrics(pd.DataFrame([record]))
    cdc = rare[(rare["aggregation"] == "per_type") & (rare["cell_type"] == "cDC")].iloc[0]
    assert cdc["predicted_positives"] == 1  # exactly 0.01 is prediction-positive
    assert cdc["tp"] == 0
    assert cdc["fp"] == 1
    assert cdc["fn"] == 2

    gd_t = rare[(rare["aggregation"] == "per_type") & (rare["cell_type"] == "gdT")].iloc[0]
    assert not bool(gd_t["auprc_defined"])
    assert gd_t["auprc_undefined_reason"] == "truth contains no positive class"
    assert set(rare[rare["aggregation"] == "per_type"]["cell_type"]) == set(
        RARE_REFERENCE_TYPES
    )

    macro = rare[rare["aggregation"] == "macro_evaluable_types"].iloc[0]
    assert macro["n_types_in_aggregate"] == 2  # cDC and Treg have both truth classes
    assert macro["precision_n_types_expected"] == 2
    assert macro["precision_n_types_defined"] == 1
    assert not bool(macro["precision_defined"])
    assert np.isnan(macro["precision"])
    assert "Treg: no predicted positives" in macro["precision_undefined_reason"]
    assert macro["f1_n_types_defined"] == 1
    assert not bool(macro["f1_defined"])
    assert np.isnan(macro["f1"])
    assert macro["auprc_n_types_defined"] == 2
    assert bool(macro["auprc_defined"])


def _write_primary_run_groups(
    tmp_path: Path,
    *,
    all_failed: bool = False,
) -> tuple[list[Path], Path]:
    registry = tmp_path / "datasets.yaml"
    registry.write_text("{}\n")
    groups = [tmp_path / "group_a", tmp_path / "group_b"]
    manifests = {group: [] for group in groups}
    comparisons = {group: [] for group in groups}
    experiment_path = tmp_path / "configs" / "experiments" / "shapemix_primary_ablation.yaml"
    protocol_path = tmp_path / "docs" / "ShapeMix" / "benchmark_protocol.md"
    code_path = tmp_path / "scripts" / "code.py"
    experiment_path.parent.mkdir(parents=True)
    protocol_path.parent.mkdir(parents=True)
    code_path.parent.mkdir(parents=True)
    experiment_path.write_text("protocol: 1\n")
    protocol_path.write_text("# frozen protocol v1\n")
    code_path.write_text("VALUE = 1\n")
    code_files = {"scripts/code.py": _sha256(code_path)}
    shared_provenance = {
        "experiment_config_path": str(experiment_path.relative_to(tmp_path)),
        "experiment_config_source_sha256": _sha256(experiment_path),
        "experiment_config_resolved_sha256": "e" * 64,
        "benchmark_protocol_path": str(protocol_path.relative_to(tmp_path)),
        "benchmark_protocol_sha256": _sha256(protocol_path),
        "registry_path": str(registry.relative_to(tmp_path)),
        "registry_source_sha256": _sha256(registry),
        "git_commit": "1" * 40,
        "git_worktree_dirty": True,
        "code_manifest_sha256": _canonical_hash(code_files),
    }
    truth = pd.DataFrame(
        [[1.0] + [0.0] * 15, [0.0] * 11 + [1.0] + [0.0] * 4],
        index=["spot_0", "spot_1"],
        columns=FROZEN_CELL_TYPES,
    )

    for outer_index, outer in enumerate(PRIMARY_OUTER_SPLIT_SEEDS):
        group = groups[0] if outer_index < 2 else groups[1]
        shard_index = 0 if outer_index < 2 else 1
        group.mkdir(parents=True, exist_ok=True)
        for condition in PRIMARY_CONDITIONS:
            for inner in PRIMARY_INNER_MIXTURE_SEEDS:
                dataset_id = f"dataset_{condition}_{outer}_{inner}"
                dataset_config = {
                    "dataset_id": dataset_id,
                    "benchmark_scope": "primary",
                    "simulation": {
                        "condition": condition,
                        "outer_split_seed": outer,
                        "inner_mixture_seed": inner,
                    },
                    "modalities": {
                        "atac": {"truth": {"cell_types": list(FROZEN_CELL_TYPES)}}
                    },
                }
                for arm in PRIMARY_ARMS:
                    run_id = f"{dataset_id}__{arm}"
                    run_dir = group / run_id
                    native = run_dir / "results" / "raw_method_output"
                    native.mkdir(parents=True)
                    correct = 0.95 if arm == "shapemix_length" else 0.8
                    prediction = truth * correct
                    prediction["CD4 Naive"] += 1.0 - correct
                    truth.to_csv(run_dir / "results" / "truth.csv")
                    prediction.to_csv(run_dir / "results" / "proportions.csv")
                    method_config = _config(arm == "shapemix_length")
                    (run_dir / "inputs.yaml").write_text(
                        yaml.safe_dump(
                            {
                                "dataset_id": dataset_id,
                                "modality": "atac",
                                "dataset_config": dataset_config,
                            }
                        )
                    )
                    (native / "signature_summary.yaml").write_text(
                        yaml.safe_dump({"signature": {"content_sha256": "a" * 64}})
                    )
                    (run_dir / "results" / "diagnostics.json").write_text(
                        json.dumps(
                            {
                                "fit": {
                                    "success": True,
                                    "selected_restart": 0,
                                    "runtime_seconds": 1.5,
                                    "objective": {
                                        "count_log_likelihood": -10.0,
                                        "shape_log_likelihood": -5.0 if arm == "shapemix_length" else 0.0,
                                        "abundance_log_prior": -1.0,
                                        "total_log_objective": -16.0,
                                    },
                                }
                            }
                        )
                    )
                    pd.DataFrame(
                        [
                            {
                                "component": "total_count",
                                "bin_name": None,
                                "root_mean_squared_error": 1.0,
                            }
                        ]
                    ).to_csv(native / "reconstruction_summary.csv", index=False)

                    output_manifest_path = run_dir / "output_sha256.yaml"
                    output_files = {
                        str(path.relative_to(run_dir)): _sha256(path)
                        for path in sorted(run_dir.rglob("*"))
                        if path.is_file()
                    }
                    output_manifest_path.write_text(
                        yaml.safe_dump(
                            {
                                "schema_version": 1,
                                "algorithm": "sha256",
                                "exclusions": ["run.yaml", "output_sha256.yaml"],
                                "files": output_files,
                            },
                            sort_keys=False,
                        )
                    )
                    output_manifest_sha256 = _sha256(output_manifest_path)
                    status = "failed" if all_failed else "success"
                    performance = {
                        "wall_runtime_seconds": 2.5,
                        "peak_rss_bytes": 100 * 1024**2,
                        "peak_rss_mb": 100.0,
                        "rss_measurement": {
                            "source": "test_peak_rss",
                            "semantics": "synthetic test measurement",
                        },
                        "scope": "run_one",
                    }
                    run_yaml = {
                        "run_id": run_id,
                        "dataset_id": dataset_id,
                        "modality": "atac",
                        "feature_set": "highly_variable",
                        "method": "shapemix",
                        "method_run_id": arm,
                        "method_config": method_config,
                        "status": status,
                        "execution_action": "executed",
                        "execution_provenance": {
                            "schema_version": 1,
                            "experiment_config": {
                                "path": shared_provenance["experiment_config_path"],
                                "source_sha256": shared_provenance[
                                    "experiment_config_source_sha256"
                                ],
                                "resolved_sha256": shared_provenance[
                                    "experiment_config_resolved_sha256"
                                ],
                            },
                            "benchmark_protocol": {
                                "path": shared_provenance["benchmark_protocol_path"],
                                "source_sha256": shared_provenance[
                                    "benchmark_protocol_sha256"
                                ],
                            },
                            "registry": {
                                "path": shared_provenance["registry_path"],
                                "source_sha256": shared_provenance[
                                    "registry_source_sha256"
                                ],
                            },
                            "code": {
                                "git_commit": shared_provenance["git_commit"],
                                "git_worktree_dirty": shared_provenance["git_worktree_dirty"],
                                "manifest_sha256": shared_provenance[
                                    "code_manifest_sha256"
                                ],
                                "files": code_files,
                            },
                            "shard": {
                                "index": shard_index,
                                "count": 2,
                                "index_base": 0,
                                "assignment": "sorted_dataset_modality_feature_set_modulo",
                                "unit_key_fields": ["dataset", "modality", "feature_set"],
                                "selected_units": [],
                                "run_group_suffix": f"__shard_0{shard_index}_of_02",
                            },
                        },
                        "performance": performance,
                        "dataset_config_sha256": "d" * 64,
                        "output_manifest": {
                            "path": "output_sha256.yaml",
                            "schema_version": 1,
                            "sha256": output_manifest_sha256,
                        },
                    }
                    metadata_without_manifest = dict(run_yaml)
                    metadata_without_manifest.pop("output_manifest")
                    output_payload = yaml.safe_load(output_manifest_path.read_text())
                    output_payload["run_metadata_sha256"] = _canonical_hash(
                        metadata_without_manifest
                    )
                    output_payload[
                        "run_metadata_hash_encoding"
                    ] = "canonical_json_sorted_keys_compact_utf8"
                    output_manifest_path.write_text(
                        yaml.safe_dump(output_payload, sort_keys=False)
                    )
                    output_manifest_sha256 = _sha256(output_manifest_path)
                    run_yaml["output_manifest"]["sha256"] = output_manifest_sha256
                    (run_dir / "run.yaml").write_text(yaml.safe_dump(run_yaml, sort_keys=False))

                    manifests[group].append(
                        {
                            "run_id": run_id,
                            "dataset_id": dataset_id,
                            "method": "shapemix",
                            "method_run_id": arm,
                            "run_dir": str(run_dir),
                            "status": status,
                            "method_config": json.dumps(
                                method_config, sort_keys=True, separators=(",", ":")
                            ),
                            **shared_provenance,
                            "dataset_config_sha256": "d" * 64,
                            "output_manifest_path": "output_sha256.yaml",
                            "output_manifest_sha256": output_manifest_sha256,
                            "wall_runtime_seconds": performance["wall_runtime_seconds"],
                            "peak_rss_bytes": performance["peak_rss_bytes"],
                            "peak_rss_mb": performance["peak_rss_mb"],
                            "rss_measurement": json.dumps(
                                performance["rss_measurement"],
                                sort_keys=True,
                                separators=(",", ":"),
                            ),
                            "performance_scope": performance["scope"],
                            "execution_action": "executed",
                            "shard_index": shard_index,
                            "shard_count": 2,
                        }
                    )
                    if all_failed:
                        comparisons[group].append(
                            {
                                "run_id": run_id,
                                "status": "failed",
                                "error": "optimizer failed",
                            }
                        )
                    else:
                        for metric in PRIMARY_METRICS:
                            evaluated = evaluate_proportion_metric(
                                metric, truth, prediction, FROZEN_CELL_TYPES
                            )
                            comparisons[group].append(
                                {
                                    "run_id": run_id,
                                    "status": "success",
                                    "metric": evaluated.metric_id,
                                    "metric_name": evaluated.metric_name,
                                    "metric_version": evaluated.metric_version,
                                    "cell_type_universe": json.dumps(
                                        list(FROZEN_CELL_TYPES), separators=(",", ":")
                                    ),
                                    "value": evaluated.value,
                                }
                            )

    for group in groups:
        pd.DataFrame(manifests[group]).to_csv(group / "runs.csv", index=False)
        if all_failed:
            pd.DataFrame(
                [
                    {
                        **row,
                        "method_run_id": next(
                            item["method_run_id"]
                            for item in manifests[group]
                            if item["run_id"] == row["run_id"]
                        ),
                    }
                    for row in comparisons[group]
                ]
            ).to_csv(group / "failures.csv", index=False)
        pd.DataFrame(comparisons[group]).to_csv(group / "comparison.csv", index=False)
    return groups, registry


def test_cli_core_combines_run_groups_and_writes_one_donor_reports(tmp_path):
    groups, registry = _write_primary_run_groups(tmp_path)
    output = tmp_path / "summary"
    manifest = summarize_benchmark(
        groups,
        output,
        registry=registry,
        project_root=tmp_path,
    )

    assert manifest["counts"] == {
        "runs": 40,
        "failed_runs": 0,
        "nested_pairs": 20,
        "unavailable_pair_metrics": 0,
        "unavailable_outer_metrics": 0,
    }
    assert not manifest["biological_replication"]
    assert manifest["donors"] == 1
    assert "cannot establish donor-level" in manifest["scientific_interpretation_limit"]
    assert all(manifest["directional_support_rule_met"].values())

    expected_outputs = {
        "run_metrics.csv",
        "paired_effects.csv",
        "outer_effects.csv",
        "primary_summary.csv",
        "cell_type_metrics.csv",
        "cell_type_paired_effects.csv",
        "rare_cell_metrics.csv",
        "rare_cell_paired_effects.csv",
        "performance.csv",
        "reconstruction.csv",
        "failures.csv",
        "provenance.csv",
        "summary.yaml",
    }
    assert {path.name for path in output.iterdir()} == expected_outputs
    assert len(pd.read_csv(output / "paired_effects.csv")) == 40
    assert len(pd.read_csv(output / "outer_effects.csv")) == 20
    cell_type_effects = pd.read_csv(output / "cell_type_paired_effects.csv")
    assert len(cell_type_effects) == 16 * 4 * (20 + 10)
    assert set(cell_type_effects["effect_level"]) == {"inner_pair", "outer_mean"}
    rare_effects = pd.read_csv(output / "rare_cell_paired_effects.csv")
    assert len(rare_effects) == 7 * 4 * (20 + 10)
    assert set(rare_effects["effect_level"]) == {"inner_pair", "outer_mean"}
    assert set(pd.read_csv(output / "primary_summary.csv")["analysis_status"]) == {"complete"}
    assert set(pd.read_csv(output / "performance.csv")["peak_memory_available"]) == {True}
    assert len(pd.read_csv(output / "provenance.csv")) == 40
    assert manifest["execution_provenance"]["all_run_output_manifests_rehashed"]


def test_all_failed_campaign_writes_explicit_incomplete_summary(tmp_path):
    groups, registry = _write_primary_run_groups(tmp_path, all_failed=True)
    output = tmp_path / "summary"
    manifest = summarize_benchmark(
        groups,
        output,
        registry=registry,
        project_root=tmp_path,
    )

    assert manifest["counts"]["failed_runs"] == 40
    assert manifest["counts"]["unavailable_pair_metrics"] == 40
    assert manifest["counts"]["unavailable_outer_metrics"] == 20
    assert len(pd.read_csv(output / "run_metrics.csv")) == 0
    assert len(pd.read_csv(output / "reconstruction.csv")) == 0
    assert len(pd.read_csv(output / "failures.csv")) == 40
    assert len(pd.read_csv(output / "provenance.csv")) == 40
    primary = pd.read_csv(output / "primary_summary.csv")
    assert set(primary["analysis_status"]) == {"incomplete"}
    assert primary["mean_outer_effect"].isna().all()
    paired = pd.read_csv(output / "paired_effects.csv")
    assert set(paired["pair_status"]) == {"unavailable"}


def test_summary_rehashes_outputs_and_rejects_post_run_tampering(tmp_path):
    groups, registry = _write_primary_run_groups(tmp_path)
    runs = pd.read_csv(groups[0] / "runs.csv")
    proportions = Path(runs.loc[0, "run_dir"]) / "results" / "proportions.csv"
    proportions.write_text(proportions.read_text() + "\n")

    with pytest.raises(ValueError, match="output hash mismatch"):
        summarize_benchmark(
            groups,
            tmp_path / "summary",
            registry=registry,
            project_root=tmp_path,
        )
