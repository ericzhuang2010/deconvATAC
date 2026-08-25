#!/usr/bin/env python
"""Summarize the frozen one-factor-at-a-time PBMC sensitivity campaign."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import average_precision_score
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXPERIMENT = ROOT / "configs/experiments/shapemix_pbmc_stress_v1.yaml"
DEFAULT_BATCH = (
    ROOT
    / "results/sensitivity/shapemix_pbmc_stress_v1"
    / "shapemix_pbmc_stress_protocol_v1_cuda"
)
METHOD_IDS = ("shapemix_length", "shapemix_count_only", "nnls")
METRICS = ("rmse_v1", "jsd_v2")
GLOBAL_CONTROL_FACTORS = {
    "depth": "keep_1p00",
    "cells": "mean_10",
    "rare_nk": "observed",
    "features": "peaks_5000",
    "bins": "three",
}


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, Mapping):
        raise ValueError(f"Expected a YAML mapping: {path}")
    return dict(value)


def dataset_design(dataset_id: str) -> dict[str, Any]:
    descriptor = read_yaml(
        ROOT / "data/processed/datasets" / dataset_id / "dataset.yaml"
    )
    sensitivity = descriptor.get("sensitivity")
    if not isinstance(sensitivity, Mapping):
        raise ValueError(f"Dataset lacks frozen sensitivity metadata: {dataset_id}")
    result = dict(sensitivity)
    if result.get("campaign") != "shapemix_pbmc_stress_v1":
        raise ValueError(f"Unexpected sensitivity campaign: {dataset_id}")
    return result


def enriched_scores(
    experiment: Mapping[str, Any],
    batch: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    comparison = pd.read_csv(batch / "comparison.csv")
    required = {"dataset_id", "method_run_id", "metric", "value", "status", "run_id"}
    missing = required.difference(comparison.columns)
    if missing:
        raise ValueError(f"comparison.csv is missing columns: {sorted(missing)}")
    failures = comparison[comparison["status"] != "success"].copy()
    scores = comparison[comparison["status"] == "success"].copy()
    scores = scores[
        scores["method_run_id"].isin(METHOD_IDS) & scores["metric"].isin(METRICS)
    ]
    designs = {
        dataset_id: dataset_design(dataset_id)
        for dataset_id in experiment["datasets"]
    }
    unknown = set(scores["dataset_id"]).difference(designs)
    if unknown:
        raise ValueError(f"Unexpected datasets in comparison.csv: {sorted(unknown)}")
    for field in (
        "factor",
        "level",
        "control_level",
        "mixture_seed",
        "reference_variant",
        "mean_cells_per_spot",
        "depth_retain_probability",
        "rare_nk_fraction",
    ):
        scores[field] = scores["dataset_id"].map(
            {dataset: design.get(field) for dataset, design in designs.items()}
        )
    expected = len(experiment["datasets"]) * len(METHOD_IDS) * len(METRICS)
    if len(scores) != expected or scores.duplicated(
        ["dataset_id", "method_run_id", "metric"]
    ).any():
        raise ValueError(f"Expected {expected} unique score rows, found {len(scores)}")
    if not np.isfinite(scores["value"].to_numpy(dtype=float)).all():
        raise ValueError("Sensitivity score values contain non-finite entries.")
    return scores, failures


def paired_shape_count(scores: pd.DataFrame) -> pd.DataFrame:
    index = [
        "dataset_id",
        "factor",
        "level",
        "control_level",
        "mixture_seed",
        "metric",
    ]
    selected = scores[
        scores["method_run_id"].isin(("shapemix_length", "shapemix_count_only"))
    ]
    wide = selected.pivot(index=index, columns="method_run_id", values="value").reset_index()
    wide["length_minus_count_only"] = (
        wide["shapemix_length"] - wide["shapemix_count_only"]
    )
    wide["length_improved"] = wide["length_minus_count_only"] < 0
    return wide


def expanded_factor_scores(scores: pd.DataFrame) -> pd.DataFrame:
    pieces = [scores[scores["factor"] != "anchor"].copy()]
    anchor = scores[scores["factor"] == "anchor"]
    for factor, level in GLOBAL_CONTROL_FACTORS.items():
        control = anchor.copy()
        control["factor"] = factor
        control["level"] = level
        control["control_level"] = level
        pieces.append(control)
    expanded = pd.concat(pieces, ignore_index=True)
    expected_levels = {
        "depth": 4,
        "cells": 4,
        "rare_nk": 4,
        "features": 3,
        "bins": 3,
        "subtype": 2,
        "reference_support": 4,
    }
    observed = expanded.groupby("factor")["level"].nunique().to_dict()
    if observed != expected_levels:
        raise ValueError(f"Expanded factor levels differ from protocol: {observed}")
    return expanded


def factor_level_scores(expanded: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        expanded.groupby(
            ["factor", "level", "control_level", "method_run_id", "metric"],
            sort=True,
        )["value"]
        .agg(["count", "mean", "std", "min", "max"])
        .reset_index()
        .rename(
            columns={
                "count": "mixture_seeds",
                "mean": "mean_value",
                "std": "seed_sd",
                "min": "minimum_value",
                "max": "maximum_value",
            }
        )
    )
    if set(grouped["mixture_seeds"]) != {2}:
        raise ValueError("Every factor level must have both evaluation mixture seeds.")
    return grouped


def control_differences(expanded: pd.DataFrame) -> pd.DataFrame:
    key = ["factor", "mixture_seed", "method_run_id", "metric"]
    controls = expanded[expanded["level"] == expanded["control_level"]]
    if controls.duplicated(key).any():
        raise ValueError("A factor has duplicate controls within seed/method/metric.")
    control_values = controls.set_index(key)["value"]
    rows = expanded.copy()
    rows["control_value"] = [
        control_values.loc[(row.factor, row.mixture_seed, row.method_run_id, row.metric)]
        for row in rows.itertuples(index=False)
    ]
    rows["value_minus_control"] = rows["value"] - rows["control_value"]
    rows["negative_favors_level"] = rows["value_minus_control"] < 0
    return rows


def per_type_scores(
    experiment: Mapping[str, Any],
    batch: Path,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for dataset_id in experiment["datasets"]:
        design = dataset_design(dataset_id)
        for method_id in METHOD_IDS:
            run_dir = batch / f"{dataset_id}__atac__highly_variable__{method_id}"
            run = read_yaml(run_dir / "run.yaml")
            if run.get("status") != "success":
                continue
            truth = pd.read_csv(run_dir / "results/truth.csv", index_col=0)
            predicted = pd.read_csv(run_dir / "results/proportions.csv", index_col=0)
            if not truth.index.equals(predicted.index) or not truth.columns.equals(predicted.columns):
                raise ValueError(f"Truth/prediction axes differ: {run_dir}")
            for cell_type in truth.columns:
                observed = truth[cell_type].to_numpy(dtype=float)
                estimate = predicted[cell_type].to_numpy(dtype=float)
                if not np.isfinite(observed).all() or not np.isfinite(estimate).all():
                    raise ValueError(f"Non-finite truth/prediction: {run_dir}")
                pearson = pearsonr(observed, estimate).statistic if np.std(observed) > 0 and np.std(estimate) > 0 else np.nan
                spearman = spearmanr(observed, estimate).statistic if np.unique(observed).size > 1 and np.unique(estimate).size > 1 else np.nan
                error = estimate - observed
                records.append(
                    {
                        "dataset_id": dataset_id,
                        "factor": design["factor"],
                        "level": design["level"],
                        "mixture_seed": design["mixture_seed"],
                        "method_id": method_id,
                        "cell_type": cell_type,
                        "spots": len(observed),
                        "truth_mean": float(observed.mean()),
                        "predicted_mean": float(estimate.mean()),
                        "mean_bias": float(error.mean()),
                        "mean_absolute_error": float(np.abs(error).mean()),
                        "root_mean_squared_error": float(np.sqrt(np.mean(error**2))),
                        "pearson": float(pearson),
                        "spearman": float(spearman),
                    }
                )
    return pd.DataFrame.from_records(records)


def nk_detection(
    experiment: Mapping[str, Any],
    batch: Path,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    eligible = [
        dataset_id
        for dataset_id in experiment["datasets"]
        if dataset_design(dataset_id)["factor"] in {"anchor", "rare_nk"}
    ]
    for dataset_id in eligible:
        design = dataset_design(dataset_id)
        for method_id in METHOD_IDS:
            run_dir = batch / f"{dataset_id}__atac__highly_variable__{method_id}"
            truth = pd.read_csv(run_dir / "results/truth.csv", index_col=0)["NK"]
            predicted = pd.read_csv(run_dir / "results/proportions.csv", index_col=0)["NK"]
            present = truth.to_numpy(dtype=float) > 0
            called = predicted.to_numpy(dtype=float) >= 0.01
            true_positive = int(np.sum(present & called))
            false_positive = int(np.sum(~present & called))
            false_negative = int(np.sum(present & ~called))
            precision = true_positive / (true_positive + false_positive) if true_positive + false_positive else np.nan
            recall = true_positive / (true_positive + false_negative) if true_positive + false_negative else np.nan
            f1 = 2 * precision * recall / (precision + recall) if np.isfinite(precision) and np.isfinite(recall) and precision + recall else np.nan
            auprc = average_precision_score(present, predicted) if np.unique(present).size == 2 else np.nan
            records.append(
                {
                    "dataset_id": dataset_id,
                    "level": "observed" if design["factor"] == "anchor" else design["level"],
                    "mixture_seed": design["mixture_seed"],
                    "method_id": method_id,
                    "truth_positive_spots": int(present.sum()),
                    "truth_negative_spots": int((~present).sum()),
                    "predicted_threshold": 0.01,
                    "true_positive": true_positive,
                    "false_positive": false_positive,
                    "false_negative": false_negative,
                    "precision": precision,
                    "recall": recall,
                    "f1": f1,
                    "auprc": auprc,
                }
            )
    return pd.DataFrame.from_records(records)


def run_resources(experiment: Mapping[str, Any], batch: Path) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for dataset_id in experiment["datasets"]:
        design = dataset_design(dataset_id)
        for method_id in METHOD_IDS:
            run_dir = batch / f"{dataset_id}__atac__highly_variable__{method_id}"
            run = read_yaml(run_dir / "run.yaml")
            performance = run.get("performance", {})
            with (run_dir / "results/diagnostics.json").open() as handle:
                diagnostics = json.load(handle)
            fit = diagnostics.get("fit", {})
            execution = fit.get("execution", {}) if isinstance(fit, Mapping) else {}
            records.append(
                {
                    "dataset_id": dataset_id,
                    "factor": design["factor"],
                    "level": design["level"],
                    "mixture_seed": design["mixture_seed"],
                    "method_id": method_id,
                    "status": run.get("status"),
                    "wall_runtime_seconds": performance.get("wall_runtime_seconds"),
                    "peak_rss_bytes": performance.get("peak_rss_bytes"),
                    "device": fit.get("device") if isinstance(fit, Mapping) else None,
                    "count_cache_mode": execution.get("count_cache_mode"),
                    "peak_device_memory_allocated_bytes": execution.get("peak_device_memory_allocated_bytes"),
                    "peak_device_memory_reserved_bytes": execution.get("peak_device_memory_reserved_bytes"),
                }
            )
    return pd.DataFrame.from_records(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-config", type=Path, default=DEFAULT_EXPERIMENT)
    parser.add_argument("--batch-dir", type=Path, default=DEFAULT_BATCH)
    parser.add_argument("--overwrite-summary", action="store_true")
    args = parser.parse_args()
    experiment = read_yaml(args.experiment_config)
    batch = args.batch_dir if args.batch_dir.is_absolute() else ROOT / args.batch_dir
    if not batch.is_dir():
        raise FileNotFoundError(batch)
    outputs = (
        "dataset_metric_scores.csv",
        "failed_metric_rows.csv",
        "paired_shape_count_dataset_effects.csv",
        "factor_level_metric_summary.csv",
        "factor_control_differences.csv",
        "per_type_scores.csv",
        "nk_detection.csv",
        "run_resources.csv",
        "evidence_summary.yaml",
    )
    existing = [name for name in outputs if (batch / name).exists()]
    if existing and not args.overwrite_summary:
        raise FileExistsError(f"Refusing to overwrite existing summaries: {existing}")

    scores, failures = enriched_scores(experiment, batch)
    paired = paired_shape_count(scores)
    expanded = expanded_factor_scores(scores)
    level_summary = factor_level_scores(expanded)
    differences = control_differences(expanded)
    per_type = per_type_scores(experiment, batch)
    detection = nk_detection(experiment, batch)
    resources = run_resources(experiment, batch)
    tables = {
        "dataset_metric_scores.csv": scores,
        "failed_metric_rows.csv": failures,
        "paired_shape_count_dataset_effects.csv": paired,
        "factor_level_metric_summary.csv": level_summary,
        "factor_control_differences.csv": differences,
        "per_type_scores.csv": per_type,
        "nk_detection.csv": detection,
        "run_resources.csv": resources,
    }
    for name, table in tables.items():
        table.to_csv(batch / name, index=False)
    report = {
        "schema_version": 1,
        "status": "complete" if failures.empty else "incomplete_with_failures",
        "campaign": str(batch.relative_to(ROOT)),
        "evidence_class": "exact_heldout_pseudospot_truth",
        "scientific_scope": "one-donor conditional diagnostic sensitivity",
        "datasets": len(experiment["datasets"]),
        "jobs": len(experiment["datasets"]) * len(METHOD_IDS),
        "mixture_seeds": [307, 401],
        "methods": list(METHOD_IDS),
        "metrics": list(METRICS),
        "failed_metric_rows": len(failures),
        "interpretation": (
            "One factor changes at a time. Metric-minus-control and "
            "length-minus-count-only values are negative when the first term is better."
        ),
    }
    _summary_path = batch / "evidence_summary.yaml"
    with _summary_path.open("w") as handle:
        yaml.safe_dump(report, handle, sort_keys=False)
    print(_summary_path)
    if not failures.empty:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
