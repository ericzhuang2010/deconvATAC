#!/usr/bin/env python
"""Summarize GSE194122 at the held-out-donor analysis level."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import yaml
from scipy.stats import t as student_t


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXPERIMENT = ROOT / "configs/experiments/shapemix_gse194122_lodo.yaml"
DEFAULT_BATCH = (
    ROOT
    / "results/external_validation/shapemix_gse194122_lodo_v1"
    / "shapemix_gse194122_broad7_lodo_protocol_v1_cuda"
)
METHOD_IDS = ("shapemix_length", "shapemix_count_only", "nnls")
METRICS = ("rmse_v1", "jsd_v2")


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, Mapping):
        raise ValueError(f"Expected a YAML mapping: {path}")
    return dict(value)


def dataset_design(dataset_id: str) -> dict[str, Any]:
    descriptor = read_yaml(
        ROOT / "data/processed/datasets" / dataset_id / "dataset.yaml"
    )
    simulation = descriptor.get("simulation")
    if not isinstance(simulation, Mapping):
        raise ValueError(f"Dataset lacks simulation metadata: {dataset_id}")
    result = {
        "donor": int(simulation["outer_split_seed"]),
        "condition": str(simulation["condition"]),
        "mixture_seed": int(simulation["inner_mixture_seed"]),
    }
    if result["donor"] not in range(1, 11):
        raise ValueError(f"Invalid held-out donor: {dataset_id}")
    if result["condition"] not in {"observed_abundance", "equal_celltype"}:
        raise ValueError(f"Invalid mixture condition: {dataset_id}")
    if result["mixture_seed"] not in {101, 211}:
        raise ValueError(f"Invalid mixture seed: {dataset_id}")
    return result


def enriched_comparison(
    experiment: Mapping[str, Any],
    batch: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    comparison = pd.read_csv(batch / "comparison.csv")
    required = {"dataset_id", "method_run_id", "metric", "value", "status", "run_id"}
    missing = required.difference(comparison.columns)
    if missing:
        raise ValueError(f"comparison.csv is missing columns: {sorted(missing)}")
    failures = comparison[comparison["status"] != "success"].copy()
    success = comparison[comparison["status"] == "success"].copy()
    designs = {
        dataset_id: dataset_design(dataset_id)
        for dataset_id in experiment["datasets"]
    }
    unknown = set(success["dataset_id"]).difference(designs)
    if unknown:
        raise ValueError(f"Unexpected datasets in comparison.csv: {sorted(unknown)}")
    for field in ("donor", "condition", "mixture_seed"):
        success[field] = success["dataset_id"].map(
            {dataset: design[field] for dataset, design in designs.items()}
        )
    success = success[success["method_run_id"].isin(METHOD_IDS)]
    success = success[success["metric"].isin(METRICS)]
    expected = len(experiment["datasets"]) * len(METHOD_IDS) * len(METRICS)
    if len(success) != expected or success.duplicated(
        ["dataset_id", "method_run_id", "metric"]
    ).any():
        raise ValueError(
            f"Expected {expected} unique core metric rows, found {len(success)}"
        )
    if not np.isfinite(success["value"].to_numpy(dtype=float)).all():
        raise ValueError("Core metric values contain non-finite entries")
    return success, failures


def paired_effects(scores: pd.DataFrame) -> pd.DataFrame:
    index = ["dataset_id", "donor", "condition", "mixture_seed", "metric"]
    selected = scores[
        scores["method_run_id"].isin(["shapemix_length", "shapemix_count_only"])
    ]
    wide = selected.pivot(index=index, columns="method_run_id", values="value")
    wide = wide.reset_index()
    wide["length_minus_count_only"] = (
        wide["shapemix_length"] - wide["shapemix_count_only"]
    )
    wide["length_improved"] = wide["length_minus_count_only"] < 0
    return wide


def donor_effects(paired: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        paired.groupby(["donor", "condition", "metric"], sort=True)[
            "length_minus_count_only"
        ]
        .agg(["count", "mean", "std"])
        .reset_index()
        .rename(
            columns={
                "count": "inner_mixtures",
                "mean": "mean_length_minus_count_only",
                "std": "inner_mixture_sd",
            }
        )
    )
    overall = (
        paired.groupby(["donor", "metric"], sort=True)["length_minus_count_only"]
        .agg(["count", "mean", "std"])
        .reset_index()
        .rename(
            columns={
                "count": "inner_mixtures",
                "mean": "mean_length_minus_count_only",
                "std": "inner_mixture_sd",
            }
        )
    )
    overall.insert(1, "condition", "all")
    return pd.concat([grouped, overall], ignore_index=True)


def donor_level_summary(donors: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (condition, metric), group in donors.groupby(
        ["condition", "metric"], sort=True
    ):
        values = group["mean_length_minus_count_only"].to_numpy(dtype=float)
        n = len(values)
        if n != 10:
            raise ValueError(
                f"Expected ten donor effects for {condition}/{metric}, found {n}"
            )
        mean = float(values.mean())
        sd = float(values.std(ddof=1))
        half_width = float(
            student_t.ppf(0.975, df=n - 1) * sd / math.sqrt(n)
        )
        rows.append(
            {
                "condition": condition,
                "metric": metric,
                "donors": n,
                "mean_length_minus_count_only": mean,
                "sample_sd": sd,
                "median_length_minus_count_only": float(np.median(values)),
                "ci95_lower": mean - half_width,
                "ci95_upper": mean + half_width,
                "donors_length_improved": int((values < 0).sum()),
                "interpretation": "negative favors shapemix_length",
            }
        )
    return pd.DataFrame.from_records(rows)


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
            predicted = pd.read_csv(
                run_dir / "results/proportions.csv",
                index_col=0,
            )
            if (
                not truth.index.equals(predicted.index)
                or not truth.columns.equals(predicted.columns)
                or truth.empty
            ):
                raise ValueError(f"Truth/prediction axes differ: {run_dir}")
            truth_values = truth.to_numpy(dtype=float)
            predicted_values = predicted.to_numpy(dtype=float)
            if not np.isfinite(truth_values).all() or not np.isfinite(
                predicted_values
            ).all():
                raise ValueError(f"Non-finite truth/prediction: {run_dir}")
            absolute = np.abs(predicted_values - truth_values)
            for index, cell_type in enumerate(truth.columns):
                records.append(
                    {
                        "dataset_id": dataset_id,
                        **design,
                        "method_id": method_id,
                        "cell_type": cell_type,
                        "spots": truth.shape[0],
                        "mean_absolute_error": float(absolute[:, index].mean()),
                        "truth_mean": float(truth_values[:, index].mean()),
                        "predicted_mean": float(predicted_values[:, index].mean()),
                        "mean_bias": float(
                            (predicted_values[:, index] - truth_values[:, index]).mean()
                        ),
                    }
                )
    frame = pd.DataFrame.from_records(records)
    expected = len(experiment["datasets"]) * len(METHOD_IDS) * 7
    if len(frame) != expected:
        raise ValueError(f"Expected {expected} per-type rows, found {len(frame)}")
    return frame


def donor_per_type_effects(per_type: pd.DataFrame) -> pd.DataFrame:
    index = ["dataset_id", "donor", "condition", "mixture_seed", "cell_type"]
    selected = per_type[
        per_type["method_id"].isin(["shapemix_length", "shapemix_count_only"])
    ]
    wide = selected.pivot(
        index=index,
        columns="method_id",
        values="mean_absolute_error",
    ).reset_index()
    wide["mae_length_minus_count_only"] = (
        wide["shapemix_length"] - wide["shapemix_count_only"]
    )
    return (
        wide.groupby(["donor", "condition", "cell_type"], sort=True)[
            "mae_length_minus_count_only"
        ]
        .mean()
        .reset_index()
    )


def run_resources(
    experiment: Mapping[str, Any],
    batch: Path,
) -> pd.DataFrame:
    rows = []
    for dataset_id in experiment["datasets"]:
        design = dataset_design(dataset_id)
        for method_id in METHOD_IDS:
            run_dir = batch / f"{dataset_id}__atac__highly_variable__{method_id}"
            run = read_yaml(run_dir / "run.yaml")
            performance = run.get("performance", {})
            diagnostics_path = run_dir / "results/diagnostics.json"
            with diagnostics_path.open() as handle:
                diagnostics = json.load(handle)
            fit = diagnostics.get("fit", {})
            execution = fit.get("execution", {}) if isinstance(fit, Mapping) else {}
            rows.append(
                {
                    "dataset_id": dataset_id,
                    **design,
                    "method_id": method_id,
                    "status": run.get("status"),
                    "wall_runtime_seconds": performance.get(
                        "wall_runtime_seconds"
                    ),
                    "peak_rss_bytes": performance.get("peak_rss_bytes"),
                    "device": fit.get("device") if isinstance(fit, Mapping) else None,
                    "count_cache_mode": execution.get("count_cache_mode"),
                    "peak_device_memory_allocated_bytes": execution.get(
                        "peak_device_memory_allocated_bytes"
                    ),
                    "peak_device_memory_reserved_bytes": execution.get(
                        "peak_device_memory_reserved_bytes"
                    ),
                    "run_dir": str(run_dir.relative_to(ROOT)),
                }
            )
    return pd.DataFrame.from_records(rows)


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
    output_names = (
        "dataset_metric_scores.csv",
        "failed_metric_rows.csv",
        "paired_shape_count_dataset_effects.csv",
        "donor_metric_effects.csv",
        "donor_level_summary.csv",
        "per_type_absolute_error.csv",
        "donor_per_type_effects.csv",
        "run_resources.csv",
        "evidence_summary.yaml",
    )
    existing = [name for name in output_names if (batch / name).exists()]
    if existing and not args.overwrite_summary:
        raise FileExistsError(f"Refusing to overwrite existing summaries: {existing}")

    scores, failures = enriched_comparison(experiment, batch)
    paired = paired_effects(scores)
    donors = donor_effects(paired)
    summary = donor_level_summary(donors)
    per_type = per_type_scores(experiment, batch)
    donor_types = donor_per_type_effects(per_type)
    resources = run_resources(experiment, batch)

    scores.to_csv(batch / "dataset_metric_scores.csv", index=False)
    failures.to_csv(batch / "failed_metric_rows.csv", index=False)
    paired.to_csv(batch / "paired_shape_count_dataset_effects.csv", index=False)
    donors.to_csv(batch / "donor_metric_effects.csv", index=False)
    summary.to_csv(batch / "donor_level_summary.csv", index=False)
    per_type.to_csv(batch / "per_type_absolute_error.csv", index=False)
    donor_types.to_csv(batch / "donor_per_type_effects.csv", index=False)
    resources.to_csv(batch / "run_resources.csv", index=False)

    report = {
        "schema_version": 1,
        "status": "complete" if failures.empty else "incomplete_with_failures",
        "campaign": str(batch.relative_to(ROOT)),
        "evidence_class": "external_donor_heldout_exact_truth",
        "analysis_unit": "donor",
        "donors": 10,
        "conditions": ["observed_abundance", "equal_celltype"],
        "inner_mixture_seeds": [101, 211],
        "methods": list(METHOD_IDS),
        "metrics": list(METRICS),
        "metric_rows": len(scores),
        "failed_metric_rows": len(failures),
        "interpretation": (
            "length-minus-count-only metric effects are averaged over inner "
            "mixtures within donor; negative values favor fragment length"
        ),
    }
    with (batch / "evidence_summary.yaml").open("w") as handle:
        yaml.safe_dump(report, handle, sort_keys=False)
    print(batch / "evidence_summary.yaml")
    if not failures.empty:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
