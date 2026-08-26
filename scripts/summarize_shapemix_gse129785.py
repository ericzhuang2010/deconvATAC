#!/usr/bin/env python
"""Summarize GSE129785 without treating nominal ratios as exact truth."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import yaml
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXPERIMENT = ROOT / "configs/experiments/shapemix_gse129785_external_v2.yaml"
DEFAULT_BATCH = (
    ROOT
    / "results/external_validation/shapemix_gse129785_v2"
    / "shapemix_gse129785_external_protocol_v2_cpu_log_abundance"
)
METHOD_IDS = ("shapemix_length", "shapemix_count_only", "nnls")
RARE_MAX_NOMINAL = 0.01
RARE_DETECTION_FRACTION = 0.5
T_CELL_TYPES = (
    "Regulatory T cells",
    "Naive CD4 T cells",
    "Memory CD4 T cells",
    "Naive CD8 T cells",
    "Memory CD8 T cells",
)


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, Mapping):
        raise ValueError(f"Expected a YAML mapping: {path}")
    return dict(value)


def run_directory(batch: Path, dataset_id: str, method_id: str) -> Path:
    path = batch / f"{dataset_id}__atac__selected_reference_peaks__{method_id}"
    required = (
        "run.yaml",
        "inputs.yaml",
        "output_sha256.yaml",
        "results/proportions.csv",
        "results/diagnostics.json",
    )
    missing = [relative for relative in required if not (path / relative).is_file()]
    if missing:
        raise FileNotFoundError(f"{path} is missing {missing}")
    run = read_yaml(path / "run.yaml")
    if run.get("status") != "success":
        raise ValueError(f"Run is not successful: {path}")
    return path


def read_prediction(path: Path, declared_types: list[str]) -> pd.DataFrame:
    frame = pd.read_csv(path / "results/proportions.csv", index_col=0)
    if frame.empty or frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError(f"Invalid prediction axes: {path}")
    if frame.columns.tolist() != declared_types:
        raise ValueError(f"Prediction cell-type axis changed: {path}")
    values = frame.to_numpy(dtype=np.float64)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.allclose(values.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-6)
    ):
        raise ValueError(f"Invalid prediction values: {path}")
    return frame


def collapse_prediction(
    row: pd.Series,
    family: str,
) -> tuple[tuple[str, str], np.ndarray]:
    if family == "cd4_memory_cd8_naive":
        names = ("CD4 Memory", "CD8 Naive")
        values = np.asarray(
            [row["Memory CD4 T cells"], row["Naive CD8 T cells"]],
            dtype=np.float64,
        )
    elif family == "monocyte_t_cell":
        names = ("Monocytes", "T cells")
        values = np.asarray(
            [row["Monocytes"], row[list(T_CELL_TYPES)].sum()],
            dtype=np.float64,
        )
    else:
        raise ValueError(f"Unknown physical-dilution family: {family}")
    off_target = max(0.0, 1.0 - float(values.sum()))
    return names, np.asarray([values[0], values[1], off_target], dtype=np.float64)


def score_nominal_datasets(
    experiment: Mapping[str, Any],
    batch: Path,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for dataset_id in experiment["datasets"]:
        descriptor_path = ROOT / "data/processed/datasets" / dataset_id / "dataset.yaml"
        descriptor = read_yaml(descriptor_path)
        physical = descriptor.get("physical_dilution")
        if not isinstance(physical, Mapping):
            continue
        validation = descriptor.get("validation", {}).get(
            "nominal_broad_proportions", {}
        )
        nominal_path = ROOT / str(validation.get("path", ""))
        if validation.get("exact_truth") is not False:
            raise ValueError(f"Nominal evidence is not explicitly non-exact: {dataset_id}")
        nominal = pd.read_csv(nominal_path, index_col=0)
        if nominal.shape != (1, 2):
            raise ValueError(f"Expected one two-component nominal row: {nominal_path}")
        target_components = nominal.iloc[0].to_numpy(dtype=np.float64)
        if (
            not np.isfinite(target_components).all()
            or (target_components < 0).any()
            or not np.isclose(target_components.sum(), 1.0, rtol=0.0, atol=1.0e-12)
        ):
            raise ValueError(f"Invalid nominal proportions: {nominal_path}")
        family = str(physical["family"])
        declared_types = list(descriptor["modalities"]["atac"]["cell_types"])
        for method_id in METHOD_IDS:
            run_dir = run_directory(batch, dataset_id, method_id)
            prediction = read_prediction(run_dir, declared_types)
            if prediction.shape[0] != 1 or prediction.index.tolist() != nominal.index.tolist():
                raise ValueError(f"Nominal and predicted sample axes differ: {dataset_id}")
            component_names, collapsed = collapse_prediction(prediction.iloc[0], family)
            if list(nominal.columns) != list(component_names):
                raise ValueError(f"Nominal component axis changed: {dataset_id}")
            target = np.asarray(
                [target_components[0], target_components[1], 0.0],
                dtype=np.float64,
            )
            rmse = float(np.sqrt(np.mean(np.square(target - collapsed))))
            jsd = float(
                np.square(jensenshannon(target, collapsed, base=2))
            )
            rare_index = int(np.argmin(target_components))
            rare_nominal = float(target_components[rare_index])
            rare_predicted = float(collapsed[rare_index])
            is_rare_level = rare_nominal <= RARE_MAX_NOMINAL
            records.append(
                {
                    "evidence_class": (
                        "nominal_sample_level"
                        if family == "cd4_memory_cd8_naive"
                        else "broad_nominal_sample_level"
                    ),
                    "dataset_id": dataset_id,
                    "sample": nominal.index[0],
                    "family": family,
                    "method_id": method_id,
                    "component_1": component_names[0],
                    "nominal_component_1": target_components[0],
                    "predicted_component_1": collapsed[0],
                    "component_2": component_names[1],
                    "nominal_component_2": target_components[1],
                    "predicted_component_2": collapsed[1],
                    "predicted_off_target_mass": collapsed[2],
                    "nominal_rmse_v1_descriptive": rmse,
                    "nominal_jsd_v2_descriptive": jsd,
                    "rare_level": is_rare_level,
                    "rare_component": component_names[rare_index] if is_rare_level else None,
                    "rare_nominal": rare_nominal if is_rare_level else np.nan,
                    "rare_predicted": rare_predicted if is_rare_level else np.nan,
                    "rare_absolute_error": (
                        abs(rare_predicted - rare_nominal) if is_rare_level else np.nan
                    ),
                    "rare_detected_at_half_nominal": (
                        rare_predicted >= RARE_DETECTION_FRACTION * rare_nominal
                        if is_rare_level
                        else None
                    ),
                    "run_dir": str(run_dir.relative_to(ROOT)),
                }
            )
    frame = pd.DataFrame.from_records(records)
    expected_rows = 14 * len(METHOD_IDS)
    if len(frame) != expected_rows:
        raise ValueError(f"Expected {expected_rows} nominal rows, found {len(frame)}")
    return frame


def nominal_series_summary(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (family, method_id), group in scores.groupby(["family", "method_id"], sort=True):
        base = {
            "family": family,
            "method_id": method_id,
            "samples": len(group),
            "mean_nominal_rmse_v1_descriptive": group[
                "nominal_rmse_v1_descriptive"
            ].mean(),
            "mean_nominal_jsd_v2_descriptive": group[
                "nominal_jsd_v2_descriptive"
            ].mean(),
            "mean_predicted_off_target_mass": group[
                "predicted_off_target_mass"
            ].mean(),
            "rare_levels": int(group["rare_level"].sum()),
            "rare_detected_at_half_nominal": int(
                group.loc[group["rare_level"], "rare_detected_at_half_nominal"].sum()
            ),
        }
        for index in (1, 2):
            nominal = group[f"nominal_component_{index}"].to_numpy(dtype=float)
            predicted = group[f"predicted_component_{index}"].to_numpy(dtype=float)
            slope, intercept = np.polyfit(nominal, predicted, 1)
            rank = spearmanr(nominal, predicted).statistic
            rows.append(
                {
                    **base,
                    "component": group[f"component_{index}"].iloc[0],
                    "calibration_slope": float(slope),
                    "calibration_intercept": float(intercept),
                    "spearman_monotonicity": float(rank),
                }
            )
    return pd.DataFrame.from_records(rows)


def paired_nominal_effects(scores: pd.DataFrame) -> pd.DataFrame:
    index = ["dataset_id", "sample", "family"]
    values = [
        "nominal_rmse_v1_descriptive",
        "nominal_jsd_v2_descriptive",
        "predicted_off_target_mass",
        "rare_absolute_error",
    ]
    selected = scores[scores["method_id"].isin(["shapemix_length", "shapemix_count_only"])]
    wide = selected.pivot(index=index, columns="method_id", values=values)
    records = []
    for key, row in wide.iterrows():
        record = dict(zip(index, key))
        for value in values:
            length = row[(value, "shapemix_length")]
            count = row[(value, "shapemix_count_only")]
            record[f"{value}_length"] = length
            record[f"{value}_count_only"] = count
            record[f"{value}_length_minus_count_only"] = length - count
        records.append(record)
    return pd.DataFrame.from_records(records)


def prediction_only_tables(
    experiment: Mapping[str, Any],
    batch: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    long_records: list[dict[str, Any]] = []
    for dataset_id in experiment["datasets"]:
        descriptor = read_yaml(
            ROOT / "data/processed/datasets" / dataset_id / "dataset.yaml"
        )
        if isinstance(descriptor.get("physical_dilution"), Mapping):
            continue
        declared_types = list(descriptor["modalities"]["atac"]["cell_types"])
        design = descriptor.get("evaluation_design", {})
        comparison_key = str(design.get("comparison_key", "sample"))
        sample_labels = {}
        for sample in design.get("samples", []):
            if isinstance(sample, Mapping) and "gsm" in sample:
                sample_labels[str(sample["gsm"])] = str(
                    sample.get(comparison_key, sample["gsm"])
                )
        for method_id in METHOD_IDS:
            run_dir = run_directory(batch, dataset_id, method_id)
            prediction = read_prediction(run_dir, declared_types)
            for sample, row in prediction.iterrows():
                for cell_type, value in row.items():
                    long_records.append(
                        {
                            "evidence_class": "prediction_only",
                            "dataset_id": dataset_id,
                            "sample": sample,
                            "sample_group": sample_labels.get(str(sample), str(sample)),
                            "comparison_key": comparison_key,
                            "method_id": method_id,
                            "cell_type": cell_type,
                            "predicted_proportion": float(value),
                            "run_dir": str(run_dir.relative_to(ROOT)),
                        }
                    )
    long = pd.DataFrame.from_records(long_records)
    if long.empty:
        raise ValueError("No prediction-only records were found")
    summary = (
        long.groupby(["dataset_id", "method_id", "cell_type"], sort=True)[
            "predicted_proportion"
        ]
        .agg(["count", "mean", "std", "min", "max"])
        .reset_index()
        .rename(
            columns={
                "count": "samples",
                "mean": "mean_prediction",
                "std": "sample_sd",
                "min": "minimum_prediction",
                "max": "maximum_prediction",
            }
        )
    )
    paired = long[long["method_id"].isin(["shapemix_length", "shapemix_count_only"])].pivot(
        index=["dataset_id", "sample", "sample_group", "cell_type"],
        columns="method_id",
        values="predicted_proportion",
    )
    paired = paired.reset_index()
    paired["length_minus_count_only"] = (
        paired["shapemix_length"] - paired["shapemix_count_only"]
    )
    paired["absolute_length_count_difference"] = paired[
        "length_minus_count_only"
    ].abs()
    return long, summary, paired


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
    outputs = {
        "nominal_sample_scores.csv": None,
        "nominal_series_summary.csv": None,
        "paired_shape_count_nominal.csv": None,
        "prediction_only_proportions.csv": None,
        "prediction_only_summary.csv": None,
        "prediction_only_shape_count_differences.csv": None,
        "evidence_summary.yaml": None,
    }
    existing = [name for name in outputs if (batch / name).exists()]
    if existing and not args.overwrite_summary:
        raise FileExistsError(f"Refusing to overwrite existing summaries: {existing}")

    nominal = score_nominal_datasets(experiment, batch)
    series = nominal_series_summary(nominal)
    paired_nominal = paired_nominal_effects(nominal)
    prediction_long, prediction_summary, prediction_paired = prediction_only_tables(
        experiment, batch
    )
    nominal.to_csv(batch / "nominal_sample_scores.csv", index=False)
    series.to_csv(batch / "nominal_series_summary.csv", index=False)
    paired_nominal.to_csv(batch / "paired_shape_count_nominal.csv", index=False)
    prediction_long.to_csv(batch / "prediction_only_proportions.csv", index=False)
    prediction_summary.to_csv(batch / "prediction_only_summary.csv", index=False)
    prediction_paired.to_csv(
        batch / "prediction_only_shape_count_differences.csv",
        index=False,
    )
    report = {
        "schema_version": 1,
        "status": "complete",
        "campaign": str(batch.relative_to(ROOT)),
        "evidence_classes": {
            "nominal_sample_level": {
                "samples": 7,
                "family": "cd4_memory_cd8_naive",
                "claim_limit": "descriptive comparison to nominal physical input ratios",
            },
            "broad_nominal_sample_level": {
                "samples": 7,
                "family": "monocyte_t_cell",
                "claim_limit": "broad monocyte-versus-total-T contrast only",
            },
            "prediction_only": {
                "datasets": 2,
                "claim_limit": "replicate/preparation robustness without accuracy claims",
            },
        },
        "rare_detection_rule": {
            "eligible_nominal_fraction_max": RARE_MAX_NOMINAL,
            "detected_if_predicted_at_least_fraction_of_nominal": (
                RARE_DETECTION_FRACTION
            ),
        },
        "methods": list(METHOD_IDS),
        "nominal_score_rows": len(nominal),
        "prediction_only_long_rows": len(prediction_long),
    }
    with (batch / "evidence_summary.yaml").open("w") as handle:
        yaml.safe_dump(report, handle, sort_keys=False)
    print(batch / "evidence_summary.yaml")


if __name__ == "__main__":
    main()
