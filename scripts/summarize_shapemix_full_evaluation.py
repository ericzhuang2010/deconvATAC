#!/usr/bin/env python3
"""Create a non-pooled cross-family index of the frozen ShapeMix evaluation."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs/experiments/shapemix_full_evaluation_v1.yaml"
EFFECT_COLUMNS = (
    "stage",
    "family",
    "evidence_class",
    "analysis_unit",
    "context",
    "endpoint",
    "estimate",
    "lower",
    "upper",
    "units",
    "n_units",
    "direction",
    "inferential_status",
)


def project_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def repository_path(path: Path) -> str:
    return path.absolute().relative_to(ROOT.absolute()).as_posix()


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, Mapping):
        raise ValueError(f"Expected a YAML mapping: {path}")
    return dict(value)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def require_file(value: str | Path, label: str) -> Path:
    path = project_path(value)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path


def complete_status(summary: Mapping[str, Any]) -> bool:
    return str(summary.get("status", "")).lower() in {"complete", "passed"}


def base_effect(source: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "stage": source["stage"],
        "family": source["family"],
        "evidence_class": source["evidence_class"],
        "analysis_unit": source["analysis_unit"],
    }


def primary_effects(source: Mapping[str, Any], table: pd.DataFrame) -> list[dict[str, Any]]:
    required = {
        "condition",
        "metric",
        "mean_outer_effect",
        "bootstrap_percentile_95_lower",
        "bootstrap_percentile_95_upper",
        "n_outer_available",
    }
    if not required.issubset(table.columns):
        raise ValueError(f"Primary effect table lacks {sorted(required.difference(table.columns))}")
    return [
        {
            **base_effect(source),
            "context": str(row.condition),
            "endpoint": str(row.metric),
            "estimate": float(row.mean_outer_effect),
            "lower": float(row.bootstrap_percentile_95_lower),
            "upper": float(row.bootstrap_percentile_95_upper),
            "units": "length_minus_count_only",
            "n_units": int(row.n_outer_available),
            "direction": "negative_favors_length",
            "inferential_status": "conditional_one_donor_bootstrap",
        }
        for row in table.itertuples(index=False)
    ]


def nominal_effects(source: Mapping[str, Any], table: pd.DataFrame) -> list[dict[str, Any]]:
    metric_columns = {
        "nominal_rmse_v1_descriptive": "nominal_rmse_v1_descriptive_length_minus_count_only",
        "nominal_jsd_v2_descriptive": "nominal_jsd_v2_descriptive_length_minus_count_only",
        "predicted_off_target_mass": "predicted_off_target_mass_length_minus_count_only",
        "rare_absolute_error": "rare_absolute_error_length_minus_count_only",
    }
    required = {"family", *metric_columns.values()}
    if not required.issubset(table.columns):
        raise ValueError(f"Nominal effect table lacks {sorted(required.difference(table.columns))}")
    rows: list[dict[str, Any]] = []
    for family, group in table.groupby("family", sort=True):
        for endpoint, column in metric_columns.items():
            values = pd.to_numeric(group[column], errors="coerce").dropna()
            if values.empty:
                continue
            rows.append(
                {
                    **base_effect(source),
                    "context": str(family),
                    "endpoint": endpoint,
                    "estimate": float(values.mean()),
                    "lower": np.nan,
                    "upper": np.nan,
                    "units": "length_minus_count_only",
                    "n_units": len(values),
                    "direction": "negative_favors_length",
                    "inferential_status": "descriptive_physical_samples",
                }
            )
    return rows


def donor_effects(source: Mapping[str, Any], table: pd.DataFrame) -> list[dict[str, Any]]:
    required = {
        "condition",
        "metric",
        "donors",
        "mean_length_minus_count_only",
        "ci95_lower",
        "ci95_upper",
    }
    if not required.issubset(table.columns):
        raise ValueError(f"Donor effect table lacks {sorted(required.difference(table.columns))}")
    return [
        {
            **base_effect(source),
            "context": str(row.condition),
            "endpoint": str(row.metric),
            "estimate": float(row.mean_length_minus_count_only),
            "lower": float(row.ci95_lower),
            "upper": float(row.ci95_upper),
            "units": "length_minus_count_only",
            "n_units": int(row.donors),
            "direction": "negative_favors_length",
            "inferential_status": "donor_level_t_interval",
        }
        for row in table.itertuples(index=False)
    ]


def diagnostic_effects(source: Mapping[str, Any], table: pd.DataFrame) -> list[dict[str, Any]]:
    required = {"factor", "level", "metric", "length_minus_count_only"}
    if not required.issubset(table.columns):
        raise ValueError(f"Diagnostic effect table lacks {sorted(required.difference(table.columns))}")
    rows: list[dict[str, Any]] = []
    for (factor, level, metric), group in table.groupby(
        ["factor", "level", "metric"], sort=True
    ):
        values = pd.to_numeric(group["length_minus_count_only"], errors="coerce").dropna()
        if values.empty:
            continue
        rows.append(
            {
                **base_effect(source),
                "context": f"{factor}={level}",
                "endpoint": str(metric),
                "estimate": float(values.mean()),
                "lower": np.nan,
                "upper": np.nan,
                "units": "length_minus_count_only",
                "n_units": len(values),
                "direction": "negative_favors_length",
                "inferential_status": "descriptive_two_seed_diagnostic",
            }
        )
    return rows


def real_spatial_effects(source: Mapping[str, Any], table: pd.DataFrame) -> list[dict[str, Any]]:
    required = {"dataset_id", "pearson_r", "spearman_r", "mean_absolute_difference"}
    if not required.issubset(table.columns):
        raise ValueError(f"Real-spatial effect table lacks {sorted(required.difference(table.columns))}")
    rows: list[dict[str, Any]] = []
    endpoints = {
        "median_cell_type_map_pearson": ("pearson_r", "higher_is_more_similar"),
        "median_cell_type_map_spearman": ("spearman_r", "higher_is_more_similar"),
        "median_cell_type_mean_absolute_difference": (
            "mean_absolute_difference",
            "lower_is_more_similar",
        ),
    }
    for dataset_id, group in table.groupby("dataset_id", sort=True):
        for endpoint, (column, direction) in endpoints.items():
            values = pd.to_numeric(group[column], errors="coerce").dropna()
            if values.empty:
                continue
            rows.append(
                {
                    **base_effect(source),
                    "context": str(dataset_id),
                    "endpoint": endpoint,
                    "estimate": float(values.median()),
                    "lower": np.nan,
                    "upper": np.nan,
                    "units": "paired_map_summary",
                    "n_units": len(values),
                    "direction": direction,
                    "inferential_status": "descriptive_no_exact_truth",
                }
            )
    return rows


ADAPTERS = {
    "primary_outer_effects": primary_effects,
    "nominal_physical": nominal_effects,
    "donor_heldout": donor_effects,
    "diagnostic_stress": diagnostic_effects,
    "real_spatial": real_spatial_effects,
}


def normalize_resources(source: Mapping[str, Any], table: pd.DataFrame) -> pd.DataFrame:
    aliases = {
        "dataset_id": ("dataset_id",),
        "method_id": ("method_id", "method_run_id"),
        "status": ("status",),
        "wall_runtime_seconds": ("wall_runtime_seconds", "runtime_seconds"),
        "peak_rss_bytes": ("peak_rss_bytes", "peak_memory_bytes"),
        "device": ("device",),
        "count_cache_mode": ("count_cache_mode",),
    }
    output = pd.DataFrame(index=table.index)
    output["stage"] = str(source["stage"])
    output["family"] = str(source["family"])
    for target, candidates in aliases.items():
        selected = next((name for name in candidates if name in table.columns), None)
        output[target] = table[selected] if selected else np.nan
    return output


def synthesize(config_path: Path) -> Path:
    config = read_yaml(config_path)
    if config.get("pooling_policy") != "never_pool_across_evidence_classes":
        raise ValueError("The frozen full-evaluation pooling policy changed")
    output = project_path(str(config["output_directory"]))
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite immutable synthesis: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    evidence_rows: list[dict[str, Any]] = []
    effect_rows: list[dict[str, Any]] = []
    resource_frames: list[pd.DataFrame] = []
    try:
        for source_value in config["sources"]:
            source = dict(source_value)
            summary_path = require_file(source["summary"], "campaign summary")
            summary = read_yaml(summary_path)
            if not complete_status(summary):
                raise ValueError(f"Campaign summary is not complete: {summary_path}")
            evidence_rows.append(
                {
                    "stage": source["stage"],
                    "family": source["family"],
                    "evidence_class": source["evidence_class"],
                    "analysis_unit": source["analysis_unit"],
                    "exact_truth": source.get("exact_truth"),
                    "status": summary.get("status"),
                    "claim_limit": source["claim_limit"],
                    "summary_path": repository_path(summary_path),
                    "summary_sha256": digest(summary_path),
                }
            )
            if source.get("effect_adapter"):
                effect_path = require_file(source["effect_table"], "effect table")
                adapter = ADAPTERS[str(source["effect_adapter"])]
                effect_rows.extend(adapter(source, pd.read_csv(effect_path)))
            if source.get("resources"):
                resource_path = require_file(source["resources"], "resource table")
                resource_frames.append(normalize_resources(source, pd.read_csv(resource_path)))

        evidence = pd.DataFrame.from_records(evidence_rows)
        effects = pd.DataFrame.from_records(effect_rows, columns=EFFECT_COLUMNS)
        resources = pd.concat(resource_frames, ignore_index=True)
        evidence.to_csv(temporary / "evidence_table.tsv", sep="\t", index=False)
        effects.to_csv(temporary / "effect_table.tsv", sep="\t", index=False)
        resources.to_csv(temporary / "resource_table.tsv", sep="\t", index=False)
        report = {
            "schema_version": 1,
            "status": "complete",
            "campaign_id": config["campaign_id"],
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "config": repository_path(config_path),
            "config_sha256": digest(config_path),
            "pooling_policy": config["pooling_policy"],
            "effect_definition": config["effect_definition"],
            "evidence_classes_are_not_pooled": True,
            "protocol_v1_result_is_preserved": True,
            "campaigns": len(evidence),
            "effect_rows": len(effects),
            "resource_rows": len(resources),
            "outputs": {
                name: {"sha256": digest(temporary / name)}
                for name in ("evidence_table.tsv", "effect_table.tsv", "resource_table.tsv")
            },
        }
        with (temporary / "evidence_summary.yaml").open("w") as handle:
            yaml.safe_dump(report, handle, sort_keys=False)
        temporary.rename(output)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run synthesis through run_shapemix_low_impact.sh")
    output = synthesize(project_path(args.config))
    print(output / "evidence_summary.yaml")


if __name__ == "__main__":
    main()
