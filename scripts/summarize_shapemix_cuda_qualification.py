#!/usr/bin/env python
"""Summarize and fail closed on ShapeMix CUDA qualification campaigns."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
VARIANTS = (
    "shapemix_length_cpu",
    "shapemix_length_cuda_cached",
    "shapemix_length_cuda_cached_repeat",
    "shapemix_length_cuda_streamed",
)
METRICS = ("rmse_v1", "jsd_v2")


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a mapping.")
    return value


def _run_directories(campaign_dir: Path) -> dict[str, Path]:
    directories: dict[str, Path] = {}
    for variant in VARIANTS:
        matches = sorted(
            path
            for path in campaign_dir.iterdir()
            if path.is_dir() and path.name.endswith(f"__{variant}")
        )
        if len(matches) != 1:
            raise ValueError(
                f"{campaign_dir} must contain exactly one run for {variant}; "
                f"found {len(matches)}."
            )
        directories[variant] = matches[0]
    return directories


def _proportions(run_dir: Path) -> pd.DataFrame:
    frame = pd.read_csv(run_dir / "results" / "proportions.csv", index_col=0)
    if frame.empty or not np.isfinite(frame.to_numpy(dtype=float)).all():
        raise ValueError(f"{run_dir} contains empty or non-finite proportions.")
    return frame


def _max_abs(left: pd.DataFrame, right: pd.DataFrame) -> float:
    if not left.index.equals(right.index) or not left.columns.equals(right.columns):
        raise ValueError("Qualification proportion axes do not match exactly.")
    return float(np.max(np.abs(left.to_numpy(dtype=float) - right.to_numpy(dtype=float))))


def _metric_values(campaign_dir: Path) -> dict[tuple[str, str], float]:
    comparison = pd.read_csv(campaign_dir / "comparison.csv")
    values: dict[tuple[str, str], float] = {}
    for variant in VARIANTS:
        selected = comparison[
            comparison["run_id"].astype(str).str.endswith(f"__{variant}")
        ]
        for metric in METRICS:
            rows = selected[selected["metric"] == metric]
            if len(rows) != 1:
                raise ValueError(
                    f"Expected one {metric} row for {variant} in {campaign_dir}."
                )
            value = float(rows.iloc[0]["value"])
            if not np.isfinite(value):
                raise ValueError(f"Non-finite {metric} for {variant}.")
            values[(variant, metric)] = value
    return values


def _run_audit(run_dir: Path) -> dict[str, Any]:
    metadata = _read_yaml(run_dir / "run.yaml")
    if metadata.get("status") != "success":
        raise ValueError(f"{run_dir} is not a successful finalized run.")
    performance = metadata.get("performance")
    if not isinstance(performance, dict):
        raise ValueError(f"{run_dir} has no performance metadata.")
    wall = float(performance.get("wall_runtime_seconds", np.nan))
    if not np.isfinite(wall) or wall <= 0:
        raise ValueError(f"{run_dir} has invalid wall runtime.")
    with (run_dir / "results" / "diagnostics.json").open() as handle:
        diagnostics = json.load(handle)
    fit = diagnostics.get("fit")
    if not isinstance(fit, dict) or fit.get("success") is not True:
        raise ValueError(f"{run_dir} has no successful ShapeMix fit diagnostics.")
    return {
        "wall_runtime_seconds": wall,
        "selected_restart": fit.get("selected_restart"),
        "stopping_reason": fit.get("stopping_reason"),
        "execution": fit.get("execution", {}),
    }


def summarize(root: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    all_pass = True
    for scale in ("smoke", "full_size"):
        campaign_dir = root / scale
        directories = _run_directories(campaign_dir)
        predictions = {
            variant: _proportions(run_dir)
            for variant, run_dir in directories.items()
        }
        audits = {
            variant: _run_audit(run_dir)
            for variant, run_dir in directories.items()
        }
        metrics = _metric_values(campaign_dir)
        cached = "shapemix_length_cuda_cached"
        repeated = "shapemix_length_cuda_cached_repeat"
        streamed = "shapemix_length_cuda_streamed"
        cpu = "shapemix_length_cpu"
        cpu_cuda_tolerance = 1.0e-5 if scale == "smoke" else 1.0e-4
        cpu_cuda_max_abs = _max_abs(predictions[cpu], predictions[cached])
        repeat_max_abs = _max_abs(predictions[cached], predictions[repeated])
        stream_max_abs = _max_abs(predictions[cached], predictions[streamed])
        metric_max_abs = max(
            abs(metrics[(cpu, metric)] - metrics[(cached, metric)])
            for metric in METRICS
        )
        speedup = (
            audits[cpu]["wall_runtime_seconds"]
            / audits[cached]["wall_runtime_seconds"]
        )
        repeat_same_status = (
            audits[cached]["selected_restart"] == audits[repeated]["selected_restart"]
            and audits[cached]["stopping_reason"] == audits[repeated]["stopping_reason"]
        )
        cache_modes = {
            variant: audits[variant]["execution"].get("count_cache_mode")
            for variant in (cached, repeated, streamed)
        }
        passed = (
            cpu_cuda_max_abs <= cpu_cuda_tolerance
            and repeat_max_abs <= 1.0e-7
            and stream_max_abs <= 1.0e-7
            and metric_max_abs <= 1.0e-5
            and repeat_same_status
            and cache_modes[cached] == "full_cuda"
            and cache_modes[repeated] == "full_cuda"
            and cache_modes[streamed] == "streamed_host_chunks"
            and (scale != "full_size" or speedup >= 2.0)
        )
        all_pass = all_pass and passed
        rows.append(
            {
                "scale": scale,
                "cpu_cuda_proportion_max_abs": cpu_cuda_max_abs,
                "cpu_cuda_proportion_tolerance": cpu_cuda_tolerance,
                "cuda_repeat_proportion_max_abs": repeat_max_abs,
                "cuda_streamed_cached_proportion_max_abs": stream_max_abs,
                "cpu_cuda_metric_max_abs": metric_max_abs,
                "cpu_wall_runtime_seconds": audits[cpu]["wall_runtime_seconds"],
                "cuda_cached_wall_runtime_seconds": audits[cached]["wall_runtime_seconds"],
                "cuda_speedup": speedup,
                "repeat_same_status": repeat_same_status,
                "cached_cache_mode": cache_modes[cached],
                "repeat_cache_mode": cache_modes[repeated],
                "streamed_cache_mode": cache_modes[streamed],
                "passed": passed,
            }
        )

    table = pd.DataFrame(rows)
    report = {
        "schema_version": 1,
        "qualification_root": str(root.relative_to(ROOT)),
        "status": "passed" if all_pass else "failed",
        "requirements": {
            "smoke_cpu_cuda_proportion_max_abs": 1.0e-5,
            "full_size_cpu_cuda_proportion_max_abs": 1.0e-4,
            "cpu_cuda_metric_max_abs": 1.0e-5,
            "cuda_repeat_proportion_max_abs": 1.0e-7,
            "cuda_streamed_cached_proportion_max_abs": 1.0e-7,
            "full_size_minimum_cuda_speedup": 2.0,
        },
        "scales": rows,
    }
    return table, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=ROOT / "results" / "development" / "shapemix_gpu_qualification_v1",
    )
    args = parser.parse_args()
    root = args.root if args.root.is_absolute() else ROOT / args.root
    table, report = summarize(root)
    table.to_csv(root / "qualification_summary.csv", index=False)
    with (root / "qualification_report.yaml").open("w") as handle:
        yaml.safe_dump(report, handle, sort_keys=False)
    print(root / "qualification_report.yaml")
    if report["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
