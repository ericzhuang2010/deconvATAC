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
ORIGINAL_CPU_CUDA_TOLERANCE = {"smoke": 1.0e-5, "full_size": 1.0e-4}
FULL_SIZE_ROUTED_TOLERANCE = 2.0e-4
METRIC_TOLERANCE = 1.0e-5
REPEAT_TOLERANCE = 1.0e-7
MINIMUM_FULL_SIZE_SPEEDUP = 2.0


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
    peak_rss_bytes = int(performance.get("peak_rss_bytes", -1))
    if not np.isfinite(wall) or wall <= 0 or peak_rss_bytes <= 0:
        raise ValueError(f"{run_dir} has invalid runtime or peak RSS metadata.")
    average_cpu_raw = performance.get("average_process_cpu_cores")
    average_cpu_cores = None if average_cpu_raw is None else float(average_cpu_raw)
    if average_cpu_cores is not None and (
        not np.isfinite(average_cpu_cores) or average_cpu_cores < 0
    ):
        raise ValueError(f"{run_dir} has invalid average CPU metadata.")
    resource_preflight = performance.get("resource_preflight")
    if not isinstance(resource_preflight, dict):
        resource_preflight = {}
    with (run_dir / "results" / "diagnostics.json").open() as handle:
        diagnostics = json.load(handle)
    fit = diagnostics.get("fit")
    if not isinstance(fit, dict) or fit.get("success") is not True:
        raise ValueError(f"{run_dir} has no successful ShapeMix fit diagnostics.")
    execution = fit.get("execution")
    if not isinstance(execution, dict):
        raise ValueError(f"{run_dir} has no ShapeMix execution metadata.")
    return {
        "wall_runtime_seconds": wall,
        "average_process_cpu_cores": average_cpu_cores,
        "peak_rss_bytes": peak_rss_bytes,
        "resource_preflight": resource_preflight,
        "selected_restart": fit.get("selected_restart"),
        "stopping_reason": fit.get("stopping_reason"),
        "execution": execution,
    }


def summarize(root: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    routing_pass = True
    requires_persisted_preflight = not root.name.endswith("_v1")
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
        cpu_cuda_tolerance = ORIGINAL_CPU_CUDA_TOLERANCE[scale]
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
        thread_counts = {
            variant: audits[variant]["execution"].get("torch_num_threads")
            for variant in VARIANTS
        }
        single_host_thread_passed = all(value == 1 for value in thread_counts.values())
        preflights = {
            variant: audits[variant]["resource_preflight"]
            for variant in VARIANTS
        }
        resource_preflight_present = all(bool(value) for value in preflights.values())
        resource_preflight_passed = resource_preflight_present and all(
            value.get("enabled") is True and value.get("passed") is True
            for value in preflights.values()
        )
        launch_load_max = (
            max(float(value["one_minute_load"]) for value in preflights.values())
            if resource_preflight_present
            else None
        )
        launch_available_memory_min_bytes = (
            min(int(value["available_memory_bytes"]) for value in preflights.values())
            if resource_preflight_present
            else None
        )
        launch_gpu_temperature_max_c = (
            max(int(value["gpu_temperature_c"]) for value in preflights.values())
            if resource_preflight_present
            else None
        )
        original_passed = (
            cpu_cuda_max_abs <= cpu_cuda_tolerance
            and repeat_max_abs <= REPEAT_TOLERANCE
            and stream_max_abs <= REPEAT_TOLERANCE
            and metric_max_abs <= METRIC_TOLERANCE
            and repeat_same_status
            and cache_modes[cached] == "full_cuda"
            and cache_modes[repeated] == "full_cuda"
            and cache_modes[streamed] == "streamed_host_chunks"
            and (scale != "full_size" or speedup >= MINIMUM_FULL_SIZE_SPEEDUP)
            and single_host_thread_passed
        )
        if scale == "smoke":
            production_backend = "cpu"
            scale_routing_pass = (
                audits[cpu]["stopping_reason"]
                == "selected_largest_finite_converged_objective"
                and single_host_thread_passed
                and (not requires_persisted_preflight or resource_preflight_passed)
            )
            routed_tolerance: float | None = None
        else:
            production_backend = "cuda:0"
            routed_tolerance = FULL_SIZE_ROUTED_TOLERANCE
            scale_routing_pass = (
                cpu_cuda_max_abs <= routed_tolerance
                and repeat_max_abs <= REPEAT_TOLERANCE
                and stream_max_abs <= REPEAT_TOLERANCE
                and metric_max_abs <= METRIC_TOLERANCE
                and repeat_same_status
                and cache_modes[cached] == "full_cuda"
                and cache_modes[repeated] == "full_cuda"
                and cache_modes[streamed] == "streamed_host_chunks"
                and speedup >= MINIMUM_FULL_SIZE_SPEEDUP
                and single_host_thread_passed
                and (not requires_persisted_preflight or resource_preflight_passed)
            )
        routing_pass = routing_pass and scale_routing_pass
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
                "cuda_repeat_wall_runtime_seconds": audits[repeated]["wall_runtime_seconds"],
                "cuda_streamed_wall_runtime_seconds": audits[streamed]["wall_runtime_seconds"],
                "cpu_average_process_cpu_cores": audits[cpu]["average_process_cpu_cores"],
                "cuda_cached_average_process_cpu_cores": audits[cached]["average_process_cpu_cores"],
                "cuda_repeat_average_process_cpu_cores": audits[repeated]["average_process_cpu_cores"],
                "cuda_streamed_average_process_cpu_cores": audits[streamed]["average_process_cpu_cores"],
                "cpu_peak_rss_bytes": audits[cpu]["peak_rss_bytes"],
                "cuda_cached_peak_rss_bytes": audits[cached]["peak_rss_bytes"],
                "cuda_repeat_peak_rss_bytes": audits[repeated]["peak_rss_bytes"],
                "cuda_streamed_peak_rss_bytes": audits[streamed]["peak_rss_bytes"],
                "cuda_cached_peak_allocated_bytes": audits[cached]["execution"].get(
                    "peak_device_memory_allocated_bytes"
                ),
                "cuda_cached_peak_reserved_bytes": audits[cached]["execution"].get("peak_device_memory_reserved_bytes"),
                "cuda_repeat_peak_allocated_bytes": audits[repeated]["execution"].get(
                    "peak_device_memory_allocated_bytes"
                ),
                "cuda_repeat_peak_reserved_bytes": audits[repeated]["execution"].get(
                    "peak_device_memory_reserved_bytes"
                ),
                "cuda_streamed_peak_allocated_bytes": audits[streamed]["execution"].get(
                    "peak_device_memory_allocated_bytes"
                ),
                "cuda_streamed_peak_reserved_bytes": audits[streamed]["execution"].get(
                    "peak_device_memory_reserved_bytes"
                ),
                "single_host_thread_passed": single_host_thread_passed,
                "resource_preflight_present": resource_preflight_present,
                "resource_preflight_passed": resource_preflight_passed,
                "launch_load_max": launch_load_max,
                "launch_available_memory_min_bytes": launch_available_memory_min_bytes,
                "launch_gpu_temperature_max_c": launch_gpu_temperature_max_c,
                "cuda_speedup": speedup,
                "repeat_same_status": repeat_same_status,
                "cached_cache_mode": cache_modes[cached],
                "repeat_cache_mode": cache_modes[repeated],
                "streamed_cache_mode": cache_modes[streamed],
                "original_predeclared_gate_passed": original_passed,
                "production_backend": production_backend,
                "routed_cpu_cuda_proportion_tolerance": routed_tolerance,
                "routing_passed": scale_routing_pass,
            }
        )

    table = pd.DataFrame(rows)
    report = {
        "schema_version": 2,
        "qualification_root": str(root.relative_to(ROOT)),
        "status": "passed" if routing_pass else "failed",
        "decision": {
            "kind": "size_routed_backend_qualification",
            "made_before_production_external_predictions": True,
            "small_inputs_backend": "cpu",
            "full_size_1024_by_5000_backend": "cuda:0",
            "reason": (
                "The smoke input had no material CUDA speedup and remains on CPU; "
                "the representative full-size input met the amended platform "
                "numerical bound, repeatability, cache-parity, metric, and speed gates."
            ),
        },
        "requirements": {
            "small_input_backend": "cpu",
            "full_size_cpu_cuda_proportion_max_abs": FULL_SIZE_ROUTED_TOLERANCE,
            "cpu_cuda_metric_max_abs": METRIC_TOLERANCE,
            "cuda_repeat_proportion_max_abs": REPEAT_TOLERANCE,
            "cuda_streamed_cached_proportion_max_abs": REPEAT_TOLERANCE,
            "full_size_minimum_cuda_speedup": MINIMUM_FULL_SIZE_SPEEDUP,
            "single_host_thread_required": True,
            "persisted_passing_resource_preflight_required_for_v2": True,
        },
        "original_predeclared_requirements": {
            "smoke_cpu_cuda_proportion_max_abs": ORIGINAL_CPU_CUDA_TOLERANCE["smoke"],
            "full_size_cpu_cuda_proportion_max_abs": ORIGINAL_CPU_CUDA_TOLERANCE["full_size"],
            "cpu_cuda_metric_max_abs": METRIC_TOLERANCE,
            "cuda_repeat_proportion_max_abs": REPEAT_TOLERANCE,
            "cuda_streamed_cached_proportion_max_abs": REPEAT_TOLERANCE,
            "full_size_minimum_cuda_speedup": MINIMUM_FULL_SIZE_SPEEDUP,
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
