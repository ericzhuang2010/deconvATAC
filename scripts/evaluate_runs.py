#!/usr/bin/env python
"""Evaluate standardized runs with the frozen proportion-metric contract."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.metrics import (
    available_proportion_metrics,
    declared_cell_types_from_metadata,
    evaluate_proportion_metric,
    resolve_proportion_metric,
)


DEFAULT_METRICS = ("rmse_v1", "jsd_v2")


def read_yaml(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise TypeError(f"{path} must contain a YAML mapping.")
    return value


def find_truth(run_dir: Path, explicit_truth: Optional[str]) -> Path:
    if explicit_truth is not None:
        return Path(explicit_truth)
    local_truth = run_dir / "results" / "truth.csv"
    if local_truth.exists():
        return local_truth
    raise FileNotFoundError(
        f"No truth file found for {run_dir}; pass --truth or include results/truth.csv."
    )


def evaluate_run_directory(
    run_dir: Path,
    metric_names: Sequence[str],
    *,
    explicit_truth: Optional[str] = None,
    cell_types: Optional[Sequence[str]] = None,
) -> list[dict]:
    """Return versioned metric rows for one standardized run directory."""

    metadata = read_yaml(run_dir / "run.yaml")
    inputs_metadata = read_yaml(run_dir / "inputs.yaml")
    declared_types = declared_cell_types_from_metadata(
        metadata,
        inputs_metadata,
        override=cell_types,
    )
    predicted = pd.read_csv(run_dir / "results" / "proportions.csv", index_col=0)
    truth = pd.read_csv(find_truth(run_dir, explicit_truth), index_col=0)
    universe_json = json.dumps(list(declared_types), separators=(",", ":"), ensure_ascii=False)

    rows: list[dict] = []
    for metric_name in metric_names:
        evaluation = evaluate_proportion_metric(metric_name, truth, predicted, declared_types)
        rows.append(
            {
                "run_id": metadata.get("run_id", run_dir.name),
                "dataset_id": metadata.get("dataset_id"),
                "modality": metadata.get("modality"),
                "feature_set": metadata.get("feature_set"),
                "method": metadata.get("method"),
                "method_run_id": metadata.get("method_run_id", metadata.get("method")),
                "status": "success",
                "metric": evaluation.metric_id,
                "metric_name": evaluation.metric_name,
                "metric_version": evaluation.metric_version,
                "evaluation_contract_version": evaluation.contract_version,
                "row_sum_atol": evaluation.row_sum_atol,
                "cell_type_universe": universe_json,
                "n_spots": evaluation.n_spots,
                "n_cell_types": evaluation.n_cell_types,
                "value": evaluation.value,
                "run_dir": str(run_dir),
                "error": None,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate standardized runs with strict, versioned proportion metrics."
    )
    parser.add_argument("--runs", nargs="+", required=True)
    parser.add_argument(
        "--metrics",
        nargs="+",
        default=list(DEFAULT_METRICS),
        choices=available_proportion_metrics(),
    )
    parser.add_argument("--truth")
    parser.add_argument(
        "--cell-types",
        nargs="+",
        help=(
            "Explicit ordered universe, applied to every run. Normally this is read from "
            "run.yaml or captured truth.cell_types in inputs.yaml."
        ),
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    canonical_ids = [resolve_proportion_metric(name).metric_id for name in args.metrics]
    if len(canonical_ids) != len(set(canonical_ids)):
        parser.error("--metrics contains aliases for the same canonical endpoint more than once.")

    rows = []
    for run in args.runs:
        rows.extend(
            evaluate_run_directory(
                Path(run),
                args.metrics,
                explicit_truth=args.truth,
                cell_types=args.cell_types,
            )
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)
    print(output)


if __name__ == "__main__":
    main()
