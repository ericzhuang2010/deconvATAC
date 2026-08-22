#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Callable, Optional

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.metrics import jsd, rmse


METRICS: dict[str, Callable[[pd.DataFrame, pd.DataFrame], float]] = {
    "rmse": rmse,
    "jsd": jsd,
}


def read_run_metadata(run_dir: Path) -> dict:
    path = run_dir / "run.yaml"
    if not path.exists():
        return {}
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def find_truth(run_dir: Path, explicit_truth: Optional[str]) -> Path:
    if explicit_truth is not None:
        return Path(explicit_truth)
    local_truth = run_dir / "results" / "truth.csv"
    if local_truth.exists():
        return local_truth
    raise FileNotFoundError(f"No truth file found for {run_dir}; pass --truth or include results/truth.csv.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate standardized deconvolution run directories.")
    parser.add_argument("--runs", nargs="+", required=True)
    parser.add_argument("--metrics", nargs="+", default=["rmse", "jsd"], choices=sorted(METRICS))
    parser.add_argument("--truth")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    for run in args.runs:
        run_dir = Path(run)
        metadata = read_run_metadata(run_dir)
        predicted = pd.read_csv(run_dir / "results" / "proportions.csv", index_col=0)
        truth = pd.read_csv(find_truth(run_dir, args.truth), index_col=0)

        for metric_name in args.metrics:
            rows.append(
                {
                    "run_id": metadata.get("run_id", run_dir.name),
                    "dataset_id": metadata.get("dataset_id"),
                    "modality": metadata.get("modality"),
                    "feature_set": metadata.get("feature_set"),
                    "method": metadata.get("method"),
                    "metric": metric_name,
                    "value": METRICS[metric_name](truth, predicted),
                }
            )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)
    print(output)


if __name__ == "__main__":
    main()
