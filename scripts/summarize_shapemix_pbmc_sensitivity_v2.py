#!/usr/bin/env python
"""Summarize the support-safe PBMC sensitivity v2 campaign."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from scripts import summarize_shapemix_pbmc_sensitivity as summary


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_ID = "shapemix_pbmc_stress_v2"
summary.DEFAULT_EXPERIMENT = ROOT / "configs/experiments/shapemix_pbmc_stress_v2.yaml"
summary.DEFAULT_BATCH = (
    ROOT
    / "results/sensitivity/shapemix_pbmc_stress_v2"
    / "shapemix_pbmc_stress_protocol_v2_cuda"
)


def dataset_design(dataset_id: str) -> dict[str, Any]:
    descriptor = summary.read_yaml(
        ROOT / "data/processed/datasets" / dataset_id / "dataset.yaml"
    )
    sensitivity = descriptor.get("sensitivity")
    if not isinstance(sensitivity, dict):
        raise ValueError(f"Dataset lacks frozen sensitivity metadata: {dataset_id}")
    result = dict(sensitivity)
    if result.get("campaign") != CAMPAIGN_ID:
        raise ValueError(f"Unexpected sensitivity campaign: {dataset_id}")
    return result


def main() -> None:
    summary.dataset_design = dataset_design
    summary.main()


if __name__ == "__main__":
    main()
