from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from scripts.evaluate_runs import evaluate_run_directory


def _write_run(tmp_path: Path, *, include_universe: bool = True) -> Path:
    run_dir = tmp_path / "run"
    results_dir = run_dir / "results"
    results_dir.mkdir(parents=True)
    truth = pd.DataFrame(
        [[0.75, 0.25], [0.1, 0.9]],
        index=["s1", "s2"],
        columns=["A", "B"],
    )
    predicted = pd.DataFrame(
        [[0.7, 0.3], [0.2, 0.8]],
        index=["s1", "s2"],
        columns=["A", "B"],
    )
    truth.to_csv(results_dir / "truth.csv")
    predicted.to_csv(results_dir / "proportions.csv")
    (run_dir / "run.yaml").write_text(
        yaml.safe_dump(
            {
                "run_id": "toy_run",
                "dataset_id": "toy",
                "modality": "atac",
                "feature_set": "all",
                "method": "nnls",
                "method_run_id": "nnls_variant",
            }
        )
    )
    truth_spec = {"path": "truth.csv"}
    if include_universe:
        truth_spec["cell_types"] = ["A", "B"]
    (run_dir / "inputs.yaml").write_text(
        yaml.safe_dump(
            {
                "dataset_config": {
                    "modalities": {"atac": {"truth": truth_spec}}
                }
            }
        )
    )
    return run_dir


def test_standalone_evaluator_uses_shared_registry_and_records_contract(tmp_path):
    rows = evaluate_run_directory(_write_run(tmp_path), ["rmse", "jsd"])

    assert {row["metric"] for row in rows} == {"rmse_v1", "jsd_v2"}
    assert {row["metric_version"] for row in rows} == {"v1", "v2"}
    assert {row["evaluation_contract_version"] for row in rows} == {"v2"}
    assert all(json.loads(row["cell_type_universe"]) == ["A", "B"] for row in rows)
    assert all(row["method_run_id"] == "nnls_variant" for row in rows)
    assert all(row["n_spots"] == 2 and row["n_cell_types"] == 2 for row in rows)


def test_standalone_evaluator_requires_declaration_or_explicit_override(tmp_path):
    run_dir = _write_run(tmp_path, include_universe=False)
    with pytest.raises(ValueError, match="No declared cell-type universe"):
        evaluate_run_directory(run_dir, ["rmse_v1"])

    rows = evaluate_run_directory(run_dir, ["rmse_v1"], cell_types=["A", "B"])
    assert rows[0]["metric"] == "rmse_v1"


def test_standalone_evaluator_propagates_strict_contract_failures(tmp_path):
    run_dir = _write_run(tmp_path)
    predicted = pd.read_csv(run_dir / "results" / "proportions.csv", index_col=0)
    predicted = predicted.drop(index="s2")
    predicted.to_csv(run_dir / "results" / "proportions.csv")

    with pytest.raises(ValueError, match="exactly the truth spot set"):
        evaluate_run_directory(run_dir, ["jsd_v2"])
