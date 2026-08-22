from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from scripts.run_deconvolution import run_experiment


def _write_toy_dataset(tmp_path: Path) -> Path:
    reference = ad.AnnData(
        X=np.array([[5, 0], [4, 1], [0, 5]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B"]}, index=["c1", "c2", "c3"]),
        var=pd.DataFrame({"highly_variable": [True, True]}, index=["f1", "f2"]),
    )
    spatial = ad.AnnData(
        X=np.array([[9, 1], [1, 8]], dtype=float),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["f1", "f2"]),
    )
    spatial.obsm["spatial"] = np.array([[0, 0], [1, 0]], dtype=float)

    truth = pd.DataFrame([[0.9, 0.1], [0.2, 0.8]], index=["s1", "s2"], columns=["A", "B"])

    reference_path = tmp_path / "reference.h5ad"
    spatial_path = tmp_path / "spatial.h5ad"
    truth_path = tmp_path / "truth.csv"
    reference.write_h5ad(reference_path)
    spatial.write_h5ad(spatial_path)
    truth.to_csv(truth_path)

    dataset_config = {
        "dataset_id": "toy",
        "labels_key": "cell_type",
        "spatial_key": "spatial",
        "modalities": {
            "atac": {
                "reference": {"path": str(reference_path)},
                "spatial": {"path": str(spatial_path)},
                "truth": {"path": str(truth_path)},
                "feature_sets": {
                    "highly_variable": {"var_column": "highly_variable"},
                    "all": {"mode": "all"},
                },
            }
        },
    }
    config_path = tmp_path / "toy.yaml"
    config_path.write_text(yaml.safe_dump(dataset_config))

    registry_path = tmp_path / "datasets.yaml"
    registry_path.write_text(yaml.safe_dump({"toy": {"config": str(config_path)}}))
    return registry_path


def test_run_experiment_writes_batch_comparison(tmp_path):
    registry_path = _write_toy_dataset(tmp_path)
    experiment_path = tmp_path / "experiment.yaml"
    experiment_path.write_text(
        yaml.safe_dump(
            {
                "run_group": "toy_batch",
                "datasets": ["toy"],
                "modalities": ["atac"],
                "feature_sets": ["highly_variable"],
                "methods": ["nnls"],
                "method_configs": {"nnls": {"method": "nnls", "params": {}}},
                "metrics": ["rmse", "jsd"],
                "skip_missing": True,
            }
        )
    )

    comparison_path = run_experiment(
        experiment_config_path=experiment_path,
        registry=registry_path,
        output_root_override=tmp_path / "results",
        overwrite=True,
    )

    assert comparison_path == tmp_path / "results" / "toy_batch" / "comparison.csv"
    comparison = pd.read_csv(comparison_path)
    assert set(comparison["metric"]) == {"rmse", "jsd"}
    assert set(comparison["status"]) == {"success"}
    assert (tmp_path / "results" / "toy_batch" / "runs.csv").exists()
    assert (
        tmp_path
        / "results"
        / "toy_batch"
        / "toy__atac__highly_variable__nnls"
        / "results"
        / "proportions.csv"
    ).exists()
