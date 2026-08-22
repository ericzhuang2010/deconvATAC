from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from deconvatac.data import load_deconvolution_input


def _write_test_dataset(tmp_path: Path) -> Path:
    reference = ad.AnnData(
        X=np.array([[5, 0, 1], [4, 1, 0], [0, 5, 1]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B"]}, index=["c1", "c2", "c3"]),
        var=pd.DataFrame({"highly_variable": [True, True, False]}, index=["f1", "f2", "f3"]),
    )
    spatial = ad.AnnData(
        X=np.array([[9, 1, 1], [1, 8, 1]], dtype=float),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["f1", "f2", "f3"]),
    )
    spatial.obsm["spatial"] = np.array([[0, 0], [1, 0]], dtype=float)
    spatial.obsm["proportions"] = pd.DataFrame(
        [[0.9, 0.1], [0.2, 0.8]],
        index=spatial.obs_names,
        columns=["A", "B"],
    )

    reference_path = tmp_path / "reference.h5ad"
    spatial_path = tmp_path / "spatial.h5ad"
    reference.write_h5ad(reference_path)
    spatial.write_h5ad(spatial_path)

    dataset_config = {
        "dataset_id": "toy",
        "labels_key": "cell_type",
        "spatial_key": "spatial",
        "modalities": {
            "atac": {
                "reference": {"path": str(reference_path)},
                "spatial": {"path": str(spatial_path)},
                "truth": {"obsm_key": "proportions"},
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


def test_load_deconvolution_input_aligns_features_and_truth(tmp_path):
    registry_path = _write_test_dataset(tmp_path)

    data = load_deconvolution_input(
        dataset_id="toy",
        modality="atac",
        feature_set="highly_variable",
        registry_path=registry_path,
    )

    assert data.dataset_id == "toy"
    assert data.modality == "atac"
    assert list(data.reference.var_names) == ["f1", "f2"]
    assert list(data.spatial.var_names) == ["f1", "f2"]
    assert data.truth is not None
    assert list(data.truth.columns) == ["A", "B"]


def test_load_deconvolution_input_uses_feature_list_file(tmp_path):
    registry_path = _write_test_dataset(tmp_path)
    feature_list = tmp_path / "features.txt"
    feature_list.write_text("f3\nf1\n")

    config_path = tmp_path / "toy.yaml"
    config = yaml.safe_load(config_path.read_text())
    config["modalities"]["atac"]["feature_sets"]["from_file"] = {"path": str(feature_list)}
    config_path.write_text(yaml.safe_dump(config))

    data = load_deconvolution_input(
        dataset_id="toy",
        modality="atac",
        feature_set="from_file",
        registry_path=registry_path,
    )

    assert list(data.reference.var_names) == ["f3", "f1"]
    assert list(data.spatial.var_names) == ["f3", "f1"]


def test_load_deconvolution_input_casts_truth_file_index_to_strings(tmp_path):
    registry_path = _write_test_dataset(tmp_path)
    truth_path = tmp_path / "truth.csv"
    pd.DataFrame([[0.9, 0.1], [0.2, 0.8]], index=[0, 1], columns=["A", "B"]).to_csv(truth_path)

    spatial_path = tmp_path / "spatial_numeric_obs.h5ad"
    spatial = ad.read_h5ad(tmp_path / "spatial.h5ad")
    spatial.obs_names = ["0", "1"]
    spatial.write_h5ad(spatial_path)

    config_path = tmp_path / "toy.yaml"
    config = yaml.safe_load(config_path.read_text())
    config["modalities"]["atac"]["spatial"] = {"path": str(spatial_path)}
    config["modalities"]["atac"]["truth"] = {"path": str(truth_path)}
    config_path.write_text(yaml.safe_dump(config))

    data = load_deconvolution_input(
        dataset_id="toy",
        modality="atac",
        feature_set="highly_variable",
        registry_path=registry_path,
    )

    assert data.truth is not None
    assert list(data.truth.index) == ["0", "1"]
