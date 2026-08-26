from pathlib import Path

import pytest
import yaml

from scripts.validate_shapemix_file_layout import (
    LayoutError,
    ROOT,
    canonical_result_root,
    project_path,
    validate_experiment,
    validate_source_manifest,
)


def test_project_paths_and_result_scopes_fail_closed():
    assert project_path("data/processed/value", "test") == Path("data/processed/value")
    with pytest.raises(LayoutError, match="project-relative"):
        project_path("/absolute/value", "test")
    with pytest.raises(LayoutError, match="escape"):
        project_path("../outside", "test")
    assert canonical_result_root(
        "results/external_validation/campaign_v1",
        "output",
    ) == Path("results/external_validation/campaign_v1")
    with pytest.raises(LayoutError, match="canonical_scope"):
        canonical_result_root("results/external/campaign_v1", "output")


def test_frozen_spatial_reference_manifests_follow_canonical_lifecycle():
    for filename in (
        "shapemix_gse216371_reference.yaml",
        "shapemix_gse246791_reference.yaml",
        "shapemix_gse244618_reference.yaml",
    ):
        validate_source_manifest(ROOT / "configs/data_sources" / filename)


def test_gse129785_campaign_passes_layout_and_truth_contract():
    batch = validate_experiment(
        ROOT / "configs/experiments/shapemix_gse129785_external_v2.yaml",
        ROOT / "data/registry/datasets.yaml",
        allow_existing_results=True,
    )
    assert batch.relative_to(ROOT) == Path(
        "results/external_validation/shapemix_gse129785_v2/"
        "shapemix_gse129785_external_protocol_v2_cpu_log_abundance"
    )


def test_noncanonical_experiment_output_is_rejected(tmp_path: Path):
    source = ROOT / "configs/experiments/shapemix_gse129785_external.yaml"
    config = yaml.safe_load(source.read_text())
    config["output_root"] = "results/external/gse129785"
    destination = tmp_path / "experiment.yaml"
    destination.write_text(yaml.safe_dump(config, sort_keys=False))

    with pytest.raises(LayoutError, match="canonical_scope"):
        validate_experiment(destination, ROOT / "data/registry/datasets.yaml")


def test_prediction_only_campaign_rejects_truth_descriptor(tmp_path: Path):
    source = ROOT / "configs/experiments/shapemix_cuda_smoke_qualification.yaml"
    config = yaml.safe_load(source.read_text())
    config["evaluation_mode"] = "prediction_only"
    config["metrics"] = []
    destination = tmp_path / "experiment.yaml"
    destination.write_text(yaml.safe_dump(config, sort_keys=False))

    with pytest.raises(LayoutError, match="declares truth"):
        validate_experiment(
            destination,
            ROOT / "data/registry/datasets.yaml",
            allow_existing_results=True,
        )
