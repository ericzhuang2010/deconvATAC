from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from scripts.materialize_shapemix_real_spatial import canonicalize_barcodes


ROOT = Path(__file__).resolve().parents[1]


def _yaml(path: Path):
    with path.open() as handle:
        return yaml.safe_load(handle)


def test_real_spatial_fragment_barcode_mapping_is_explicit_and_injective() -> None:
    suffix_policy = {
        "fragment_terminal_suffix_to_strip": "-1",
        "require_suffix_on_every_fragment_barcode": True,
    }
    assert canonicalize_barcodes(["pixel-a-1", "pixel-b-1"], suffix_policy) == [
        "pixel-a",
        "pixel-b",
    ]
    assert canonicalize_barcodes(["pixel-a", "pixel-b"], {"canonical_form": "identity"}) == [
        "pixel-a",
        "pixel-b",
    ]
    with pytest.raises(ValueError, match="required terminal suffix"):
        canonicalize_barcodes(["pixel-a-1", "pixel-b"], suffix_policy)
    with pytest.raises(ValueError, match="created duplicates"):
        canonicalize_barcodes(["pixel-a", "pixel-a"], {"canonical_form": "identity"})


@pytest.mark.parametrize(
    ("family", "expected_sections"),
    (("gse205055", 6), ("gse263333", 2)),
)
def test_real_spatial_templates_and_experiments_have_one_frozen_job_matrix(
    family: str,
    expected_sections: int,
) -> None:
    template = _yaml(ROOT / f"configs/datasets/shapemix_{family}_real_spatial_v1.yaml")
    experiment = _yaml(ROOT / f"configs/experiments/shapemix_{family}_real_spatial_v1.yaml")
    sections = template["sections"]
    dataset_ids = [section["dataset_id"] for section in sections]

    assert template["status"] == "frozen_before_predictions"
    assert template["truth_policy"] == "orthogonal_validation_only"
    assert len(sections) == expected_sections
    assert len(set(dataset_ids)) == expected_sections
    assert experiment["datasets"] == dataset_ids
    assert experiment["modalities"] == ["atac"]
    assert experiment["feature_sets"] == {"atac": ["all"]}
    assert [run["id"] for run in experiment["method_runs"]] == [
        "shapemix_length",
        "shapemix_count_only",
        "nnls",
    ]
    assert experiment["evaluation_mode"] == "prediction_only"
    assert experiment["metrics"] == []
    assert experiment["overwrite"] is False


def test_real_spatial_templates_keep_validation_modalities_out_of_shapemix_inputs() -> None:
    for family in ("gse205055", "gse263333"):
        template = _yaml(ROOT / f"configs/datasets/shapemix_{family}_real_spatial_v1.yaml")
        for section in template["sections"]:
            assert section["atac_gsm"].startswith("GSM")
            assert section["reference_id"].startswith("gse")
            assert "truth" not in section
            assert "validation_epigenome_gsms" in section
            assert "rna_gsm" in section
