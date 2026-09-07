from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from deconvatac.pp.fragment_shapes import (
    count_fragment_shapes,
    count_fragment_shapes_from_records,
)
import scripts.materialize_shapemix_real_spatial as materializer
from scripts.materialize_shapemix_real_spatial import (
    canonicalize_barcodes,
    merge_fragment_shape_results,
    spatial_count_workers,
)


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
    assert template["materialization"]["fragment_count_partition"] == "tabix_contig"
    assert (
        template["materialization"]["fragment_count_merge"]
        == "exact_canonical_integer_csr"
    )
    assert template["materialization"]["fragment_count_workers_max"] == {
        "co_tenant": 2,
        "exclusive": 8,
    }


def test_real_spatial_templates_keep_validation_modalities_out_of_shapemix_inputs() -> None:
    for family in ("gse205055", "gse263333"):
        template = _yaml(ROOT / f"configs/datasets/shapemix_{family}_real_spatial_v1.yaml")
        for section in template["sections"]:
            assert section["atac_gsm"].startswith("GSM")
            assert section["reference_id"].startswith("gse")
            assert "truth" not in section
            assert "validation_epigenome_gsms" in section
            assert "rna_gsm" in section


def test_real_spatial_parallel_worker_limit_tracks_resource_profile(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("DECONVATAC_RESOURCE_PROFILE", raising=False)
    assert spatial_count_workers(20) == 2
    assert spatial_count_workers(1) == 1
    monkeypatch.setenv("DECONVATAC_RESOURCE_PROFILE", "exclusive")
    assert spatial_count_workers(20) == 8
    assert spatial_count_workers(4) == 4
    monkeypatch.setenv("DECONVATAC_RESOURCE_PROFILE", "invalid")
    with pytest.raises(ValueError, match="Unsupported"):
        spatial_count_workers(20)


def test_disjoint_contig_fragment_shards_merge_exactly() -> None:
    barcodes = ["cellA", "cellB"]
    peaks = [
        ("chr1", 0, 100, "chr1:0-100"),
        ("chr1", 100, 200, "chr1:100-200"),
        ("chr2", 0, 100, "chr2:0-100"),
    ]
    chr1 = [
        "chr1\t0\t50\tcellA\t1",
        "chr1\t50\t150\tcellB\t2",
        "chr1\t100\t180\tunknown\t1",
    ]
    chr2 = [
        "chr2\t0\t90\tcellA\t3",
        "chr2\t10\t300\tcellB\t1",
    ]
    shards = [
        count_fragment_shapes_from_records(
            records,
            barcodes,
            peaks,
            right_cut_offset=0,
            chunk_size=2,
        )
        for records in (chr1, chr2)
    ]
    merged = merge_fragment_shape_results(reversed(shards))
    serial = count_fragment_shapes_from_records(
        chr1 + chr2,
        barcodes,
        peaks,
        right_cut_offset=0,
        chunk_size=2,
    )

    assert merged.barcodes == serial.barcodes
    assert merged.peaks == serial.peaks
    assert merged.bins == serial.bins
    assert merged.qc.to_dict() == serial.qc.to_dict()
    for layer in (bin_.layer for bin_ in serial.bins):
        np.testing.assert_array_equal(
            merged.layers[layer].toarray(),
            serial.layers[layer].toarray(),
        )
    np.testing.assert_array_equal(merged.X.toarray(), serial.X.toarray())


def test_contig_parallel_counter_matches_serial_bgzf(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pysam = pytest.importorskip("pysam")
    source = tmp_path / "fragments.tsv"
    source.write_text(
        "chr1\t0\t50\tcellA\t1\n"
        "chr1\t50\t150\tcellB\t2\n"
        "chr2\t0\t90\tcellA\t3\n"
        "chr2\t10\t300\tcellB\t1\n"
    )
    compressed = tmp_path / "fragments.tsv.gz"
    pysam.tabix_compress(str(source), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="bed", force=True)

    barcodes = ["cellA", "cellB"]
    peaks = [
        ("chr1", 0, 100, "chr1:0-100"),
        ("chr1", 100, 200, "chr1:100-200"),
        ("chr2", 0, 100, "chr2:0-100"),
    ]
    bins = [
        {
            "name": value.name,
            "min_inclusive": value.min_inclusive,
            "max_exclusive": value.max_exclusive,
            "layer": value.layer,
        }
        for value in count_fragment_shapes_from_records(
            [],
            barcodes,
            peaks,
            right_cut_offset=0,
        ).bins
    ]
    monkeypatch.setenv("DECONVATAC_RESOURCE_PROFILE", "exclusive")
    monkeypatch.setattr(materializer, "repository_path", lambda value: str(value))

    parallel = materializer.count_spatial_fragment_shapes(
        compressed,
        barcodes,
        peaks,
        right_cut_offset=0,
        bins=bins,
        chunk_size=2,
    )
    serial = count_fragment_shapes(
        compressed,
        barcodes,
        peaks,
        right_cut_offset=0,
        bins=bins,
        chunk_size=2,
    )

    assert parallel.qc.to_dict() == serial.qc.to_dict()
    for layer in (bin_.layer for bin_ in serial.bins):
        np.testing.assert_array_equal(
            parallel.layers[layer].toarray(),
            serial.layers[layer].toarray(),
        )
