import gzip

import pysam
from pathlib import Path

from scripts.download_gse129785 import load_config, resources_from_config
from scripts.preprocess_gse129785 import (
    DEFAULT_CONFIG,
    REFERENCE_TYPES,
    ROOT,
    is_bgzf,
    merged_query_regions,
    parse_author_peaks,
    physical_dilution_descriptor,
)
from deconvatac.pp import PeakInterval


def test_frozen_gse129785_scope_resolves_unique_resources():
    config = load_config(DEFAULT_CONFIG)
    resources = resources_from_config(config, "all")

    assert len(resources) == 43
    assert sum(resource.gsm is not None for resource in resources) == 30
    assert len({resource.name for resource in resources}) == len(resources)
    assert len([sample for sample in config["samples"] if sample["role"] == "physical_dilution"]) == 14
    assert tuple(
        sample["cell_type"]
        for sample in config["samples"]
        if sample["role"] == "sorted_reference"
    ) == REFERENCE_TYPES


def test_bgzf_detection_distinguishes_plain_gzip(tmp_path: Path):
    plain = tmp_path / "plain.tsv.gz"
    blocked = tmp_path / "blocked.tsv.gz"
    with gzip.open(plain, "wt") as handle:
        handle.write("chr1\t0\t10\tcell\t1\n")
    with pysam.BGZFile(str(blocked), "wb") as handle:
        handle.write(b"chr1\t0\t10\tcell\t1\n")

    assert not is_bgzf(plain)
    assert is_bgzf(blocked)


def test_author_peak_parser_preserves_source_order(tmp_path: Path):
    path = tmp_path / "peaks.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("Feature\nchr2_100_600\nchr1_25_525\n")

    peaks = parse_author_peaks(path)

    assert peaks == (
        PeakInterval("chr2", 99, 600, "chr2:99-600"),
        PeakInterval("chr1", 24, 525, "chr1:24-525"),
    )


def test_peak_queries_include_end_at_peak_start_without_duplicate_windows():
    peaks = (
        PeakInterval("chr1", 100, 200, "p1"),
        PeakInterval("chr1", 300, 400, "p2"),
        PeakInterval("chr1", 501, 600, "p3"),
        PeakInterval("chr2", 0, 50, "p4"),
    )

    regions = merged_query_regions(peaks, max_fragment_length=100)

    assert regions == (
        ("chr1", 99, 400),
        ("chr1", 500, 600),
        ("chr2", 0, 50),
    )


def test_physical_dilution_descriptor_keeps_nominal_evidence_out_of_truth():
    dataset_id = "gse129785_shapemix_physical_dilution_test"
    dataset_root = ROOT / "data/processed/datasets" / dataset_id
    descriptor = physical_dilution_descriptor(
        {
            "gsm": "GSM_TEST",
            "family": "cd4_memory_cd8_naive",
            "components": ["CD4 Memory", "CD8 Naive"],
            "fractions": [0.001, 0.999],
        },
        dataset_id,
        ROOT / "data/processed/references/gse129785_immune/atac/reference.h5ad",
        dataset_root / "atac/spatial.h5ad",
        dataset_root / "atac/features/selected_reference_peaks.txt",
        dataset_root / "validation/nominal_broad_proportions.csv",
    )

    assert "truth" not in descriptor["modalities"]["atac"]
    nominal = descriptor["validation"]["nominal_broad_proportions"]
    assert nominal == {
        "path": (
            "data/processed/datasets/"
            "gse129785_shapemix_physical_dilution_test/"
            "validation/nominal_broad_proportions.csv"
        ),
        "evidence_class": "nominal_sample_level",
        "exact_truth": False,
    }
    assert descriptor["benchmark_scope"] == "external_validation_nominal"
    assert (
        descriptor["physical_dilution"]["evidence_limitation"]
        == "Nominal sample-level input proportions; sorting and cell-recovery "
        "uncertainty remain."
    )
