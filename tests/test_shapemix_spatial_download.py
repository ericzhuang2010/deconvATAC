import gzip
from pathlib import Path

import pytest

from scripts.download_shapemix_spatial import (
    load_config,
    md5_file,
    normalize_etag,
    resources_from_config,
    safe_tar_member_name,
    validate_gzip_schema,
)


ROOT = Path(__file__).resolve().parents[1]


def test_gse205055_scope_includes_complete_super_series_metadata_and_archive():
    path = ROOT / "configs/data_sources/shapemix_gse205055.yaml"
    resources = resources_from_config(load_config(path))

    assert len(resources) == 9
    assert resources[-1].name == "GSE205055_RAW.tar"
    assert resources[-1].expected_bytes == 8177080320
    assert {resource.accession for resource in resources[:-1]} == {
        "GSE205055",
        "GSE205051",
        "GSE205052",
        "GSE205054",
        "GSE205180",
        "GSE205181",
        "GSE217091",
        "GSE218593",
    }


def test_gse263333_scope_includes_complete_series_archive():
    path = ROOT / "configs/data_sources/shapemix_gse263333.yaml"
    resources = resources_from_config(load_config(path))

    assert len(resources) == 2
    assert resources[-1].name == "GSE263333_RAW.tar"
    assert resources[-1].expected_bytes == 3306485760


def test_tar_member_safety_rejects_escapes_and_absolute_paths():
    assert safe_tar_member_name("GSM1_file.tsv.gz")
    assert safe_tar_member_name("nested/GSM1_file.tsv.gz")
    assert not safe_tar_member_name("../outside")
    assert not safe_tar_member_name("nested/../../outside")
    assert not safe_tar_member_name("/absolute")


def test_etag_normalization_accepts_http_quoting_and_weak_prefix():
    assert normalize_etag('"abc-123"') == "abc-123"
    assert normalize_etag('W/"abc-123"') == "abc-123"
    assert normalize_etag(None) is None


def test_reference_manifests_resolve_exact_frozen_scopes():
    expected = {
        "shapemix_gse216371_reference.yaml": (3, 76261606450),
        "shapemix_gse246791_reference.yaml": (17, 8018764799),
        "shapemix_gse244618_reference.yaml": (16, 7559398338),
        "shapemix_gse246791_fragment_reads.yaml": (24, 149964894951),
    }
    for filename, (count, total_bytes) in expected.items():
        resources = resources_from_config(load_config(ROOT / "configs/data_sources" / filename))
        assert len(resources) == count
        assert sum(resource.expected_bytes or 0 for resource in resources) == total_bytes
        assert all(str(resource.destination).startswith(str(ROOT / "data/raw/sources")) for resource in resources)
        assert all(str(resource.staging_path).startswith(str(ROOT / "data/work/downloads")) for resource in resources)


def test_ucsc_mm10_reference_manifest_is_pinned_and_organized():
    path = ROOT / "configs/data_sources/shapemix_ucsc_mm10_initial.yaml"
    resources = resources_from_config(load_config(path))
    assert len(resources) == 2
    assert sum(resource.expected_bytes or 0 for resource in resources) == 870142765
    assert resources[0].expected_md5 == "db005b65828db31735f384e4c5787be5"
    assert resources[1].expected_md5 == "5a103c9a15dd660c295a089ef5035672"
    assert all(
        str(resource.destination).startswith(
            str(ROOT / "data/raw/sources/ucsc/mm10_initial")
        )
        for resource in resources
    )
    assert all(
        str(resource.staging_path).startswith(
            str(ROOT / "data/work/downloads/ucsc_mm10_initial")
        )
        for resource in resources
    )


def test_gzip_schema_accepts_h5ad_magic_and_bedpe_rows(tmp_path: Path):
    h5ad = tmp_path / "sample.h5ad.gz"
    with gzip.open(h5ad, "wb") as handle:
        handle.write(b"\x89HDF\r\n\x1a\nrest")
    validate_gzip_schema(h5ad, h5ad.name)

    bedpe = tmp_path / "sample.bedpe.gz"
    with gzip.open(bedpe, "wt") as handle:
        handle.write("chr1\t1\t2\tchr1\t10\t11\tbarcode\n")
    validate_gzip_schema(bedpe, bedpe.name)


def test_gzip_schema_rejects_non_hdf5_h5ad(tmp_path: Path):
    path = tmp_path / "bad.h5ad.gz"
    with gzip.open(path, "wb") as handle:
        handle.write(b"not hdf5")
    with pytest.raises(ValueError, match="HDF5"):
        validate_gzip_schema(path, path.name)


def test_md5_file_matches_provider_style_digest(tmp_path: Path):
    path = tmp_path / "payload"
    path.write_bytes(b"abc")
    assert md5_file(path) == "900150983cd24fb0d6963f7d28e17f72"
