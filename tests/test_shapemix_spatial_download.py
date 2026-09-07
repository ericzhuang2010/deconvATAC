import gzip
import threading
import urllib.error
from pathlib import Path

import pytest

import scripts.download_shapemix_spatial as spatial_download
from scripts.download_shapemix_spatial import (
    Resource,
    acquire,
    download_once,
    load_config,
    md5_file,
    normalize_etag,
    prioritize_incomplete_resources,
    remote_metadata,
    resources_from_config,
    safe_tar_member_name,
    select_resources,
    validate_gzip_schema,
    worker_limits,
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


def test_remote_metadata_retries_transient_http_error(monkeypatch: pytest.MonkeyPatch):
    calls = 0
    delays = []

    class Response:
        headers = {"Content-Length": "123", "ETag": '"abc"'}

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

    def urlopen(_request, timeout):
        nonlocal calls
        assert timeout == 7
        calls += 1
        if calls == 1:
            raise urllib.error.HTTPError(
                "https://example.invalid/payload", 403, "Forbidden", {}, None
            )
        return Response()

    monkeypatch.setattr(spatial_download.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(spatial_download.time, "sleep", delays.append)

    assert remote_metadata(
        "https://example.invalid/payload",
        timeout=7,
        attempts=2,
        retry_delay_seconds=0.25,
    ) == (123, "abc")
    assert calls == 2
    assert delays == [0.25]


def test_remote_metadata_does_not_retry_permanent_http_error(
    monkeypatch: pytest.MonkeyPatch,
):
    calls = 0

    def urlopen(_request, timeout):
        nonlocal calls
        assert timeout == 7
        calls += 1
        raise urllib.error.HTTPError(
            "https://example.invalid/missing", 404, "Not Found", {}, None
        )

    monkeypatch.setattr(spatial_download.urllib.request, "urlopen", urlopen)

    with pytest.raises(urllib.error.HTTPError):
        remote_metadata(
            "https://example.invalid/missing",
            timeout=7,
            attempts=5,
            retry_delay_seconds=0,
        )
    assert calls == 1


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


def test_resource_selector_retains_only_requested_representative_pair():
    resources = resources_from_config(
        load_config(ROOT / "configs/data_sources/shapemix_gse246791_fragment_reads.yaml")
    )
    selected = select_resources(resources, ["GSM7877011"])
    assert [resource.name for resource in selected] == [
        "SRR26585986_1.fastq.gz",
        "SRR26585986_2.fastq.gz",
    ]


def test_resource_selector_rejects_unknown_or_repeated_accessions():
    resources = resources_from_config(
        load_config(ROOT / "configs/data_sources/shapemix_gse246791_fragment_reads.yaml")
    )
    with pytest.raises(ValueError, match="absent"):
        select_resources(resources, ["GSM_missing"])
    with pytest.raises(ValueError, match="unique"):
        select_resources(resources, ["GSM7877011", "GSM7877011"])


def test_adult_fragment_download_separates_transfer_and_cpu_worker_limits():
    config = load_config(
        ROOT / "configs/data_sources/shapemix_gse246791_fragment_reads.yaml"
    )
    assert worker_limits(config, 8) == (8, 4)
    with pytest.raises(ValueError, match="must be <= 8"):
        worker_limits(config, 9)


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


def test_download_once_reconnects_and_resumes_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    destination = tmp_path / "raw" / "payload.gz"
    staging = tmp_path / "work" / "payload.gz.part"
    staging.parent.mkdir(parents=True)
    staging.write_bytes(b"a")
    resource = Resource(
        name="payload.gz",
        role="test",
        accession="GSM_test",
        url="https://example.invalid/payload.gz",
        destination=destination,
        staging_path=staging,
    )
    calls = []

    def run(command, check):
        assert check
        calls.append(command)
        with staging.open("ab") as handle:
            handle.write(b"b" if len(calls) == 1 else b"c")
        if len(calls) == 1:
            raise spatial_download.subprocess.CalledProcessError(1, command)

    delays = []
    monkeypatch.setattr(spatial_download.subprocess, "run", run)
    monkeypatch.setattr(spatial_download.time, "sleep", delays.append)

    download_once(
        resource,
        expected_bytes=3,
        timeout=7,
        attempts=2,
        retry_delay_seconds=0.25,
    )

    assert staging.read_bytes() == b"abc"
    assert len(calls) == 2
    assert "--continue" in calls[0]
    assert delays == [0.25]


def test_acquire_removes_checksum_invalid_staging_file_for_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    destination = tmp_path / "raw" / "payload.gz"
    staging = tmp_path / "work" / "payload.gz.part"
    staging.parent.mkdir(parents=True)
    staging.write_bytes(b"bad")
    resource = Resource(
        name="payload.gz",
        role="test",
        accession="GSM_test",
        url="https://example.invalid/payload.gz",
        destination=destination,
        staging_path=staging,
        expected_bytes=3,
        expected_md5="900150983cd24fb0d6963f7d28e17f72",
    )
    monkeypatch.setattr(
        spatial_download, "remote_metadata", lambda _url, _timeout: (3, None)
    )

    with pytest.raises(ValueError, match="removed invalid staging file"):
        acquire(resource, timeout=1, validation_slots=threading.BoundedSemaphore(1))

    assert not staging.exists()
    assert not destination.exists()


def test_incomplete_resources_are_scheduled_first_in_stable_order(tmp_path: Path):
    def make_resource(name: str) -> Resource:
        return Resource(
            name=name,
            role="test",
            accession="GSM_test",
            url=f"https://example.invalid/{name}",
            destination=tmp_path / "raw" / name,
            staging_path=tmp_path / "work" / f"{name}.part",
        )

    complete = make_resource("complete.gz")
    complete.destination.parent.mkdir(parents=True)
    complete.destination.write_bytes(b"complete")
    incomplete_a = make_resource("incomplete-a.gz")
    incomplete_b = make_resource("incomplete-b.gz")

    scheduled = prioritize_incomplete_resources(
        (complete, incomplete_a, incomplete_b)
    )

    assert [resource.name for resource in scheduled] == [
        "incomplete-a.gz",
        "incomplete-b.gz",
        "complete.gz",
    ]
