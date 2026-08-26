from __future__ import annotations

import gzip
import hashlib
import shutil
from pathlib import Path

import h5py

import scripts.prepare_shapemix_gse246791 as preparation


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _minimal_h5ad(path: Path) -> None:
    with h5py.File(path, "w") as handle:
        handle.attrs["encoding-type"] = "anndata"
        handle.create_dataset("X", data=[[1, 0], [0, 1]])
        handle.create_group("obs")
        handle.create_group("var")


def test_selected_sources_match_frozen_major_region_scope() -> None:
    sources = preparation.selected_sources()
    assert len(sources) == 12
    assert [source.major_region for source in sources] == [
        "AMY",
        "CB",
        "HPF",
        "HY",
        "Isocortex",
        "MB",
        "MY",
        "OLF",
        "PAL",
        "Pons",
        "STR",
        "TH",
    ]
    assert len({source.gsm for source in sources}) == 12
    assert all(source.source_path.is_file() for source in sources)
    assert all(
        source.output_path.parent == preparation.SOURCE_OBJECT_ROOT
        for source in sources
    )


def test_repository_path_preserves_symlinked_data_path(
    tmp_path: Path,
    monkeypatch,
) -> None:
    repository = tmp_path / "repository"
    external = tmp_path / "external"
    repository.mkdir()
    external.mkdir()
    (repository / "data").symlink_to(external, target_is_directory=True)
    monkeypatch.setattr(preparation, "ROOT", repository)
    path = repository / "data" / "processed" / "object.h5ad"
    assert preparation._repository_path(path) == "data/processed/object.h5ad"


def test_decompress_source_is_atomic_and_reuses_valid_output(
    tmp_path: Path,
    monkeypatch,
) -> None:
    plain = tmp_path / "author.h5ad"
    source_path = tmp_path / "author.h5ad.gz"
    output_path = tmp_path / "processed" / "author.h5ad"
    _minimal_h5ad(plain)
    expected_output_sha256 = _sha256(plain)
    with plain.open("rb") as input_handle, gzip.open(source_path, "wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle)
    plain.unlink()
    monkeypatch.setattr(preparation, "ROOT", tmp_path)
    source = preparation.SelectedSource(
        gsm="GSM_TEST",
        sample="sample",
        major_region="region",
        name=source_path.name,
        source_path=source_path,
        expected_source_bytes=source_path.stat().st_size,
        source_sha256=_sha256(source_path),
        output_path=output_path,
    )

    first = preparation.decompress_source(source)
    assert output_path.is_file()
    assert not Path(f"{output_path}.partial").exists()
    assert first["output_sha256"] == expected_output_sha256
    assert first["validation"]["required_h5ad_root_objects"] == "passed"

    second = preparation.decompress_source(source)
    assert second["output_sha256"] == first["output_sha256"]
    assert second["output_bytes"] == first["output_bytes"]


def test_schema_snapshot_records_shapes_without_loading_values(
    tmp_path: Path,
    monkeypatch,
) -> None:
    output_path = tmp_path / "processed" / "author.h5ad"
    output_path.parent.mkdir(parents=True)
    _minimal_h5ad(output_path)
    monkeypatch.setattr(preparation, "ROOT", tmp_path)
    source = preparation.SelectedSource(
        gsm="GSM_TEST",
        sample="sample",
        major_region="region",
        name="author.h5ad.gz",
        source_path=tmp_path / "author.h5ad.gz",
        expected_source_bytes=0,
        source_sha256="unused",
        output_path=output_path,
    )

    snapshot = preparation.schema_snapshot(source)
    objects = {record["path"]: record for record in snapshot["objects"]}
    assert objects["X"]["shape"] == [2, 2]
    assert objects["X"]["kind"] == "dataset"
    assert objects["obs"]["kind"] == "group"
    assert snapshot["root_attributes"]["encoding-type"] == "anndata"
