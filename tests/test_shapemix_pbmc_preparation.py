import hashlib
from dataclasses import dataclass

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse

from deconvatac.data import ordered_feature_sha256
from deconvatac.pp import FRAGMENT_SHAPE_LAYER_NAMES, PeakInterval
from scripts.prepare_shapemix_pbmc import (
    FROZEN_CELL_TYPE_COUNTS,
    FROZEN_CELL_TYPES,
    atomic_output_directory,
    build_cell_universe,
    build_shape_cache_from_records,
    canonical_split_csv,
    derive_shape_pool,
    matrices_equal,
    official_order_peak_union,
    prepare_split_from_objects,
    select_reference_only_peaks,
    sha256_file,
    split_membership_sha256,
    stratified_reference_heldout_split,
    universe_sha256,
    validate_pinned_files,
    validate_shape_cache,
    validate_shape_cache_against_official,
    validate_split_membership,
    write_split_outputs,
)


SOURCE_HASHES = {"fragments": "a" * 64, "tabix_index": "b" * 64}


@dataclass
class FakeSelection:
    peak_ids: tuple[str, ...]
    original_peak_ids: tuple[str, ...]

    @property
    def indices(self):
        lookup = {name: index for index, name in enumerate(self.original_peak_ids)}
        return np.asarray([lookup[name] for name in self.peak_ids], dtype=np.int64)

    @property
    def candidate_feature_sha256(self):
        return ordered_feature_sha256(self.original_peak_ids)

    @property
    def selected_feature_sha256(self):
        return ordered_feature_sha256(self.peak_ids)

    def to_frame(self):
        return pd.DataFrame(
            {
                "rank": np.arange(1, len(self.peak_ids) + 1),
                "peak_id": self.peak_ids,
                "original_index": self.indices,
                "score": np.linspace(2.0, 1.0, len(self.peak_ids)),
                "nonzero_reference_cells": 2,
                "total_reference_count": 4,
            }
        )


def toy_universe(cell_types=("A", "B"), cells_per_type=4):
    return pd.DataFrame(
        [
            {"barcode": f"{cell_type}_{index:02d}", "cell_type": cell_type}
            for cell_type in cell_types
            for index in range(cells_per_type)
        ]
    )


def toy_cache_and_reference():
    cell_types = ("A", "B")
    universe = toy_universe(cell_types, cells_per_type=4)
    peaks = (
        PeakInterval("chr1", 0, 100, "chr1:0-100"),
        PeakInterval("chr1", 100, 200, "chr1:100-200"),
        PeakInterval("chr1", 200, 300, "chr1:200-300"),
    )
    records = []
    for barcode in universe["barcode"]:
        records.extend(
            [
                f"chr1\t10\t109\t{barcode}\t1",  # short: first and second peaks
                f"chr1\t20\t120\t{barcode}\t2",  # mono: first and second peaks
                f"chr1\t30\t280\t{barcode}\t3",  # long: first and third peaks
            ]
        )
    cache = build_shape_cache_from_records(
        records,
        universe,
        peaks,
        source_sha256=SOURCE_HASHES,
        cell_types=cell_types,
        chunk_size=3,
    )
    reference = ad.AnnData(
        X=cache.X.copy(),
        obs=cache.obs.copy(),
        var=cache.var.copy(),
    )
    return universe, cell_types, peaks, cache, reference


def test_pinned_file_validation_is_complete_and_fails_closed(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.write_bytes(b"alpha")
    second.write_bytes(b"beta")
    paths = {"first": first, "second": second}
    hashes = {
        "first": hashlib.sha256(b"alpha").hexdigest(),
        "second": hashlib.sha256(b"beta").hexdigest(),
    }

    assert validate_pinned_files(paths, hashes) == hashes
    with pytest.raises(ValueError, match="missing=.*second"):
        validate_pinned_files({"first": first}, hashes)
    with pytest.raises(ValueError, match="expected"):
        validate_pinned_files(paths, {**hashes, "second": "0" * 64})
    second.unlink()
    with pytest.raises(FileNotFoundError, match="second"):
        validate_pinned_files(paths, hashes)


def test_exact_frozen_16_type_universe_and_order_are_enforced():
    mapping_rows = []
    for cell_type, count in FROZEN_CELL_TYPE_COUNTS.items():
        mapping_rows.extend(
            {
                "barcode": f"{FROZEN_CELL_TYPES.index(cell_type):02d}_{index:04d}",
                "cell_type": cell_type,
            }
            for index in range(count)
        )
    mapping = pd.DataFrame(mapping_rows).sample(frac=1.0, random_state=3)
    reference_barcodes = mapping["barcode"].tolist() + ["excluded_low_count"]
    mapping = pd.concat(
        [
            mapping,
            pd.DataFrame([{"barcode": "excluded_low_count", "cell_type": "pDC"}]),
        ],
        ignore_index=True,
    )
    label_lookup = mapping.set_index("barcode")["cell_type"]

    universe = build_cell_universe(
        reference_barcodes,
        mapping,
        cell_types=FROZEN_CELL_TYPES,
        expected_counts=FROZEN_CELL_TYPE_COUNTS,
        reference_labels=label_lookup.loc[reference_barcodes].tolist(),
    )

    assert len(universe) == 9_500
    assert tuple(universe["cell_type"].drop_duplicates()) == FROZEN_CELL_TYPES
    assert universe.groupby("cell_type", sort=False).size().to_dict() == FROZEN_CELL_TYPE_COUNTS
    assert universe_sha256(universe, FROZEN_CELL_TYPES) == universe_sha256(
        universe.sample(frac=1.0, random_state=11), FROZEN_CELL_TYPES
    )

    reduced = mapping[mapping["barcode"] != "00_0000"]
    with pytest.raises(ValueError, match="CD14 Mono.*2550.*2551"):
        build_cell_universe(
            reference_barcodes,
            reduced,
            cell_types=FROZEN_CELL_TYPES,
            expected_counts=FROZEN_CELL_TYPE_COUNTS,
        )


def test_stratified_split_is_seeded_disjoint_canonical_and_input_order_invariant():
    cell_types = ("A", "B", "C")
    universe = toy_universe(cell_types, cells_per_type=10)
    first = stratified_reference_heldout_split(
        universe, 1103, cell_types=cell_types
    )
    reordered = stratified_reference_heldout_split(
        universe.sample(frac=1.0, random_state=22), 1103, cell_types=cell_types
    )
    changed_seed = stratified_reference_heldout_split(
        universe, 2203, cell_types=cell_types
    )

    pd.testing.assert_frame_equal(first, reordered)
    assert not first.equals(changed_seed)
    assert first.groupby(["cell_type", "pool"]).size().to_dict() == {
        (cell_type, pool): count
        for cell_type in cell_types
        for pool, count in (("heldout", 3), ("reference", 7))
    }
    reference = set(first.loc[first["pool"] == "reference", "barcode"])
    heldout = set(first.loc[first["pool"] == "heldout", "barcode"])
    assert reference.isdisjoint(heldout)
    assert reference | heldout == set(universe["barcode"])

    expected_hash = hashlib.sha256(canonical_split_csv(first, cell_types)).hexdigest()
    assert split_membership_sha256(first, cell_types) == expected_hash
    assert split_membership_sha256(first.sample(frac=1), cell_types) == expected_hash
    assert expected_hash == "1b3a760bb9fadfe2a15b7a9a57d0d6de63827925ef1a7de33820c1451532a525"


def test_frozen_split_has_exact_protocol_pool_sizes():
    universe = pd.DataFrame(
        [
            {"barcode": f"{order:02d}_{index:04d}", "cell_type": cell_type}
            for order, (cell_type, count) in enumerate(FROZEN_CELL_TYPE_COUNTS.items())
            for index in range(count)
        ]
    )
    membership = stratified_reference_heldout_split(
        universe, 0, cell_types=FROZEN_CELL_TYPES
    )
    assert (membership["pool"] == "reference").sum() == 6_644
    assert (membership["pool"] == "heldout").sum() == 2_856
    assert all(
        set(group["pool"]) == {"reference", "heldout"}
        for _, group in membership.groupby("cell_type", sort=False)
    )


def test_split_validator_rejects_overlap_label_changes_and_wrong_sizes():
    cell_types = ("A", "B")
    universe = toy_universe(cell_types, cells_per_type=4)
    membership = stratified_reference_heldout_split(
        universe, 0, cell_types=cell_types
    )

    duplicated = pd.concat([membership, membership.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="cannot occur in both"):
        validate_split_membership(duplicated, universe, cell_types)
    changed_label = membership.copy()
    changed_label.loc[0, "cell_type"] = "B"
    with pytest.raises(ValueError, match="changes one or more frozen"):
        validate_split_membership(changed_label, universe, cell_types)
    wrong_size = membership.copy()
    row = wrong_size.index[
        (wrong_size["cell_type"] == "A") & (wrong_size["pool"] == "heldout")
    ][0]
    wrong_size.loc[row, "pool"] = "reference"
    with pytest.raises(ValueError, match="Reference pool for 'A'"):
        validate_split_membership(wrong_size, universe, cell_types)


def test_reference_only_selection_hook_never_receives_heldout_rows():
    universe = toy_universe(("A", "B"), cells_per_type=5)
    membership = stratified_reference_heldout_split(
        universe, 9, cell_types=("A", "B")
    )
    reference = ad.AnnData(
        X=sparse.csr_matrix(np.ones((10, 3), dtype=np.int64)),
        obs=universe.set_index("barcode"),
        var=pd.DataFrame(index=["p0", "p1", "p2"]),
    )
    seen = {}

    def selector(training, cell_type_key, **kwargs):
        seen["barcodes"] = tuple(training.obs_names)
        seen["cell_type_key"] = cell_type_key
        seen["kwargs"] = kwargs
        return FakeSelection(("p2", "p0"), tuple(reference.var_names))

    result = select_reference_only_peaks(
        reference,
        membership,
        cell_types=("A", "B"),
        n_top_peaks=2,
        min_reference_cells=1,
        selector=selector,
    )
    expected = tuple(membership.loc[membership["pool"] == "reference", "barcode"])
    heldout = set(membership.loc[membership["pool"] == "heldout", "barcode"])
    assert seen["barcodes"] == expected
    assert heldout.isdisjoint(seen["barcodes"])
    assert seen["cell_type_key"] == "cell_type"
    assert seen["kwargs"]["cell_types"] == ["A", "B"]
    assert seen["kwargs"]["scale"] == 1.0e4
    assert result.peak_ids == ("p2", "p0")


def test_official_order_union_does_not_change_per_split_rankings():
    peaks = tuple(
        PeakInterval("chr1", index * 10, (index + 1) * 10, f"p{index}")
        for index in range(5)
    )
    first = FakeSelection(("p3", "p1"), tuple(peak.name for peak in peaks))
    second = FakeSelection(("p2", "p3"), tuple(peak.name for peak in peaks))
    union = official_order_peak_union([first, second], peaks)
    assert tuple(peak.name for peak in union) == ("p1", "p2", "p3")
    assert first.peak_ids == ("p3", "p1")
    assert second.peak_ids == ("p2", "p3")


def test_shape_cache_slice_refreshes_matrix_metadata_but_preserves_source_qc():
    universe, cell_types, peaks, cache, reference = toy_cache_and_reference()
    validate_shape_cache(
        cache,
        universe,
        [peak.name for peak in peaks],
        source_sha256=SOURCE_HASHES,
        cell_types=cell_types,
    )
    validate_shape_cache_against_official(cache, reference, universe)
    original_qc = dict(cache.uns["fragment_shape"]["preprocessing_counters"])
    selected = derive_shape_pool(
        cache,
        ["B_01", "A_00"],
        ["chr1:200-300", "chr1:0-100"],
        split_sha256="c" * 64,
        pool="reference",
    )

    assert tuple(selected.obs_names) == ("B_01", "A_00")
    assert tuple(selected.var_names) == ("chr1:200-300", "chr1:0-100")
    assert selected.uns["fragment_shape"]["feature_sha256"] == ordered_feature_sha256(
        selected.var_names
    )
    assert selected.uns["fragment_shape"]["split_sha256"] == "c" * 64
    assert selected.uns["fragment_shape"]["preprocessing_counters"] == original_qc
    total = sum(
        (selected.layers[layer] for layer in FRAGMENT_SHAPE_LAYER_NAMES),
        sparse.csr_matrix(selected.shape, dtype=np.int64),
    )
    assert matrices_equal(selected.X, total)
    counters = selected.uns["fragment_shape"]["matrix_counters"]
    assert counters["assigned_cut_sites"] == int(selected.X.sum())
    for layer in FRAGMENT_SHAPE_LAYER_NAMES:
        assert counters[f"cut_sites_per_bin.{layer}"] == int(selected.layers[layer].sum())

    broken = reference.copy()
    broken.X[0, 0] += 1
    with pytest.raises(ValueError, match="does not reconstruct"):
        validate_shape_cache_against_official(cache, broken, universe)


def test_prepare_split_orders_ranked_peaks_and_matches_official_counts():
    universe, cell_types, _, cache, reference = toy_cache_and_reference()

    def selector(training, cell_type_key, **kwargs):
        return FakeSelection(
            ("chr1:200-300", "chr1:0-100"), tuple(reference.var_names)
        )

    prepared = prepare_split_from_objects(
        reference,
        cache,
        universe,
        17,
        cell_types=cell_types,
        n_top_peaks=2,
        min_reference_cells=1,
        selector=selector,
    )
    reference_barcodes = prepared.membership.loc[
        prepared.membership["pool"] == "reference", "barcode"
    ].tolist()
    heldout_barcodes = prepared.membership.loc[
        prepared.membership["pool"] == "heldout", "barcode"
    ].tolist()
    assert set(reference_barcodes).isdisjoint(heldout_barcodes)
    assert tuple(prepared.reference_cells.var_names) == (
        "chr1:200-300",
        "chr1:0-100",
    )
    assert matrices_equal(
        prepared.reference_cells.X,
        reference[reference_barcodes, list(prepared.reference_cells.var_names)].X,
    )
    assert matrices_equal(
        prepared.heldout_test_cells.X,
        reference[heldout_barcodes, list(prepared.heldout_test_cells.var_names)].X,
    )


def test_split_outputs_and_manifest_are_published_atomically(tmp_path):
    universe, cell_types, _, cache, reference = toy_cache_and_reference()
    cache_roundtrip_path = tmp_path / "shape_cache_roundtrip.h5ad"
    cache.write_h5ad(cache_roundtrip_path)
    cache = ad.read_h5ad(cache_roundtrip_path)
    membership = stratified_reference_heldout_split(
        universe, 0, cell_types=cell_types
    )
    selection = FakeSelection(
        ("chr1:200-300", "chr1:0-100"), tuple(reference.var_names)
    )
    prepared = prepare_split_from_objects(
        reference,
        cache,
        universe,
        0,
        cell_types=cell_types,
        n_top_peaks=2,
        min_reference_cells=1,
        membership=membership,
        peak_selection=selection,
    )
    output_dir = tmp_path / "split_000"
    inputs = {
        "fragments": tmp_path / "fragments.tsv.gz",
        "tabix_index": tmp_path / "fragments.tsv.gz.tbi",
    }
    write_split_outputs(
        output_dir,
        prepared,
        outer_split_seed=0,
        cell_types=cell_types,
        input_paths=inputs,
        source_sha256=SOURCE_HASHES,
        n_top_peaks=2,
        min_reference_cells=1,
        shape_cache_path=tmp_path / "cache.h5ad",
    )

    assert sorted(path.name for path in output_dir.iterdir()) == [
        "heldout_test_cells.h5ad",
        "manifest.yaml",
        "peak_selection.csv",
        "reference_cells.h5ad",
        "selected_peaks.txt",
        "signal_audit.csv",
        "signal_audit_summary.yaml",
        "split.csv",
    ]
    manifest = yaml.safe_load((output_dir / "manifest.yaml").read_text())
    assert manifest["benchmark_scope"] == "development_smoke"
    assert manifest["limitation"].startswith("Conditional resampling from one donor")
    assert manifest["split"]["sha256"] == sha256_file(output_dir / "split.csv")
    assert manifest["rng"] == {
        "numpy_version": np.__version__,
        "bit_generator": "PCG64",
        "seed_sequence": "SeedSequence([20260822, outer_split_seed, one_based_cell_type_order])",
        "namespace": 20260822,
        "outer_split_seed": 0,
    }
    assert [record["path"] for record in manifest["outputs"]] == [
        "reference_cells.h5ad",
        "heldout_test_cells.h5ad",
        "split.csv",
        "selected_peaks.txt",
        "peak_selection.csv",
        "signal_audit.csv",
        "signal_audit_summary.yaml",
    ]
    assert manifest["signal_audit"] == {
        "scope": "training_reference_only",
        "heldout_counts_or_labels_read": False,
        "review_required_before_step_4": True,
        "summary_path": "signal_audit_summary.yaml",
    }
    assert (output_dir / "selected_peaks.txt").read_text().splitlines() == list(
        selection.peak_ids
    )
    restored_reference = ad.read_h5ad(output_dir / "reference_cells.h5ad")
    restored_metadata = restored_reference.uns["fragment_shape"]
    assert restored_metadata["source_sha256"] == SOURCE_HASHES
    assert restored_metadata["split_sha256"] == prepared.split_sha256
    assert restored_metadata["feature_sha256"] == selection.selected_feature_sha256
    assert restored_metadata["preprocessing_counters"] == cache.uns["fragment_shape"][
        "preprocessing_counters"
    ]
    assert restored_metadata["matrix_counters"]["assigned_cut_sites"] == int(
        restored_reference.X.sum()
    )
    with pytest.raises(FileExistsError, match="will not be overwritten"):
        write_split_outputs(
            output_dir,
            prepared,
            outer_split_seed=0,
            cell_types=cell_types,
            input_paths=inputs,
            source_sha256=SOURCE_HASHES,
            n_top_peaks=2,
            min_reference_cells=1,
            shape_cache_path=tmp_path / "cache.h5ad",
        )


def test_atomic_directory_removes_staging_tree_after_failure(tmp_path):
    destination = tmp_path / "not_published"
    with pytest.raises(RuntimeError, match="injected"):
        with atomic_output_directory(destination) as staging:
            (staging / "partial.txt").write_text("partial")
            raise RuntimeError("injected failure")
    assert not destination.exists()
    assert not list(tmp_path.glob(".not_published.staging.*"))
