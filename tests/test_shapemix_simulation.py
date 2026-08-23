from __future__ import annotations

import copy
import hashlib
import json
from collections import Counter
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse

from deconvatac.data import (
    DeconvolutionInput,
    FragmentShapeSpec,
    ordered_feature_sha256,
    validate_deconvolution_input,
)
from scripts.regenerate_shapemix_pbmc_simulations import (
    PBMC_CELL_TYPES,
    SEED_NAMESPACE,
    ShapeMixtureSimulation,
    condition_probabilities,
    dataset_id_for_simulation,
    is_primary_simulation,
    main,
    simulate_shapemix_spots,
    subset_shape_cells,
    write_simulation_dataset,
)


CELL_TYPES = ("A", "B", "C")
LAYERS = (
    "fragment_length_lt_100",
    "fragment_length_100_249",
    "fragment_length_ge_250",
)


def _fragment_metadata(adata: ad.AnnData) -> dict:
    layer_totals = {layer: int(adata.layers[layer].sum()) for layer in LAYERS}
    assigned = sum(layer_totals.values())
    retained = max((assigned + 1) // 2, 1)
    assigned_fragments = max((assigned + 1) // 2, 1)
    return {
        "schema_version": 1,
        "axis": "parent_fragment_length_bp",
        "count_unit": "deduplicated_cut_sites",
        "read_support_policy": "ignore",
        "peak_assignment": "containing_nonoverlapping_peak",
        "left_cut_offset": 0,
        "right_cut_offset": 0,
        "bins": {
            "0": {
                "name": "short",
                "order": 0,
                "min_inclusive": 0,
                "max_exclusive": 100,
                "layer": LAYERS[0],
            },
            "1": {
                "name": "mono",
                "order": 1,
                "min_inclusive": 100,
                "max_exclusive": 250,
                "layer": LAYERS[1],
            },
            "2": {
                "name": "long",
                "order": 2,
                "min_inclusive": 250,
                "layer": LAYERS[2],
            },
        },
        "source_sha256": {"fragments": "a" * 64, "tabix_index": "b" * 64},
        "feature_sha256": ordered_feature_sha256(adata.var_names),
        "split_sha256": "c" * 64,
        "coordinate_validation": {
            "selected_right_cut_offset": 0,
            "mismatched_entries": 0,
            "absolute_error": 0,
            "matrix_match": "exact",
        },
        "software_versions": {"deconvatac": "0.0.1", "pysam": "0.24.0"},
        "preprocessing_counters": {
            "total_rows": retained,
            "header_rows": 0,
            "invalid_schema_rows": 0,
            "invalid_coordinate_rows": 0,
            "unknown_barcodes": 0,
            "filtered_contigs": 0,
            "valid_rows": retained,
            "retained_fragments": retained,
            "fragments_with_assigned_cut_sites": assigned_fragments,
            "cut_sites_outside_peaks": 2 * retained - assigned,
            "assigned_cut_sites": assigned,
            "read_support_total": retained,
            **{f"cut_sites_per_bin.{layer}": count for layer, count in layer_totals.items()},
        },
        "matrix_counters": {
            "assigned_cut_sites": assigned,
            **{f"cut_sites_per_bin.{layer}": count for layer, count in layer_totals.items()},
        },
    }


def _shape_cells(prefix: str = "h", n_features: int = 4) -> ad.AnnData:
    names = [f"{prefix}_{cell_type.lower()}_{index}" for cell_type in CELL_TYPES for index in range(4)]
    labels = [cell_type for cell_type in CELL_TYPES for _ in range(4)]
    features = pd.Index(
        [f"chr{1 + index // 50}:{20 * index}-{20 * index + 10}" for index in range(n_features)]
    )
    short_base = np.asarray(
        [[1 + row % 3, row % 2, 0, 1] for row in range(len(names))],
        dtype=np.int64,
    )
    mono_base = np.asarray(
        [[0, 1, 1 + row % 2, row % 3] for row in range(len(names))],
        dtype=np.int64,
    )
    long_base = np.asarray(
        [[row % 2, 0, 1, 1 + (row + 1) % 2] for row in range(len(names))],
        dtype=np.int64,
    )
    repeats = int(np.ceil(n_features / 4))
    short = np.tile(short_base, (1, repeats))[:, :n_features]
    mono = np.tile(mono_base, (1, repeats))[:, :n_features]
    long = np.tile(long_base, (1, repeats))[:, :n_features]
    adata = ad.AnnData(
        X=sparse.csr_matrix((len(names), len(features)), dtype=np.int64),
        obs=pd.DataFrame({"cell_type": labels}, index=names),
        var=pd.DataFrame(index=features),
    )
    for layer, values in zip(LAYERS, (short, mono, long)):
        adata.layers[layer] = sparse.csr_matrix(values)
    adata.X = sum(
        (adata.layers[layer] for layer in LAYERS),
        sparse.csr_matrix(adata.shape, dtype=np.int64),
    ).tocsr()
    adata.uns["fragment_shape"] = _fragment_metadata(adata)
    return adata


def _simulate(
    heldout: ad.AnnData,
    *,
    condition: str = "equal_celltype",
    outer_seed: int = 7,
    mixture_seed: int = 11,
    depth: float | None = None,
):
    observed = {"A": 9, "B": 3, "C": 2}
    probabilities = condition_probabilities(condition, CELL_TYPES, observed)
    return simulate_shapemix_spots(
        heldout,
        cell_types=CELL_TYPES,
        sampling_probabilities=probabilities,
        condition=condition,
        outer_split_seed=outer_seed,
        inner_mixture_seed=mixture_seed,
        num_spots=16,
        mean_cells_per_spot=4,
        grid_shape=(4, 4),
        reference_barcodes=["reference_only_1", "reference_only_2"],
        depth_retain_probability=depth,
    )


def _assert_sparse_equal(left, right) -> None:
    difference = left.tocsr() - right.tocsr()
    difference.eliminate_zeros()
    assert difference.nnz == 0


def test_exact_layer_aggregation_truth_grid_and_immutable_provenance() -> None:
    heldout = _shape_cells()
    source_metadata = copy.deepcopy(heldout.uns["fragment_shape"])

    result = _simulate(heldout)

    assert result.spatial.obs_names.tolist() == [f"spot_{index:04d}" for index in range(16)]
    assert result.spatial.obsm["spatial"].tolist() == [
        [column, row] for row in range(4) for column in range(4)
    ]
    for spot_index, row in enumerate(result.source_cells_by_spot):
        assert row["spot_id"] == f"spot_{spot_index:04d}"
        assert row["cell_count"] == len(row["cell_id"]) == len(row["cell_type"])
        source_positions = heldout.obs_names.get_indexer(row["cell_id"])
        assert (source_positions >= 0).all()
        for layer in LAYERS:
            expected = heldout.layers[layer][source_positions].sum(axis=0)
            _assert_sparse_equal(result.spatial.layers[layer][spot_index], sparse.csr_matrix(expected))

        counts = Counter(row["cell_type"])
        expected_truth = np.asarray([counts[cell_type] / row["cell_count"] for cell_type in CELL_TYPES])
        np.testing.assert_allclose(result.truth.iloc[spot_index].to_numpy(), expected_truth)
        for cell_type, count in counts.items():
            if count <= int((heldout.obs["cell_type"] == cell_type).sum()):
                ids = [
                    barcode
                    for barcode, sampled_type in zip(row["cell_id"], row["cell_type"])
                    if sampled_type == cell_type
                ]
                assert len(ids) == len(set(ids))

    layer_sum = sum(
        (result.spatial.layers[layer] for layer in LAYERS),
        sparse.csr_matrix(result.spatial.shape, dtype=np.int64),
    )
    _assert_sparse_equal(result.spatial.X, layer_sum)
    np.testing.assert_allclose(result.truth.sum(axis=1), 1.0)
    assert heldout.uns["fragment_shape"] == source_metadata
    spatial_metadata = result.spatial.uns["fragment_shape"]
    assert spatial_metadata["preprocessing_counters"] == source_metadata["preprocessing_counters"]
    assert spatial_metadata["source_sha256"] == source_metadata["source_sha256"]
    assert spatial_metadata["split_sha256"] == source_metadata["split_sha256"]
    assert spatial_metadata["feature_sha256"] == ordered_feature_sha256(result.spatial.var_names)
    assert spatial_metadata["matrix_counters"]["assigned_cut_sites"] == int(result.spatial.X.sum())

    reference = _shape_cells(prefix="r")
    declared = FragmentShapeSpec.from_mapping(
        {
            key: source_metadata[key]
            for key in (
                "schema_version",
                "axis",
                "count_unit",
                "read_support_policy",
                "peak_assignment",
                "bins",
            )
        }
    )
    validate_deconvolution_input(
        DeconvolutionInput(
            dataset_id="toy_shapemix",
            modality="atac",
            feature_set="all",
            reference=reference,
            spatial=result.spatial,
            labels_key="cell_type",
            truth=result.truth,
            fragment_shape=declared,
            cell_types=list(CELL_TYPES),
        )
    )


def test_sources_are_heldout_only_and_reference_overlap_is_rejected() -> None:
    heldout = _shape_cells()
    result = _simulate(heldout)
    heldout_barcodes = set(heldout.obs_names)
    assert all(
        barcode in heldout_barcodes
        for row in result.source_cells_by_spot
        for barcode in row["cell_id"]
    )

    with pytest.raises(ValueError, match="must be disjoint"):
        simulate_shapemix_spots(
            heldout,
            cell_types=CELL_TYPES,
            sampling_probabilities=condition_probabilities("equal_celltype", CELL_TYPES),
            condition="equal_celltype",
            outer_split_seed=7,
            inner_mixture_seed=11,
            num_spots=4,
            mean_cells_per_spot=2,
            reference_barcodes=[heldout.obs_names[0]],
        )


def test_seed_streams_are_deterministic_independent_and_input_order_invariant() -> None:
    heldout = _shape_cells()
    first = _simulate(heldout)
    # Exercise another condition between repeats: no process-global stream is shared.
    observed = _simulate(heldout, condition="observed_abundance")
    repeat = _simulate(heldout)
    reordered = _simulate(heldout[np.arange(heldout.n_obs)[::-1]].copy())

    assert first.source_cells_by_spot == repeat.source_cells_by_spot
    assert first.source_cells_by_spot == reordered.source_cells_by_spot
    pd.testing.assert_frame_equal(first.truth, repeat.truth)
    pd.testing.assert_frame_equal(first.truth, reordered.truth)
    for layer in LAYERS:
        _assert_sparse_equal(first.spatial.layers[layer], repeat.spatial.layers[layer])
        _assert_sparse_equal(first.spatial.layers[layer], reordered.spatial.layers[layer])
    assert first.source_cells_by_spot != observed.source_cells_by_spot
    assert first.seed_streams == {
        "cell_counts_and_types": (SEED_NAMESPACE, 7, 11, 0),
        "source_cells": (SEED_NAMESPACE, 7, 11, 0, 1),
    }
    changed_outer = _simulate(heldout, outer_seed=8)
    changed_inner = _simulate(heldout, mixture_seed=12)
    assert first.source_cells_by_spot != changed_outer.source_cells_by_spot
    assert first.source_cells_by_spot != changed_inner.source_cells_by_spot


def test_secondary_depth_thinning_conserves_layers_and_does_not_resample_cells() -> None:
    heldout = _shape_cells()
    primary = _simulate(heldout)
    thinned = _simulate(heldout, depth=0.35)
    thinned_repeat = _simulate(heldout, depth=0.35)

    assert primary.source_cells_by_spot == thinned.source_cells_by_spot
    pd.testing.assert_frame_equal(primary.truth, thinned.truth)
    assert int(thinned.spatial.X.sum()) < int(primary.spatial.X.sum())
    for layer in LAYERS:
        assert (thinned.spatial.layers[layer].toarray() <= primary.spatial.layers[layer].toarray()).all()
        _assert_sparse_equal(thinned.spatial.layers[layer], thinned_repeat.spatial.layers[layer])
    layer_sum = sum(
        (thinned.spatial.layers[layer] for layer in LAYERS),
        sparse.csr_matrix(thinned.spatial.shape, dtype=np.int64),
    )
    _assert_sparse_equal(thinned.spatial.X, layer_sum)
    assert thinned.spatial.uns["simulation"]["depth_thinning_enabled"]
    assert thinned.depth_retain_probability == 0.35
    assert set(thinned.seed_streams).issuperset(
        {f"depth_thinning.{layer}" for layer in LAYERS}
    )
    assert dataset_id_for_simulation("equal_celltype", 7, 11).endswith("split_007_mix_011")
    assert dataset_id_for_simulation("equal_celltype", 7, 11, 0.35).endswith(
        "split_007_mix_011_depth_keep_0p35"
    )
    assert dataset_id_for_simulation("equal_celltype", 7, 11, smoke=True).endswith(
        "split_007_mix_011_smoke"
    )


def test_shape_subsetting_preserves_source_qc_and_refreshes_matrix_metadata() -> None:
    cells = _shape_cells()
    source_metadata = copy.deepcopy(cells.uns["fragment_shape"])
    selected = subset_shape_cells(
        cells,
        observation_mask=cells.obs["cell_type"].isin(["A", "C"]),
        feature_names=[cells.var_names[2], cells.var_names[0]],
    )

    assert selected.obs["cell_type"].tolist() == ["A"] * 4 + ["C"] * 4
    assert selected.var_names.tolist() == [cells.var_names[2], cells.var_names[0]]
    metadata = selected.uns["fragment_shape"]
    assert metadata["preprocessing_counters"] == source_metadata["preprocessing_counters"]
    assert metadata["split_sha256"] == source_metadata["split_sha256"]
    assert metadata["feature_sha256"] == ordered_feature_sha256(selected.var_names)
    assert metadata["matrix_counters"]["assigned_cut_sites"] == int(selected.X.sum())
    assert cells.uns["fragment_shape"] == source_metadata


def test_atomic_writer_records_hashes_and_refuses_overwrite(tmp_path: Path) -> None:
    heldout = _shape_cells()
    reference = _shape_cells(prefix="r")
    result = _simulate(heldout)
    source_dir = tmp_path / "split_007"
    source_dir.mkdir()
    reference_path = source_dir / "reference_cells.h5ad"
    heldout_path = source_dir / "heldout_test_cells.h5ad"
    split_manifest_path = source_dir / "manifest.yaml"
    reference.write_h5ad(reference_path)
    heldout.write_h5ad(heldout_path)
    split_manifest_path.write_text("schema_version: 1\nstatus: complete\n")

    dataset_id = dataset_id_for_simulation("equal_celltype", 7, 11)
    output_root = tmp_path / "datasets"
    config_path = write_simulation_dataset(
        result,
        output_root=output_root,
        dataset_id=dataset_id,
        reference_path=reference_path,
        heldout_path=heldout_path,
        split_manifest_path=split_manifest_path,
    )

    target = output_root / dataset_id
    assert config_path == target / "dataset.yaml"
    assert sorted(path.relative_to(target).as_posix() for path in target.rglob("*") if path.is_file()) == [
        "atac/features/highly_variable.txt",
        "atac/spatial.h5ad",
        "dataset.yaml",
        "simulation/manifest.yaml",
        "simulation/source_cells_by_spot.jsonl",
        "truth/proportions.csv",
    ]
    config = yaml.safe_load(config_path.read_text())
    assert config["modalities"]["atac"]["truth"]["cell_types"] == list(CELL_TYPES)
    assert config["simulation"]["primary_dataset"] is False
    assert config["benchmark_scope"] == "development"
    manifest = yaml.safe_load((target / "simulation" / "manifest.yaml").read_text())
    assert manifest["scientific_scope"].startswith("Conditional resampling variability")
    assert manifest["simulation"]["source_sampling"]["heldout_only"] is True
    for output in manifest["outputs"].values():
        output_path = Path(output["path"])
        if not output_path.is_absolute():
            output_path = Path.cwd() / output_path
        assert output_path.stat().st_size == output["bytes"]
        assert hashlib.sha256(output_path.read_bytes()).hexdigest() == output["sha256"]

    provenance = [
        json.loads(line)
        for line in (target / "simulation" / "source_cells_by_spot.jsonl").read_text().splitlines()
    ]
    assert provenance == list(result.source_cells_by_spot)
    roundtrip = ad.read_h5ad(target / "atac" / "spatial.h5ad")
    for layer in LAYERS:
        assert sparse.isspmatrix_csr(roundtrip.layers[layer])
        _assert_sparse_equal(roundtrip.layers[layer], result.spatial.layers[layer])

    original_hash = hashlib.sha256(config_path.read_bytes()).hexdigest()
    with pytest.raises(FileExistsError, match="never overwritten"):
        write_simulation_dataset(
            result,
            output_root=output_root,
            dataset_id=dataset_id,
            reference_path=reference_path,
            heldout_path=heldout_path,
            split_manifest_path=split_manifest_path,
        )
    assert hashlib.sha256(config_path.read_bytes()).hexdigest() == original_hash
    assert not list(output_root.glob(f".{dataset_id}.*"))


def test_probability_contracts_are_explicit() -> None:
    equal = condition_probabilities("equal_celltype", CELL_TYPES)
    assert equal == {"A": 1 / 3, "B": 1 / 3, "C": 1 / 3}
    observed = condition_probabilities("observed_abundance", CELL_TYPES, {"A": 5, "B": 3, "C": 2})
    assert observed == {"A": 0.5, "B": 0.3, "C": 0.2}
    with pytest.raises(ValueError, match="exactly"):
        condition_probabilities("observed_abundance", CELL_TYPES, {"A": 1, "B": 1})


@pytest.mark.parametrize("condition", ["equal_celltype", "observed_abundance"])
def test_primary_classification_requires_all_frozen_fields(condition: str) -> None:
    spatial = ad.AnnData(X=sparse.csr_matrix((1024, 5000), dtype=np.int64))
    spatial.uns["simulation"] = {
        "grid_rows": 32,
        "grid_columns": 32,
        "mean_cells_per_spot": 10.0,
    }
    probabilities = condition_probabilities(condition, PBMC_CELL_TYPES)
    normalized = np.asarray(list(probabilities.values()), dtype=np.float64)
    normalized /= normalized.sum()
    stored_probabilities = {
        cell_type: float(normalized[index])
        for index, cell_type in enumerate(PBMC_CELL_TYPES)
    }
    simulation = ShapeMixtureSimulation(
        spatial=spatial,
        truth=pd.DataFrame(),
        source_cells_by_spot=(),
        cell_types=PBMC_CELL_TYPES,
        sampling_probabilities=stored_probabilities,
        condition=condition,
        outer_split_seed=1103,
        inner_mixture_seed=101,
        seed_streams={},
        depth_retain_probability=None,
    )
    assert is_primary_simulation(simulation)

    narrowed = replace(simulation, spatial=simulation.spatial[:, :4999].copy())
    assert not is_primary_simulation(narrowed)


def test_smoke_cli_validates_split_manifest_and_writes_sliced_reference(tmp_path: Path) -> None:
    reference = _shape_cells(prefix="r", n_features=100)
    heldout = _shape_cells(n_features=100)
    split_dir = tmp_path / "split_000"
    split_dir.mkdir()
    membership = pd.DataFrame(
        [
            {
                "barcode": str(barcode),
                "cell_type": str(adata.obs.loc[barcode, "cell_type"]),
                "pool": pool,
            }
            for pool, adata in (("reference", reference), ("heldout", heldout))
            for barcode in adata.obs_names
        ]
    )
    split_path = split_dir / "split.csv"
    membership.to_csv(split_path, index=False, lineterminator="\n")
    split_sha256 = hashlib.sha256(split_path.read_bytes()).hexdigest()
    for pool, adata in (("reference", reference), ("heldout", heldout)):
        adata.obs["split_pool"] = pool
        adata.uns["fragment_shape"]["split_sha256"] = split_sha256
        adata.uns["shapemix_preparation"] = {
            "pool": pool,
            "split_sha256": split_sha256,
        }

    reference_path = split_dir / "reference_cells.h5ad"
    heldout_path = split_dir / "heldout_test_cells.h5ad"
    selected_path = split_dir / "selected_peaks.txt"
    reference.write_h5ad(reference_path)
    heldout.write_h5ad(heldout_path)
    selected_path.write_text("\n".join(heldout.var_names) + "\n")

    def output_record(path: Path, role: str) -> dict:
        return {
            "path": path.name,
            "role": role,
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }

    with (split_dir / "manifest.yaml").open("w") as handle:
        yaml.safe_dump(
            {
                "schema_version": 1,
                "benchmark_scope": "development_smoke",
                "cell_types": list(CELL_TYPES),
                "rng": {"outer_split_seed": 0},
                "split": {
                    "sha256": split_sha256,
                    "reference_cells": reference.n_obs,
                    "heldout_cells": heldout.n_obs,
                    "counts_by_cell_type": {
                        cell_type: {"reference": 4, "heldout": 4}
                        for cell_type in CELL_TYPES
                    },
                },
                "outputs": [
                    output_record(reference_path, "training_reference_fragment_shapes"),
                    output_record(heldout_path, "heldout_source_fragment_shapes"),
                    output_record(split_path, "canonical_split_membership"),
                    output_record(selected_path, "ranked_selected_peak_ids"),
                ],
            },
            handle,
            sort_keys=False,
        )
    output_root = tmp_path / "datasets"

    written = main(
        [
            "--split-dir",
            str(split_dir),
            "--outer-split-seed",
            "0",
            "--inner-mixture-seed",
            "0",
            "--conditions",
            "equal_celltype",
            "--output-root",
            str(output_root),
            "--num-spots",
            "32",
            "--mean-cells-per-spot",
            "2",
            "--grid-shape",
            "4",
            "8",
            "--cell-types",
            *CELL_TYPES,
            "--smoke",
        ]
    )

    assert len(written) == 1
    assert written[0].parent.name.endswith("_smoke")
    config = yaml.safe_load(written[0].read_text())
    assert config["benchmark_scope"] == "development_smoke"
    assert config["simulation"]["primary_dataset"] is False
    local_reference = written[0].parent / "atac" / "reference_cells.h5ad"
    assert local_reference.is_file()
    assert ad.read_h5ad(local_reference).var_names.tolist() == heldout.var_names.tolist()
    assert ad.read_h5ad(written[0].parent / "atac" / "spatial.h5ad").shape == (32, 100)

    bad_manifest = yaml.safe_load((split_dir / "manifest.yaml").read_text())
    bad_manifest["rng"]["outer_split_seed"] = 1
    with (split_dir / "manifest.yaml").open("w") as handle:
        yaml.safe_dump(bad_manifest, handle)
    with pytest.raises(ValueError, match="outer_split_seed"):
        main(
            [
                "--split-dir",
                str(split_dir),
                "--outer-split-seed",
                "0",
                "--inner-mixture-seed",
                "0",
                "--conditions",
                "equal_celltype",
                "--output-root",
                str(output_root),
                "--num-spots",
                "32",
                "--cell-types",
                *CELL_TYPES,
                "--smoke",
            ]
        )

    bad_manifest["rng"]["outer_split_seed"] = 0
    with (split_dir / "manifest.yaml").open("w") as handle:
        yaml.safe_dump(bad_manifest, handle, sort_keys=False)
    with pytest.raises(ValueError, match="exactly 32 spots"):
        main(
            [
                "--split-dir",
                str(split_dir),
                "--outer-split-seed",
                "0",
                "--inner-mixture-seed",
                "0",
                "--conditions",
                "equal_celltype",
                "--output-root",
                str(output_root),
                "--num-spots",
                "4",
                "--cell-types",
                *CELL_TYPES,
                "--smoke",
            ]
        )

    heldout_output = next(
        record
        for record in bad_manifest["outputs"]
        if record["role"] == "heldout_source_fragment_shapes"
    )
    heldout_output["sha256"] = "0" * 64
    with (split_dir / "manifest.yaml").open("w") as handle:
        yaml.safe_dump(bad_manifest, handle, sort_keys=False)
    with pytest.raises(ValueError, match="SHA-256"):
        main(
            [
                "--split-dir",
                str(split_dir),
                "--outer-split-seed",
                "0",
                "--inner-mixture-seed",
                "0",
                "--conditions",
                "equal_celltype",
                "--output-root",
                str(output_root),
                "--num-spots",
                "32",
                "--cell-types",
                *CELL_TYPES,
                "--smoke",
            ]
        )
