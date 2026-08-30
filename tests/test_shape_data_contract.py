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
    load_deconvolution_input,
    ordered_feature_sha256,
    validate_deconvolution_input,
)


DECLARED_FRAGMENT_SHAPE = {
    "schema_version": 1,
    "axis": "parent_fragment_length_bp",
    "count_unit": "deduplicated_cut_sites",
    "read_support_policy": "ignore",
    "peak_assignment": "containing_nonoverlapping_peak",
    "bins": [
        {
            "name": "short",
            "min_inclusive": 0,
            "max_exclusive": 100,
            "layer": "fragment_length_lt_100",
        },
        {
            "name": "mono",
            "min_inclusive": 100,
            "max_exclusive": 250,
            "layer": "fragment_length_100_249",
        },
        {
            "name": "long",
            "min_inclusive": 250,
            "max_exclusive": None,
            "layer": "fragment_length_ge_250",
        },
    ],
}


def _declared_spec() -> FragmentShapeSpec:
    return FragmentShapeSpec.from_mapping(DECLARED_FRAGMENT_SHAPE)


def _resolved_spec(
    rows: int,
    feature_names: pd.Index,
    layer_totals: dict[str, int],
) -> FragmentShapeSpec:
    assigned_cut_sites = sum(layer_totals.values())
    return replace(
        _declared_spec(),
        left_cut_offset=0,
        right_cut_offset=0,
        source_sha256={"fragments": "a" * 64, "tabix_index": "b" * 64},
        feature_sha256=ordered_feature_sha256(feature_names),
        split_sha256="d" * 64,
        coordinate_validation={
            "selected_right_cut_offset": 0,
            "mismatched_entries": 0,
            "absolute_error": 0,
            "matrix_match": "exact",
        },
        software_versions={"deconvatac": "0.0.1", "pysam": "0.24.0"},
        preprocessing_counters={
            "total_rows": rows,
            "header_rows": 0,
            "invalid_schema_rows": 0,
            "invalid_coordinate_rows": 0,
            "unknown_barcodes": 0,
            "filtered_contigs": 0,
            "valid_rows": rows,
            "retained_fragments": rows,
            "fragments_with_assigned_cut_sites": rows,
            "cut_sites_outside_peaks": 2 * rows - assigned_cut_sites,
            "assigned_cut_sites": assigned_cut_sites,
            "read_support_total": rows,
            "cut_sites_per_bin": layer_totals,
        },
        matrix_counters={
            "assigned_cut_sites": assigned_cut_sites,
            **{
                f"cut_sites_per_bin.{layer_name}": count
                for layer_name, count in layer_totals.items()
            },
        },
    )


def _add_shape_layers(adata: ad.AnnData, values: tuple[np.ndarray, np.ndarray, np.ndarray], rows: int) -> None:
    for layer_name, layer_values in zip(_declared_spec().layer_names, values):
        adata.layers[layer_name] = sparse.csr_matrix(layer_values)
    adata.X = sum((adata.layers[name] for name in _declared_spec().layer_names), sparse.csr_matrix(adata.shape))
    layer_totals = {name: int(adata.layers[name].sum()) for name in _declared_spec().layer_names}
    adata.uns["fragment_shape"] = _resolved_spec(rows, adata.var_names, layer_totals).to_uns()


def _refresh_stored_matrix_metadata(adata: ad.AnnData) -> None:
    metadata = adata.uns["fragment_shape"]
    metadata["feature_sha256"] = ordered_feature_sha256(adata.var_names)
    layer_totals = {name: int(adata.layers[name].sum()) for name in _declared_spec().layer_names}
    metadata["matrix_counters"] = {
        "assigned_cut_sites": sum(layer_totals.values()),
        **{
            f"cut_sites_per_bin.{layer_name}": count
            for layer_name, count in layer_totals.items()
        },
    }


def _shape_input() -> DeconvolutionInput:
    reference = ad.AnnData(
        X=sparse.csr_matrix((3, 2)),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B"]}, index=["c1", "c2", "c3"]),
        var=pd.DataFrame(index=["p1", "p2"]),
    )
    _add_shape_layers(
        reference,
        (
            np.array([[1, 0], [2, 0], [0, 1]]),
            np.array([[0, 1], [0, 1], [2, 0]]),
            np.array([[0, 0], [1, 0], [0, 1]]),
        ),
        rows=6,
    )

    spatial = ad.AnnData(
        X=sparse.csr_matrix((2, 2)),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["p1", "p2"]),
    )
    _add_shape_layers(
        spatial,
        (
            np.array([[3, 0], [0, 2]]),
            np.array([[0, 2], [1, 0]]),
            np.array([[1, 0], [0, 1]]),
        ),
        rows=5,
    )
    spatial.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 0.0]])
    truth = pd.DataFrame([[0.75, 0.25], [0.1, 0.9]], index=spatial.obs_names, columns=["A", "B"])

    return DeconvolutionInput(
        dataset_id="shape_toy",
        modality="atac",
        feature_set="all",
        spatial=spatial,
        reference=reference,
        labels_key="cell_type",
        truth=truth,
        fragment_shape=_declared_spec(),
        cell_types=["A", "B"],
    )


def test_valid_fragment_shape_contract() -> None:
    data = _shape_input()

    validate_deconvolution_input(data)

    assert data.fragment_shape is not None
    assert data.fragment_shape.layer_names == (
        "fragment_length_lt_100",
        "fragment_length_100_249",
        "fragment_length_ge_250",
    )


def test_ordered_feature_sha256_has_a_stable_unambiguous_encoding() -> None:
    assert ordered_feature_sha256(["p1", "p2"]) == (
        "854dfd24c1d6de85ae577ff7248a58f4320df1d580719945fd8868c2d3185b68"
    )
    assert ordered_feature_sha256(["ab", "c"]) != ordered_feature_sha256(["a", "bc"])
    assert ordered_feature_sha256(["p1", "p2"]) != ordered_feature_sha256(["p2", "p1"])


def test_fragment_shape_is_optional_for_legacy_inputs() -> None:
    data = _shape_input()
    data.fragment_shape = None
    data.cell_types = None
    data.reference.layers.clear()
    data.spatial.layers.clear()
    data.reference.uns.clear()
    data.spatial.uns.clear()
    data.reference.X = np.ones(data.reference.shape)
    data.spatial.X = np.ones(data.spatial.shape)

    validate_deconvolution_input(data)


def test_fragment_shape_requires_declared_cell_types() -> None:
    data = _shape_input()
    data.cell_types = None

    with pytest.raises(ValueError, match="cell_types must be declared"):
        validate_deconvolution_input(data)


def test_declared_cell_types_without_truth_are_valid_for_real_data() -> None:
    data = _shape_input()
    data.truth = None

    validate_deconvolution_input(data)


def test_missing_fragment_shape_layer_is_rejected() -> None:
    data = _shape_input()
    del data.spatial.layers["fragment_length_100_249"]

    with pytest.raises(ValueError, match="missing declared layer"):
        validate_deconvolution_input(data)


@pytest.mark.parametrize("matrix_kind", ["dense", "csc"])
def test_fragment_shape_layers_must_be_csr(matrix_kind: str) -> None:
    data = _shape_input()
    layer_name = "fragment_length_lt_100"
    values = data.reference.layers[layer_name].toarray()
    data.reference.layers[layer_name] = values if matrix_kind == "dense" else sparse.csc_matrix(values)

    with pytest.raises(ValueError, match="scipy CSR"):
        validate_deconvolution_input(data)


def test_collapsed_matrix_must_equal_exact_layer_sum() -> None:
    data = _shape_input()
    data.spatial.X = data.spatial.X.copy()
    data.spatial.X[0, 0] += 1

    with pytest.raises(ValueError, match="exact elementwise sum"):
        validate_deconvolution_input(data)


@pytest.mark.parametrize(
    ("value", "message"),
    [(-1.0, "negative"), (1.5, "non-integer"), (np.nan, "non-finite")],
)
def test_invalid_fragment_shape_counts_are_rejected(value: float, message: str) -> None:
    data = _shape_input()
    layer = data.reference.layers["fragment_length_lt_100"].astype(float)
    layer.data[0] = value
    data.reference.layers["fragment_length_lt_100"] = layer

    with pytest.raises(ValueError, match=message):
        validate_deconvolution_input(data)


def test_shape_metadata_allows_independent_reference_and_spatial_provenance() -> None:
    data = _shape_input()
    data.spatial.uns["fragment_shape"]["split_sha256"] = "e" * 64
    data.spatial.uns["fragment_shape"]["source_sha256"] = {
        "fragments": "f" * 64,
        "tabix_index": "0" * 64,
    }
    data.spatial.uns["fragment_shape"]["software_versions"] = {
        "deconvatac": "0.0.2",
        "pysam": "0.24.1",
    }

    validate_deconvolution_input(data)


def test_shape_metadata_must_match_yaml_and_both_feature_axes() -> None:
    data = _shape_input()
    data.spatial.uns["fragment_shape"]["feature_sha256"] = "e" * 64

    with pytest.raises(ValueError, match="current ordered feature axis"):
        validate_deconvolution_input(data)


def test_shape_metadata_requires_resolved_provenance() -> None:
    data = _shape_input()
    del data.reference.uns["fragment_shape"]["split_sha256"]

    with pytest.raises(ValueError, match="split_sha256"):
        validate_deconvolution_input(data)


def test_feature_hash_must_match_the_current_ordered_axis() -> None:
    data = _shape_input()
    data.reference.uns["fragment_shape"]["feature_sha256"] = ordered_feature_sha256(["p2", "p1"])

    with pytest.raises(ValueError, match="does not match the current ordered feature axis"):
        validate_deconvolution_input(data)


@pytest.mark.parametrize("missing_source", ["fragments", "tabix_index"])
def test_source_hashes_require_fragments_and_tabix_index(missing_source: str) -> None:
    data = _shape_input()
    del data.reference.uns["fragment_shape"]["source_sha256"][missing_source]

    with pytest.raises(ValueError, match=f"missing required sources: {missing_source}"):
        validate_deconvolution_input(data)


@pytest.mark.parametrize(
    ("field_name", "value"),
    [
        ("selected_right_cut_offset", -1),
        ("mismatched_entries", 1),
        ("absolute_error", 1),
    ],
)
def test_coordinate_validation_must_record_an_exact_offset_zero_match(
    field_name: str,
    value: int,
) -> None:
    data = _shape_input()
    data.reference.uns["fragment_shape"]["coordinate_validation"][field_name] = value

    with pytest.raises(ValueError, match=field_name):
        validate_deconvolution_input(data)


def test_coordinate_validation_accepts_exact_offset_minus_one() -> None:
    data = _shape_input()
    for adata in (data.reference, data.spatial):
        adata.uns["fragment_shape"]["right_cut_offset"] = -1
        adata.uns["fragment_shape"]["coordinate_validation"][
            "selected_right_cut_offset"
        ] = -1

    validate_deconvolution_input(data)


def test_coordinate_validation_requires_matrix_match_field() -> None:
    data = _shape_input()
    del data.reference.uns["fragment_shape"]["coordinate_validation"]["matrix_match"]

    with pytest.raises(ValueError, match="missing required fields: matrix_match"):
        validate_deconvolution_input(data)

def test_coordinate_validation_accepts_exact_strand_aware_bedpe_semantics() -> None:
    data = _shape_input()
    semantic_validation = {
        "selected_right_cut_offset": 0,
        "matrix_match": "not_available",
        "validation_method": "deposited_strand_aware_bedpe_5prime",
        "semantic_match": "exact",
    }
    for adata in (data.reference, data.spatial):
        adata.uns["fragment_shape"]["coordinate_validation"] = semantic_validation.copy()

    validate_deconvolution_input(data)


def test_coordinate_validation_accepts_exact_bed_parent_fragment_semantics() -> None:
    data = _shape_input()
    semantic_validation = {
        "selected_right_cut_offset": 0,
        "matrix_match": "not_available",
        "validation_method": "deposited_bed_half_open_parent_fragments",
        "semantic_match": "exact",
        "fragment_total_match": "exact",
        "audit": "data/processed/shapemix/example/manifests/coordinate_audit.yaml",
    }
    for adata in (data.reference, data.spatial):
        adata.uns["fragment_shape"]["coordinate_validation"] = semantic_validation.copy()

    validate_deconvolution_input(data)


def test_bed_parent_fragment_source_basename_is_accepted() -> None:
    data = _shape_input()
    semantic_validation = {
        "selected_right_cut_offset": 0,
        "matrix_match": "not_available",
        "validation_method": "deposited_bed_half_open_parent_fragments",
        "semantic_match": "exact",
        "fragment_total_match": "exact",
        "audit": "data/processed/shapemix/example/manifests/coordinate_audit.yaml",
    }
    for adata in (data.reference, data.spatial):
        adata.uns["fragment_shape"]["source_sha256"] = {
            "GSM6671097_E11A.bed.gz": "a" * 64
        }
        adata.uns["fragment_shape"]["coordinate_validation"] = (
            semantic_validation.copy()
        )

    validate_deconvolution_input(data)


def test_coordinate_validation_rejects_bed_parent_fragments_without_total_match() -> None:
    data = _shape_input()
    semantic_validation = {
        "selected_right_cut_offset": 0,
        "matrix_match": "not_available",
        "validation_method": "deposited_bed_half_open_parent_fragments",
        "semantic_match": "exact",
        "fragment_total_match": "mismatch",
        "audit": "data/processed/shapemix/example/manifests/coordinate_audit.yaml",
    }
    data.reference.uns["fragment_shape"]["coordinate_validation"] = semantic_validation
    with pytest.raises(ValueError, match="fragment_total_match must be 'exact'"):
        validate_deconvolution_input(data)


def test_coordinate_validation_accepts_passed_representative_matrix_concordance() -> None:
    data = _shape_input()
    representative_validation = {
        "selected_right_cut_offset": 0,
        "matrix_match": "representative_recovered_read_concordance",
        "validation_method": "deposited_matrix_concordance_thresholds",
        "passed": True,
        "audit": "data/processed/shapemix/example/manifests/coordinate_audit.yaml",
    }
    for adata in (data.reference, data.spatial):
        adata.uns["fragment_shape"]["coordinate_validation"] = (
            representative_validation.copy()
        )

    validate_deconvolution_input(data)


def test_coordinate_validation_rejects_failed_representative_matrix_concordance() -> None:
    data = _shape_input()
    representative_validation = {
        "selected_right_cut_offset": 0,
        "matrix_match": "representative_recovered_read_concordance",
        "validation_method": "deposited_matrix_concordance_thresholds",
        "passed": False,
        "audit": "data/processed/shapemix/example/manifests/coordinate_audit.yaml",
    }
    for adata in (data.reference, data.spatial):
        adata.uns["fragment_shape"]["coordinate_validation"] = (
            representative_validation.copy()
        )

    with pytest.raises(ValueError, match="passed must be true"):
        validate_deconvolution_input(data)


def test_multisample_fragment_source_hashes_preserve_source_names() -> None:
    data = _shape_input()
    multisample_sources = {
        "GSM1_sample.fragments.tsv.gz": "0" * 64,
        "GSM1_sample.fragments.tsv.gz.tbi": "1" * 64,
        "GSM2_sample.bedpe.gz": "2" * 64,
        "cell_annotations.tsv.gz": "3" * 64,
    }
    for adata in (data.reference, data.spatial):
        adata.uns["fragment_shape"]["source_sha256"] = multisample_sources.copy()

    validate_deconvolution_input(data)


def test_preprocessing_counters_require_every_qc_field() -> None:
    data = _shape_input()
    del data.reference.uns["fragment_shape"]["preprocessing_counters"]["header_rows"]

    with pytest.raises(ValueError, match="missing required fields: header_rows"):
        validate_deconvolution_input(data)


def test_read_support_total_must_cover_every_valid_row() -> None:
    data = _shape_input()
    counters = data.reference.uns["fragment_shape"]["preprocessing_counters"]
    counters["read_support_total"] = counters["valid_rows"] - 1

    with pytest.raises(ValueError, match="read_support_total must be at least valid_rows"):
        validate_deconvolution_input(data)


@pytest.mark.parametrize(
    ("field_name", "increment", "message"),
    [
        ("valid_rows", -1, "valid_rows.*total_rows"),
        ("cut_sites_outside_peaks", 1, "cut_sites_outside_peaks.*retained_fragments"),
        ("unknown_barcodes", 7, "cannot exceed valid_rows"),
    ],
)
def test_preprocessing_counter_identities_and_bounds_are_enforced(
    field_name: str,
    increment: int,
    message: str,
) -> None:
    data = _shape_input()
    counters = data.reference.uns["fragment_shape"]["preprocessing_counters"]
    counters[field_name] += increment

    with pytest.raises(ValueError, match=message):
        validate_deconvolution_input(data)


def test_per_bin_counters_must_sum_to_assigned_cut_sites() -> None:
    data = _shape_input()
    counters = data.reference.uns["fragment_shape"]["preprocessing_counters"]["cut_sites_per_bin"]
    counters["fragment_length_100_249"] += 1

    with pytest.raises(ValueError, match="sum to assigned_cut_sites"):
        validate_deconvolution_input(data)


def test_flattened_per_bin_counters_are_accepted() -> None:
    data = _shape_input()
    for adata in (data.reference, data.spatial):
        counters = adata.uns["fragment_shape"]["preprocessing_counters"]
        nested = counters.pop("cut_sites_per_bin")
        counters.update({f"cut_sites_per_bin.{layer}": count for layer, count in nested.items()})

    validate_deconvolution_input(data)


def test_matrix_counters_are_required() -> None:
    data = _shape_input()
    del data.reference.uns["fragment_shape"]["matrix_counters"]

    with pytest.raises(ValueError, match="matrix_counters must be a non-empty mapping"):
        validate_deconvolution_input(data)


@pytest.mark.parametrize("invalid_total", [0, False])
def test_matrix_counter_total_must_be_a_positive_integer(invalid_total) -> None:
    data = _shape_input()
    data.reference.uns["fragment_shape"]["matrix_counters"]["assigned_cut_sites"] = invalid_total

    with pytest.raises(ValueError, match="assigned_cut_sites must be a positive integer"):
        validate_deconvolution_input(data)


def test_matrix_per_bin_counters_must_match_current_layers() -> None:
    data = _shape_input()
    counters = data.reference.uns["fragment_shape"]["matrix_counters"]
    counters["cut_sites_per_bin.fragment_length_lt_100"] -= 1
    counters["cut_sites_per_bin.fragment_length_100_249"] += 1

    with pytest.raises(ValueError, match="does not match stored layer"):
        validate_deconvolution_input(data)


def test_matrix_counters_reject_an_empty_stored_matrix() -> None:
    data = _shape_input()
    for layer_name in data.fragment_shape.layer_names:
        data.reference.layers[layer_name] = sparse.csr_matrix(data.reference.shape, dtype=np.int64)
    data.reference.X = sparse.csr_matrix(data.reference.shape, dtype=np.int64)
    data.reference.uns["fragment_shape"]["matrix_counters"] = {
        "assigned_cut_sites": 0,
        **{
            f"cut_sites_per_bin.{layer_name}": 0
            for layer_name in data.fragment_shape.layer_names
        },
    }

    with pytest.raises(ValueError, match="assigned_cut_sites must be a positive integer"):
        validate_deconvolution_input(data)


def test_nested_matrix_bin_counters_are_accepted() -> None:
    data = _shape_input()
    for adata in (data.reference, data.spatial):
        counters = adata.uns["fragment_shape"]["matrix_counters"]
        counters["cut_sites_per_bin"] = {
            layer_name: counters.pop(f"cut_sites_per_bin.{layer_name}")
            for layer_name in data.fragment_shape.layer_names
        }

    validate_deconvolution_input(data)


def test_bin_order_and_semantics_are_versioned() -> None:
    data = _shape_input()
    reversed_bins = tuple(reversed(data.fragment_shape.bins))
    data.fragment_shape = replace(data.fragment_shape, bins=reversed_bins)

    with pytest.raises(ValueError, match="ordered schema-version-1 bins"):
        validate_deconvolution_input(data)


def test_declared_cell_types_must_match_reference_universe() -> None:
    data = _shape_input()
    data.cell_types = ["A", "C"]
    data.truth.columns = ["A", "C"]

    with pytest.raises(ValueError, match="reference label universe"):
        validate_deconvolution_input(data)


def test_declared_cell_type_order_must_match_truth_columns() -> None:
    data = _shape_input()
    data.cell_types = ["B", "A"]

    with pytest.raises(ValueError, match="truth columns"):
        validate_deconvolution_input(data)


def test_feature_slicing_preserves_shape_contract() -> None:
    data = _shape_input()
    data.reference = data.reference[:, ["p2"]].copy()
    data.spatial = data.spatial[:, ["p2"]].copy()

    with pytest.raises(ValueError, match="does not match the current ordered feature axis"):
        validate_deconvolution_input(data)

    _refresh_stored_matrix_metadata(data.reference)
    _refresh_stored_matrix_metadata(data.spatial)
    validate_deconvolution_input(data)

    for adata in (data.reference, data.spatial):
        expected = sum(
            (adata.layers[layer] for layer in data.fragment_shape.layer_names),
            sparse.csr_matrix(adata.shape),
        )
        assert (adata.X != expected).nnz == 0


def _write_shape_dataset(tmp_path: Path, *, with_truth: bool = True) -> Path:
    data = _shape_input()
    reference_path = tmp_path / "reference.h5ad"
    spatial_path = tmp_path / "spatial.h5ad"
    data.reference.write_h5ad(reference_path)
    if with_truth:
        data.spatial.obsm["truth"] = data.truth
    data.spatial.write_h5ad(spatial_path)

    config = {
        "dataset_id": "shape_toy",
        "modalities": {
            "atac": {
                "reference": {"path": str(reference_path)},
                "spatial": {"path": str(spatial_path)},
                "labels_key": "cell_type",
                "spatial_key": "spatial",
                "truth": {"obsm_key": "truth", "cell_types": ["A", "B"]},
                "fragment_shape": DECLARED_FRAGMENT_SHAPE,
                "feature_sets": {
                    "all": {"mode": "all"},
                    "single_peak": {"features": ["p2"]},
                },
            }
        },
    }
    if not with_truth:
        del config["modalities"]["atac"]["truth"]
        config["modalities"]["atac"]["cell_types"] = ["A", "B"]

    config_path = tmp_path / "shape_toy.yaml"
    config_path.write_text(yaml.safe_dump(config))
    registry_path = tmp_path / "datasets.yaml"
    registry_path.write_text(yaml.safe_dump({"shape_toy": {"config": str(config_path)}}))
    return registry_path


def test_loader_parses_shape_spec_and_ordered_cell_types(tmp_path: Path) -> None:
    registry_path = _write_shape_dataset(tmp_path)

    data = load_deconvolution_input(
        dataset_id="shape_toy",
        modality="atac",
        registry_path=registry_path,
    )

    assert isinstance(data.fragment_shape, FragmentShapeSpec)
    assert data.cell_types == ["A", "B"]
    assert list(data.truth.columns) == ["A", "B"]
    assert all(sparse.isspmatrix_csr(data.reference.layers[name]) for name in data.fragment_shape.layer_names)


def test_loader_accepts_declared_cell_types_without_truth(tmp_path: Path) -> None:
    registry_path = _write_shape_dataset(tmp_path, with_truth=False)

    data = load_deconvolution_input(
        dataset_id="shape_toy",
        modality="atac",
        registry_path=registry_path,
    )

    assert data.cell_types == ["A", "B"]
    assert data.truth is None


def test_loader_feature_slicing_preserves_shape_layers(tmp_path: Path) -> None:
    registry_path = _write_shape_dataset(tmp_path)

    data = load_deconvolution_input(
        dataset_id="shape_toy",
        modality="atac",
        feature_set="single_peak",
        registry_path=registry_path,
    )

    assert list(data.reference.var_names) == ["p2"]
    assert list(data.spatial.var_names) == ["p2"]
    original_hash = ordered_feature_sha256(["p1", "p2"])
    sliced_hash = ordered_feature_sha256(["p2"])
    assert original_hash != sliced_hash
    expected_matrix_totals = (4, 5)
    expected_source_rows = (6, 5)
    for adata, matrix_total, source_rows in zip(
        (data.reference, data.spatial),
        expected_matrix_totals,
        expected_source_rows,
    ):
        assert adata.uns["fragment_shape"]["feature_sha256"] == sliced_hash
        metadata = adata.uns["fragment_shape"]
        assert metadata["preprocessing_counters"]["total_rows"] == source_rows
        assert metadata["preprocessing_counters"]["assigned_cut_sites"] == 10
        assert metadata["matrix_counters"]["assigned_cut_sites"] == matrix_total
        for layer_name in data.fragment_shape.layer_names:
            assert metadata["matrix_counters"][f"cut_sites_per_bin.{layer_name}"] == int(
                adata.layers[layer_name].sum()
            )
        layer_sum = sum(
            (adata.layers[layer] for layer in data.fragment_shape.layer_names),
            sparse.csr_matrix(adata.shape),
        )
        assert (adata.X != layer_sum).nnz == 0


def test_loader_rejects_a_stale_original_hash_before_slicing(tmp_path: Path) -> None:
    registry_path = _write_shape_dataset(tmp_path)
    reference_path = tmp_path / "reference.h5ad"
    reference = ad.read_h5ad(reference_path)
    reference.uns["fragment_shape"]["feature_sha256"] = "f" * 64
    reference.write_h5ad(reference_path)

    with pytest.raises(ValueError, match="original reference.*current ordered feature axis"):
        load_deconvolution_input(
            dataset_id="shape_toy",
            modality="atac",
            feature_set="single_peak",
            registry_path=registry_path,
        )
