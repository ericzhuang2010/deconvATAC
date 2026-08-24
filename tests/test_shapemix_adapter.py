import json
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse

from deconvatac.data import (
    DeconvolutionInput,
    DeconvolutionResult,
    FragmentShapeSpec,
    ordered_feature_sha256,
)
from deconvatac.methods.shapemix import ShapeMixDeconvolver


FRAGMENT_SHAPE = {
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

FAST_CONFIG = {
    "max_steps": 3,
    "patience": 2,
    "tolerance": 1.0e6,
    "restarts": 1,
    "spot_batch_size": 2,
    "peak_chunk_size": 2,
}


def _declared_spec() -> FragmentShapeSpec:
    return FragmentShapeSpec.from_mapping(FRAGMENT_SHAPE)


def _resolved_spec(
    feature_names: pd.Index,
    layer_totals: dict[str, int],
) -> FragmentShapeSpec:
    assigned = sum(layer_totals.values())
    return replace(
        _declared_spec(),
        left_cut_offset=0,
        right_cut_offset=0,
        source_sha256={"fragments": "a" * 64, "tabix_index": "b" * 64},
        feature_sha256=ordered_feature_sha256(feature_names),
        split_sha256="c" * 64,
        coordinate_validation={
            "selected_right_cut_offset": 0,
            "mismatched_entries": 0,
            "absolute_error": 0,
            "matrix_match": "exact",
        },
        software_versions={"deconvatac": "0.0.1", "pysam": "0.24.0"},
        preprocessing_counters={
            "total_rows": assigned,
            "header_rows": 0,
            "invalid_schema_rows": 0,
            "invalid_coordinate_rows": 0,
            "unknown_barcodes": 0,
            "filtered_contigs": 0,
            "valid_rows": assigned,
            "retained_fragments": assigned,
            "fragments_with_assigned_cut_sites": assigned,
            "cut_sites_outside_peaks": assigned,
            "assigned_cut_sites": assigned,
            "read_support_total": assigned,
            "cut_sites_per_bin": layer_totals,
        },
        matrix_counters={
            "assigned_cut_sites": assigned,
            **{
                f"cut_sites_per_bin.{layer_name}": count
                for layer_name, count in layer_totals.items()
            },
        },
    )


def _add_shape_layers(adata: ad.AnnData, values: tuple[np.ndarray, ...]) -> None:
    for layer_name, layer_values in zip(_declared_spec().layer_names, values):
        adata.layers[layer_name] = sparse.csr_matrix(layer_values)
    adata.X = sum(
        (adata.layers[name] for name in _declared_spec().layer_names),
        sparse.csr_matrix(adata.shape, dtype=np.int64),
    )
    layer_totals = {
        name: int(adata.layers[name].sum()) for name in _declared_spec().layer_names
    }
    adata.uns["fragment_shape"] = _resolved_spec(
        adata.var_names, layer_totals
    ).to_uns()


def _shape_input(output_dir=None) -> DeconvolutionInput:
    reference = ad.AnnData(
        X=sparse.csr_matrix((4, 3), dtype=np.int64),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B"]},
            index=["cell_a1", "cell_a2", "cell_b1", "cell_b2"],
        ),
        var=pd.DataFrame(index=["peak_1", "peak_2", "peak_3"]),
    )
    _add_shape_layers(
        reference,
        (
            np.array([[4, 1, 0], [3, 0, 1], [0, 2, 3], [1, 1, 4]]),
            np.array([[1, 3, 1], [1, 2, 0], [2, 1, 1], [1, 2, 1]]),
            np.array([[1, 0, 1], [0, 1, 1], [1, 1, 2], [0, 1, 1]]),
        ),
    )

    spatial = ad.AnnData(
        X=sparse.csr_matrix((3, 3), dtype=np.int64),
        obs=pd.DataFrame(index=["spot_b", "spot_a", "spot_c"]),
        var=pd.DataFrame(index=["peak_1", "peak_2", "peak_3"]),
    )
    _add_shape_layers(
        spatial,
        (
            np.array([[7, 1, 1], [1, 3, 7], [4, 2, 4]]),
            np.array([[2, 5, 1], [3, 3, 2], [2, 4, 1]]),
            np.array([[1, 1, 2], [1, 2, 3], [1, 2, 2]]),
        ),
    )
    spatial.obsm["spatial"] = np.array(
        [[1.0, 0.0], [0.0, 0.0], [2.0, 0.0]], dtype=float
    )
    truth = pd.DataFrame(
        [[0.8, 0.2], [0.2, 0.8], [0.5, 0.5]],
        index=spatial.obs_names,
        columns=["A", "B"],
    )
    return DeconvolutionInput(
        dataset_id="shapemix_adapter_toy",
        modality="atac",
        feature_set="all",
        spatial=spatial,
        reference=reference,
        labels_key="cell_type",
        truth=truth,
        output_dir=output_dir,
        metadata={
            "dataset_config": {
                "dataset_id": "shapemix_adapter_toy",
                "simulation": {
                    "outer_split_seed": 0,
                    "inner_mixture_seed": 0,
                },
            }
        },
        fragment_shape=_declared_spec(),
        cell_types=["A", "B"],
    )


@pytest.mark.parametrize("use_shape", [False, True])
def test_adapter_runs_both_nested_arms_and_writes_compact_outputs(
    tmp_path, use_shape
):
    output_dir = tmp_path / ("shape" if use_shape else "count")
    data = _shape_input(output_dir=output_dir)
    result = ShapeMixDeconvolver(use_shape=use_shape, **FAST_CONFIG).run(data)

    assert isinstance(result, DeconvolutionResult)
    assert result.method == "shapemix"
    assert result.proportions.index.equals(data.spatial.obs_names)
    assert result.abundance.index.equals(data.spatial.obs_names)
    assert list(result.proportions.columns) == data.cell_types
    assert list(result.abundance.columns) == data.cell_types
    assert np.isfinite(result.proportions.to_numpy()).all()
    assert np.isfinite(result.abundance.to_numpy()).all()
    np.testing.assert_allclose(
        result.proportions.sum(axis=1).to_numpy(), 1.0, rtol=0.0, atol=1.0e-6
    )
    assert result.diagnostics["use_shape"] is use_shape
    assert result.diagnostics["fit"]["success"] is True
    assert "history" not in result.diagnostics["fit"]["restarts"][0]
    assert result.diagnostics["native_outputs"] is not None
    json.dumps(result.diagnostics, allow_nan=False)

    native_dir = output_dir / "results" / "raw_method_output"
    expected_files = {
        "training_history.csv",
        "restart_summary.csv",
        "reconstruction_summary.csv",
        "residuals_by_peak_and_bin.csv.gz",
        "signature_summary.yaml",
    }
    assert {path.name for path in native_dir.iterdir()} == expected_files
    history = pd.read_csv(native_dir / "training_history.csv")
    restarts = pd.read_csv(native_dir / "restart_summary.csv")
    reconstruction = pd.read_csv(native_dir / "reconstruction_summary.csv")
    residuals = pd.read_csv(native_dir / "residuals_by_peak_and_bin.csv.gz")
    with (native_dir / "signature_summary.yaml").open() as handle:
        signature_summary = yaml.safe_load(handle)
    assert history["step"].tolist() == [0, 1]
    assert restarts.loc[0, "selected"]
    assert reconstruction["component"].tolist() == [
        "total_count",
        "shape_bin",
        "shape_bin",
        "shape_bin",
    ]
    assert residuals.shape[0] == data.spatial.n_vars * 3
    assert residuals["peak_id"].tolist() == [
        peak
        for peak in data.spatial.var_names
        for _ in range(3)
    ]
    assert signature_summary["source"] == "training_reference_only"
    assert signature_summary["outer_split_seed"] == 0
    assert signature_summary["signature"]["cell_types"] == ["A", "B"]
    assert not any(
        path.suffix in {".npy", ".npz", ".pt", ".pkl", ".h5", ".h5ad"}
        for path in output_dir.rglob("*")
        if path.is_file()
    )

    result.write(output_dir)
    assert (output_dir / "results" / "proportions.csv").is_file()
    assert (output_dir / "results" / "abundance.csv").is_file()
    assert (output_dir / "results" / "diagnostics.json").is_file()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("missing_declaration", "requires declared fragment_shape"),
        ("missing_object_metadata", "missing 'fragment_shape' metadata"),
        ("missing_layer", "missing declared layer"),
        ("rna", "only ATAC"),
        ("missing_cell_types", "ordered cell_types"),
        ("missing_seed", "inner_mixture_seed must be an integer"),
        ("axis_mismatch", "features must be aligned"),
    ],
)
def test_adapter_rejects_invalid_contracts(mutation, message):
    data = _shape_input()
    if mutation == "missing_declaration":
        data.fragment_shape = None
    elif mutation == "missing_object_metadata":
        del data.spatial.uns["fragment_shape"]
    elif mutation == "missing_layer":
        del data.spatial.layers[_declared_spec().layer_names[0]]
    elif mutation == "rna":
        data.modality = "rna"
    elif mutation == "missing_cell_types":
        data.cell_types = None
    elif mutation == "missing_seed":
        del data.metadata["dataset_config"]["simulation"]["inner_mixture_seed"]
    elif mutation == "axis_mismatch":
        data.spatial = data.spatial[:, ["peak_2", "peak_1", "peak_3"]].copy()
    else:
        raise AssertionError(mutation)

    with pytest.raises((TypeError, ValueError), match=message):
        ShapeMixDeconvolver(**FAST_CONFIG).run(data)


@pytest.mark.parametrize(("value", "message"), [(np.nan, "non-finite"), (1.5, "non-integer")])
def test_adapter_rejects_invalid_count_values(value, message):
    data = _shape_input()
    layer_name = _declared_spec().layer_names[0]
    layer = data.reference.layers[layer_name].astype(float).copy()
    layer.data[0] = value
    data.reference.layers[layer_name] = layer

    with pytest.raises(ValueError, match=message):
        ShapeMixDeconvolver(**FAST_CONFIG).run(data)


def test_adapter_is_deterministic_and_output_dir_is_optional():
    adapter = ShapeMixDeconvolver(use_shape=True, **FAST_CONFIG)
    first = adapter.run(_shape_input())
    second = adapter.run(_shape_input())

    np.testing.assert_array_equal(first.abundance.to_numpy(), second.abundance.to_numpy())
    np.testing.assert_array_equal(
        first.proportions.to_numpy(), second.proportions.to_numpy()
    )
    assert first.diagnostics["signature"]["content_sha256"] == second.diagnostics[
        "signature"
    ]["content_sha256"]
    assert first.diagnostics["native_outputs"] is None
    assert second.diagnostics["native_outputs"] is None


def test_adapter_accepts_explicit_external_dataset_seeds():
    data = _shape_input()
    dataset_config = data.metadata["dataset_config"]
    del dataset_config["simulation"]
    dataset_config["shapemix_seeds"] = {
        "outer_split_seed": 17,
        "inner_mixture_seed": 23,
    }

    result = ShapeMixDeconvolver(use_shape=True, **FAST_CONFIG).run(data)

    assert result.diagnostics["seeds"]["outer_split_seed"] == 17
    assert result.diagnostics["seeds"]["inner_mixture_seed"] == 23


def test_adapter_configuration_is_strict():
    with pytest.raises(ValueError, match="Unknown ShapeMix parameter"):
        ShapeMixDeconvolver(use_shapes=True)
