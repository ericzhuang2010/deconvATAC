from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import anndata as ad
import pandas as pd

from .registry import get_dataset_config, resolve_project_path
from .schemas import DeconvolutionInput, FragmentShapeSpec
from .validators import (
    ordered_feature_sha256,
    validate_deconvolution_input,
    validate_fragment_shape_feature_axis,
)


def _slice_adata_features(
    adata: ad.AnnData,
    selected_features: pd.Index,
    *,
    validate_fragment_shape: bool,
    role: str,
) -> ad.AnnData:
    """Validate the source axis, slice it, and bind metadata to the new axis."""
    if validate_fragment_shape:
        validate_fragment_shape_feature_axis(adata, f"original {role}")
    selected = adata[:, selected_features].to_memory()
    if validate_fragment_shape:
        metadata = dict(selected.uns["fragment_shape"])
        metadata["feature_sha256"] = ordered_feature_sha256(selected.var_names)
        stored_spec = FragmentShapeSpec.from_mapping(metadata)
        layer_totals = {
            layer_name: int(selected.layers[layer_name].sum())
            for layer_name in stored_spec.layer_names
        }
        metadata["matrix_counters"] = {
            "assigned_cut_sites": sum(layer_totals.values()),
            **{
                f"cut_sites_per_bin.{layer_name}": count
                for layer_name, count in layer_totals.items()
            },
        }
        selected.uns["fragment_shape"] = metadata
    return selected


def _read_adata(
    spec: dict[str, Any],
    modality: str,
    project_root: Optional[Union[str, Path]],
    selected_features: Optional[pd.Index] = None,
    validate_fragment_shape: bool = False,
    role: str = "input",
) -> ad.AnnData:
    path = resolve_project_path(spec["path"], project_root=project_root)
    selected_modality = spec.get("modality", modality)

    if path.suffix == ".h5mu":
        import muon as mu

        mdata = mu.read_h5mu(path, backed="r" if selected_features is not None else None)
        if selected_modality not in mdata.mod:
            raise KeyError(f"Modality '{selected_modality}' is missing from {path}.")
        adata = mdata.mod[selected_modality]
        if selected_features is not None:
            return _slice_adata_features(
                adata,
                selected_features,
                validate_fragment_shape=validate_fragment_shape,
                role=role,
            )
        return adata.copy()

    if path.suffix == ".h5ad":
        if selected_features is None:
            return ad.read_h5ad(path)
        adata = ad.read_h5ad(path, backed="r")
        return _slice_adata_features(
            adata,
            selected_features,
            validate_fragment_shape=validate_fragment_shape,
            role=role,
        )

    raise ValueError(f"Unsupported input file type for {path}. Expected .h5ad or .h5mu.")


def _read_var_names_and_var(
    spec: dict[str, Any],
    modality: str,
    project_root: Optional[Union[str, Path]],
) -> tuple[pd.Index, pd.DataFrame]:
    path = resolve_project_path(spec["path"], project_root=project_root)
    selected_modality = spec.get("modality", modality)

    if path.suffix == ".h5mu":
        import muon as mu

        mdata = mu.read_h5mu(path, backed="r")
        if selected_modality not in mdata.mod:
            raise KeyError(f"Modality '{selected_modality}' is missing from {path}.")
        adata = mdata.mod[selected_modality]
        return adata.var_names.copy(), adata.var.copy()

    if path.suffix == ".h5ad":
        adata = ad.read_h5ad(path, backed="r")
        return adata.var_names.copy(), adata.var.copy()

    raise ValueError(f"Unsupported input file type for {path}. Expected .h5ad or .h5mu.")


def _features_from_file(path: Union[str, Path], project_root: Optional[Union[str, Path]]) -> list[str]:
    resolved = resolve_project_path(path, project_root=project_root)
    with resolved.open() as handle:
        return [line.strip() for line in handle if line.strip()]


def _resolve_selected_features(
    reference_spec: dict[str, Any],
    spatial_spec: dict[str, Any],
    modality: str,
    feature_set: str,
    feature_spec: Optional[dict[str, Any]],
    project_root: Optional[Union[str, Path]],
) -> Optional[pd.Index]:
    if feature_spec is None:
        if feature_set == "all":
            return None
        else:
            raise KeyError(f"Feature set '{feature_set}' is not configured for this dataset/modality.")
    elif feature_spec.get("mode") == "all":
        return None
    elif "features" in feature_spec:
        selected = list(feature_spec["features"])
    elif "path" in feature_spec:
        selected = _features_from_file(feature_spec["path"], project_root=project_root)
    elif "var_column" in feature_spec:
        _, reference_var = _read_var_names_and_var(reference_spec, modality=modality, project_root=project_root)
        reference_var_names = reference_var.index
        column = feature_spec["var_column"]
        if column in reference_var:
            selected = list(reference_var_names[reference_var[column].astype(bool).values])
        else:
            _, spatial_var = _read_var_names_and_var(spatial_spec, modality=modality, project_root=project_root)
            spatial_var_names = spatial_var.index
            if column in spatial_var:
                selected = list(spatial_var_names[spatial_var[column].astype(bool).values])
            else:
                raise KeyError(f"Feature column '{column}' is missing from both reference.var and spatial.var.")
    else:
        raise ValueError(f"Unsupported feature spec for feature set '{feature_set}': {feature_spec}")

    reference_var_names, _ = _read_var_names_and_var(reference_spec, modality=modality, project_root=project_root)
    spatial_var_names, _ = _read_var_names_and_var(spatial_spec, modality=modality, project_root=project_root)
    common = reference_var_names.intersection(spatial_var_names)
    selected_index = pd.Index(selected).intersection(common)
    if len(selected_index) == 0:
        raise ValueError(f"Feature set '{feature_set}' has no overlap between reference and spatial data.")

    return selected_index


def _truth_from_spatial(
    spatial: ad.AnnData,
    truth_spec: Optional[dict[str, Any]],
    project_root: Optional[Union[str, Path]],
) -> Optional[pd.DataFrame]:
    if not truth_spec:
        return None

    if "path" in truth_spec:
        path = resolve_project_path(truth_spec["path"], project_root=project_root)
        truth = pd.read_csv(path, index_col=0)
        truth.index = truth.index.astype(str)
        return truth

    obsm_key = truth_spec.get("obsm_key")
    if obsm_key is None:
        return None
    if obsm_key not in spatial.obsm:
        raise KeyError(f"Truth obsm key '{obsm_key}' is missing from spatial.obsm.")

    truth = spatial.obsm[obsm_key]
    if isinstance(truth, pd.DataFrame):
        return truth.copy()

    columns = truth_spec.get("columns") or truth_spec.get("cell_types")
    if columns is None and "proportion_names" in spatial.uns:
        columns = list(spatial.uns["proportion_names"])
    if columns is None:
        columns = [f"cell_type_{idx}" for idx in range(truth.shape[1])]

    return pd.DataFrame(truth, index=spatial.obs_names, columns=columns)


def _ordered_cell_types(truth_spec: Optional[dict[str, Any]]) -> Optional[list[str]]:
    """Read an optional, ordered truth cell-type universe."""
    if not truth_spec or "cell_types" not in truth_spec:
        return None

    cell_types = truth_spec["cell_types"]
    if not isinstance(cell_types, list):
        raise TypeError("truth.cell_types must be an ordered YAML list.")
    return cell_types.copy()


def _fragment_shape_spec(
    modality_config: dict[str, Any],
    config: dict[str, Any],
) -> Optional[FragmentShapeSpec]:
    """Read opt-in shape metadata, preferring modality-specific metadata."""
    if "fragment_shape" in modality_config:
        raw_spec = modality_config["fragment_shape"]
    else:
        raw_spec = config.get("fragment_shape")
    if raw_spec is None:
        return None
    return FragmentShapeSpec.from_mapping(raw_spec)


def load_deconvolution_input(
    dataset_id: str,
    modality: str,
    feature_set: str = "all",
    registry_path: Optional[Union[str, Path]] = None,
    project_root: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
) -> DeconvolutionInput:
    """Load and validate standardized deconvolution input data."""
    config = get_dataset_config(dataset_id, registry_path=registry_path, project_root=project_root)
    modality_configs = config.get("modalities", {})
    if modality not in modality_configs:
        available = ", ".join(sorted(modality_configs))
        raise KeyError(f"Dataset '{dataset_id}' does not define modality '{modality}'. Available: {available}")

    modality_config = modality_configs[modality]
    truth_spec = modality_config.get("truth") or config.get("truth")
    fragment_shape = _fragment_shape_spec(modality_config, config)
    cell_types = _ordered_cell_types(truth_spec)
    feature_sets = modality_config.get("feature_sets", {})
    selected_features = _resolve_selected_features(
        reference_spec=modality_config["reference"],
        spatial_spec=modality_config["spatial"],
        modality=modality,
        feature_set=feature_set,
        feature_spec=feature_sets.get(feature_set),
        project_root=project_root,
    )
    reference = _read_adata(
        modality_config["reference"],
        modality=modality,
        project_root=project_root,
        selected_features=selected_features,
        validate_fragment_shape=fragment_shape is not None,
        role="reference",
    )
    spatial = _read_adata(
        modality_config["spatial"],
        modality=modality,
        project_root=project_root,
        selected_features=selected_features,
        validate_fragment_shape=fragment_shape is not None,
        role="spatial",
    )

    truth = _truth_from_spatial(
        spatial=spatial,
        truth_spec=truth_spec,
        project_root=project_root,
    )

    data = DeconvolutionInput(
        dataset_id=dataset_id,
        modality=modality,
        feature_set=feature_set,
        spatial=spatial,
        reference=reference,
        labels_key=modality_config.get("labels_key", config.get("labels_key", "cell_type")),
        spatial_key=modality_config.get("spatial_key", config.get("spatial_key", "spatial")),
        truth=truth,
        output_dir=Path(output_dir) if output_dir is not None else None,
        metadata={"dataset_config": config, "modality_config": modality_config},
        fragment_shape=fragment_shape,
        cell_types=cell_types,
    )
    validate_deconvolution_input(data)
    return data
