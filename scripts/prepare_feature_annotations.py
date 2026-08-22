#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Optional

import anndata as ad
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data.registry import get_dataset_config, resolve_project_path
from deconvatac.pp import highly_accessible_peaks, highly_variable_peaks


def read_adata(spec: dict[str, Any], modality: str) -> ad.AnnData:
    path = resolve_project_path(spec["path"], project_root=ROOT)
    selected_modality = spec.get("modality", modality)
    if path.suffix == ".h5mu":
        import muon as mu

        mdata = mu.read_h5mu(path)
        if selected_modality not in mdata.mod:
            raise KeyError(f"Modality '{selected_modality}' is missing from {path}.")
        return mdata.mod[selected_modality].copy()
    if path.suffix == ".h5ad":
        return ad.read_h5ad(path)
    raise ValueError(f"Unsupported input file type: {path}")


def relpath(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def copy_feature_mask(reference: ad.AnnData, spatial: ad.AnnData, column: str) -> None:
    selected = reference.var_names[reference.var[column].astype(bool).values]
    spatial.var[column] = spatial.var_names.isin(selected)


def compute_feature_columns(
    reference: ad.AnnData,
    modality: str,
    labels_key: str,
    n_top_features: int,
    layer: Optional[str],
) -> list[str]:
    columns = []

    highly_variable_peaks(
        adata=reference,
        cluster_key=labels_key,
        layer=layer,
        n_top_features=n_top_features,
    )
    columns.append("highly_variable")

    if modality == "atac":
        highly_accessible_peaks(
            adata=reference,
            layer=layer,
            n_top_features=n_top_features,
        )
        columns.append("highly_accessible")

    return columns


def build_processed_config(
    dataset_id: str,
    source_config: dict[str, Any],
    processed_modalities: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    return {
        "dataset_id": dataset_id,
        "source": "processed_feature_annotations",
        "description": source_config.get("description"),
        "labels_key": source_config.get("labels_key", "cell_type"),
        "spatial_key": source_config.get("spatial_key", "spatial"),
        "modalities": processed_modalities,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute feature-selection annotations and write processed inputs.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--registry", default=str(ROOT / "data" / "registry" / "datasets.yaml"))
    parser.add_argument("--modalities", nargs="+")
    parser.add_argument("--output-root", default=str(ROOT / "data" / "processed" / "datasets"))
    parser.add_argument("--n-top-features", type=int, default=20000)
    parser.add_argument("--layer")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    source_config = get_dataset_config(args.dataset, registry_path=args.registry, project_root=ROOT)
    available_modalities = source_config.get("modalities", {})
    modalities = args.modalities or list(available_modalities)
    output_dataset_dir = Path(args.output_root) / args.dataset
    output_dataset_dir.mkdir(parents=True, exist_ok=True)

    processed_modalities: dict[str, dict[str, Any]] = {}
    summary_rows = []

    for modality in modalities:
        if modality not in available_modalities:
            raise KeyError(f"Dataset '{args.dataset}' does not define modality '{modality}'.")

        modality_config = available_modalities[modality]
        labels_key = modality_config.get("labels_key", source_config.get("labels_key", "cell_type"))
        spatial_key = modality_config.get("spatial_key", source_config.get("spatial_key", "spatial"))

        reference = read_adata(modality_config["reference"], modality=modality)
        spatial = read_adata(modality_config["spatial"], modality=modality)

        columns = compute_feature_columns(
            reference=reference,
            modality=modality,
            labels_key=labels_key,
            n_top_features=args.n_top_features,
            layer=args.layer,
        )

        for column in columns:
            copy_feature_mask(reference, spatial, column)
            summary_rows.append(
                {
                    "dataset_id": args.dataset,
                    "modality": modality,
                    "feature_set": column,
                    "reference_selected": int(reference.var[column].sum()),
                    "spatial_selected": int(spatial.var[column].sum()),
                    "n_reference_features": int(reference.n_vars),
                    "n_spatial_features": int(spatial.n_vars),
                }
            )

        modality_dir = output_dataset_dir / modality
        modality_dir.mkdir(parents=True, exist_ok=True)
        reference_path = modality_dir / "reference.h5ad"
        spatial_path = modality_dir / "spatial.h5ad"
        if not args.overwrite and (reference_path.exists() or spatial_path.exists()):
            raise FileExistsError(f"Processed files already exist under {modality_dir}. Use --overwrite.")

        reference.write_h5ad(reference_path)
        spatial.write_h5ad(spatial_path)

        feature_sets = {column: {"var_column": column} for column in columns}
        feature_sets["all"] = {"mode": "all"}

        processed_modalities[modality] = {
            "reference": {"path": relpath(reference_path)},
            "spatial": {"path": relpath(spatial_path)},
            "labels_key": labels_key,
            "spatial_key": spatial_key,
            "truth": modality_config.get("truth"),
            "feature_sets": feature_sets,
        }

    processed_config = build_processed_config(args.dataset, source_config, processed_modalities)
    config_path = output_dataset_dir / "dataset.yaml"
    if config_path.exists() and not args.overwrite:
        raise FileExistsError(f"{config_path} already exists. Use --overwrite.")
    with config_path.open("w") as handle:
        yaml.safe_dump(processed_config, handle, sort_keys=False)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output_dataset_dir / "feature_annotation_summary.csv", index=False)
    print(config_path)


if __name__ == "__main__":
    main()
