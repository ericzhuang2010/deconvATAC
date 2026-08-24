#!/usr/bin/env python
"""Convert the 10x PBMC multiome matrix into labeled ATAC/RNA references."""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path

import pandas as pd
import scanpy as sc
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATRIX = (
    ROOT
    / "data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/"
    / "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
)
DEFAULT_QC = (
    ROOT
    / "data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/"
    / "pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv"
)
DEFAULT_MAPPING = ROOT / "data/raw/sources/snapatac2/pbmc10k_multiome/cell_type_mapping.csv"
DEFAULT_OUTPUT_DIR = ROOT / "data/processed/references/pbmc_granulocyte_sorted_10k_multiome"


def parse_peak_coordinates(var: pd.DataFrame) -> pd.DataFrame:
    peak_names = var["gene_ids"].fillna(var.index.to_series()).astype(str)
    coords = peak_names.str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    if coords.isna().any().any():
        bad = peak_names[coords.isna().any(axis=1)].head().tolist()
        raise ValueError(f"Could not parse peak coordinates for examples: {bad}")
    coords["start"] = coords["start"].astype(int)
    coords["end"] = coords["end"].astype(int)
    return coords


def add_common_obs(adata, mapping: pd.DataFrame):
    required_columns = ["cell_type", "cell_type_source"]
    missing = [column for column in required_columns if column not in mapping.columns]
    if missing:
        raise KeyError(f"cell_type_mapping.csv is missing required columns: {missing}")
    for column in required_columns:
        adata.obs[column] = mapping.loc[adata.obs_names, column].values
    for column in [
        "cell_type_confidence",
        "cell_type_broad",
        "cell_type_fine",
        "cell_type_broad_confidence",
        "cell_type_fine_confidence",
    ]:
        if column in mapping.columns:
            adata.obs[column] = mapping.loc[adata.obs_names, column].values
    adata.obs["donor_id"] = "healthy_donor_10x_pbmc_granulocyte_sorted_10k"
    adata.obs["organism"] = "human"
    adata.obs["tissue"] = "peripheral_blood_mononuclear_cells"
    adata.obs["assay"] = "10x_multiome_cellranger_arc_2.0.0"
    adata.obs["source_dataset_id"] = "pbmc_granulocyte_sorted_10k"


def attach_qc_metrics(adata, qc_path: Path):
    if not qc_path.exists():
        return
    qc = pd.read_csv(qc_path)
    if "barcode" not in qc.columns:
        raise KeyError(f"{qc_path} does not contain a barcode column.")
    qc = qc.set_index("barcode")
    keep_columns = [
        "is_cell",
        "gex_raw_reads",
        "gex_mapped_reads",
        "gex_umis_count",
        "gex_genes_count",
        "atac_raw_reads",
        "atac_fragments",
        "atac_TSS_fragments",
        "atac_peak_region_fragments",
        "atac_peak_region_cutsites",
    ]
    available = [column for column in keep_columns if column in qc.columns]
    for column in available:
        adata.obs[f"qc_{column}"] = qc.loc[adata.obs_names, column].values


def write_reference_manifest(output_dir: Path, mapping_path: Path, n_input_cells: int, n_labeled_cells: int, n_atac_features: int, n_rna_features: int):
    manifest = {
        "reference_id": "pbmc_granulocyte_sorted_10k_multiome",
        "description": "10x PBMC granulocyte-sorted 10k Multiome references labeled with SnapATAC2/GET prepare_pbmc annotations.",
        "labels_key": "cell_type",
        "prepared_at": date.today().isoformat(),
        "source_dataset_id": "pbmc_granulocyte_sorted_10k",
        "label_source": str(mapping_path.relative_to(ROOT)),
        "fallback_policy": "Do not use CellTypist or any other fallback label source without explicit user approval.",
        "counts": {
            "input_cells": int(n_input_cells),
            "labeled_cells": int(n_labeled_cells),
            "dropped_cells_missing_label": int(n_input_cells - n_labeled_cells),
            "atac_features": int(n_atac_features),
            "rna_features": int(n_rna_features),
        },
        "modalities": {
            "atac": {
                "path": "data/processed/references/pbmc_granulocyte_sorted_10k_multiome/atac/reference.h5ad",
                "feature_type": "Peaks",
                "source_filename": "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
            },
            "rna": {
                "path": "data/processed/references/pbmc_granulocyte_sorted_10k_multiome/rna/reference.h5ad",
                "feature_type": "Gene Expression",
                "source_filename": "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
            },
        },
    }
    with (output_dir / "reference.yaml").open("w") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    parser.add_argument("--qc", type=Path, default=DEFAULT_QC)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    atac_path = args.output_dir / "atac" / "reference.h5ad"
    rna_path = args.output_dir / "rna" / "reference.h5ad"
    if not args.overwrite and (atac_path.exists() or rna_path.exists() or (args.output_dir / "reference.yaml").exists()):
        raise FileExistsError(f"Refusing to overwrite existing outputs in {args.output_dir}; pass --overwrite.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    atac_path.parent.mkdir(parents=True, exist_ok=True)
    rna_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_10x_h5(args.matrix, gex_only=False)
    mapping = pd.read_csv(args.mapping).set_index("barcode")
    common_mask = adata.obs_names.isin(mapping.index)
    labeled = adata[common_mask].copy()
    mapping = mapping.loc[labeled.obs_names]

    if mapping["cell_type"].isna().any():
        raise ValueError("cell_type_mapping.csv contains missing cell_type values.")

    add_common_obs(labeled, mapping)
    attach_qc_metrics(labeled, args.qc)

    atac = labeled[:, labeled.var["feature_types"] == "Peaks"].copy()
    rna = labeled[:, labeled.var["feature_types"] == "Gene Expression"].copy()

    atac.var["feature_type"] = "Peaks"
    atac.var[["chrom", "start", "end"]] = parse_peak_coordinates(atac.var)
    atac.var_names_make_unique()

    rna.var["feature_type"] = "Gene Expression"
    rna.var["gene_id"] = rna.var["gene_ids"].astype(str)
    rna.var["gene_name"] = rna.var_names.astype(str)
    rna.var_names_make_unique()

    atac.write_h5ad(atac_path)
    rna.write_h5ad(rna_path)
    write_reference_manifest(args.output_dir, args.mapping, adata.n_obs, labeled.n_obs, atac.n_vars, rna.n_vars)

    print(f"wrote {atac_path}")
    print(f"wrote {rna_path}")
    print(f"wrote {args.output_dir / 'reference.yaml'}")
    print(f"input_cells {adata.n_obs}")
    print(f"labeled_cells {labeled.n_obs}")
    print(f"atac_features {atac.n_vars}")
    print(f"rna_features {rna.n_vars}")


if __name__ == "__main__":
    main()
