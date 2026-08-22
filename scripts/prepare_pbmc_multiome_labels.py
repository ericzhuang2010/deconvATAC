#!/usr/bin/env python
"""Prepare PBMC multiome barcode-to-cell-type labels."""

from __future__ import annotations

import argparse
import hashlib
import os
from datetime import date
from pathlib import Path

import h5py
import pandas as pd
import scanpy as sc
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATRIX = (
    ROOT
    / "data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/"
    / "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
)
DEFAULT_RNA_H5AD = ROOT / "data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad"
DEFAULT_SNAPATAC2_OUTPUT_DIR = ROOT / "data/raw/sources/snapatac2/pbmc10k_multiome"
DEFAULT_CELLTYPIST_OUTPUT_DIR = ROOT / "data/raw/sources/celltypist/pbmc_granulocyte_sorted_10k"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_rna_counts(matrix_path: Path):
    adata = sc.read_10x_h5(matrix_path, gex_only=False)
    rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
    rna.var_names_make_unique()
    return rna


def read_10x_barcodes(matrix_path: Path) -> list[str]:
    with h5py.File(matrix_path, "r") as handle:
        return decode_h5_strings(handle["matrix/barcodes"][:])


def decode_h5_strings(values) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def read_snapatac2_mapping(rna_h5ad: Path, label_key: str) -> tuple[pd.DataFrame, dict[str, int]]:
    with h5py.File(rna_h5ad, "r") as handle:
        obs = handle["obs"]
        index_key = obs.attrs.get("_index", "cells")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        if index_key not in obs:
            raise KeyError(f"{rna_h5ad} obs does not contain index field {index_key!r}.")
        if label_key not in obs:
            raise KeyError(f"{rna_h5ad} obs does not contain label field {label_key!r}.")

        barcodes = decode_h5_strings(obs[index_key][:])
        label_obj = obs[label_key]
        if isinstance(label_obj, h5py.Group) and {"categories", "codes"}.issubset(label_obj.keys()):
            categories = decode_h5_strings(label_obj["categories"][:])
            codes = label_obj["codes"][:]
            labels = [categories[int(code)] if int(code) >= 0 else pd.NA for code in codes]
        else:
            labels = decode_h5_strings(label_obj[:])

    mapping = pd.DataFrame({"barcode": barcodes, "cell_type": labels})
    total_rows = int(mapping.shape[0])
    mapping = mapping.dropna(subset=["cell_type"]).copy()
    mapping = mapping[mapping["cell_type"].astype(str).str.len() > 0].copy()
    mapping["cell_type"] = mapping["cell_type"].astype(str)
    mapping["cell_type_source"] = "snapatac2_pbmc10k_multiome_get_prepare_pbmc"
    counts = {
        "rna_h5ad_rows": total_rows,
        "labeled_rows": int(mapping.shape[0]),
        "missing_cell_type_rows": int(total_rows - mapping.shape[0]),
        "cell_type_categories": int(mapping["cell_type"].nunique()),
    }
    return mapping, counts


def normalize_for_celltypist(rna):
    rna = rna.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    return rna


def annotate_with_celltypist(rna, output_dir: Path, broad_model: str, fine_model: str) -> pd.DataFrame:
    os.environ.setdefault("CELLTYPIST_FOLDER", str(output_dir))

    import celltypist
    from celltypist import models

    model_files = [Path(models.models_path) / broad_model, Path(models.models_path) / fine_model]
    if any(not path.exists() for path in model_files):
        models.download_models(model=[broad_model, fine_model])

    broad = celltypist.annotate(rna, model=broad_model, mode="best match", majority_voting=False)
    fine = celltypist.annotate(rna, model=fine_model, mode="best match", majority_voting=False)

    mapping = pd.DataFrame(index=rna.obs_names)
    mapping.index.name = "barcode"
    mapping["cell_type_broad"] = broad.predicted_labels["predicted_labels"].astype(str)
    mapping["cell_type_fine"] = fine.predicted_labels["predicted_labels"].astype(str)
    mapping["cell_type"] = mapping["cell_type_broad"]
    mapping["cell_type_source"] = "celltypist_immune_all_high_low"
    mapping["cell_type_confidence"] = broad.probability_matrix.max(axis=1).astype(float)
    mapping["cell_type_broad_confidence"] = broad.probability_matrix.max(axis=1).astype(float)
    mapping["cell_type_fine_confidence"] = fine.probability_matrix.max(axis=1).astype(float)
    return mapping.reset_index()


def write_snapatac2_manifest(
    output_dir: Path,
    mapping_path: Path,
    summary_path: Path,
    rna_h5ad: Path,
    matrix_path: Path,
    counts: dict[str, int],
) -> None:
    matrix_barcodes = pd.Index(read_10x_barcodes(matrix_path))
    mapping_barcodes = pd.Index(pd.read_csv(mapping_path, usecols=["barcode"])["barcode"].astype(str))
    retained_in_matrix = int(matrix_barcodes.isin(mapping_barcodes).sum())
    labels_not_in_matrix = int((~mapping_barcodes.isin(matrix_barcodes)).sum())
    manifest = {
        "source": "snapatac2_pbmc10k_multiome_get_prepare_pbmc",
        "source_file": str(rna_h5ad.relative_to(ROOT)),
        "prepared_at": date.today().isoformat(),
        "role": "canonical_pbmc_barcode_to_cell_type_label_source",
        "fallback_policy": "Do not use CellTypist or any other fallback label source without explicit user approval.",
        "input_files": [
            {
                "path": str(matrix_path.relative_to(ROOT)),
                "role": "10x_filtered_feature_barcode_matrix",
                "barcodes": int(matrix_barcodes.shape[0]),
                "sha256": sha256_file(matrix_path),
            },
            {
                "path": str(rna_h5ad.relative_to(ROOT)),
                "role": "rna_h5ad_with_cell_type_labels",
                "bytes": rna_h5ad.stat().st_size,
                "sha256": sha256_file(rna_h5ad),
            },
        ],
        "output_files": [
            {
                "path": str(mapping_path.relative_to(ROOT)),
                "role": "barcode_to_cell_type_mapping",
                "bytes": mapping_path.stat().st_size,
                "sha256": sha256_file(mapping_path),
            },
            {
                "path": str(summary_path.relative_to(ROOT)),
                "role": "cell_type_count_summary",
                "bytes": summary_path.stat().st_size,
                "sha256": sha256_file(summary_path),
            },
        ],
        "obs": {
            "barcode_field": "cells",
            "label_field": "cell_type",
            **counts,
        },
        "original_10x_matrix": {
            "cells": int(matrix_barcodes.shape[0]),
            "retained_labeled_cells_expected": retained_in_matrix,
            "dropped_unmatched_or_unlabeled_cells_expected": int(matrix_barcodes.shape[0] - retained_in_matrix),
            "labels_not_found_in_10x_matrix": labels_not_in_matrix,
        },
        "package_versions": {
            "h5py": h5py.__version__,
            "pandas": pd.__version__,
        },
        "references": [
            "https://github.com/GET-Foundation/get_model/blob/master/tutorials/prepare_pbmc.ipynb",
            "https://github.com/scverse/SnapATAC2/blob/main/snapatac2/datasets.py",
        ],
    }
    with (output_dir / "manifest.yaml").open("w") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False)


def write_celltypist_manifest(output_dir: Path, mapping_path: Path, summary_path: Path, matrix_path: Path, broad_model: str, fine_model: str):
    import celltypist

    manifest = {
        "source": "celltypist",
        "source_reason": (
            "Fallback label source used because the SnapATAC2 PBMC10k Multiome RNA h5ad URL referenced "
            "by the GET prepare_pbmc.ipynb workflow was inaccessible from this environment."
        ),
        "prepared_at": date.today().isoformat(),
        "input_matrix": str(matrix_path.relative_to(ROOT)),
        "output_files": [
            {
                "path": str(mapping_path.relative_to(ROOT)),
                "role": "barcode_to_cell_type_mapping",
                "bytes": mapping_path.stat().st_size,
                "sha256": sha256_file(mapping_path),
            },
            {
                "path": str(summary_path.relative_to(ROOT)),
                "role": "cell_type_count_summary",
                "bytes": summary_path.stat().st_size,
                "sha256": sha256_file(summary_path),
            },
        ],
        "models": [
            {"name": broad_model, "role": "cell_type_broad"},
            {"name": fine_model, "role": "cell_type_fine"},
        ],
        "package_versions": {
            "celltypist": celltypist.__version__,
            "scanpy": sc.__version__,
            "pandas": pd.__version__,
        },
        "references": [
            "https://github.com/GET-Foundation/get_model/blob/master/tutorials/prepare_pbmc.ipynb",
            "https://github.com/scverse/SnapATAC2/blob/main/snapatac2/datasets.py",
            "https://www.celltypist.org/",
        ],
    }
    with (output_dir / "manifest.yaml").open("w") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--rna-h5ad", type=Path, default=DEFAULT_RNA_H5AD)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--broad-model", default="Immune_All_High.pkl")
    parser.add_argument("--fine-model", default="Immune_All_Low.pkl")
    parser.add_argument(
        "--allow-celltypist-fallback",
        action="store_true",
        help="Required explicit acknowledgement before using CellTypist fallback labels.",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.allow_celltypist_fallback:
        output_dir = args.output_dir or DEFAULT_CELLTYPIST_OUTPUT_DIR
    else:
        output_dir = args.output_dir or DEFAULT_SNAPATAC2_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    mapping_path = output_dir / "cell_type_mapping.csv"
    summary_path = output_dir / "cell_type_summary.csv"
    if not args.overwrite and (mapping_path.exists() or summary_path.exists()):
        raise FileExistsError(f"Refusing to overwrite existing outputs in {output_dir}; pass --overwrite.")

    if args.allow_celltypist_fallback:
        os.environ.setdefault("CELLTYPIST_FOLDER", str(output_dir))
        rna_counts = read_rna_counts(args.matrix)
        rna = normalize_for_celltypist(rna_counts)
        mapping = annotate_with_celltypist(rna, output_dir, args.broad_model, args.fine_model)
        mapping.to_csv(mapping_path, index=False)
        summary = (
            mapping.groupby(["cell_type", "cell_type_fine"], observed=False)
            .size()
            .reset_index(name="n_cells")
            .sort_values(["cell_type", "cell_type_fine"])
        )
        summary.to_csv(summary_path, index=False)
        write_celltypist_manifest(output_dir, mapping_path, summary_path, args.matrix, args.broad_model, args.fine_model)
    else:
        mapping, counts = read_snapatac2_mapping(args.rna_h5ad, label_key="cell_type")
        mapping.to_csv(mapping_path, index=False)
        summary = (
            mapping.groupby(["cell_type"], observed=False)
            .size()
            .reset_index(name="n_cells")
            .sort_values(["n_cells", "cell_type"], ascending=[False, True])
        )
        summary.to_csv(summary_path, index=False)
        write_snapatac2_manifest(output_dir, mapping_path, summary_path, args.rna_h5ad, args.matrix, counts)

    print(f"wrote {mapping_path}")
    print(f"wrote {summary_path}")
    print(f"labeled_cells {mapping.shape[0]}")


if __name__ == "__main__":
    main()
