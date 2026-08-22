#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Optional

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data.registry import get_dataset_config, resolve_project_path


def read_reference(spec: dict[str, Any], modality: str) -> ad.AnnData:
    path = resolve_project_path(spec["path"], project_root=ROOT)
    selected_modality = spec.get("modality", modality)
    if path.suffix == ".h5mu":
        import muon as mu

        mdata = mu.read_h5mu(path, backed="r")
        if selected_modality not in mdata.mod:
            raise KeyError(f"Modality '{selected_modality}' is missing from {path}.")
        return mdata.mod[selected_modality]
    if path.suffix == ".h5ad":
        return ad.read_h5ad(path, backed="r")
    raise ValueError(f"Unsupported reference file type: {path}")


def matrix_from_layer(adata: ad.AnnData, layer: Optional[str]):
    if layer is None:
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"layer '{layer}' is missing from reference.layers.")
    return adata.layers[layer]


def as_1d(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def top_feature_names(var_names: pd.Index, scores: np.ndarray, n_top_features: int) -> pd.Index:
    n_top = min(int(n_top_features), scores.shape[0])
    if n_top <= 0:
        raise ValueError("n_top_features must be positive.")
    indices = np.argpartition(scores, -n_top)[-n_top:]
    ordered = indices[np.argsort(scores[indices])[::-1]]
    return var_names[ordered]


def write_feature_list(path: Path, features: pd.Index) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(map(str, features)) + "\n")


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def compute_scores(
    reference: ad.AnnData,
    labels_key: str,
    modality: str,
    n_top_features: int,
    layer: Optional[str],
    chunk_size: int,
) -> dict[str, pd.Index]:
    if labels_key not in reference.obs:
        raise KeyError(f"labels_key '{labels_key}' is missing from reference.obs.")

    matrix = matrix_from_layer(reference, layer)
    labels = reference.obs[labels_key].astype(str).to_numpy()
    clusters = pd.Index(pd.unique(labels))
    cluster_to_index = {cluster: idx for idx, cluster in enumerate(clusters)}

    group_sums = np.zeros((len(clusters), reference.n_vars), dtype=np.float64)
    accessibility_counts = np.zeros(reference.n_vars, dtype=np.float64) if modality == "atac" else None

    for start in range(0, reference.n_obs, chunk_size):
        end = min(start + chunk_size, reference.n_obs)
        chunk = matrix[start:end]
        chunk_labels = labels[start:end]

        if accessibility_counts is not None:
            if sparse.issparse(chunk):
                accessibility_counts += np.asarray(chunk.getnnz(axis=0)).ravel()
            else:
                accessibility_counts += np.asarray(chunk > 0).sum(axis=0).ravel()

        for cluster in pd.unique(chunk_labels):
            row_mask = chunk_labels == cluster
            group_sums[cluster_to_index[cluster]] += as_1d(chunk[row_mask].sum(axis=0))

        print(f"processed rows {end}/{reference.n_obs}", flush=True)

    row_sums = group_sums.sum(axis=1, keepdims=True)
    normalized = np.divide(group_sums, row_sums, out=np.zeros_like(group_sums), where=row_sums != 0)
    normalized = np.log2(normalized + 1)
    variances = np.var(normalized, axis=0)

    feature_sets = {
        "highly_variable": top_feature_names(reference.var_names, variances, n_top_features),
    }
    if accessibility_counts is not None:
        feature_sets["highly_accessible"] = top_feature_names(
            reference.var_names,
            accessibility_counts,
            n_top_features,
        )
    return feature_sets


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute reusable feature-set lists from shared references.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--registry", default=str(ROOT / "data" / "registry" / "datasets.yaml"))
    parser.add_argument("--modalities", nargs="+")
    parser.add_argument("--feature-set-id", default=None)
    parser.add_argument("--output-root", default=str(ROOT / "data" / "processed" / "feature_sets"))
    parser.add_argument("--n-top-features", type=int, default=20000)
    parser.add_argument("--layer")
    parser.add_argument("--chunk-size", type=int, default=2048)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    config = get_dataset_config(args.dataset, registry_path=args.registry, project_root=ROOT)
    modalities = args.modalities or list(config.get("modalities", {}))
    feature_set_id = args.feature_set_id or args.dataset
    output_root = Path(args.output_root) / feature_set_id
    summary_rows = []

    for modality in modalities:
        modality_config = config["modalities"][modality]
        labels_key = modality_config.get("labels_key", config.get("labels_key", "cell_type"))
        reference = read_reference(modality_config["reference"], modality)

        feature_sets = compute_scores(
            reference=reference,
            labels_key=labels_key,
            modality=modality,
            n_top_features=args.n_top_features,
            layer=args.layer,
            chunk_size=args.chunk_size,
        )

        for feature_set, features in feature_sets.items():
            output_path = output_root / modality / f"{feature_set}.txt"
            if output_path.exists() and not args.overwrite:
                raise FileExistsError(f"{output_path} already exists. Use --overwrite.")
            write_feature_list(output_path, features)
            summary_rows.append(
                {
                    "feature_set_id": feature_set_id,
                    "dataset_source": args.dataset,
                    "modality": modality,
                    "feature_set": feature_set,
                    "n_features": len(features),
                    "path": display_path(output_path),
                }
            )

    summary = pd.DataFrame(summary_rows)
    summary_path = output_root / "feature_set_summary.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    if summary_path.exists():
        existing = pd.read_csv(summary_path)
        if not existing.empty:
            keys = ["feature_set_id", "modality", "feature_set"]
            incoming_keys = set(map(tuple, summary[keys].to_numpy()))
            keep_existing = [
                tuple(row) not in incoming_keys
                for row in existing[keys].to_numpy()
            ]
            summary = pd.concat([existing.loc[keep_existing], summary], ignore_index=True)
    summary.to_csv(summary_path, index=False)
    print(summary_path)


if __name__ == "__main__":
    main()
