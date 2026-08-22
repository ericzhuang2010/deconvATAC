#!/usr/bin/env python
from __future__ import annotations

import argparse
import gc
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


PBMC_DATASETS: dict[str, dict[str, Any]] = {
    "pbmc_granulocyte_sorted_10k_sim_equal_celltype": {
        "sampling_design": "equal_celltype",
        "description": "PBMC synthetic spatial spots with selected SnapATAC2/GET cell types sampled at equal probability.",
    },
    "pbmc_granulocyte_sorted_10k_sim_observed_abundance": {
        "sampling_design": "observed_abundance",
        "description": "PBMC synthetic spatial spots with selected SnapATAC2/GET cell types sampled by observed reference abundance.",
    },
}


def relpath(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_reference(path: Path) -> ad.AnnData:
    return ad.read_h5ad(path)


def validate_references(atac: ad.AnnData, rna: ad.AnnData, labels_key: str) -> None:
    if not atac.obs_names.equals(rna.obs_names):
        raise ValueError("PBMC ATAC and RNA references must have identical cells in identical order.")
    for modality, adata in {"atac": atac, "rna": rna}.items():
        if labels_key not in adata.obs:
            raise KeyError(f"{modality} reference is missing obs[{labels_key!r}].")
        if adata.obs[labels_key].isna().any():
            raise ValueError(f"{modality} reference has missing labels in obs[{labels_key!r}].")


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


def compute_feature_sets(
    reference: ad.AnnData,
    labels_key: str,
    selected_cell_types: list[str],
    modality: str,
    n_top_features: int,
    chunk_size: int,
) -> dict[str, pd.Index]:
    matrix = reference.X.tocsr() if sparse.issparse(reference.X) else np.asarray(reference.X)
    labels = reference.obs[labels_key].astype(str).to_numpy()
    selected = pd.Index(selected_cell_types)
    cluster_to_index = {cluster: idx for idx, cluster in enumerate(selected)}

    group_sums = np.zeros((len(selected), reference.n_vars), dtype=np.float64)
    accessibility_counts = np.zeros(reference.n_vars, dtype=np.float64) if modality == "atac" else None

    for start in range(0, reference.n_obs, chunk_size):
        end = min(start + chunk_size, reference.n_obs)
        chunk = matrix[start:end]
        chunk_labels = labels[start:end]
        selected_mask = np.isin(chunk_labels, selected)
        if not selected_mask.any():
            continue
        selected_chunk = chunk[selected_mask]
        selected_labels = chunk_labels[selected_mask]

        if accessibility_counts is not None:
            if sparse.issparse(selected_chunk):
                accessibility_counts += np.asarray(selected_chunk.getnnz(axis=0)).ravel()
            else:
                accessibility_counts += np.asarray(selected_chunk > 0).sum(axis=0).ravel()

        for cluster in pd.unique(selected_labels):
            row_mask = selected_labels == cluster
            group_sums[cluster_to_index[cluster]] += as_1d(selected_chunk[row_mask].sum(axis=0))

        print(f"{modality}: processed reference rows {end}/{reference.n_obs}", flush=True)

    row_sums = group_sums.sum(axis=1, keepdims=True)
    normalized = np.divide(group_sums, row_sums, out=np.zeros_like(group_sums), where=row_sums != 0)
    normalized = np.log2(normalized + 1)
    variances = np.var(normalized, axis=0)

    feature_sets = {"highly_variable": top_feature_names(reference.var_names, variances, n_top_features)}
    if accessibility_counts is not None:
        feature_sets["highly_accessible"] = top_feature_names(
            reference.var_names,
            accessibility_counts,
            n_top_features,
        )
    return feature_sets


def write_feature_lists(dataset_dir: Path, feature_sets: dict[str, dict[str, pd.Index]]) -> None:
    for modality, modality_feature_sets in feature_sets.items():
        feature_dir = dataset_dir / modality / "features"
        feature_dir.mkdir(parents=True, exist_ok=True)
        for name, features in modality_feature_sets.items():
            (feature_dir / f"{name}.txt").write_text("\n".join(map(str, features)) + "\n")


def sample_source_positions(
    rng: np.random.Generator,
    cell_types: list[str],
    cell_type_probabilities: np.ndarray,
    positions_by_cell_type: dict[str, np.ndarray],
    n_cells: int,
) -> tuple[np.ndarray, list[str]]:
    sampled_cell_types = rng.choice(cell_types, size=n_cells, replace=True, p=cell_type_probabilities)
    positions = [
        int(rng.choice(positions_by_cell_type[cell_type]))
        for cell_type in sampled_cell_types
    ]
    return np.asarray(positions, dtype=int), sampled_cell_types.astype(str).tolist()


def sparse_sum_row(matrix, positions: np.ndarray) -> sparse.csr_matrix:
    row = matrix[positions].sum(axis=0)
    if sparse.issparse(row):
        return row.tocsr()
    return sparse.csr_matrix(np.asarray(row))


def simulate_dataset(
    dataset_id: str,
    dataset_spec: dict[str, Any],
    references: dict[str, ad.AnnData],
    selected_cell_types: list[str],
    observed_probabilities: pd.Series,
    output_root: Path,
    num_spots: int,
    mean_cells_per_spot: float,
    min_cell_type_cells: int,
    seed: int,
    labels_key: str,
    overwrite: bool,
) -> Path:
    output_dir = output_root / dataset_id
    if output_dir.exists() and not overwrite:
        raise FileExistsError(f"{output_dir} already exists. Use --overwrite.")

    for modality in ["atac", "rna"]:
        (output_dir / modality / "features").mkdir(parents=True, exist_ok=True)
    (output_dir / "truth").mkdir(parents=True, exist_ok=True)
    (output_dir / "simulation").mkdir(parents=True, exist_ok=True)

    labels = references["atac"].obs[labels_key].astype(str)
    all_cell_types = pd.Index(pd.unique(labels))
    positions_by_cell_type = {
        cell_type: np.flatnonzero(labels.to_numpy() == cell_type)
        for cell_type in selected_cell_types
    }

    if dataset_spec["sampling_design"] == "equal_celltype":
        probabilities = pd.Series(1.0 / len(selected_cell_types), index=selected_cell_types)
    elif dataset_spec["sampling_design"] == "observed_abundance":
        probabilities = observed_probabilities.loc[selected_cell_types]
        probabilities = probabilities / probabilities.sum()
    else:
        raise ValueError(f"Unsupported sampling design: {dataset_spec['sampling_design']}")

    grid_size = int(np.sqrt(num_spots))
    if grid_size * grid_size != num_spots:
        raise ValueError("num_spots must be a perfect square for the synthetic grid.")
    x, y = np.meshgrid(np.arange(grid_size), np.arange(grid_size))
    coords = np.vstack((x.flatten(), y.flatten())).T * 10

    rng = np.random.default_rng(seed)
    cell_counts = rng.poisson(lam=mean_cells_per_spot, size=num_spots)
    cell_counts[cell_counts == 0] = 1

    matrices = {
        modality: adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
        for modality, adata in references.items()
    }
    rows = {modality: [] for modality in references}
    density = np.zeros((num_spots, len(all_cell_types)), dtype=np.float64)
    cell_type_to_density_index = {cell_type: idx for idx, cell_type in enumerate(all_cell_types)}
    sampled_cells_rows = []

    for spot_idx, n_cells in enumerate(cell_counts):
        positions, sampled_cell_types = sample_source_positions(
            rng=rng,
            cell_types=selected_cell_types,
            cell_type_probabilities=probabilities.to_numpy(dtype=float),
            positions_by_cell_type=positions_by_cell_type,
            n_cells=int(n_cells),
        )
        sampled_cell_ids = references["atac"].obs_names[positions].astype(str).tolist()

        for modality, matrix in matrices.items():
            rows[modality].append(sparse_sum_row(matrix, positions))

        for cell_type, count in pd.Series(sampled_cell_types).value_counts().items():
            density[spot_idx, cell_type_to_density_index[str(cell_type)]] = float(count)

        sampled_cells_rows.append(
            {
                "spot_id": f"spot_{spot_idx:04d}",
                "region": [0] * int(n_cells),
                "cell_id": sampled_cell_ids,
                "cell_type": sampled_cell_types,
            }
        )
        if (spot_idx + 1) % 100 == 0:
            print(f"{dataset_id}: sampled spots {spot_idx + 1}/{num_spots}", flush=True)

    spatial_paths = {}
    proportions = density / density.sum(axis=1, keepdims=True)
    for modality, adata in references.items():
        x_matrix = sparse.vstack(rows[modality], format="csr")
        spatial = ad.AnnData(X=x_matrix, obs=pd.DataFrame(index=[f"spot_{idx:04d}" for idx in range(num_spots)]))
        spatial.var = adata.var.copy()
        spatial.obs["cell_count"] = cell_counts
        spatial.obs["synthetic_region"] = "0"
        spatial.obs["simulation_design"] = dataset_spec["sampling_design"]
        spatial.obs["leiden_lsi"] = "0"
        spatial.obs["leiden_pca"] = "0"
        spatial.uns["density"] = pd.DataFrame(density, columns=all_cell_types, index=spatial.obs_names)
        spatial.uns["proportion_names"] = all_cell_types.astype(str).to_numpy()
        spatial.obsm["proportions"] = proportions
        spatial.obsm["spatial"] = coords.astype(int)

        path = output_dir / modality / "spatial.h5ad"
        print(f"{dataset_id}: writing {relpath(path)}", flush=True)
        spatial.write_h5ad(path, compression="gzip")
        spatial_paths[modality] = path
        del spatial

    truth_path = output_dir / "truth" / "proportions.csv"
    pd.DataFrame(
        proportions,
        index=[f"spot_{idx:04d}" for idx in range(num_spots)],
        columns=all_cell_types,
    ).to_csv(truth_path)

    sample_path = output_dir / "simulation" / "source_cells_by_spot.jsonl"
    with sample_path.open("w") as handle:
        for row in sampled_cells_rows:
            handle.write(json.dumps(row) + "\n")

    config = build_dataset_config(
        dataset_id=dataset_id,
        dataset_spec=dataset_spec,
        output_dir=output_dir,
        spatial_paths=spatial_paths,
        truth_path=truth_path,
        sample_path=sample_path,
        selected_cell_types=selected_cell_types,
        observed_probabilities=observed_probabilities,
        probabilities=probabilities,
        num_spots=num_spots,
        mean_cells_per_spot=mean_cells_per_spot,
        min_cell_type_cells=min_cell_type_cells,
        seed=seed,
    )
    config_path = output_dir / "dataset.yaml"
    with config_path.open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return config_path


def build_dataset_config(
    dataset_id: str,
    dataset_spec: dict[str, Any],
    output_dir: Path,
    spatial_paths: dict[str, Path],
    truth_path: Path,
    sample_path: Path,
    selected_cell_types: list[str],
    observed_probabilities: pd.Series,
    probabilities: pd.Series,
    num_spots: int,
    mean_cells_per_spot: float,
    min_cell_type_cells: int,
    seed: int,
) -> dict[str, Any]:
    observed_selected = observed_probabilities.loc[selected_cell_types]
    observed_selected = observed_selected / observed_selected.sum()
    return {
        "dataset_id": dataset_id,
        "source": "pbmc_multiome_simulation",
        "description": dataset_spec["description"],
        "reference_id": "pbmc_granulocyte_sorted_10k_multiome",
        "labels_key": "cell_type",
        "spatial_key": "spatial",
        "simulation": {
            "parameters": {
                "sampling_design": dataset_spec["sampling_design"],
                "num_spots": num_spots,
                "grid_shape": [int(np.sqrt(num_spots)), int(np.sqrt(num_spots))],
                "mean_cells_per_spot": float(mean_cells_per_spot),
                "cell_count_distribution": "poisson_min_1",
                "selected_cell_types": selected_cell_types,
                "selected_cell_type_min_source_cells": int(min_cell_type_cells),
                "observed_abundance_among_selected_cell_types": {
                    str(key): float(value) for key, value in observed_selected.items()
                },
                "sampling_probabilities": {str(key): float(value) for key, value in probabilities.items()},
                "seed": seed,
            },
            "source_cells_by_spot": relpath(sample_path),
        },
        "modalities": {
            "atac": {
                "reference": {
                    "path": "data/raw/references/pbmc_granulocyte_sorted_10k_multiome/atac/reference.h5ad"
                },
                "spatial": {"path": relpath(spatial_paths["atac"])},
                "labels_key": "cell_type",
                "spatial_key": "spatial",
                "truth": {"path": relpath(truth_path)},
                "feature_sets": {
                    "all": {"mode": "all"},
                    "highly_variable": {"path": relpath(output_dir / "atac" / "features" / "highly_variable.txt")},
                    "highly_accessible": {"path": relpath(output_dir / "atac" / "features" / "highly_accessible.txt")},
                },
            },
            "rna": {
                "reference": {
                    "path": "data/raw/references/pbmc_granulocyte_sorted_10k_multiome/rna/reference.h5ad"
                },
                "spatial": {"path": relpath(spatial_paths["rna"])},
                "labels_key": "cell_type",
                "spatial_key": "spatial",
                "truth": {"path": relpath(truth_path)},
                "feature_sets": {
                    "all": {"mode": "all"},
                    "highly_variable": {"path": relpath(output_dir / "rna" / "features" / "highly_variable.txt")},
                },
            },
        },
    }


def update_registry(registry_path: Path, generated_configs: dict[str, Path]) -> None:
    with registry_path.open() as handle:
        registry = yaml.safe_load(handle) or {}
    for dataset_id, config_path in generated_configs.items():
        registry[dataset_id] = {"config": relpath(config_path)}
    with registry_path.open("w") as handle:
        yaml.safe_dump(registry, handle, sort_keys=False)


def normalize_dataset_ids(dataset_ids: list[str]) -> list[str]:
    aliases = {
        "equal": "pbmc_granulocyte_sorted_10k_sim_equal_celltype",
        "equal_celltype": "pbmc_granulocyte_sorted_10k_sim_equal_celltype",
        "observed": "pbmc_granulocyte_sorted_10k_sim_observed_abundance",
        "observed_abundance": "pbmc_granulocyte_sorted_10k_sim_observed_abundance",
    }
    normalized = [aliases.get(dataset_id, dataset_id) for dataset_id in dataset_ids]
    unknown = sorted(set(normalized).difference(PBMC_DATASETS))
    if unknown:
        raise KeyError(f"Unknown PBMC dataset(s): {', '.join(unknown)}")
    return normalized


def main() -> None:
    parser = argparse.ArgumentParser(description="Regenerate PBMC multiome synthetic spatial datasets.")
    parser.add_argument("--datasets", nargs="+", default=sorted(PBMC_DATASETS))
    parser.add_argument("--output-root", default=str(ROOT / "data" / "processed" / "datasets"))
    parser.add_argument(
        "--atac-reference",
        default=str(ROOT / "data" / "raw" / "references" / "pbmc_granulocyte_sorted_10k_multiome" / "atac" / "reference.h5ad"),
    )
    parser.add_argument(
        "--rna-reference",
        default=str(ROOT / "data" / "raw" / "references" / "pbmc_granulocyte_sorted_10k_multiome" / "rna" / "reference.h5ad"),
    )
    parser.add_argument("--registry", default=str(ROOT / "data" / "registry" / "datasets.yaml"))
    parser.add_argument("--num-spots", type=int, default=1024)
    parser.add_argument("--mean-cells-per-spot", type=float, default=10.0)
    parser.add_argument("--min-cell-type-cells", type=int, default=100)
    parser.add_argument("--n-top-features", type=int, default=20000)
    parser.add_argument("--chunk-size", type=int, default=2048)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-registry-update", action="store_true")
    args = parser.parse_args()

    dataset_ids = normalize_dataset_ids(args.datasets)
    output_root = Path(args.output_root)

    print("reading PBMC references", flush=True)
    references = {
        "atac": read_reference(Path(args.atac_reference)),
        "rna": read_reference(Path(args.rna_reference)),
    }
    validate_references(references["atac"], references["rna"], labels_key="cell_type")

    counts = references["atac"].obs["cell_type"].astype(str).value_counts()
    selected_cell_types = counts[counts >= args.min_cell_type_cells].index.astype(str).tolist()
    if len(selected_cell_types) < 2:
        raise ValueError("Need at least two selected cell types for PBMC simulation.")
    observed_probabilities = counts / counts.sum()
    observed_probabilities = observed_probabilities.astype(float)
    print("selected cell types:", ", ".join(selected_cell_types), flush=True)

    print("computing PBMC feature sets", flush=True)
    feature_sets = {
        "atac": compute_feature_sets(
            references["atac"],
            labels_key="cell_type",
            selected_cell_types=selected_cell_types,
            modality="atac",
            n_top_features=args.n_top_features,
            chunk_size=args.chunk_size,
        ),
        "rna": compute_feature_sets(
            references["rna"],
            labels_key="cell_type",
            selected_cell_types=selected_cell_types,
            modality="rna",
            n_top_features=args.n_top_features,
            chunk_size=args.chunk_size,
        ),
    }

    generated_configs = {}
    for offset, dataset_id in enumerate(dataset_ids):
        config_path = simulate_dataset(
            dataset_id=dataset_id,
            dataset_spec=PBMC_DATASETS[dataset_id],
            references=references,
            selected_cell_types=selected_cell_types,
            observed_probabilities=observed_probabilities,
            output_root=output_root,
            num_spots=args.num_spots,
            mean_cells_per_spot=args.mean_cells_per_spot,
            min_cell_type_cells=args.min_cell_type_cells,
            seed=args.seed + offset,
            labels_key="cell_type",
            overwrite=args.overwrite,
        )
        write_feature_lists(output_root / dataset_id, feature_sets)
        generated_configs[dataset_id] = config_path
        gc.collect()

    if not args.no_registry_update:
        update_registry(Path(args.registry), generated_configs)
        print(f"updated {relpath(Path(args.registry))}", flush=True)

    print("generated configs:", flush=True)
    for dataset_id, config_path in generated_configs.items():
        print(f"  {dataset_id}: {relpath(config_path)}", flush=True)


if __name__ == "__main__":
    main()
