#!/usr/bin/env python
from __future__ import annotations

import argparse
import gc
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.tl.simulate import Sampler, conway_maxwell_poisson


HEART_DATASET_ALIASES = {
    "1": "human_cardiac_niches_sim_1zone_3ct_low_density",
    "heart1": "human_cardiac_niches_sim_1zone_3ct_low_density",
    "heart_1": "human_cardiac_niches_sim_1zone_3ct_low_density",
    "2": "human_cardiac_niches_sim_1zone_10ct",
    "heart2": "human_cardiac_niches_sim_1zone_10ct",
    "heart_2": "human_cardiac_niches_sim_1zone_10ct",
    "3": "human_cardiac_niches_sim_4zone_stripes",
    "heart3": "human_cardiac_niches_sim_4zone_stripes",
    "heart_3": "human_cardiac_niches_sim_4zone_stripes",
    "4": "human_cardiac_niches_sim_4zone_circles",
    "heart4": "human_cardiac_niches_sim_4zone_circles",
    "heart_4": "human_cardiac_niches_sim_4zone_circles",
}


HEART_PARAMS: dict[str, dict[str, Any]] = {
    "human_cardiac_niches_sim_1zone_3ct_low_density": {
        "n_regions": 1,
        "cell_type_number": [3],
        "cell_number_nu": [20],
        "cell_number_mean": [5],
        "region_type": "stripes",
        "num_spots": 1000,
        "balance": "balanced",
    },
    "human_cardiac_niches_sim_1zone_10ct": {
        "n_regions": 1,
        "cell_type_number": [10],
        "cell_number_nu": [20],
        "cell_number_mean": [15],
        "region_type": "stripes",
        "num_spots": 1000,
        "balance": "balanced",
    },
    "human_cardiac_niches_sim_4zone_stripes": {
        "n_regions": 4,
        "cell_type_number": [10, 5, 10, 5],
        "cell_number_nu": [20, 20, 20, 20],
        "cell_number_mean": [15, 10, 15, 5],
        "region_type": "stripes",
        "num_spots": 1000,
        "balance": "balanced",
    },
    "human_cardiac_niches_sim_4zone_circles": {
        "n_regions": 4,
        "cell_type_number": [3, 5, 3, 5],
        "cell_number_nu": [20, 20, 20, 20],
        "cell_number_mean": [15, 10, 15, 5],
        "region_type": "circles",
        "num_spots": 1000,
        "balance": "balanced",
    },
}


def relpath(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def normalize_dataset_ids(dataset_ids: list[str]) -> list[str]:
    normalized = []
    for dataset_id in dataset_ids:
        lowered = dataset_id.lower()
        normalized.append(HEART_DATASET_ALIASES.get(lowered, lowered))
    unknown = sorted(set(normalized).difference(HEART_PARAMS))
    if unknown:
        available = ", ".join(sorted(HEART_PARAMS))
        raise KeyError(f"Unknown Heart dataset(s): {', '.join(unknown)}. Available: {available}")
    return normalized


def read_reference(path: Path) -> ad.AnnData:
    return ad.read_h5ad(path, backed="r")


def validate_references(atac: ad.AnnData, rna: ad.AnnData, labels_key: str) -> None:
    for modality, adata in {"atac": atac, "rna": rna}.items():
        if labels_key not in adata.obs:
            raise KeyError(f"{modality} reference is missing obs[{labels_key!r}].")
        if not isinstance(adata.obs[labels_key].dtype, pd.CategoricalDtype):
            adata.obs[labels_key] = adata.obs[labels_key].astype("category")

    missing_rna_cells = atac.obs_names.difference(rna.obs_names)
    if len(missing_rna_cells) > 0:
        raise ValueError(
            "RNA reference is missing cell IDs sampled from the ATAC reference. "
            f"First missing ID: {missing_rna_cells[0]}"
        )


def sparse_sum_row(adata: ad.AnnData, obs_names: list[str]) -> sparse.csr_matrix:
    row = adata[obs_names, :].X.sum(axis=0)
    if sparse.issparse(row):
        return row.tocsr()
    return sparse.csr_matrix(np.asarray(row))


def generate_sparse_spatial(
    references: dict[str, ad.AnnData],
    cell_type_key: str,
    **params: Any,
) -> tuple[mu.MuData, pd.DataFrame]:
    """Generate Heart spatial data without materializing dense feature matrices."""
    np.random.seed(0)
    sampler = Sampler(
        reference=references["atac"],
        cell_type_key=cell_type_key,
        num_spots=params["num_spots"],
        n_regions=params["n_regions"],
        cell_number_mean=params["cell_number_mean"],
        cell_number_nu=params["cell_number_nu"],
        cell_type_number=params["cell_type_number"],
        balance=params["balance"],
        region_type=params["region_type"],
    )

    used_clusters = {
        region: np.random.choice(sampler.cluster_p.index, size=cell_type_number, p=sampler.cluster_p, replace=False)
        for region, cell_type_number in enumerate(sampler.cell_type_number)
    }
    sampler.define_regions(used_clusters)

    cell_count = np.array(
        [
            conway_maxwell_poisson(sampler.cell_number_mean[region], sampler.cell_number_nu[region])
            for region in sampler.regions
        ]
    )
    cell_count[cell_count == 0] = 1

    for region, _ in enumerate(sampler.cell_type_number):
        print(
            f"Region {region}: {used_clusters[region]} "
            f"(mean cells per cluster: {cell_count[sampler.regions == region].mean()})",
            flush=True,
        )

    if "gradient" in sampler.region_type:
        spot_clusters = [used_clusters[0] for _ in sampler.regions]
    else:
        spot_clusters = [used_clusters[region] for region in sampler.regions]

    rows = {modality: [] for modality in references}
    density = np.zeros((len(sampler.regions), len(sampler.clusters)))
    sampled_cells_rows = []
    spot_params = list(zip(cell_count, spot_clusters, sampler.regions))

    for spot_idx, (num_cell, clusters, region) in enumerate(spot_params):
        cluster_mask = sampler.obs[cell_type_key].isin(clusters).values
        if sampler.region_type == "gradient_celltype":
            p = sampler.cell_p[region][cluster_mask] / sampler.cell_p[region][cluster_mask].sum()
        else:
            p = sampler.cell_p[cluster_mask] / sampler.cell_p[cluster_mask].sum()

        sampled_cells = np.random.choice(sampler.obs.index[cluster_mask], size=num_cell, p=p)
        sampled_cells_list = sampled_cells.tolist()

        for modality, adata in references.items():
            rows[modality].append(sparse_sum_row(adata, sampled_cells_list))

        density[spot_idx, :] = sampler.obs.loc[sampled_cells, cell_type_key].value_counts().loc[sampler.clusters].values
        sampled_cells_rows.append(
            {
                "spot_id": str(spot_idx),
                "region": [int(region)] * int(num_cell),
                "cell_id": sampled_cells_list,
                "cell_type": sampler.obs.loc[sampled_cells, cell_type_key].astype(str).values.tolist(),
            }
        )

        if (spot_idx + 1) % 100 == 0:
            print(f"sampled spots {spot_idx + 1}/{len(spot_params)}", flush=True)

    x_coords, y_coords = sampler.get_coords()
    coords = np.vstack((x_coords.flatten(), y_coords.flatten())).T * 10
    spatial_mod = {}

    for modality, adata in references.items():
        x = sparse.vstack(rows[modality], format="csr")
        spatial_ann = ad.AnnData(x)
        spatial_ann.var.index = adata.var_names.copy()
        spatial_ann.obs["cell_count"] = density.sum(axis=1)
        spatial_ann.uns["density"] = pd.DataFrame(density, columns=sampler.clusters, index=spatial_ann.obs_names)
        spatial_ann.obsm["proportions"] = density / density.sum(axis=1)[:, None]
        spatial_ann.uns["proportion_names"] = sampler.clusters.astype(str).values
        spatial_ann.obsm["spatial"] = coords.astype("int")
        spatial_mod[modality] = spatial_ann

    return mu.MuData(spatial_mod), pd.DataFrame(sampled_cells_rows)


def maybe_leiden(adata: ad.AnnData, key_added: str, skip_leiden: bool) -> None:
    if skip_leiden:
        return
    try:
        sc.tl.leiden(adata, key_added=key_added)
    except ImportError as error:
        raise ImportError(
            "Leiden clustering requires the optional packages 'leidenalg' and 'igraph'. "
            "Install them or rerun with --skip-leiden to regenerate files without Leiden labels."
        ) from error


def process_spatial(mudata: mu.MuData, modality: str, skip_leiden: bool) -> ad.AnnData:
    adata = mudata.mod[modality]
    adata.layers["counts"] = adata.X.copy()

    if modality == "rna":
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, use_rep="X_pca")
        maybe_leiden(adata, key_added="leiden_pca", skip_leiden=skip_leiden)
        adata.layers["log_norm"] = adata.X.copy()
        adata.X = adata.layers["counts"].copy()
        adata.layers.pop("counts")
        return adata

    if modality == "atac":
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.layers["log_norm"] = adata.X.copy()
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, use_rep="X_pca")
        maybe_leiden(adata, key_added="leiden_pca", skip_leiden=skip_leiden)

        adata.X = adata.layers["counts"].copy()
        mu.atac.pp.tfidf(adata)
        mu.atac.tl.lsi(adata)
        sc.pp.neighbors(adata, use_rep="X_lsi")
        maybe_leiden(adata, key_added="leiden_lsi", skip_leiden=skip_leiden)

        adata.layers["tfidf_normalized"] = adata.X.copy()
        adata.X = adata.layers["counts"].copy()
        adata.layers.pop("counts")
        return adata

    raise ValueError(f"Unsupported modality: {modality}")


def write_samples(output_dir: Path, dataset_id: str, samples: pd.DataFrame) -> Path:
    jsonl_path = output_dir / "source_cells_by_spot.jsonl"
    with jsonl_path.open("w") as handle:
        for row in samples.to_dict(orient="records"):
            handle.write(json.dumps(row) + "\n")
    print(f"{dataset_id}: wrote source-cell provenance to {relpath(jsonl_path)}", flush=True)
    return jsonl_path


def build_dataset_config(
    dataset_id: str,
    spatial_paths: dict[str, Path],
    atac_reference_path: Path,
    rna_reference_path: Path,
    truth_path: Path,
    params: dict[str, Any],
    sample_path: Path,
) -> dict[str, Any]:
    return {
        "dataset_id": dataset_id,
        "source": "regenerated_heart_simulation",
        "description": f"Regenerated Human cardiac niches simulation {dataset_id}.",
        "labels_key": "cell_type",
        "spatial_key": "spatial",
        "simulation": {
            "parameters": params,
            "source_cells_by_spot": relpath(sample_path),
        },
        "modalities": {
            "atac": {
                "reference": {"path": relpath(atac_reference_path)},
                "spatial": {"path": relpath(spatial_paths["atac"])},
                "labels_key": "cell_type",
                "spatial_key": "spatial",
                "truth": {"path": relpath(truth_path)},
                "feature_sets": {
                    "all": {"mode": "all"},
                    "highly_variable": {
                        "path": relpath(spatial_paths["atac"].parent / "features" / "highly_variable.txt")
                    },
                    "highly_accessible": {
                        "path": relpath(spatial_paths["atac"].parent / "features" / "highly_accessible.txt")
                    },
                },
            },
            "rna": {
                "reference": {"path": relpath(rna_reference_path)},
                "spatial": {"path": relpath(spatial_paths["rna"])},
                "labels_key": "cell_type",
                "spatial_key": "spatial",
                "truth": {"path": relpath(truth_path)},
                "feature_sets": {
                    "all": {"mode": "all"},
                    "highly_variable": {
                        "path": relpath(spatial_paths["rna"].parent / "features" / "highly_variable.txt")
                    },
                },
            },
        },
    }


def write_dataset_config(output_dir: Path, config: dict[str, Any]) -> Path:
    config_path = output_dir / "dataset.yaml"
    with config_path.open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return config_path


def update_registry(registry_path: Path, generated_configs: dict[str, Path]) -> None:
    with registry_path.open() as handle:
        registry = yaml.safe_load(handle) or {}

    for legacy_id in ("heart_1", "heart_2", "heart_3", "heart_4"):
        registry.pop(legacy_id, None)
    for dataset_id, config_path in generated_configs.items():
        registry[dataset_id] = {"config": relpath(config_path)}

    with registry_path.open("w") as handle:
        yaml.safe_dump(registry, handle, sort_keys=False)


def regenerate_dataset(
    dataset_id: str,
    references: dict[str, ad.AnnData],
    output_root: Path,
    atac_reference_path: Path,
    rna_reference_path: Path,
    overwrite: bool,
    skip_leiden: bool,
) -> Path:
    output_dir = output_root / dataset_id
    spatial_paths = {
        "atac": output_dir / "atac" / "spatial.h5ad",
        "rna": output_dir / "rna" / "spatial.h5ad",
    }
    if output_dir.exists() and any(path.exists() for path in spatial_paths.values()) and not overwrite:
        raise FileExistsError(f"{output_dir} already contains spatial files. Use --overwrite.")

    for modality in ["atac", "rna"]:
        (output_dir / modality / "features").mkdir(parents=True, exist_ok=True)
    (output_dir / "truth").mkdir(parents=True, exist_ok=True)
    simulation_dir = output_dir / "simulation"
    simulation_dir.mkdir(parents=True, exist_ok=True)

    print(f"{dataset_id}: generating spatial data", flush=True)
    mudata, samples = generate_sparse_spatial(references, cell_type_key="cell_type", **HEART_PARAMS[dataset_id])

    print(f"{dataset_id}: processing ATAC", flush=True)
    atac = process_spatial(mudata, "atac", skip_leiden=skip_leiden)
    print(f"{dataset_id}: processing RNA", flush=True)
    rna = process_spatial(mudata, "rna", skip_leiden=skip_leiden)

    print(f"{dataset_id}: writing {relpath(spatial_paths['atac'])}", flush=True)
    atac.write_h5ad(spatial_paths["atac"])
    print(f"{dataset_id}: writing {relpath(spatial_paths['rna'])}", flush=True)
    rna.write_h5ad(spatial_paths["rna"])

    truth_dir = output_dir / "truth"
    truth_path = truth_dir / "proportions.csv"
    pd.DataFrame(
        atac.obsm["proportions"],
        index=atac.obs_names,
        columns=atac.uns["proportion_names"],
    ).to_csv(truth_path)

    sample_path = write_samples(simulation_dir, dataset_id, samples)
    config = build_dataset_config(
        dataset_id=dataset_id,
        spatial_paths=spatial_paths,
        atac_reference_path=atac_reference_path,
        rna_reference_path=rna_reference_path,
        truth_path=truth_path,
        params=HEART_PARAMS[dataset_id],
        sample_path=sample_path,
    )
    config_path = write_dataset_config(output_dir, config)

    del mudata, atac, rna, samples
    gc.collect()
    return config_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Regenerate Heart spatial ATAC/RNA simulation files.")
    parser.add_argument("--datasets", nargs="+", default=sorted(HEART_PARAMS))
    parser.add_argument("--output-root", default=str(ROOT / "data" / "processed" / "datasets"))
    parser.add_argument(
        "--atac-reference",
        default=str(ROOT / "data" / "raw" / "references" / "human_cardiac_niches" / "atac" / "reference.h5ad"),
    )
    parser.add_argument(
        "--rna-reference",
        default=str(ROOT / "data" / "raw" / "references" / "human_cardiac_niches" / "rna" / "reference.h5ad"),
    )
    parser.add_argument("--registry", default=str(ROOT / "data" / "registry" / "datasets.yaml"))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-leiden", action="store_true")
    parser.add_argument("--no-registry-update", action="store_true")
    args = parser.parse_args()

    dataset_ids = normalize_dataset_ids(args.datasets)
    output_root = Path(args.output_root)
    atac_reference_path = Path(args.atac_reference)
    rna_reference_path = Path(args.rna_reference)

    print("reading references in backed mode", flush=True)
    references = {
        "atac": read_reference(atac_reference_path),
        "rna": read_reference(rna_reference_path),
    }
    validate_references(references["atac"], references["rna"], labels_key="cell_type")

    generated_configs = {}
    for dataset_id in dataset_ids:
        generated_configs[dataset_id] = regenerate_dataset(
            dataset_id=dataset_id,
            references=references,
            output_root=output_root,
            atac_reference_path=atac_reference_path,
            rna_reference_path=rna_reference_path,
            overwrite=args.overwrite,
            skip_leiden=args.skip_leiden,
        )

    if not args.no_registry_update:
        update_registry(Path(args.registry), generated_configs)
        print(f"updated {relpath(Path(args.registry))}", flush=True)

    print("generated configs:", flush=True)
    for dataset_id, config_path in generated_configs.items():
        print(f"  {dataset_id}: {relpath(config_path)}", flush=True)


if __name__ == "__main__":
    main()
