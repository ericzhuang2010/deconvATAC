#!/usr/bin/env python3
"""Summarize real-spatial ShapeMix runs without manufacturing composition truth."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_VALIDATION_CONFIG = (
    ROOT / "configs/experiments/shapemix_real_spatial_validation_v1.yaml"
)
METHOD_IDS = ("shapemix_length", "shapemix_count_only", "nnls")
SHAPE_METHOD_IDS = ("shapemix_length", "shapemix_count_only")
OUTPUT_NAMES = (
    "map_concordance.csv",
    "spatial_continuity.csv",
    "boundary_agreement.csv",
    "cross_modality_concordance.csv",
    "replicate_consistency.csv",
    "reconstruction_warnings.csv",
    "run_resources.csv",
    "evidence_summary.yaml",
)


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, Mapping):
        raise ValueError(f"Expected a YAML mapping: {path}")
    return dict(value)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def project_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def descriptor(dataset_id: str) -> dict[str, Any]:
    value = read_yaml(ROOT / "data/processed/datasets" / dataset_id / "dataset.yaml")
    if value.get("benchmark_scope") != "real_spatial_orthogonal_validation":
        raise ValueError(f"Dataset is not a real-spatial validation section: {dataset_id}")
    validation = value.get("validation")
    if not isinstance(validation, Mapping) or validation.get("exact_truth") is not False:
        raise ValueError(f"Real-spatial validation contract is absent: {dataset_id}")
    if "truth" in value or "truth" in value["modalities"]["atac"]:
        raise ValueError(f"Real-spatial descriptor improperly declares truth: {dataset_id}")
    return value


def run_directory(batch: Path, dataset_id: str, method_id: str) -> Path:
    path = batch / f"{dataset_id}__atac__all__{method_id}"
    required = (
        "run.yaml",
        "inputs.yaml",
        "output_sha256.yaml",
        "results/proportions.csv",
        "results/diagnostics.json",
    )
    missing = [name for name in required if not (path / name).is_file()]
    if missing:
        raise FileNotFoundError(f"{path} is missing {missing}")
    run = read_yaml(path / "run.yaml")
    if run.get("status") != "success":
        raise ValueError(f"Run is not successful: {path}")
    return path


def read_prediction(run_dir: Path, cell_types: Sequence[str]) -> pd.DataFrame:
    value = pd.read_csv(run_dir / "results/proportions.csv", index_col=0)
    if value.empty or value.index.has_duplicates or value.columns.has_duplicates:
        raise ValueError(f"Invalid prediction axes: {run_dir}")
    if value.columns.tolist() != list(cell_types):
        raise ValueError(f"Prediction cell-type universe changed: {run_dir}")
    matrix = value.to_numpy(dtype=np.float64)
    if (
        not np.isfinite(matrix).all()
        or (matrix < 0).any()
        or not np.allclose(matrix.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-6)
    ):
        raise ValueError(f"Invalid prediction proportions: {run_dir}")
    return value


def load_spatial(dataset: Mapping[str, Any]) -> ad.AnnData:
    path = project_path(dataset["modalities"]["atac"]["spatial"]["path"])
    value = ad.read_h5ad(path)
    if "spatial" not in value.obsm:
        raise KeyError(f"Spatial coordinates are absent: {path}")
    coordinates = np.asarray(value.obsm["spatial"], dtype=np.float64)
    if coordinates.shape != (value.n_obs, 2) or not np.isfinite(coordinates).all():
        raise ValueError(f"Invalid spatial coordinates: {path}")
    return value


def safe_correlation(left: np.ndarray, right: np.ndarray, method: str) -> tuple[float, bool]:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    if left.shape != right.shape or left.ndim != 1 or len(left) < 2:
        raise ValueError("Correlation inputs must be aligned vectors of length at least two")
    constant = bool(np.ptp(left) == 0 or np.ptp(right) == 0)
    if constant:
        return float("nan"), True
    if method == "pearson":
        return float(pearsonr(left, right).statistic), False
    if method == "spearman":
        return float(spearmanr(left, right).statistic), False
    raise ValueError(method)


def neighbor_edges(coordinates: np.ndarray, neighbors: int) -> np.ndarray:
    coordinates = np.asarray(coordinates, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2 or len(coordinates) < 2:
        raise ValueError("A spatial graph requires at least two two-dimensional coordinates")
    k = min(int(neighbors) + 1, len(coordinates))
    if k < 2:
        raise ValueError("neighbors must be positive")
    indices = cKDTree(coordinates).query(coordinates, k=k)[1]
    if indices.ndim == 1:
        indices = indices[:, None]
    edges = {
        tuple(sorted((row, int(column))))
        for row in range(len(coordinates))
        for column in indices[row, 1:]
        if row != int(column)
    }
    if not edges:
        raise ValueError("The spatial graph has no edges")
    return np.asarray(sorted(edges), dtype=np.int64)


def morans_i(values: np.ndarray, edges: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    centered = values - values.mean()
    denominator = float(np.square(centered).sum())
    if denominator == 0:
        return float("nan")
    numerator = float(
        (centered[edges[:, 0]] * centered[edges[:, 1]]).sum() * 2.0
    )
    total_weight = float(len(edges) * 2)
    return float(len(values) / total_weight * numerator / denominator)


def map_and_spatial_tables(
    experiment: Mapping[str, Any],
    batch: Path,
    validation_config: Mapping[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    map_rows: list[dict[str, Any]] = []
    continuity_rows: list[dict[str, Any]] = []
    boundary_rows: list[dict[str, Any]] = []
    neighbors = int(validation_config["spatial_graph"]["nearest_neighbors"])
    boundary_quantile = float(validation_config["spatial_graph"]["boundary_edge_quantile"])
    for dataset_id in experiment["datasets"]:
        dataset = descriptor(str(dataset_id))
        cell_types = list(dataset["modalities"]["atac"]["cell_types"])
        spatial = load_spatial(dataset)
        edges = neighbor_edges(np.asarray(spatial.obsm["spatial"]), neighbors)
        predictions = {
            method: read_prediction(run_directory(batch, str(dataset_id), method), cell_types)
            for method in METHOD_IDS
        }
        for method, value in predictions.items():
            if value.index.tolist() != spatial.obs_names.astype(str).tolist():
                raise ValueError(f"Prediction/spatial spot axes differ: {dataset_id}/{method}")
            for cell_type in cell_types:
                vector = value[cell_type].to_numpy(dtype=np.float64)
                gradients = np.abs(vector[edges[:, 0]] - vector[edges[:, 1]])
                continuity_rows.append(
                    {
                        "dataset_id": dataset_id,
                        "method_id": method,
                        "cell_type": cell_type,
                        "spots": len(vector),
                        "undirected_edges": len(edges),
                        "mean_neighbor_absolute_difference": float(gradients.mean()),
                        "median_neighbor_absolute_difference": float(np.median(gradients)),
                        "morans_i_knn": morans_i(vector, edges),
                        "interpretation": "descriptive_no_spatial_prior",
                    }
                )
        length = predictions["shapemix_length"]
        count = predictions["shapemix_count_only"]
        for cell_type in cell_types:
            left = length[cell_type].to_numpy(dtype=np.float64)
            right = count[cell_type].to_numpy(dtype=np.float64)
            pearson, pearson_constant = safe_correlation(left, right, "pearson")
            spearman, spearman_constant = safe_correlation(left, right, "spearman")
            difference = left - right
            map_rows.append(
                {
                    "dataset_id": dataset_id,
                    "cell_type": cell_type,
                    "spots": len(left),
                    "pearson_r": pearson,
                    "spearman_r": spearman,
                    "constant_map": pearson_constant or spearman_constant,
                    "mean_absolute_difference": float(np.abs(difference).mean()),
                    "median_absolute_difference": float(np.median(np.abs(difference))),
                    "maximum_absolute_difference": float(np.abs(difference).max()),
                    "mean_length_minus_count_only": float(difference.mean()),
                    "evidence_class": "paired_prediction_map_stability",
                }
            )
            length_gradient = np.abs(left[edges[:, 0]] - left[edges[:, 1]])
            count_gradient = np.abs(right[edges[:, 0]] - right[edges[:, 1]])
            length_threshold = float(np.quantile(length_gradient, boundary_quantile))
            count_threshold = float(np.quantile(count_gradient, boundary_quantile))
            length_boundary = set(np.flatnonzero(length_gradient >= length_threshold).tolist())
            count_boundary = set(np.flatnonzero(count_gradient >= count_threshold).tolist())
            union = length_boundary.union(count_boundary)
            boundary_rows.append(
                {
                    "dataset_id": dataset_id,
                    "cell_type": cell_type,
                    "undirected_edges": len(edges),
                    "boundary_edge_quantile": boundary_quantile,
                    "length_boundary_edges": len(length_boundary),
                    "count_boundary_edges": len(count_boundary),
                    "boundary_edge_jaccard": (
                        len(length_boundary.intersection(count_boundary)) / len(union)
                        if union
                        else float("nan")
                    ),
                    "interpretation": "paired_method_boundary_agreement_descriptive",
                }
            )
    return (
        pd.DataFrame.from_records(map_rows),
        pd.DataFrame.from_records(continuity_rows),
        pd.DataFrame.from_records(boundary_rows),
    )


def feature_lookup(names: Iterable[str]) -> dict[str, int]:
    result: dict[str, int] = {}
    for index, value in enumerate(names):
        key = str(value).casefold()
        result.setdefault(key, index)
    return result


def marker_score(
    matrix: ad.AnnData,
    markers: Sequence[str],
    *,
    scale_factor: float,
    minimum_features: int,
) -> tuple[pd.Series, list[str], list[str]]:
    lookup = feature_lookup(matrix.var_names.astype(str))
    present = [str(marker) for marker in markers if str(marker).casefold() in lookup]
    missing = [str(marker) for marker in markers if str(marker).casefold() not in lookup]
    indices = [lookup[marker.casefold()] for marker in present]
    if len(indices) < minimum_features:
        raise ValueError(
            f"Only {len(indices)} marker features are present; need {minimum_features}; "
            f"available={present} missing={missing}"
        )
    x = sparse.csr_matrix(matrix.X)
    library = np.asarray(x.sum(axis=1)).ravel().astype(np.float64)
    selected = x[:, indices].toarray().astype(np.float64)
    nonzero = library > 0
    selected[nonzero] = np.log1p(
        selected[nonzero] / library[nonzero, None] * scale_factor
    )
    selected[~nonzero] = 0.0
    standard_deviation = selected.std(axis=0, ddof=0)
    variable = standard_deviation > 0
    if int(variable.sum()) < minimum_features:
        raise ValueError("Too few nonconstant marker features remain after normalization")
    standardized = (
        selected[:, variable] - selected[:, variable].mean(axis=0)
    ) / standard_deviation[variable]
    score = standardized.mean(axis=1)
    return (
        pd.Series(score, index=matrix.obs_names.astype(str), name="marker_score"),
        [present[index] for index in np.flatnonzero(variable)],
        missing,
    )


def marker_peak_panel(reference_id: str, validation_config: Mapping[str, Any]) -> dict[str, Any]:
    panel_path = project_path(
        validation_config["reference_panels"][reference_id]["marker_features"]
    )
    panel = read_yaml(panel_path)
    if panel.get("reference_id") != reference_id or panel.get("outcome_data_used") is not False:
        raise ValueError(f"Invalid reference-only marker panel: {panel_path}")
    return panel


def correlate_score(
    prediction: pd.DataFrame,
    cell_type: str,
    score: pd.Series,
) -> tuple[int, float, bool]:
    shared = prediction.index.intersection(score.index)
    if len(shared) < 3:
        raise ValueError(f"Fewer than three aligned spots for {cell_type}")
    correlation, constant = safe_correlation(
        prediction.loc[shared, cell_type].to_numpy(dtype=np.float64),
        score.loc[shared].to_numpy(dtype=np.float64),
        "spearman",
    )
    return len(shared), correlation, constant


def cross_modality_table(
    experiment: Mapping[str, Any],
    batch: Path,
    validation_config: Mapping[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    score_config = validation_config["score_construction"]
    scale_factor = float(score_config["library_scale_factor"])
    minimum_features = int(score_config["minimum_present_features"])
    for dataset_id in experiment["datasets"]:
        dataset = descriptor(str(dataset_id))
        modality = dataset["modalities"]["atac"]
        cell_types = list(modality["cell_types"])
        reference_id = str(dataset["evaluation_design"].get("reference_id", ""))
        if not reference_id:
            reference_path = Path(modality["reference"]["path"])
            reference_id = reference_path.parents[1].name
        panel_config = validation_config["reference_panels"][reference_id]
        predictions = {
            method: read_prediction(run_directory(batch, str(dataset_id), method), cell_types)
            for method in METHOD_IDS
        }
        spatial = load_spatial(dataset)
        peak_panel = marker_peak_panel(reference_id, validation_config)
        matrices: list[tuple[str, str | None, ad.AnnData, Mapping[str, Any], str]] = []
        validation = dataset["validation"]
        if "rna" in validation:
            matrices.append(
                (
                    "rna_marker_score",
                    str(validation["rna"]["gsm"]),
                    ad.read_h5ad(project_path(validation["rna"]["path"])),
                    panel_config["rna_markers"],
                    "independent_orthogonal_validation",
                )
            )
        if "protein" in validation and panel_config.get("protein_markers"):
            matrices.append(
                (
                    "protein_marker_score",
                    str(validation["protein"]["gsm"]),
                    ad.read_h5ad(project_path(validation["protein"]["path"])),
                    panel_config["protein_markers"],
                    "independent_orthogonal_validation",
                )
            )
        for evidence_type, gsm, matrix, markers_by_type, evidence_class in matrices:
            for cell_type, markers in markers_by_type.items():
                if cell_type not in cell_types:
                    raise ValueError(f"Marker panel type is absent from reference: {reference_id}/{cell_type}")
                score, present, missing = marker_score(
                    matrix,
                    list(markers),
                    scale_factor=scale_factor,
                    minimum_features=minimum_features,
                )
                for method, prediction in predictions.items():
                    spots, correlation, constant = correlate_score(prediction, cell_type, score)
                    rows.append(
                        {
                            "dataset_id": dataset_id,
                            "method_id": method,
                            "cell_type": cell_type,
                            "evidence_type": evidence_type,
                            "evidence_gsm": gsm,
                            "assay": "rna" if evidence_type.startswith("rna") else "protein",
                            "aligned_spots": spots,
                            "spearman_r": correlation,
                            "constant_map_or_score": constant,
                            "present_markers": ";".join(present),
                            "missing_markers": ";".join(missing),
                            "evidence_class": evidence_class,
                        }
                    )
        marker_mapping = {
            cell_type: list(peak_panel["markers"][cell_type]["features"])
            for cell_type in cell_types
        }
        atac_matrices: list[tuple[str, str | None, str, ad.AnnData, str]] = [
            ("atac_reference_marker_accessibility", None, "atac", spatial, "within_assay_descriptive")
        ]
        for epigenome in validation.get("epigenome", []):
            atac_matrices.append(
                (
                    "epigenome_reference_marker_accessibility",
                    str(epigenome["gsm"]),
                    str(epigenome["assay"]),
                    ad.read_h5ad(project_path(epigenome["path"])),
                    "independent_orthogonal_validation",
                )
            )
        for evidence_type, gsm, assay, matrix, evidence_class in atac_matrices:
            for cell_type, markers in marker_mapping.items():
                score, present, missing = marker_score(
                    matrix,
                    markers,
                    scale_factor=scale_factor,
                    minimum_features=minimum_features,
                )
                for method, prediction in predictions.items():
                    spots, correlation, constant = correlate_score(prediction, cell_type, score)
                    rows.append(
                        {
                            "dataset_id": dataset_id,
                            "method_id": method,
                            "cell_type": cell_type,
                            "evidence_type": evidence_type,
                            "evidence_gsm": gsm,
                            "assay": assay,
                            "aligned_spots": spots,
                            "spearman_r": correlation,
                            "constant_map_or_score": constant,
                            "present_markers": ";".join(present),
                            "missing_markers": ";".join(missing),
                            "evidence_class": evidence_class,
                        }
                    )
    return pd.DataFrame.from_records(rows)


def replicate_table(
    experiments: Sequence[Mapping[str, Any]],
    batches: Sequence[Path],
    validation_config: Mapping[str, Any],
) -> pd.DataFrame:
    dataset_batch = {
        str(dataset_id): batch
        for experiment, batch in zip(experiments, batches, strict=True)
        for dataset_id in experiment["datasets"]
    }
    rows: list[dict[str, Any]] = []
    for pair in validation_config.get("replicate_pairs", []):
        left_id, right_id = str(pair["left"]), str(pair["right"])
        left_descriptor, right_descriptor = descriptor(left_id), descriptor(right_id)
        left_types = list(left_descriptor["modalities"]["atac"]["cell_types"])
        right_types = list(right_descriptor["modalities"]["atac"]["cell_types"])
        if left_types != right_types:
            raise ValueError(f"Replicate reference axes differ: {left_id}/{right_id}")
        for method in METHOD_IDS:
            left = read_prediction(run_directory(dataset_batch[left_id], left_id, method), left_types)
            right = read_prediction(run_directory(dataset_batch[right_id], right_id, method), right_types)
            shared = left.index.intersection(right.index)
            if len(shared) < int(pair["minimum_shared_spots"]):
                raise ValueError(f"Insufficient replicate spot overlap: {left_id}/{right_id}")
            for cell_type in left_types:
                left_values = left.loc[shared, cell_type].to_numpy(dtype=np.float64)
                right_values = right.loc[shared, cell_type].to_numpy(dtype=np.float64)
                correlation, constant = safe_correlation(left_values, right_values, "spearman")
                rows.append(
                    {
                        "pair_id": str(pair["id"]),
                        "left_dataset_id": left_id,
                        "right_dataset_id": right_id,
                        "method_id": method,
                        "cell_type": cell_type,
                        "shared_spots": len(shared),
                        "spearman_r": correlation,
                        "constant_map": constant,
                        "mean_absolute_difference": float(np.abs(left_values - right_values).mean()),
                        "interpretation": "descriptive_independent_section_consistency",
                    }
                )
    return pd.DataFrame.from_records(rows)


def diagnostics_and_resources(
    experiment: Mapping[str, Any],
    batch: Path,
    validation_config: Mapping[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    warning_threshold = float(
        validation_config["residual_warning"]["normalized_absolute_error_threshold"]
    )
    warnings: list[dict[str, Any]] = []
    resources: list[dict[str, Any]] = []
    for dataset_id in experiment["datasets"]:
        dataset = descriptor(str(dataset_id))
        spatial = load_spatial(dataset)
        spatial_x = sparse.csr_matrix(spatial.X, dtype=np.float64)
        mean_row_l2 = float(np.sqrt(spatial_x.multiply(spatial_x).sum(axis=1)).mean())
        for method in METHOD_IDS:
            run_dir = run_directory(batch, str(dataset_id), method)
            run = read_yaml(run_dir / "run.yaml")
            with (run_dir / "results/diagnostics.json").open() as handle:
                diagnostics = json.load(handle)
            performance = run.get("performance", {})
            fit = diagnostics.get("fit", {}) if isinstance(diagnostics, Mapping) else {}
            execution = fit.get("execution", {}) if isinstance(fit, Mapping) else {}
            resources.append(
                {
                    "dataset_id": dataset_id,
                    "method_id": method,
                    "wall_runtime_seconds": performance.get("wall_runtime_seconds"),
                    "peak_rss_bytes": performance.get("peak_rss_bytes"),
                    "device": fit.get("device") if isinstance(fit, Mapping) else None,
                    "count_cache_mode": execution.get("count_cache_mode"),
                    "peak_device_memory_allocated_bytes": execution.get(
                        "peak_device_memory_allocated_bytes"
                    ),
                    "peak_device_memory_reserved_bytes": execution.get(
                        "peak_device_memory_reserved_bytes"
                    ),
                    "run_dir": str(run_dir.relative_to(ROOT)),
                }
            )
            if method in SHAPE_METHOD_IDS:
                count = diagnostics["reconstruction"]["count"]
                observed_total = float(count["observed_total"])
                normalized = (
                    float(count["absolute_error_sum"]) / observed_total
                    if observed_total > 0
                    else float("inf")
                )
                residual_fraction = (
                    abs(float(count["residual_total"])) / observed_total
                    if observed_total > 0
                    else float("inf")
                )
                converged = bool(fit.get("success", False))
                background_mode = diagnostics.get("config", {}).get("background_mode")
            else:
                normalized = (
                    float(diagnostics["mean_residual"]) / mean_row_l2
                    if mean_row_l2 > 0
                    else float("inf")
                )
                residual_fraction = float("nan")
                converged = True
                background_mode = "not_available_nnls"
            warnings.append(
                {
                    "dataset_id": dataset_id,
                    "method_id": method,
                    "normalized_reconstruction_error_proxy": normalized,
                    "absolute_total_residual_fraction": residual_fraction,
                    "warning_threshold": warning_threshold,
                    "off_reference_or_mismatch_warning": (
                        not converged or not np.isfinite(normalized) or normalized > warning_threshold
                    ),
                    "converged": converged,
                    "background_mode": background_mode,
                    "interpretation": (
                        "proxy_warning_not_identified_off_reference_mass; high error may reflect "
                        "protocol, disease, reference, or model mismatch"
                    ),
                }
            )
    return pd.DataFrame.from_records(warnings), pd.DataFrame.from_records(resources)


def write_summary(
    batch: Path,
    experiment: Mapping[str, Any],
    validation_config_path: Path,
    tables: Mapping[str, pd.DataFrame],
    overwrite: bool,
) -> None:
    existing = [name for name in OUTPUT_NAMES if (batch / name).exists()]
    if existing and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing summaries: {existing}")
    for filename, frame in tables.items():
        frame.to_csv(batch / filename, index=False)
    summary = {
        "schema_version": 1,
        "status": "complete",
        "evidence_class": "real_spatial_orthogonal_validation",
        "exact_composition_truth": False,
        "experiment_run_group": experiment["run_group"],
        "datasets": len(experiment["datasets"]),
        "methods": list(METHOD_IDS),
        "validation_config": str(validation_config_path.relative_to(ROOT)),
        "validation_config_sha256": digest(validation_config_path),
        "tables": {
            filename: {
                "rows": len(frame),
                "sha256": digest(batch / filename),
            }
            for filename, frame in tables.items()
        },
        "claims": {
            "map_and_cross_modality_results_are_descriptive": True,
            "spatial_sections_are_not_treated_as_independent_spots_for_inference": True,
            "no_exact_truth_metrics_or_p_values": True,
            "histone_rna_and_protein_not_used_as_shapemix_inputs": True,
            "off_reference_values_are_warning_proxies_not_identified_mass": True,
            "anatomical_region_labels_available": False,
        },
    }
    with (batch / "evidence_summary.yaml").open("w") as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-config", action="append", type=Path, required=True)
    parser.add_argument("--batch-dir", action="append", type=Path, required=True)
    parser.add_argument("--validation-config", type=Path, default=DEFAULT_VALIDATION_CONFIG)
    parser.add_argument("--overwrite-summary", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if len(args.experiment_config) != len(args.batch_dir):
        raise ValueError("--experiment-config and --batch-dir counts must match")
    validation_path = project_path(args.validation_config)
    validation_config = read_yaml(validation_path)
    experiments = [read_yaml(project_path(path)) for path in args.experiment_config]
    batches = [project_path(path) for path in args.batch_dir]
    for batch in batches:
        if not batch.is_dir():
            raise FileNotFoundError(batch)
    replicate = replicate_table(experiments, batches, validation_config)
    for experiment, batch in zip(experiments, batches, strict=True):
        maps, continuity, boundaries = map_and_spatial_tables(
            experiment, batch, validation_config
        )
        cross_modality = cross_modality_table(experiment, batch, validation_config)
        warnings, resources = diagnostics_and_resources(
            experiment, batch, validation_config
        )
        dataset_ids = set(str(value) for value in experiment["datasets"])
        replicate_subset = replicate[
            replicate["left_dataset_id"].isin(dataset_ids)
            | replicate["right_dataset_id"].isin(dataset_ids)
        ].copy()
        write_summary(
            batch,
            experiment,
            validation_path,
            {
                "map_concordance.csv": maps,
                "spatial_continuity.csv": continuity,
                "boundary_agreement.csv": boundaries,
                "cross_modality_concordance.csv": cross_modality,
                "replicate_consistency.csv": replicate_subset,
                "reconstruction_warnings.csv": warnings,
                "run_resources.csv": resources,
            },
            args.overwrite_summary,
        )


if __name__ == "__main__":
    main()
