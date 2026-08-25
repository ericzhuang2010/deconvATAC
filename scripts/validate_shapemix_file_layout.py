#!/usr/bin/env python
"""Fail-closed validator for the canonical ShapeMix file organization."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml


ROOT = Path(__file__).resolve().parents[1]
ALLOWED_RESULT_SCOPES = {
    "development",
    "primary",
    "sensitivity",
    "external_validation",
    "real_spatial",
}
REQUIRED_RUN_FILES = {
    "run.yaml",
    "inputs.yaml",
    "environment.txt",
    "output_sha256.yaml",
    "results/proportions.csv",
    "results/diagnostics.json",
}


class LayoutError(ValueError):
    """Raised when an authoritative artifact violates the frozen layout."""


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, Mapping):
        raise LayoutError(f"Expected a YAML mapping: {path}")
    return dict(value)


def project_path(value: str, context: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        raise LayoutError(f"{context} must be project-relative: {value}")
    if ".." in path.parts:
        raise LayoutError(f"{context} cannot escape the project: {value}")
    return path


def require_prefix(path: Path, prefix: Path, context: str) -> None:
    try:
        path.relative_to(prefix)
    except ValueError as exc:
        raise LayoutError(f"{context} must be under {prefix}, found {path}") from exc


def iter_declared_paths(
    value: Any,
    key_path: tuple[str, ...] = (),
) -> Iterable[tuple[tuple[str, ...], str]]:
    if isinstance(value, Mapping):
        for key, child in value.items():
            child_path = (*key_path, str(key))
            if str(key) == "path" and isinstance(child, str):
                yield child_path, child
            else:
                yield from iter_declared_paths(child, child_path)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from iter_declared_paths(child, (*key_path, str(index)))


def canonical_result_root(value: str, context: str) -> Path:
    path = project_path(value, context)
    if (
        len(path.parts) < 3
        or path.parts[0] != "results"
        or path.parts[1] not in ALLOWED_RESULT_SCOPES
    ):
        raise LayoutError(
            f"{context} must be results/<canonical_scope>/<campaign_id>, found {path}"
        )
    return path


def load_registry(path: Path) -> dict[str, Any]:
    registry = load_yaml(path)
    for dataset_id, record in registry.items():
        if not isinstance(record, Mapping) or not isinstance(record.get("config"), str):
            raise LayoutError(f"Invalid registry entry for {dataset_id}")
        config = project_path(str(record["config"]), f"registry.{dataset_id}.config")
        expected = Path("data/processed/datasets") / str(dataset_id) / "dataset.yaml"
        if config != expected:
            raise LayoutError(
                f"registry.{dataset_id}.config must be {expected}, found {config}"
            )
    return registry


def validate_dataset(
    dataset_id: str,
    registry: Mapping[str, Any],
    evaluation_mode: str,
    require_standard_reference: bool,
) -> None:
    if dataset_id not in registry:
        raise LayoutError(f"Dataset is not registered: {dataset_id}")
    config_relative = project_path(
        str(registry[dataset_id]["config"]),
        f"registry.{dataset_id}.config",
    )
    config_path = ROOT / config_relative
    if not config_path.is_file():
        raise LayoutError(f"Registered dataset descriptor is missing: {config_relative}")
    config = load_yaml(config_path)
    if config.get("dataset_id") != dataset_id:
        raise LayoutError(f"Dataset ID mismatch in {config_relative}")

    dataset_root = Path("data/processed/datasets") / dataset_id
    for key_path, declared in iter_declared_paths(config):
        relative = project_path(declared, f"{dataset_id}.{'.'.join(key_path)}")
        if not (ROOT / relative).exists():
            raise LayoutError(f"Declared dataset artifact is missing: {relative}")

    modalities = config.get("modalities")
    if not isinstance(modalities, Mapping) or not modalities:
        raise LayoutError(f"{dataset_id} has no modality mappings")
    for modality, record in modalities.items():
        if not isinstance(record, Mapping):
            raise LayoutError(f"{dataset_id}.{modality} must be a mapping")
        reference = record.get("reference")
        spatial = record.get("spatial")
        if not isinstance(reference, Mapping) or not isinstance(reference.get("path"), str):
            raise LayoutError(f"{dataset_id}.{modality} has no reference path")
        if not isinstance(spatial, Mapping) or not isinstance(spatial.get("path"), str):
            raise LayoutError(f"{dataset_id}.{modality} has no spatial path")
        reference_path = project_path(
            str(reference["path"]), f"{dataset_id}.{modality}.reference.path"
        )
        spatial_path = project_path(
            str(spatial["path"]), f"{dataset_id}.{modality}.spatial.path"
        )
        if require_standard_reference:
            require_prefix(
                reference_path,
                Path("data/processed/references"),
                f"{dataset_id}.{modality}.reference.path",
            )
        require_prefix(
            spatial_path,
            dataset_root,
            f"{dataset_id}.{modality}.spatial.path",
        )
        truth = record.get("truth")
        if evaluation_mode == "prediction_only" and truth is not None:
            raise LayoutError(
                f"{dataset_id}.{modality} declares truth in a prediction-only campaign"
            )
        if evaluation_mode == "exact_truth":
            if not isinstance(truth, Mapping) or not isinstance(truth.get("path"), str):
                raise LayoutError(
                    f"{dataset_id}.{modality} lacks exact truth for an exact-truth campaign"
                )
            truth_path = project_path(
                str(truth["path"]), f"{dataset_id}.{modality}.truth.path"
            )
            require_prefix(
                truth_path,
                dataset_root / "truth",
                f"{dataset_id}.{modality}.truth.path",
            )

    validation = config.get("validation", {})
    if not isinstance(validation, Mapping):
        raise LayoutError(f"{dataset_id}.validation must be a mapping")
    for key_path, declared in iter_declared_paths(validation, ("validation",)):
        relative = project_path(declared, f"{dataset_id}.{'.'.join(key_path)}")
        require_prefix(
            relative,
            dataset_root / "validation",
            f"{dataset_id}.{'.'.join(key_path)}",
        )


def validate_experiment(
    path: Path,
    registry_path: Path,
    allow_existing_results: bool = False,
) -> Path:
    config = load_yaml(path)
    output_root = canonical_result_root(
        str(config.get("output_root", "")),
        f"{path}.output_root",
    )
    run_group = str(config.get("run_group", ""))
    if not run_group or Path(run_group).name != run_group:
        raise LayoutError(f"{path}.run_group must be one path-safe campaign name")
    batch_dir = ROOT / output_root / run_group
    if batch_dir.exists() and any(batch_dir.iterdir()) and not allow_existing_results:
        raise LayoutError(
            f"Campaign destination already contains files: {batch_dir.relative_to(ROOT)}"
        )

    mode = str(config.get("evaluation_mode", "exact_truth"))
    if mode not in {"exact_truth", "prediction_only"}:
        raise LayoutError(f"{path}.evaluation_mode is invalid: {mode}")
    metrics = config.get("metrics", [])
    if not isinstance(metrics, list):
        raise LayoutError(f"{path}.metrics must be a list")
    if mode == "prediction_only" and metrics:
        raise LayoutError(f"{path} requests exact-truth metrics in prediction-only mode")
    if mode == "exact_truth" and not metrics:
        raise LayoutError(f"{path} has no metrics in exact-truth mode")

    method_runs = config.get("method_runs")
    if not isinstance(method_runs, list) or not method_runs:
        raise LayoutError(f"{path} has no method_runs")
    for index, record in enumerate(method_runs):
        if not isinstance(record, Mapping):
            raise LayoutError(f"{path}.method_runs[{index}] is not a mapping")
        relative = project_path(
            str(record.get("config", "")),
            f"{path}.method_runs[{index}].config",
        )
        require_prefix(relative, Path("configs/methods"), "method config")
        if not (ROOT / relative).is_file():
            raise LayoutError(f"Method config is missing: {relative}")

    registry = load_registry(registry_path)
    datasets = config.get("datasets")
    if not isinstance(datasets, list) or not datasets:
        raise LayoutError(f"{path} has no datasets")
    standardized_reference_scope = output_root.parts[1] in {
        "sensitivity",
        "external_validation",
        "real_spatial",
    }
    for dataset_id in datasets:
        if not isinstance(dataset_id, str) or Path(dataset_id).name != dataset_id:
            raise LayoutError(f"Invalid dataset ID in {path}: {dataset_id!r}")
        validate_dataset(dataset_id, registry, mode, standardized_reference_scope)
    return batch_dir


def validate_source_manifest(path: Path) -> None:
    config = load_yaml(path)
    raw = project_path(str(config.get("raw_directory", "")), f"{path}.raw_directory")
    work = project_path(str(config.get("work_directory", "")), f"{path}.work_directory")
    processed = project_path(
        str(config.get("processed_directory", "")),
        f"{path}.processed_directory",
    )
    require_prefix(raw, Path("data/raw/sources"), f"{path}.raw_directory")
    require_prefix(work, Path("data/work/downloads"), f"{path}.work_directory")
    require_prefix(processed, Path("data/processed/shapemix"), f"{path}.processed_directory")
    resources = config.get("resources")
    if resources is not None:
        if not isinstance(resources, list) or not resources:
            raise LayoutError(f"{path}.resources must be a nonempty list")
        destinations: set[Path] = set()
        for index, record in enumerate(resources):
            if not isinstance(record, Mapping):
                raise LayoutError(f"{path}.resources[{index}] is not a mapping")
            destination = project_path(
                str(record.get("destination", "")),
                f"{path}.resources[{index}].destination",
            )
            if destination in destinations:
                raise LayoutError(f"Duplicate resource destination in {path}: {destination}")
            destinations.add(destination)
            expected = record.get("expected_bytes")
            if not isinstance(expected, int) or isinstance(expected, bool) or expected <= 0:
                raise LayoutError(f"Invalid expected_bytes for {path}.resources[{index}]")


def validate_completed_campaign(path: Path) -> None:
    relative = path.resolve().relative_to(ROOT.resolve())
    if len(relative.parts) < 3 or relative.parts[:1] != ("results",):
        raise LayoutError(f"Campaign result path is noncanonical: {relative}")
    if relative.parts[1] not in ALLOWED_RESULT_SCOPES:
        raise LayoutError(f"Campaign result scope is noncanonical: {relative}")
    run_dirs = [child for child in path.iterdir() if (child / "run.yaml").is_file()]
    for run_dir in run_dirs:
        missing = sorted(item for item in REQUIRED_RUN_FILES if not (run_dir / item).is_file())
        if missing:
            raise LayoutError(f"{run_dir.relative_to(ROOT)} is missing {missing}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-config", type=Path, action="append", default=[])
    parser.add_argument("--source-config", type=Path, action="append", default=[])
    parser.add_argument("--campaign-results", type=Path, action="append", default=[])
    parser.add_argument("--registry", type=Path, default=ROOT / "data/registry/datasets.yaml")
    parser.add_argument("--allow-existing-results", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.experiment_config and not args.source_config and not args.campaign_results:
        raise LayoutError("Declare at least one experiment, source, or result path")
    for path in args.source_config:
        validate_source_manifest(path)
        print(f"source layout passed: {path}")
    for path in args.experiment_config:
        batch_dir = validate_experiment(
            path,
            args.registry,
            allow_existing_results=args.allow_existing_results,
        )
        print(f"experiment layout passed: {path} -> {batch_dir.relative_to(ROOT)}")
    for path in args.campaign_results:
        validate_completed_campaign(path)
        print(f"result layout passed: {path}")


if __name__ == "__main__":
    main()
