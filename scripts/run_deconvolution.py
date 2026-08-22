#!/usr/bin/env python
from __future__ import annotations

import argparse
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Optional, Union

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data import load_deconvolution_input
from deconvatac.data.registry import get_dataset_config
from deconvatac.metrics import jsd, rmse
from deconvatac.methods import get_method


METRICS = {
    "rmse": rmse,
    "jsd": jsd,
}
DEFAULT_METHOD_CONFIG_DIR = ROOT / "configs" / "methods"


def read_yaml(path: Optional[Union[str, Path]]) -> dict[str, Any]:
    if path is None:
        return {}
    with Path(path).open() as handle:
        return yaml.safe_load(handle) or {}


def build_run_id(dataset: str, modality: str, feature_set: str, method: str) -> str:
    date = datetime.now().strftime("%Y-%m-%d")
    return f"{date}_{dataset}_{modality}_{feature_set}_{method}"


def format_run_id(template: str, dataset: str, modality: str, feature_set: str, method: str) -> str:
    return template.format(
        dataset=dataset,
        modality=modality,
        feature_set=feature_set,
        method=method,
    )


def write_environment(output_dir: Path) -> None:
    lines = [
        f"python={sys.version.replace(chr(10), ' ')}",
        f"platform={platform.platform()}",
    ]
    packages = {
        "deconvatac": "deconvatac",
        "anndata": "anndata",
        "muon": "muon",
        "numpy": "numpy",
        "pandas": "pandas",
        "scipy": "scipy",
        "yaml": "yaml",
    }
    for package, module_name in packages.items():
        try:
            module = __import__(module_name)
            version = getattr(module, "__version__", "unknown")
        except Exception as exc:
            version = f"unavailable: {exc}"
        lines.append(f"{package}={version}")
    (output_dir / "environment.txt").write_text("\n".join(lines) + "\n")


def validate_method_config(method: str, method_config: dict[str, Any]) -> None:
    configured_method = method_config.get("method")
    if configured_method is not None and configured_method != method:
        raise ValueError(f"Config method '{configured_method}' does not match method '{method}'.")


def run_one(
    dataset: str,
    modality: str,
    feature_set: str,
    method: str,
    method_config: Optional[dict[str, Any]] = None,
    registry: Union[str, Path] = ROOT / "data" / "registry" / "datasets.yaml",
    output_root: Union[str, Path] = ROOT / "results",
    run_id: Optional[str] = None,
    overwrite: bool = False,
) -> Path:
    method_config = method_config or {}
    validate_method_config(method, method_config)

    resolved_run_id = run_id or build_run_id(dataset, modality, feature_set, method)
    output_dir = Path(output_root or ROOT / "results") / resolved_run_id
    if output_dir.exists() and not overwrite:
        raise FileExistsError(f"{output_dir} already exists. Use --overwrite to reuse it.")
    output_dir.mkdir(parents=True, exist_ok=True)

    data = load_deconvolution_input(
        dataset_id=dataset,
        modality=modality,
        feature_set=feature_set,
        registry_path=registry,
        project_root=ROOT,
        output_dir=output_dir,
    )

    method_cls = get_method(method)
    method = method_cls(**method_config.get("params", {}))
    result = method.run(data)

    extra_metadata = {
        "run_id": resolved_run_id,
        "method_config": method_config,
        "inputs": {
            "dataset_id": dataset,
            "modality": modality,
            "feature_set": feature_set,
            "registry": str(registry),
        },
    }
    result.write(output_dir, extra_metadata=extra_metadata)

    inputs_metadata = {
        "dataset_id": dataset,
        "modality": modality,
        "feature_set": feature_set,
        "registry": str(registry),
        "dataset_config": data.metadata.get("dataset_config", {}),
    }
    if data.truth is not None:
        truth_path = output_dir / "results" / "truth.csv"
        data.truth.to_csv(truth_path)
        inputs_metadata["truth"] = "results/truth.csv"

    with (output_dir / "inputs.yaml").open("w") as handle:
        yaml.safe_dump(inputs_metadata, handle, sort_keys=False)
    write_environment(output_dir)

    print(output_dir)
    return output_dir


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def feature_sets_for_modality(feature_sets: Any, modality: str, modality_config: dict[str, Any]) -> list[str]:
    if feature_sets is None:
        return list(modality_config.get("feature_sets", {}))
    if isinstance(feature_sets, dict):
        return _as_list(feature_sets.get(modality))
    return _as_list(feature_sets)


def resolve_method_config(method: str, experiment_config: dict[str, Any]) -> dict[str, Any]:
    method_configs = experiment_config.get("method_configs", {})
    method_spec = method_configs.get(method)

    if isinstance(method_spec, dict):
        return method_spec
    if isinstance(method_spec, str):
        return read_yaml(resolve_config_path(method_spec))

    default_config = DEFAULT_METHOD_CONFIG_DIR / f"{method}.yaml"
    if default_config.exists():
        return read_yaml(default_config)
    return {}


def resolve_config_path(path: Union[str, Path]) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return ROOT / path


def iter_experiment_jobs(
    experiment_config: dict[str, Any],
    registry: Union[str, Path],
) -> Iterable[dict[str, str]]:
    datasets = _as_list(experiment_config.get("datasets"))
    methods = _as_list(experiment_config.get("methods"))
    configured_modalities = experiment_config.get("modalities")
    configured_feature_sets = experiment_config.get("feature_sets")
    skip_missing = experiment_config.get("skip_missing", True)

    if not datasets:
        raise ValueError("Experiment config must define at least one dataset.")
    if not methods:
        raise ValueError("Experiment config must define at least one method.")

    for dataset in datasets:
        dataset_config = get_dataset_config(dataset, registry_path=registry, project_root=ROOT)
        modality_configs = dataset_config.get("modalities", {})
        modalities = _as_list(configured_modalities) or list(modality_configs)

        for modality in modalities:
            if modality not in modality_configs:
                if skip_missing:
                    continue
                raise KeyError(f"Dataset '{dataset}' does not define modality '{modality}'.")

            modality_config = modality_configs[modality]
            feature_sets = feature_sets_for_modality(configured_feature_sets, modality, modality_config)
            if not feature_sets:
                raise ValueError(f"No feature sets selected for dataset '{dataset}' modality '{modality}'.")

            available_feature_sets = set(modality_config.get("feature_sets", {}))
            for feature_set in feature_sets:
                if feature_set not in available_feature_sets:
                    if skip_missing:
                        continue
                    raise KeyError(
                        f"Dataset '{dataset}' modality '{modality}' does not define feature set '{feature_set}'."
                    )

                for method in methods:
                    yield {
                        "dataset": str(dataset),
                        "modality": str(modality),
                        "feature_set": str(feature_set),
                        "method": str(method),
                    }


def read_run_metadata(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "run.yaml"
    if not path.exists():
        return {}
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def evaluate_run(run_dir: Path, metrics: list[str]) -> list[dict[str, Any]]:
    metadata = read_run_metadata(run_dir)
    predicted = pd.read_csv(run_dir / "results" / "proportions.csv", index_col=0)
    truth = pd.read_csv(run_dir / "results" / "truth.csv", index_col=0)

    rows = []
    for metric_name in metrics:
        rows.append(
            {
                "run_id": metadata.get("run_id", run_dir.name),
                "dataset_id": metadata.get("dataset_id"),
                "modality": metadata.get("modality"),
                "feature_set": metadata.get("feature_set"),
                "method": metadata.get("method"),
                "status": "success",
                "metric": metric_name,
                "value": METRICS[metric_name](truth, predicted),
                "run_dir": str(run_dir),
                "error": None,
            }
        )
    return rows


def write_comparison(
    output_path: Path,
    successful_runs: list[Path],
    failures: list[dict[str, Any]],
    metrics: list[str],
) -> None:
    rows = []
    for run_dir in successful_runs:
        rows.extend(evaluate_run(run_dir, metrics))
    rows.extend(failures)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_path, index=False)


def run_experiment(
    experiment_config_path: Union[str, Path],
    registry: Union[str, Path],
    output_root_override: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
) -> Path:
    experiment_config_path = resolve_config_path(experiment_config_path)
    experiment_config = read_yaml(experiment_config_path)

    run_group = experiment_config.get("run_group")
    if run_group is None:
        timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        run_group = f"{timestamp}_{experiment_config_path.stem}"

    output_root = Path(output_root_override or experiment_config.get("output_root", ROOT / "results"))
    if not output_root.is_absolute():
        output_root = ROOT / output_root
    batch_dir = output_root / run_group
    batch_dir.mkdir(parents=True, exist_ok=True)

    metrics = _as_list(experiment_config.get("metrics")) or ["rmse", "jsd"]
    unknown_metrics = sorted(set(metrics).difference(METRICS))
    if unknown_metrics:
        raise KeyError(f"Unknown metric(s): {', '.join(unknown_metrics)}. Available: {', '.join(sorted(METRICS))}")

    run_id_template = experiment_config.get(
        "run_id_template",
        "{dataset}__{modality}__{feature_set}__{method}",
    )
    continue_on_error = experiment_config.get("continue_on_error", False)
    jobs = list(iter_experiment_jobs(experiment_config, registry=registry))
    if not jobs:
        raise ValueError("Experiment config did not produce any runnable jobs.")

    successful_runs: list[Path] = []
    manifest_rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []

    for job in jobs:
        run_id = format_run_id(run_id_template, **job)
        method_config = resolve_method_config(job["method"], experiment_config)
        try:
            run_dir = run_one(
                dataset=job["dataset"],
                modality=job["modality"],
                feature_set=job["feature_set"],
                method=job["method"],
                method_config=method_config,
                registry=registry,
                output_root=batch_dir,
                run_id=run_id,
                overwrite=overwrite or experiment_config.get("overwrite", False),
            )
            successful_runs.append(run_dir)
            manifest_rows.append({**job, "run_id": run_id, "run_dir": str(run_dir), "status": "success"})
        except Exception as exc:
            failure = {
                "run_id": run_id,
                "dataset_id": job["dataset"],
                "modality": job["modality"],
                "feature_set": job["feature_set"],
                "method": job["method"],
                "run_dir": str(batch_dir / run_id),
                "status": "failed",
                "metric": None,
                "value": None,
                "error": str(exc),
            }
            failures.append(failure)
            manifest_rows.append({**job, "run_id": run_id, "run_dir": str(batch_dir / run_id), "status": "failed"})
            if not continue_on_error:
                pd.DataFrame(manifest_rows).to_csv(batch_dir / "runs.csv", index=False)
                pd.DataFrame(failures).to_csv(batch_dir / "failures.csv", index=False)
                raise

    pd.DataFrame(manifest_rows).to_csv(batch_dir / "runs.csv", index=False)
    if failures:
        pd.DataFrame(failures).to_csv(batch_dir / "failures.csv", index=False)

    comparison_path = experiment_config.get("comparison_output")
    if comparison_path is None:
        comparison_path = batch_dir / "comparison.csv"
    else:
        comparison_path = Path(comparison_path)
        if not comparison_path.is_absolute():
            comparison_path = ROOT / comparison_path

    write_comparison(comparison_path, successful_runs=successful_runs, failures=failures, metrics=metrics)
    print(comparison_path)
    return comparison_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run one deconvolution method or an experiment config through the unified interface."
    )
    parser.add_argument("--experiment-config", "--experiment", dest="experiment_config")
    parser.add_argument("--dataset")
    parser.add_argument("--modality", choices=["atac", "rna"])
    parser.add_argument("--feature-set", default="all")
    parser.add_argument("--method")
    parser.add_argument("--config", help="Method config for single-run mode.")
    parser.add_argument("--registry", default=str(ROOT / "data" / "registry" / "datasets.yaml"))
    parser.add_argument("--output-root")
    parser.add_argument("--run-id")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.experiment_config is not None:
        run_experiment(
            experiment_config_path=args.experiment_config,
            registry=args.registry,
            output_root_override=args.output_root,
            overwrite=args.overwrite,
        )
        return

    missing = [name for name in ("dataset", "modality", "method") if getattr(args, name) is None]
    if missing:
        parser.error(
            "single-run mode requires --dataset, --modality, and --method; "
            "batch mode requires --experiment-config."
        )

    method_config = read_yaml(args.config)
    run_one(
        dataset=args.dataset,
        modality=args.modality,
        feature_set=args.feature_set,
        method=args.method,
        method_config=method_config,
        registry=args.registry,
        output_root=args.output_root or ROOT / "results",
        run_id=args.run_id,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
