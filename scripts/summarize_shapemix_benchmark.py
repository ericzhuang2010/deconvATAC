#!/usr/bin/env python
"""Summarize the frozen ShapeMix paired, nested-seed benchmark.

This is intentionally a protocol-specific reporter, not a generic leaderboard.
Its resampling units are the five outer reference/test splits.  The two inner
mixtures are averaged inside each outer split, and spots are never treated as
replicates.  Consequently, all intervals and tests emitted here describe
conditional resampling variability in one PBMC donor; they are not donor-level
or population-level uncertainty.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data.registry import get_dataset_config
from deconvatac.metrics import align_proportions, evaluate_proportion_metric


PROTOCOL_VERSION = 1
PRIMARY_ARMS = ("shapemix_length", "shapemix_count_only")
PRIMARY_METRICS = ("rmse_v1", "jsd_v2")
PRIMARY_CONDITIONS = ("observed_abundance", "equal_celltype")
PRIMARY_OUTER_SPLIT_SEEDS = (1103, 2203, 3301, 4409, 5501)
PRIMARY_INNER_MIXTURE_SEEDS = (101, 211)
FROZEN_CELL_TYPES = (
    "CD14 Mono",
    "CD4 Naive",
    "CD8 Naive",
    "CD4 TCM",
    "CD16 Mono",
    "NK",
    "CD8 TEM_1",
    "CD8 TEM_2",
    "Intermediate B",
    "Memory B",
    "CD4 TEM",
    "cDC",
    "Treg",
    "gdT",
    "MAIT",
    "Naive B",
)
RARE_REFERENCE_TYPES = ("cDC", "Treg", "gdT", "MAIT", "Naive B")
TRUE_PRESENCE_THRESHOLD = 0.0
PREDICTED_PRESENCE_THRESHOLD = 0.01
ROW_SUM_ATOL = 1e-6
BOOTSTRAP_SEED_TUPLE = (20260822, 9001)
BOOTSTRAP_REPLICATES = 10_000
REPORTING_SCOPE = (
    "one-donor conditional resampling variability; not biological replication, "
    "donor-level uncertainty, or population-level generalization"
)
FROZEN_SHAPEMIX_PARAMS = {
    "total_likelihood": "negative_binomial",
    "conditional_shape_likelihood": "multinomial",
    "signature_rate_pseudocount": 0.5,
    "signature_shape_concentration": 1.0,
    "exposure_mode": "absorbed_in_abundance",
    "dispersion_mode": "reference_crossfit_global_scaled_by_abundance",
    "dispersion_crossfit_folds": 2,
    "dispersion_alpha_floor": 1.0e-8,
    "background_mode": "none",
    "abundance_prior": "gamma",
    "abundance_prior_shape": 2.0,
    "abundance_prior_rate": 1.0,
    "optimizer": "adam",
    "learning_rate": 0.03,
    "max_steps": 2000,
    "patience": 100,
    "tolerance": 1.0e-5,
    "restarts": 3,
    "spot_batch_size": 64,
    "peak_chunk_size": 512,
    "seed": 0,
    "device": "cpu",
    "dtype": "float32",
}


PAIR_KEY = ("outer_split_seed", "inner_mixture_seed", "condition")
RUN_ID_COLUMNS = (
    "source_run_group",
    "run_id",
    "dataset_id",
    "outer_split_seed",
    "inner_mixture_seed",
    "condition",
    "method_run_id",
)
PROVENANCE_COLUMNS = (
    "experiment_config_path",
    "experiment_config_source_sha256",
    "experiment_config_resolved_sha256",
    "benchmark_protocol_path",
    "benchmark_protocol_sha256",
    "registry_path",
    "registry_source_sha256",
    "dataset_config_sha256",
    "git_commit",
    "git_worktree_dirty",
    "code_manifest_sha256",
    "output_manifest_path",
    "output_manifest_sha256",
    "wall_runtime_seconds",
    "peak_rss_bytes",
    "peak_rss_mb",
    "rss_measurement",
    "performance_scope",
    "execution_action",
    "shard_index",
    "shard_count",
)


@dataclass(frozen=True)
class Design:
    """Frozen design metadata resolved from a dataset configuration."""

    dataset_id: str
    benchmark_scope: str
    condition: str
    outer_split_seed: int
    inner_mixture_seed: int
    cell_types: tuple[str, ...]


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, Mapping):
        raise TypeError(f"{path} must contain a YAML mapping.")
    return dict(value)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _require_sha256(value: Any, field: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError(f"{field} must be a 64-character SHA-256 digest.")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ValueError(f"{field} must be hexadecimal.") from exc
    return value.lower()


def _canonical_mapping_sha256(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _resolve_contained_path(root: Path, raw_path: Any, field: str) -> tuple[Path, str]:
    text = _require_string(raw_path, field)
    relative = Path(text)
    if relative.is_absolute() or relative in {Path("."), Path("..")}:
        raise ValueError(f"{field} must be a relative path inside {root}.")
    resolved = (root / relative).resolve()
    try:
        canonical = str(resolved.relative_to(root.resolve()))
    except ValueError as exc:
        raise ValueError(f"{field} escapes {root}: {text!r}.") from exc
    return resolved, canonical


def _resolve_project_path(project_root: Path, raw_path: Any, field: str) -> Path:
    text = _require_string(raw_path, field)
    path = Path(text)
    if not path.is_absolute():
        path = project_root / path
    resolved = path.resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{field} does not resolve to a regular file: {resolved}")
    return resolved


def _manifest_scalar_equal(observed: Any, expected: Any) -> bool:
    if isinstance(expected, (float, np.floating)) and pd.isna(expected):
        return observed is None
    if isinstance(observed, (float, np.floating)) and pd.isna(observed):
        return expected is None
    if (
        isinstance(observed, (int, float, np.integer, np.floating))
        and not isinstance(observed, (bool, np.bool_))
        and isinstance(expected, (int, float, np.integer, np.floating))
        and not isinstance(expected, (bool, np.bool_))
    ):
        return bool(np.isclose(float(observed), float(expected), rtol=1e-12, atol=0.0))
    return observed == expected


def _validate_output_hash_manifest(
    run_dir: Path,
    output_spec: Mapping[str, Any],
    run_metadata: Mapping[str, Any],
) -> tuple[str, str, int]:
    if output_spec.get("schema_version") != 1:
        raise ValueError(f"Run {run_dir.name!r} output_manifest.schema_version must be 1.")
    manifest_path, relative_manifest_path = _resolve_contained_path(
        run_dir,
        output_spec.get("path"),
        "output_manifest.path",
    )
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Run {run_dir.name!r} is missing {manifest_path}.")
    recorded_manifest_hash = _require_sha256(
        output_spec.get("sha256"), "output_manifest.sha256"
    )
    observed_manifest_hash = _sha256_file(manifest_path)
    if observed_manifest_hash != recorded_manifest_hash:
        raise ValueError(
            f"Run {run_dir.name!r} output manifest hash is stale: "
            f"recorded={recorded_manifest_hash}, observed={observed_manifest_hash}."
        )

    manifest = _read_yaml(manifest_path)
    if manifest.get("schema_version") != 1 or manifest.get("algorithm") != "sha256":
        raise ValueError(f"{manifest_path} must declare schema_version 1 and algorithm sha256.")
    if manifest.get("run_metadata_hash_encoding") != "canonical_json_sorted_keys_compact_utf8":
        raise ValueError(f"{manifest_path} has an unsupported run metadata hash encoding.")
    recorded_run_hash = _require_sha256(
        manifest.get("run_metadata_sha256"), "output manifest run_metadata_sha256"
    )
    metadata_without_manifest = copy.deepcopy(dict(run_metadata))
    metadata_without_manifest.pop("output_manifest", None)
    if _canonical_mapping_sha256(metadata_without_manifest) != recorded_run_hash:
        raise ValueError(f"{manifest_path} does not authenticate the current run.yaml metadata.")
    exclusions = manifest.get("exclusions")
    if not isinstance(exclusions, list) or set(exclusions) != {
        "run.yaml",
        relative_manifest_path,
    }:
        raise ValueError(
            f"{manifest_path} exclusions must be run.yaml and {relative_manifest_path}."
        )
    files = manifest.get("files")
    if not isinstance(files, Mapping):
        raise ValueError(f"{manifest_path} files must be a path-to-SHA mapping.")
    if list(files) != sorted(files):
        raise ValueError(f"{manifest_path} file keys must be sorted deterministically.")

    expected_paths: set[str] = set()
    for path in run_dir.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Run output trees must not contain symlinks: {path}")
        if path.is_file():
            relative = str(path.relative_to(run_dir))
            if relative not in exclusions:
                expected_paths.add(relative)
    if set(files) != expected_paths:
        missing = sorted(expected_paths.difference(files))
        extra = sorted(set(files).difference(expected_paths))
        raise ValueError(
            f"{manifest_path} does not cover the exact output file set; "
            f"missing={missing!r}, extra={extra!r}."
        )
    for relative, expected_hash in files.items():
        expected_hash = _require_sha256(expected_hash, f"output hash for {relative}")
        path, canonical = _resolve_contained_path(run_dir, relative, "output manifest file")
        if canonical != relative or not path.is_file():
            raise ValueError(f"Output manifest path is not canonical: {relative!r}.")
        observed_hash = _sha256_file(path)
        if observed_hash != expected_hash:
            raise ValueError(
                f"Run {run_dir.name!r} output hash mismatch for {relative!r}: "
                f"recorded={expected_hash}, observed={observed_hash}."
            )
    return relative_manifest_path, recorded_manifest_hash, len(files)


def validate_run_provenance(
    record: Mapping[str, Any],
    *,
    project_root: Path,
) -> dict[str, Any]:
    """Validate one run's nested provenance and every hashed output byte."""

    run_dir = Path(str(record["run_dir"]))
    run_yaml_path = run_dir / "run.yaml"
    if not run_yaml_path.is_file():
        raise FileNotFoundError(f"Run is missing required metadata: {run_yaml_path}")
    metadata = _read_yaml(run_yaml_path)
    if metadata.get("status") != record["status"]:
        raise ValueError(
            f"{run_yaml_path} status disagrees with runs.csv: "
            f"{metadata.get('status')!r} != {record['status']!r}."
        )
    provenance = metadata.get("execution_provenance")
    if not isinstance(provenance, Mapping) or provenance.get("schema_version") != 1:
        raise ValueError(f"{run_yaml_path} requires execution_provenance schema version 1.")

    experiment = provenance.get("experiment_config")
    protocol = provenance.get("benchmark_protocol")
    registry = provenance.get("registry")
    code = provenance.get("code")
    shard = provenance.get("shard")
    if not all(
        isinstance(value, Mapping) for value in (experiment, protocol, registry, code)
    ):
        raise ValueError(
            f"{run_yaml_path} requires experiment_config, benchmark_protocol, registry, and code provenance."
        )
    experiment = dict(experiment)
    protocol = dict(protocol)
    registry = dict(registry)
    code = dict(code)
    if shard is not None and not isinstance(shard, Mapping):
        raise ValueError(f"{run_yaml_path} execution_provenance.shard must be a mapping.")
    shard = None if shard is None else dict(shard)

    experiment_path = _require_string(experiment.get("path"), "experiment_config.path")
    experiment_source_sha256 = _require_sha256(
        experiment.get("source_sha256"), "experiment_config.source_sha256"
    )
    experiment_resolved_sha256 = _require_sha256(
        experiment.get("resolved_sha256"), "experiment_config.resolved_sha256"
    )
    protocol_path = _require_string(protocol.get("path"), "benchmark_protocol.path")
    protocol_sha256 = _require_sha256(
        protocol.get("source_sha256"), "benchmark_protocol.source_sha256"
    )
    registry_path = _require_string(registry.get("path"), "registry.path")
    registry_sha256 = _require_sha256(registry.get("source_sha256"), "registry.source_sha256")
    dataset_config_sha256 = _require_sha256(
        metadata.get("dataset_config_sha256"), "run.yaml dataset_config_sha256"
    )
    git_commit = _require_string(code.get("git_commit"), "code.git_commit")
    if len(git_commit) != 40:
        raise ValueError("code.git_commit must be a full 40-character commit hash.")
    try:
        int(git_commit, 16)
    except ValueError as exc:
        raise ValueError("code.git_commit must be hexadecimal.") from exc
    worktree_dirty = code.get("git_worktree_dirty")
    if not isinstance(worktree_dirty, (bool, np.bool_)):
        raise TypeError("code.worktree_dirty must be boolean.")
    code_files = code.get("files")
    if not isinstance(code_files, Mapping) or not code_files:
        raise ValueError("code.files must be a nonempty path-to-SHA mapping.")
    normalized_code_files: dict[str, str] = {}
    for path, digest in code_files.items():
        normalized_code_files[_require_string(path, "code.files path")] = _require_sha256(
            digest, f"code.files[{path!r}]"
        )
    code_manifest_sha256 = _require_sha256(
        code.get("manifest_sha256"), "code.manifest_sha256"
    )
    if _canonical_mapping_sha256(normalized_code_files) != code_manifest_sha256:
        raise ValueError("code.manifest_sha256 does not match the canonical code.files mapping.")

    if shard is None:
        shard_index = 0
        shard_count = 1
        flattened_shard_index = None
        flattened_shard_count = None
    else:
        shard_index = _as_nonnegative_int(
            shard.get("index"), "execution_provenance.shard.index"
        )
        shard_count = _as_nonnegative_int(
            shard.get("count"), "execution_provenance.shard.count"
        )
        if shard_count < 1 or shard_index >= shard_count:
            raise ValueError("Shard provenance requires 0 <= index < count and count >= 1.")
        if shard.get("index_base") != 0:
            raise ValueError("Shard provenance must declare zero-based indexing.")
        if shard.get("assignment") != "sorted_dataset_modality_feature_set_modulo":
            raise ValueError("Shard provenance has an unsupported assignment rule.")
        if shard.get("unit_key_fields") != ["dataset", "modality", "feature_set"]:
            raise ValueError("Shard provenance has unsupported unit key fields.")
        flattened_shard_index = shard_index
        flattened_shard_count = shard_count

    output_spec = metadata.get("output_manifest")
    if not isinstance(output_spec, Mapping):
        raise ValueError(f"{run_yaml_path} requires output_manifest metadata.")
    output_path, output_sha256, output_file_count = _validate_output_hash_manifest(
        run_dir, output_spec, metadata
    )

    performance = metadata.get("performance")
    if not isinstance(performance, Mapping):
        raise ValueError(f"{run_yaml_path} requires performance metadata.")
    wall_runtime = performance.get("wall_runtime_seconds")
    peak_rss_bytes = performance.get("peak_rss_bytes")
    peak_rss_mb = performance.get("peak_rss_mb")
    if not isinstance(wall_runtime, (int, float)) or not np.isfinite(wall_runtime):
        raise ValueError("performance.wall_runtime_seconds must be finite.")
    if not isinstance(peak_rss_bytes, int) or isinstance(peak_rss_bytes, bool):
        raise ValueError("performance.peak_rss_bytes must be an integer.")
    if not isinstance(peak_rss_mb, (int, float)) or not np.isfinite(peak_rss_mb):
        raise ValueError("performance.peak_rss_mb must be finite.")
    if wall_runtime < 0 or peak_rss_bytes < 0 or peak_rss_mb < 0:
        raise ValueError("Performance measurements must be nonnegative.")
    if not np.isclose(peak_rss_mb, peak_rss_bytes / (1024.0**2), rtol=0.0, atol=1e-9):
        raise ValueError("performance.peak_rss_mb is inconsistent with peak_rss_bytes.")
    rss_payload = performance.get("rss_measurement")
    if not isinstance(rss_payload, Mapping):
        raise ValueError("performance.rss_measurement must be a mapping.")
    rss_measurement = json.dumps(
        rss_payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )
    performance_scope = _require_string(performance.get("scope"), "performance.scope")
    execution_action = _require_string(metadata.get("execution_action"), "execution_action")
    if execution_action not in {"executed", "resumed"}:
        raise ValueError(f"Unsupported execution_action {execution_action!r}.")

    values = {
        "experiment_config_path": experiment_path,
        "experiment_config_source_sha256": experiment_source_sha256,
        "experiment_config_resolved_sha256": experiment_resolved_sha256,
        "benchmark_protocol_path": protocol_path,
        "benchmark_protocol_sha256": protocol_sha256,
        "registry_path": registry_path,
        "registry_source_sha256": registry_sha256,
        "dataset_config_sha256": dataset_config_sha256,
        "git_commit": git_commit,
        "git_worktree_dirty": bool(worktree_dirty),
        "code_manifest_sha256": code_manifest_sha256,
        "output_manifest_path": output_path,
        "output_manifest_sha256": output_sha256,
        "output_manifest_file_count": output_file_count,
        "wall_runtime_seconds": float(wall_runtime),
        "peak_rss_bytes": peak_rss_bytes,
        "peak_rss_mb": float(peak_rss_mb),
        "rss_measurement": rss_measurement,
        "performance_scope": performance_scope,
        "execution_action": execution_action,
        "shard_index": shard_index,
        "shard_count": shard_count,
        "_code_files": normalized_code_files,
    }
    for field in PROVENANCE_COLUMNS:
        if field not in record:
            raise ValueError(f"runs.csv is missing flattened provenance field {field!r}.")
        observed = record[field]
        expected = values[field]
        if field == "shard_index":
            expected = flattened_shard_index
        elif field == "shard_count":
            expected = flattened_shard_count
        if field == "git_worktree_dirty" and isinstance(observed, str):
            observed = observed.strip().lower() == "true"
        if field in {"shard_index", "shard_count"} and pd.isna(observed):
            observed = None
        elif field in {"shard_index", "shard_count", "peak_rss_bytes"}:
            observed = int(observed)
        if field in {"wall_runtime_seconds", "peak_rss_mb"}:
            observed = float(observed)
        if not _manifest_scalar_equal(observed, expected):
            raise ValueError(
                f"Run {record['run_id']!r} flattened {field} disagrees with run.yaml: "
                f"{observed!r} != {expected!r}."
            )
    return values


def collect_and_validate_provenance(
    records: pd.DataFrame,
    *,
    project_root: Path,
) -> pd.DataFrame:
    """Validate per-run hashes and one shared campaign provenance."""

    rows: list[dict[str, Any]] = []
    for record in records.to_dict("records"):
        values = validate_run_provenance(record, project_root=project_root)
        rows.append(
            {
                **{field: record[field] for field in RUN_ID_COLUMNS},
                "status": record["status"],
                **{key: value for key, value in values.items() if not key.startswith("_")},
                "_code_files": values["_code_files"],
            }
        )
    frame = pd.DataFrame(rows)
    shared_fields = [
        "experiment_config_path",
        "experiment_config_source_sha256",
        "experiment_config_resolved_sha256",
        "benchmark_protocol_path",
        "benchmark_protocol_sha256",
        "registry_path",
        "registry_source_sha256",
        "git_commit",
        "git_worktree_dirty",
        "code_manifest_sha256",
        "shard_count",
    ]
    mismatched = [field for field in shared_fields if frame[field].nunique(dropna=False) != 1]
    if mismatched:
        raise ValueError(
            "Combined run groups do not share one experiment/protocol/code provenance: "
            + ", ".join(mismatched)
        )

    shard_count = int(frame["shard_count"].iloc[0])
    shard_indices = set(frame["shard_index"].astype(int))
    if shard_indices != set(range(shard_count)):
        raise ValueError(
            f"Combined campaign does not contain every declared shard; "
            f"expected={list(range(shard_count))!r}, observed={sorted(shard_indices)!r}."
        )
    pairing_records = records.drop(columns=["shard_index"], errors="ignore").merge(
        frame[["run_id", "shard_index"]],
        on="run_id",
        how="left",
        validate="one_to_one",
    )
    pair_shards = pairing_records.groupby(list(PAIR_KEY), sort=False).agg(
        source_groups=("source_run_group", "nunique"),
        shard_indices=("shard_index", "nunique"),
    )
    if (pair_shards != 1).any(axis=None):
        raise ValueError("Both arms of every pair must remain together in one shard/run group.")

    first = frame.iloc[0]
    experiment_path = _resolve_project_path(
        project_root, first["experiment_config_path"], "experiment_config.path"
    )
    if _sha256_file(experiment_path) != first["experiment_config_source_sha256"]:
        raise ValueError("The recorded experiment config source hash no longer matches its file.")
    protocol_path = _resolve_project_path(
        project_root, first["benchmark_protocol_path"], "benchmark_protocol.path"
    )
    if _sha256_file(protocol_path) != first["benchmark_protocol_sha256"]:
        raise ValueError("The recorded benchmark protocol hash no longer matches its file.")
    registry_path = _resolve_project_path(
        project_root, first["registry_path"], "registry.path"
    )
    if _sha256_file(registry_path) != first["registry_source_sha256"]:
        raise ValueError("The recorded dataset registry hash no longer matches its file.")
    if protocol_path.name != "benchmark_protocol.md":
        raise ValueError("Primary ShapeMix provenance must reference benchmark_protocol.md.")
    if experiment_path.name != "shapemix_primary_ablation.yaml":
        raise ValueError("Primary ShapeMix provenance must reference shapemix_primary_ablation.yaml.")

    code_files = first["_code_files"]
    for relative, expected_hash in code_files.items():
        source_path = _resolve_project_path(project_root, relative, f"code file {relative!r}")
        if _sha256_file(source_path) != expected_hash:
            raise ValueError(f"Recorded code hash no longer matches {relative!r}.")
    return frame


def _as_nonnegative_int(value: Any, field: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{field} must be an integer, got {value!r}.")
    result = int(value)
    if result < 0:
        raise ValueError(f"{field} must be nonnegative, got {result}.")
    return result


def _require_string(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a nonempty string.")
    return value.strip()


def _json_list(values: Sequence[str]) -> str:
    return json.dumps(list(values), ensure_ascii=False, separators=(",", ":"))


def _safe_error(value: Any) -> Optional[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    return text or None


def _run_directory(raw_path: Any, run_group: Path, run_id: str) -> Path:
    if raw_path is None or (isinstance(raw_path, float) and math.isnan(raw_path)):
        return (run_group / run_id).resolve()
    path = Path(str(raw_path))
    if not path.is_absolute():
        path = run_group / path
    return path.resolve()


def _config_from_inputs(run_dir: Path) -> Optional[dict[str, Any]]:
    path = run_dir / "inputs.yaml"
    if not path.exists():
        return None
    payload = _read_yaml(path)
    config = payload.get("dataset_config")
    if config is None:
        return None
    if not isinstance(config, Mapping):
        raise TypeError(f"inputs.yaml dataset_config must be a mapping for {run_dir}.")
    return copy.deepcopy(dict(config))


def resolve_design(
    dataset_id: str,
    run_dir: Path,
    *,
    registry: Path,
    project_root: Path = ROOT,
) -> Design:
    """Resolve seeds/condition from metadata, never from a dataset filename."""

    config = _config_from_inputs(run_dir)
    if config is None:
        config = get_dataset_config(
            dataset_id,
            registry_path=registry,
            project_root=project_root,
        )
    configured_id = config.get("dataset_id", dataset_id)
    if configured_id != dataset_id:
        raise ValueError(
            f"Run dataset_id {dataset_id!r} disagrees with dataset config {configured_id!r}."
        )

    simulation = config.get("simulation")
    if not isinstance(simulation, Mapping):
        raise ValueError(f"Dataset {dataset_id!r} lacks a simulation metadata mapping.")
    modality = (config.get("modalities") or {}).get("atac")
    if not isinstance(modality, Mapping):
        raise ValueError(f"Dataset {dataset_id!r} lacks an ATAC modality configuration.")
    truth_spec = modality.get("truth")
    if not isinstance(truth_spec, Mapping):
        raise ValueError(f"Dataset {dataset_id!r} lacks modalities.atac.truth metadata.")
    raw_cell_types = truth_spec.get("cell_types")
    if not isinstance(raw_cell_types, list) or not raw_cell_types:
        raise ValueError(f"Dataset {dataset_id!r} must declare a nonempty truth.cell_types list.")
    cell_types = tuple(_require_string(value, "truth.cell_types entry") for value in raw_cell_types)
    if len(set(cell_types)) != len(cell_types):
        raise ValueError(f"Dataset {dataset_id!r} declares duplicate truth cell types.")

    return Design(
        dataset_id=dataset_id,
        benchmark_scope=_require_string(config.get("benchmark_scope"), "benchmark_scope"),
        condition=_require_string(simulation.get("condition"), "simulation.condition"),
        outer_split_seed=_as_nonnegative_int(
            simulation.get("outer_split_seed"), "simulation.outer_split_seed"
        ),
        inner_mixture_seed=_as_nonnegative_int(
            simulation.get("inner_mixture_seed"), "simulation.inner_mixture_seed"
        ),
        cell_types=cell_types,
    )


def _manifest_failure_lookup(run_group: Path) -> dict[tuple[str, str], str]:
    path = run_group / "failures.csv"
    if not path.exists():
        return {}
    failures = pd.read_csv(path)
    required = {"run_id", "method_run_id", "error"}
    missing = required.difference(failures.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}.")
    result: dict[tuple[str, str], str] = {}
    for row in failures.to_dict("records"):
        key = (str(row["run_id"]), str(row["method_run_id"]))
        if key in result:
            raise ValueError(f"Duplicate failure entry {key!r} in {path}.")
        result[key] = _safe_error(row.get("error")) or "unspecified run failure"
    return result


def load_run_groups(
    run_groups: Sequence[Path],
    *,
    registry: Path,
    project_root: Path = ROOT,
) -> pd.DataFrame:
    """Combine run manifests and attach the frozen nested design metadata."""

    if not run_groups:
        raise ValueError("At least one run group is required.")
    rows: list[dict[str, Any]] = []
    seen_run_ids: set[str] = set()
    for raw_group in run_groups:
        group = Path(raw_group).resolve()
        manifest_path = group / "runs.csv"
        if not manifest_path.exists():
            raise FileNotFoundError(f"Missing run manifest: {manifest_path}")
        manifest = pd.read_csv(manifest_path)
        required = {
            "run_id",
            "dataset_id",
            "method",
            "method_run_id",
            "run_dir",
            "status",
            *PROVENANCE_COLUMNS,
        }
        missing = required.difference(manifest.columns)
        if missing:
            raise ValueError(
                f"{manifest_path} is missing columns: {', '.join(sorted(missing))}."
            )
        failures = _manifest_failure_lookup(group)
        for raw in manifest.to_dict("records"):
            run_id = _require_string(raw.get("run_id"), "run_id")
            if run_id in seen_run_ids:
                raise ValueError(f"Duplicate run_id across run groups: {run_id!r}.")
            seen_run_ids.add(run_id)
            method = _require_string(raw.get("method"), "method")
            method_run_id = _require_string(raw.get("method_run_id"), "method_run_id")
            if method != "shapemix" or method_run_id not in PRIMARY_ARMS:
                raise ValueError(
                    f"Primary ShapeMix summary received unsupported run {run_id!r}: "
                    f"method={method!r}, method_run_id={method_run_id!r}."
                )
            status = _require_string(raw.get("status"), "status").lower()
            if status not in {"success", "failed"}:
                raise ValueError(f"Run {run_id!r} has unsupported status {status!r}.")
            dataset_id = _require_string(raw.get("dataset_id"), "dataset_id")
            run_dir = _run_directory(raw.get("run_dir"), group, run_id)
            design = resolve_design(
                dataset_id,
                run_dir,
                registry=registry,
                project_root=project_root,
            )
            error = _safe_error(raw.get("error"))
            if status == "failed":
                error = error or failures.get((run_id, method_run_id)) or "unspecified run failure"
            rows.append(
                {
                    **raw,
                    "source_run_group": str(group),
                    "run_id": run_id,
                    "dataset_id": dataset_id,
                    "run_dir": str(run_dir),
                    "method": method,
                    "method_run_id": method_run_id,
                    "status": status,
                    "error": error,
                    "benchmark_scope": design.benchmark_scope,
                    "condition": design.condition,
                    "outer_split_seed": design.outer_split_seed,
                    "inner_mixture_seed": design.inner_mixture_seed,
                    "cell_type_universe": _json_list(design.cell_types),
                    "_cell_types": design.cell_types,
                }
            )

    records = pd.DataFrame(rows)
    duplicate_pair_arms = records.duplicated(list(PAIR_KEY) + ["method_run_id"], keep=False)
    if duplicate_pair_arms.any():
        offending = records.loc[
            duplicate_pair_arms, list(PAIR_KEY) + ["method_run_id", "run_id"]
        ]
        raise ValueError(
            "More than one run represents the same nested-seed pair and arm: "
            f"{offending.to_dict('records')!r}."
        )
    return records.sort_values(
        ["condition", "outer_split_seed", "inner_mixture_seed", "method_run_id"]
    ).reset_index(drop=True)


def expected_primary_pairs() -> set[tuple[int, int, str]]:
    return {
        (outer, inner, condition)
        for outer, inner, condition in itertools.product(
            PRIMARY_OUTER_SPLIT_SEEDS,
            PRIMARY_INNER_MIXTURE_SEEDS,
            PRIMARY_CONDITIONS,
        )
    }


def validate_primary_design(records: pd.DataFrame) -> None:
    """Fail closed if the combined manifests are not the frozen primary grid."""

    if records.empty:
        raise ValueError("No ShapeMix runs were found.")
    scopes = set(records["benchmark_scope"])
    if scopes != {"primary"}:
        raise ValueError(f"Primary summary requires benchmark_scope='primary', found {sorted(scopes)!r}.")
    conditions = set(records["condition"])
    if conditions != set(PRIMARY_CONDITIONS):
        raise ValueError(
            f"Primary conditions must be {list(PRIMARY_CONDITIONS)!r}, found {sorted(conditions)!r}."
        )
    universes = {tuple(value) for value in records["_cell_types"]}
    if universes != {FROZEN_CELL_TYPES}:
        raise ValueError("Every primary dataset must declare the frozen ordered 16-type universe.")

    actual_pairs = set(map(tuple, records.loc[:, PAIR_KEY].drop_duplicates().to_numpy()))
    expected = expected_primary_pairs()
    if actual_pairs != expected:
        missing = sorted(expected.difference(actual_pairs))
        extra = sorted(actual_pairs.difference(expected))
        raise ValueError(f"Primary nested-seed grid mismatch; missing={missing!r}, extra={extra!r}.")

    grouped_arms = records.groupby(list(PAIR_KEY), sort=False)["method_run_id"].agg(
        arm_set=lambda values: set(values),
        row_count="size",
    )
    bad = grouped_arms[
        grouped_arms.apply(
            lambda row: row["arm_set"] != set(PRIMARY_ARMS)
            or int(row["row_count"]) != len(PRIMARY_ARMS),
            axis=1,
        )
    ]
    if not bad.empty:
        raise ValueError(
            "Every nested-seed dataset must contain exactly both ShapeMix arms; "
            f"mismatches={bad.to_dict()!r}."
        )


def _parse_method_config(value: Any, *, run_id: str) -> dict[str, Any]:
    if isinstance(value, Mapping):
        config = dict(value)
    elif isinstance(value, str):
        try:
            config = json.loads(value)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Run {run_id!r} has invalid method_config JSON.") from exc
    else:
        raise ValueError(f"Run {run_id!r} has no resolved method_config.")
    if not isinstance(config, Mapping):
        raise ValueError(f"Run {run_id!r} method_config must resolve to a mapping.")
    return copy.deepcopy(dict(config))


def _config_without_use_shape(config: Mapping[str, Any]) -> tuple[dict[str, Any], Any]:
    result = copy.deepcopy(dict(config))
    params = result.get("params")
    if not isinstance(params, Mapping):
        raise ValueError("ShapeMix method config must contain a params mapping.")
    params = dict(params)
    if "use_shape" not in params:
        raise ValueError("ShapeMix method config must explicitly define params.use_shape.")
    use_shape = params.pop("use_shape")
    result["params"] = params
    return result, use_shape


def _signature_hash(run_dir: Path) -> Optional[str]:
    path = run_dir / "results" / "raw_method_output" / "signature_summary.yaml"
    if not path.exists():
        return None
    payload = _read_yaml(path)
    signature = payload.get("signature")
    if not isinstance(signature, Mapping):
        raise ValueError(f"Missing signature mapping in {path}.")
    value = signature.get("content_sha256")
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError(f"Missing or invalid signature.content_sha256 in {path}.")
    return value


def validate_paired_run_contract(records: pd.DataFrame) -> None:
    """Check that paired arms differ only by ``use_shape`` and share signatures."""

    if "method_config" not in records.columns:
        raise ValueError("runs.csv must record the fully resolved method_config.")
    for key, group in records.groupby(list(PAIR_KEY), sort=False):
        for field in ("dataset_id", "modality", "feature_set"):
            if field in group.columns and group[field].nunique(dropna=False) != 1:
                raise ValueError(
                    f"Pair {key!r} does not use one identical {field} across both arms."
                )
        by_arm = group.set_index("method_run_id", verify_integrity=True)
        if set(by_arm.index) != set(PRIMARY_ARMS):
            raise ValueError(f"Pair {key!r} does not contain exactly both ShapeMix arms.")
        normalized: dict[str, dict[str, Any]] = {}
        for arm in PRIMARY_ARMS:
            row = by_arm.loc[arm]
            config = _parse_method_config(row["method_config"], run_id=str(row["run_id"]))
            normalized[arm], use_shape = _config_without_use_shape(config)
            expected = arm == "shapemix_length"
            if use_shape is not expected:
                raise ValueError(
                    f"Pair {key!r} arm {arm!r} must set use_shape={expected}, got {use_shape!r}."
                )
        if normalized[PRIMARY_ARMS[0]] != normalized[PRIMARY_ARMS[1]]:
            raise ValueError(f"Pair {key!r} ShapeMix arms differ by more than params.use_shape.")
        expected_config = {
            "method": "shapemix",
            "params": copy.deepcopy(FROZEN_SHAPEMIX_PARAMS),
        }
        if normalized[PRIMARY_ARMS[0]] != expected_config:
            raise ValueError(
                f"Pair {key!r} does not use the exact frozen protocol-v1 ShapeMix config."
            )

        successful = group[group["status"] == "success"]
        if len(successful) == 2:
            hashes = {
                arm: _signature_hash(Path(str(row["run_dir"])))
                for arm, row in successful.set_index("method_run_id").iterrows()
            }
            if None in hashes.values():
                raise ValueError(f"Pair {key!r} is missing a required signature summary: {hashes!r}.")
            if len(set(hashes.values())) != 1:
                raise ValueError(f"Pair {key!r} uses different fixed signatures: {hashes!r}.")

    # A split has one training-reference pool. Its fixed signatures must be
    # identical across both inner mixtures and both simulation conditions, not
    # merely between the two arms of one dataset.
    for outer_seed, group in records.groupby("outer_split_seed", sort=False):
        hashes: dict[str, str] = {}
        for record in group[group["status"] == "success"].to_dict("records"):
            signature_hash = _signature_hash(Path(str(record["run_dir"])))
            if signature_hash is None:
                raise ValueError(
                    f"Successful run {record['run_id']!r} is missing its signature summary."
                )
            hashes[str(record["run_id"])] = signature_hash
        if hashes and len(set(hashes.values())) != 1:
            raise ValueError(
                f"Outer split {outer_seed} uses different signatures across inner seeds or "
                f"conditions: {hashes!r}."
            )


def _canonical_metric_id(row: Mapping[str, Any]) -> str:
    metric = _require_string(row.get("metric"), "comparison metric")
    expected_versions = {"rmse_v1": ("rmse", "v1"), "jsd_v2": ("jsd", "v2")}
    if metric not in expected_versions:
        raise ValueError(
            f"Noncanonical or unsupported primary metric {metric!r}; "
            f"required={list(PRIMARY_METRICS)!r}."
        )
    expected_name, expected_version = expected_versions[metric]
    if "metric_name" in row and _safe_error(row.get("metric_name")) is not None:
        if str(row["metric_name"]) != expected_name:
            raise ValueError(
                f"Metric {metric!r} must declare metric_name={expected_name!r}, "
                f"got {row['metric_name']!r}."
            )
    if "metric_version" in row and not pd.isna(row.get("metric_version")):
        raw_version = str(row["metric_version"])
        if raw_version != expected_version:
            raise ValueError(
                f"Metric {metric!r} must declare version {expected_version}, got {row['metric_version']!r}."
            )
    return metric


def validate_metric_versions(
    comparison: pd.DataFrame,
    *,
    expected_cell_types_by_run: Optional[Mapping[str, Sequence[str]]] = None,
) -> None:
    """Reject historical aliases and mismatched metric metadata."""

    required = {"run_id", "status"}
    missing = required.difference(comparison.columns)
    if missing:
        raise ValueError(f"Comparison is missing columns: {', '.join(sorted(missing))}.")
    success = comparison[comparison["status"] == "success"]
    if success.empty:
        return
    success_required = {"metric", "metric_name", "metric_version", "value"}
    missing = success_required.difference(comparison.columns)
    if missing:
        raise ValueError(
            "Successful comparison rows require columns: " + ", ".join(sorted(missing)) + "."
        )
    for run_id, group in success.groupby("run_id", sort=False):
        metric_ids = [_canonical_metric_id(row) for row in group.to_dict("records")]
        if len(metric_ids) != len(set(metric_ids)):
            raise ValueError(f"Run {run_id!r} has duplicate primary metric rows.")
        if set(metric_ids) != set(PRIMARY_METRICS):
            raise ValueError(
                f"Run {run_id!r} must report exactly {list(PRIMARY_METRICS)!r}, "
                f"found {sorted(metric_ids)!r}."
            )
        values = pd.to_numeric(group["value"], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Run {run_id!r} has a non-finite primary metric value.")
        if expected_cell_types_by_run is not None:
            if "cell_type_universe" not in group.columns:
                raise ValueError("Comparison must record cell_type_universe for each metric.")
            expected = tuple(expected_cell_types_by_run[str(run_id)])
            universes: set[tuple[str, ...]] = set()
            for raw_universe in group["cell_type_universe"]:
                try:
                    parsed = json.loads(str(raw_universe))
                except json.JSONDecodeError as exc:
                    raise ValueError(
                        f"Run {run_id!r} comparison has invalid cell_type_universe JSON."
                    ) from exc
                if not isinstance(parsed, list) or not all(
                    isinstance(value, str) for value in parsed
                ):
                    raise ValueError(
                        f"Run {run_id!r} comparison cell_type_universe must be a JSON string list."
                    )
                universes.add(tuple(parsed))
            if universes != {expected}:
                raise ValueError(
                    f"Run {run_id!r} comparison has a stale cell-type universe: {universes!r}."
                )


def _load_and_validate_comparisons(
    run_groups: Sequence[Path], records: pd.DataFrame
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for raw_group in run_groups:
        path = Path(raw_group).resolve() / "comparison.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing versioned comparison output: {path}")
        frame = pd.read_csv(path)
        frame["source_run_group"] = str(Path(raw_group).resolve())
        frames.append(frame)
    comparison = pd.concat(frames, ignore_index=True)
    expected = {
        str(row["run_id"]): tuple(row["_cell_types"])
        for row in records.to_dict("records")
    }
    validate_metric_versions(comparison, expected_cell_types_by_run=expected)
    known = set(records["run_id"])
    observed_ids = set(comparison["run_id"])
    if observed_ids != known:
        missing = sorted(known.difference(observed_ids))
        unknown = sorted(observed_ids.difference(known))
        raise ValueError(
            f"Comparison/run manifest run sets differ; missing={missing!r}, unknown={unknown!r}."
        )
    statuses = records.set_index("run_id")["status"].to_dict()
    for run_id, group in comparison.groupby("run_id", sort=False):
        observed_statuses = set(group["status"])
        if observed_statuses != {statuses[run_id]}:
            raise ValueError(
                f"Comparison status for run {run_id!r} disagrees with runs.csv: "
                f"{observed_statuses!r} != {statuses[run_id]!r}."
            )
    return comparison


def _run_metadata_matches(record: Mapping[str, Any]) -> None:
    run_dir = Path(str(record["run_dir"]))
    path = run_dir / "run.yaml"
    if not path.exists():
        raise FileNotFoundError(f"Successful run is missing {path}.")
    metadata = _read_yaml(path)
    for field in ("run_id", "dataset_id", "method", "method_run_id"):
        if metadata.get(field) != record[field]:
            raise ValueError(
                f"{path} field {field!r}={metadata.get(field)!r} disagrees with "
                f"runs.csv value {record[field]!r}."
            )


def _aligned_outputs(
    record: Mapping[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, tuple[str, ...]]:
    run_dir = Path(str(record["run_dir"]))
    truth_path = run_dir / "results" / "truth.csv"
    prediction_path = run_dir / "results" / "proportions.csv"
    if not truth_path.exists() or not prediction_path.exists():
        raise FileNotFoundError(
            f"Successful run {record['run_id']!r} lacks truth/proportion output files."
        )
    truth = pd.read_csv(truth_path, index_col=0)
    prediction = pd.read_csv(prediction_path, index_col=0)
    cell_types = tuple(record["_cell_types"])
    return (*align_proportions(truth, prediction, cell_types, row_sum_atol=ROW_SUM_ATOL), cell_types)


def _evaluation_fields(evaluation: Any) -> dict[str, Any]:
    if hasattr(evaluation, "to_record"):
        record = evaluation.to_record()
        if isinstance(record, Mapping):
            return dict(record)
    fields = {}
    for name in (
        "metric_id",
        "metric_name",
        "metric_version",
        "value",
        "n_spots",
        "n_cell_types",
        "cell_types",
    ):
        if hasattr(evaluation, name):
            fields[name] = getattr(evaluation, name)
    if "metric_id" not in fields and hasattr(evaluation, "metric"):
        fields["metric_id"] = getattr(evaluation, "metric")
    return fields


def compute_primary_metrics(records: pd.DataFrame) -> pd.DataFrame:
    """Recompute canonical primary metrics from standardized run outputs."""

    rows: list[dict[str, Any]] = []
    for record in records.to_dict("records"):
        if record["status"] != "success":
            continue
        _run_metadata_matches(record)
        run_dir = Path(str(record["run_dir"]))
        truth = pd.read_csv(run_dir / "results" / "truth.csv", index_col=0)
        prediction = pd.read_csv(run_dir / "results" / "proportions.csv", index_col=0)
        cell_types = tuple(record["_cell_types"])
        for metric_id in PRIMARY_METRICS:
            evaluation = evaluate_proportion_metric(
                metric_id,
                truth,
                prediction,
                cell_types,
                row_sum_atol=ROW_SUM_ATOL,
            )
            fields = _evaluation_fields(evaluation)
            resolved_id = fields.pop("metric_id", metric_id)
            if resolved_id != metric_id:
                raise ValueError(
                    f"Metric registry resolved {metric_id!r} to unexpected {resolved_id!r}."
                )
            value = float(fields.pop("value"))
            if not np.isfinite(value):
                raise ValueError(f"Run {record['run_id']!r} produced non-finite {metric_id}.")
            row = {field: record[field] for field in RUN_ID_COLUMNS}
            row.update(
                {
                    "metric": metric_id,
                    "value": value,
                    "metric_name": fields.pop(
                        "metric_name", "rmse" if metric_id == "rmse_v1" else "jsd"
                    ),
                    "metric_version": fields.pop(
                        "metric_version", "v1" if metric_id == "rmse_v1" else "v2"
                    ),
                    "n_spots": int(fields.pop("n_spots", len(truth))),
                    "n_cell_types": int(fields.pop("n_cell_types", len(cell_types))),
                    "cell_type_universe": _json_list(cell_types),
                }
            )
            rows.append(row)
    columns = list(RUN_ID_COLUMNS) + [
        "metric",
        "value",
        "metric_name",
        "metric_version",
        "n_spots",
        "n_cell_types",
        "cell_type_universe",
    ]
    return pd.DataFrame(rows, columns=columns)


def verify_comparison_values(
    comparison: pd.DataFrame,
    computed: pd.DataFrame,
    *,
    atol: float = 1e-12,
) -> None:
    observed = comparison[comparison["status"] == "success"].copy()
    if computed.empty and observed.empty:
        return
    required = {"run_id", "metric", "value"}
    missing_observed = required.difference(observed.columns)
    missing_computed = required.difference(computed.columns)
    if missing_observed or missing_computed:
        raise ValueError(
            "Computed/stored comparison schemas are incomplete; "
            f"stored_missing={sorted(missing_observed)!r}, "
            f"computed_missing={sorted(missing_computed)!r}."
        )
    observed = observed[["run_id", "metric", "value"]]
    observed["value"] = pd.to_numeric(observed["value"], errors="raise")
    merged = computed[["run_id", "metric", "value"]].merge(
        observed,
        on=["run_id", "metric"],
        how="outer",
        suffixes=("_computed", "_comparison"),
        indicator=True,
        validate="one_to_one",
    )
    if set(merged["_merge"]) != {"both"}:
        raise ValueError("Computed and comparison primary metric rows do not match exactly.")
    if not np.allclose(
        merged["value_computed"], merged["value_comparison"], rtol=0.0, atol=atol
    ):
        bad = merged.loc[
            ~np.isclose(
                merged["value_computed"], merged["value_comparison"], rtol=0.0, atol=atol
            )
        ]
        raise ValueError(f"Stored comparison values disagree with strict recomputation: {bad!r}.")


def build_paired_effects(records: pd.DataFrame, run_metrics: pd.DataFrame) -> pd.DataFrame:
    """Build one explicit result per nested pair and primary endpoint."""

    metric_lookup = {
        (str(row.run_id), str(row.metric)): float(row.value)
        for row in run_metrics.itertuples()
    }
    rows: list[dict[str, Any]] = []
    for key, group in records.groupby(list(PAIR_KEY), sort=True):
        by_arm = group.set_index("method_run_id", verify_integrity=True)
        for metric in PRIMARY_METRICS:
            row: dict[str, Any] = dict(zip(PAIR_KEY, key))
            row["metric"] = metric
            arm_values: dict[str, float] = {}
            reasons: list[str] = []
            for arm in PRIMARY_ARMS:
                if arm not in by_arm.index:
                    row[f"{arm}_status"] = "missing"
                    row[f"{arm}_run_id"] = None
                    row[f"{arm}_value"] = np.nan
                    row[f"{arm}_error"] = "arm absent from run manifests"
                    reasons.append(f"{arm}: missing")
                    continue
                record = by_arm.loc[arm]
                status = str(record["status"])
                run_id = str(record["run_id"])
                row[f"{arm}_status"] = status
                row[f"{arm}_run_id"] = run_id
                row[f"{arm}_error"] = _safe_error(record.get("error"))
                value = metric_lookup.get((run_id, metric))
                if status != "success" or value is None:
                    row[f"{arm}_value"] = np.nan
                    reasons.append(
                        f"{arm}: {_safe_error(record.get('error')) or 'metric unavailable'}"
                    )
                else:
                    row[f"{arm}_value"] = value
                    arm_values[arm] = value
            if len(arm_values) == 2:
                row["pair_status"] = "complete"
                row["unavailable_reason"] = None
                row["delta_length_minus_count_only"] = (
                    arm_values["shapemix_length"] - arm_values["shapemix_count_only"]
                )
            else:
                row["pair_status"] = "unavailable"
                row["unavailable_reason"] = "; ".join(reasons)
                row["delta_length_minus_count_only"] = np.nan
            rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["condition", "metric", "outer_split_seed", "inner_mixture_seed"]
    ).reset_index(drop=True)


def average_inner_effects(paired_effects: pd.DataFrame) -> pd.DataFrame:
    """Average exactly two inner-mixture effects inside each outer split."""

    rows: list[dict[str, Any]] = []
    groups = paired_effects.groupby(["condition", "metric", "outer_split_seed"], sort=True)
    for (condition, metric, outer), group in groups:
        actual_inner = set(group["inner_mixture_seed"])
        complete = group[group["pair_status"] == "complete"]
        valid = (
            actual_inner == set(PRIMARY_INNER_MIXTURE_SEEDS)
            and len(group) == len(PRIMARY_INNER_MIXTURE_SEEDS)
            and len(complete) == len(PRIMARY_INNER_MIXTURE_SEEDS)
        )
        unavailable = group[group["pair_status"] != "complete"]
        rows.append(
            {
                "condition": condition,
                "metric": metric,
                "outer_split_seed": int(outer),
                "outer_status": "complete" if valid else "unavailable",
                "n_inner_expected": len(PRIMARY_INNER_MIXTURE_SEEDS),
                "n_inner_available": len(complete),
                "inner_mixture_seeds": json.dumps(sorted(actual_inner)),
                "delta_outer": (
                    float(complete["delta_length_minus_count_only"].mean())
                    if valid
                    else np.nan
                ),
                "unavailable_reason": (
                    None
                    if valid
                    else "; ".join(
                        unavailable["unavailable_reason"].dropna().astype(str).tolist()
                    )
                    or "not all frozen inner-mixture pairs are complete"
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["condition", "metric", "outer_split_seed"]
    ).reset_index(drop=True)


def bootstrap_mean_interval(
    values: Sequence[float],
    *,
    replicates: int = BOOTSTRAP_REPLICATES,
    seed_tuple: Sequence[int] = BOOTSTRAP_SEED_TUPLE,
) -> tuple[float, float]:
    """Percentile interval from deterministic outer-split resampling."""

    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or len(array) == 0 or not np.isfinite(array).all():
        raise ValueError("Bootstrap input must be a nonempty finite one-dimensional sequence.")
    if isinstance(replicates, bool) or not isinstance(replicates, int) or replicates <= 0:
        raise ValueError("replicates must be a positive integer.")
    seed = np.random.SeedSequence([int(value) for value in seed_tuple])
    rng = np.random.Generator(np.random.PCG64(seed))
    indices = rng.integers(0, len(array), size=(replicates, len(array)))
    bootstrap_means = array[indices].mean(axis=1)
    lower, upper = np.percentile(bootstrap_means, [2.5, 97.5])
    return float(lower), float(upper)


def exact_two_sided_sign_flip_pvalue(values: Sequence[float]) -> float:
    """Exact two-sided paired sign-flip p-value for the mean effect."""

    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or len(array) == 0 or not np.isfinite(array).all():
        raise ValueError("Sign-flip input must be a nonempty finite one-dimensional sequence.")
    observed = abs(float(array.mean()))
    absolute_null_means = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(array)):
        absolute_null_means.append(abs(float(np.mean(array * np.asarray(signs)))))
    tolerance = np.finfo(float).eps * max(1.0, observed) * 8.0
    extreme = np.count_nonzero(np.asarray(absolute_null_means) >= observed - tolerance)
    return float(extreme / len(absolute_null_means))


def summarize_outer_effects(
    outer_effects: pd.DataFrame,
    *,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
) -> pd.DataFrame:
    """Summarize complete five-outer-split effects without partial-case analysis."""

    rows: list[dict[str, Any]] = []
    for condition in PRIMARY_CONDITIONS:
        for metric in PRIMARY_METRICS:
            group = outer_effects[
                (outer_effects["condition"] == condition) & (outer_effects["metric"] == metric)
            ]
            complete = group[group["outer_status"] == "complete"]
            seeds = set(group["outer_split_seed"])
            valid = (
                seeds == set(PRIMARY_OUTER_SPLIT_SEEDS)
                and len(group) == len(PRIMARY_OUTER_SPLIT_SEEDS)
                and len(complete) == len(PRIMARY_OUTER_SPLIT_SEEDS)
            )
            row: dict[str, Any] = {
                "condition": condition,
                "metric": metric,
                "analysis_status": "complete" if valid else "incomplete",
                "reporting_scope": REPORTING_SCOPE,
                "n_outer_expected": len(PRIMARY_OUTER_SPLIT_SEEDS),
                "n_outer_available": len(complete),
                "bootstrap_seed_tuple": json.dumps(list(BOOTSTRAP_SEED_TUPLE)),
                "bootstrap_replicates": bootstrap_replicates,
                "sign_flip_label": "exploratory_exact_two_sided_paired_sign_flip",
                "unavailable_reason": (
                    None
                    if valid
                    else "formal summary requires all five frozen outer-split effects"
                ),
            }
            if valid:
                values = complete.sort_values("outer_split_seed")["delta_outer"].to_numpy(float)
                lower, upper = bootstrap_mean_interval(
                    values,
                    replicates=bootstrap_replicates,
                )
                row.update(
                    {
                        "mean_outer_effect": float(np.mean(values)),
                        "median_outer_effect": float(np.median(values)),
                        "sd_outer_effect": float(np.std(values, ddof=1)),
                        "min_outer_effect": float(np.min(values)),
                        "max_outer_effect": float(np.max(values)),
                        "bootstrap_percentile_95_lower": lower,
                        "bootstrap_percentile_95_upper": upper,
                        "n_outer_improved": int(np.count_nonzero(values < 0.0)),
                        "exploratory_sign_flip_pvalue": exact_two_sided_sign_flip_pvalue(values),
                    }
                )
            else:
                row.update(
                    {
                        "mean_outer_effect": np.nan,
                        "median_outer_effect": np.nan,
                        "sd_outer_effect": np.nan,
                        "min_outer_effect": np.nan,
                        "max_outer_effect": np.nan,
                        "bootstrap_percentile_95_lower": np.nan,
                        "bootstrap_percentile_95_upper": np.nan,
                        "n_outer_improved": np.nan,
                        "exploratory_sign_flip_pvalue": np.nan,
                    }
                )
            rows.append(row)
    summary = pd.DataFrame(rows)
    for condition in PRIMARY_CONDITIONS:
        selected = summary[summary["condition"] == condition]
        support = bool(
            len(selected) == len(PRIMARY_METRICS)
            and (selected["analysis_status"] == "complete").all()
            and (selected["mean_outer_effect"] < 0.0).all()
            and (selected["n_outer_improved"] >= 4).all()
        )
        summary.loc[summary["condition"] == condition, "directional_support_rule_met"] = support
    return summary


def _correlation(
    truth: np.ndarray, prediction: np.ndarray, *, method: str
) -> tuple[float, bool, Optional[str]]:
    if len(truth) < 2:
        return np.nan, False, "fewer than two spots"
    if np.all(truth == truth[0]):
        return np.nan, False, "truth is constant across spots"
    if np.all(prediction == prediction[0]):
        return np.nan, False, "prediction is constant across spots"
    if method == "spearman":
        truth = rankdata(truth, method="average")
        prediction = rankdata(prediction, method="average")
    value = float(np.corrcoef(truth, prediction)[0, 1])
    if not np.isfinite(value):
        return np.nan, False, "correlation is non-finite"
    return value, True, None


def compute_cell_type_metrics(records: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for record in records.to_dict("records"):
        if record["status"] != "success":
            continue
        truth, prediction, cell_types = _aligned_outputs(record)
        for cell_type in cell_types:
            truth_values = truth[cell_type].to_numpy(float)
            prediction_values = prediction[cell_type].to_numpy(float)
            correlation_values = {
                "pearson_r_v1": _correlation(truth_values, prediction_values, method="pearson"),
                "spearman_rho_v1": _correlation(
                    truth_values, prediction_values, method="spearman"
                ),
            }
            base = {field: record[field] for field in RUN_ID_COLUMNS}
            for metric, value in (
                ("mae_v1", float(np.mean(np.abs(truth_values - prediction_values)))),
                ("rmse_v1", float(np.sqrt(np.mean((truth_values - prediction_values) ** 2)))),
            ):
                rows.append(
                    {
                        **base,
                        "cell_type": cell_type,
                        "metric": metric,
                        "value": value,
                        "defined": True,
                        "undefined_reason": None,
                        "n_spots": len(truth_values),
                    }
                )
            for metric, (value, defined, reason) in correlation_values.items():
                rows.append(
                    {
                        **base,
                        "cell_type": cell_type,
                        "metric": metric,
                        "value": value,
                        "defined": defined,
                        "undefined_reason": reason,
                        "n_spots": len(truth_values),
                    }
                )
    columns = list(RUN_ID_COLUMNS) + [
        "cell_type",
        "metric",
        "value",
        "defined",
        "undefined_reason",
        "n_spots",
    ]
    return pd.DataFrame(rows, columns=columns)


def _binary_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, Any]:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    if labels.shape != scores.shape or labels.ndim != 1:
        raise ValueError("Rare-cell labels and scores must be one-dimensional with the same shape.")
    predicted = scores >= PREDICTED_PRESENCE_THRESHOLD
    tp = int(np.count_nonzero(labels & predicted))
    fp = int(np.count_nonzero(~labels & predicted))
    fn = int(np.count_nonzero(labels & ~predicted))
    tn = int(np.count_nonzero(~labels & ~predicted))
    positives = tp + fn
    negatives = fp + tn
    precision = float(tp / (tp + fp)) if tp + fp else np.nan
    recall = float(tp / positives) if positives else np.nan
    if np.isfinite(precision) and np.isfinite(recall):
        f1 = float(2.0 * precision * recall / (precision + recall)) if precision + recall else 0.0
    else:
        f1 = np.nan

    auprc = np.nan
    auprc_reason: Optional[str] = None
    if positives == 0:
        auprc_reason = "truth contains no positive class"
    elif negatives == 0:
        auprc_reason = "truth contains no negative class"
    else:
        order = np.argsort(-scores, kind="mergesort")
        ordered_scores = scores[order]
        ordered_labels = labels[order]
        tp_cumulative = np.cumsum(ordered_labels)
        fp_cumulative = np.cumsum(~ordered_labels)
        group_ends = np.flatnonzero(
            np.r_[ordered_scores[1:] != ordered_scores[:-1], True]
        )
        previous_tp = 0
        area = 0.0
        for index in group_ends:
            current_tp = int(tp_cumulative[index])
            current_fp = int(fp_cumulative[index])
            recall_increment = (current_tp - previous_tp) / positives
            area += recall_increment * current_tp / (current_tp + current_fp)
            previous_tp = current_tp
        auprc = float(area)

    return {
        "n_pairs": len(labels),
        "truth_positives": positives,
        "truth_negatives": negatives,
        "predicted_positives": tp + fp,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "precision": precision,
        "precision_defined": bool(tp + fp),
        "precision_undefined_reason": None if tp + fp else "no predicted positives",
        "recall": recall,
        "recall_defined": bool(positives),
        "recall_undefined_reason": None if positives else "truth contains no positive class",
        "f1": f1,
        "f1_defined": bool(np.isfinite(f1)),
        "f1_undefined_reason": (
            None if np.isfinite(f1) else "precision or recall is undefined"
        ),
        "auprc": auprc,
        "auprc_defined": bool(np.isfinite(auprc)),
        "auprc_undefined_reason": auprc_reason,
        "has_both_truth_classes": bool(positives and negatives),
    }


def compute_rare_cell_metrics(records: pd.DataFrame) -> pd.DataFrame:
    """Compute frozen-threshold per-type, pooled-micro, and macro metrics."""

    rows: list[dict[str, Any]] = []
    for record in records.to_dict("records"):
        if record["status"] != "success":
            continue
        truth, prediction, cell_types = _aligned_outputs(record)
        missing_rare = set(RARE_REFERENCE_TYPES).difference(cell_types)
        if missing_rare:
            raise ValueError(
                f"Run {record['run_id']!r} lacks frozen rare reference types {sorted(missing_rare)!r}."
            )
        base = {field: record[field] for field in RUN_ID_COLUMNS}
        per_type: list[dict[str, Any]] = []
        pooled_labels: list[np.ndarray] = []
        pooled_scores: list[np.ndarray] = []
        for cell_type in RARE_REFERENCE_TYPES:
            labels = truth[cell_type].to_numpy(float) > TRUE_PRESENCE_THRESHOLD
            scores = prediction[cell_type].to_numpy(float)
            metrics = _binary_metrics(labels, scores)
            row = {
                **base,
                "aggregation": "per_type",
                "cell_type": cell_type,
                "n_types_in_aggregate": 1,
                "true_presence_rule": "truth_proportion > 0",
                "predicted_presence_threshold": PREDICTED_PRESENCE_THRESHOLD,
                **metrics,
            }
            rows.append(row)
            per_type.append(row)
            pooled_labels.append(labels)
            pooled_scores.append(scores)

        micro = _binary_metrics(np.concatenate(pooled_labels), np.concatenate(pooled_scores))
        rows.append(
            {
                **base,
                "aggregation": "pooled_micro",
                "cell_type": None,
                "n_types_in_aggregate": len(RARE_REFERENCE_TYPES),
                "true_presence_rule": "truth_proportion > 0",
                "predicted_presence_threshold": PREDICTED_PRESENCE_THRESHOLD,
                **micro,
            }
        )

        eligible = [row for row in per_type if row["has_both_truth_classes"]]
        macro: dict[str, Any] = {
            **base,
            "aggregation": "macro_evaluable_types",
            "cell_type": None,
            "n_types_in_aggregate": len(eligible),
            "n_pairs": np.nan,
            "truth_positives": np.nan,
            "truth_negatives": np.nan,
            "predicted_positives": np.nan,
            "tp": np.nan,
            "fp": np.nan,
            "fn": np.nan,
            "tn": np.nan,
            "true_presence_rule": "truth_proportion > 0",
            "predicted_presence_threshold": PREDICTED_PRESENCE_THRESHOLD,
            "has_both_truth_classes": bool(eligible),
        }
        for metric in ("precision", "recall", "f1", "auprc"):
            defined_rows = [row for row in eligible if row[f"{metric}_defined"]]
            undefined_rows = [row for row in eligible if not row[f"{metric}_defined"]]
            macro[f"{metric}_n_types_expected"] = len(eligible)
            macro[f"{metric}_n_types_defined"] = len(defined_rows)
            macro_is_defined = bool(eligible) and not undefined_rows
            macro[metric] = (
                float(np.mean([row[metric] for row in defined_rows]))
                if macro_is_defined
                else np.nan
            )
            macro[f"{metric}_defined"] = macro_is_defined
            if not eligible:
                reason = "no rare type has both positive and negative truth classes"
            elif undefined_rows:
                details = "; ".join(
                    f"{row['cell_type']}: {row[f'{metric}_undefined_reason']}"
                    for row in undefined_rows
                )
                reason = (
                    "macro requires every truth-evaluable rare type; undefined constituents: "
                    + details
                )
            else:
                reason = None
            macro[f"{metric}_undefined_reason"] = reason
        rows.append(macro)
    if rows:
        return pd.DataFrame(rows)
    columns = list(RUN_ID_COLUMNS) + [
        "aggregation",
        "cell_type",
        "n_types_in_aggregate",
        "true_presence_rule",
        "predicted_presence_threshold",
        "n_pairs",
        "truth_positives",
        "truth_negatives",
        "predicted_positives",
        "tp",
        "fp",
        "fn",
        "tn",
        "precision",
        "precision_defined",
        "precision_undefined_reason",
        "recall",
        "recall_defined",
        "recall_undefined_reason",
        "f1",
        "f1_defined",
        "f1_undefined_reason",
        "auprc",
        "auprc_defined",
        "auprc_undefined_reason",
        "has_both_truth_classes",
    ]
    return pd.DataFrame(columns=columns)


def _dimension_key(values: Sequence[Any]) -> tuple[Any, ...]:
    return tuple(None if pd.isna(value) else value for value in values)


def build_secondary_nested_effects(
    records: pd.DataFrame,
    measurements: pd.DataFrame,
    *,
    dimension_columns: Sequence[str],
    dimension_values: Sequence[Sequence[Any]],
    metric_ids: Sequence[str],
) -> pd.DataFrame:
    """Pair secondary metrics and average them over the two inner mixtures."""

    lookup: dict[tuple[Any, ...], dict[str, Any]] = {}
    for row in measurements.to_dict("records"):
        key = (
            str(row["run_id"]),
            *_dimension_key([row[column] for column in dimension_columns]),
            str(row["metric"]),
        )
        if key in lookup:
            raise ValueError(f"Duplicate secondary metric row for key {key!r}.")
        lookup[key] = row

    pair_rows: list[dict[str, Any]] = []
    for pair_key, group in records.groupby(list(PAIR_KEY), sort=True):
        by_arm = group.set_index("method_run_id", verify_integrity=True)
        for raw_dimensions in dimension_values:
            dimensions = _dimension_key(raw_dimensions)
            dimension_fields = dict(zip(dimension_columns, dimensions))
            for metric in metric_ids:
                row: dict[str, Any] = {
                    "effect_level": "inner_pair",
                    **dict(zip(PAIR_KEY, pair_key)),
                    **dimension_fields,
                    "metric": metric,
                    "n_inner_expected": 1,
                    "n_inner_available": np.nan,
                }
                values: dict[str, float] = {}
                reasons: list[str] = []
                for arm in PRIMARY_ARMS:
                    if arm not in by_arm.index:
                        row[f"{arm}_status"] = "missing"
                        row[f"{arm}_run_id"] = None
                        row[f"{arm}_value"] = np.nan
                        reasons.append(f"{arm}: arm absent from run manifests")
                        continue
                    run = by_arm.loc[arm]
                    run_id = str(run["run_id"])
                    row[f"{arm}_run_id"] = run_id
                    if run["status"] != "success":
                        row[f"{arm}_status"] = "failed"
                        row[f"{arm}_value"] = np.nan
                        reasons.append(
                            f"{arm}: {_safe_error(run.get('error')) or 'run failed'}"
                        )
                        continue
                    measured = lookup.get((run_id, *dimensions, metric))
                    if measured is None:
                        row[f"{arm}_status"] = "metric_missing"
                        row[f"{arm}_value"] = np.nan
                        reasons.append(f"{arm}: secondary metric row is missing")
                        continue
                    if not bool(measured["defined"]):
                        row[f"{arm}_status"] = "metric_undefined"
                        row[f"{arm}_value"] = np.nan
                        reasons.append(
                            f"{arm}: {measured.get('undefined_reason') or 'metric undefined'}"
                        )
                        continue
                    value = float(measured["value"])
                    if not np.isfinite(value):
                        raise ValueError(f"Defined secondary metric {metric!r} is non-finite.")
                    row[f"{arm}_status"] = "success"
                    row[f"{arm}_value"] = value
                    values[arm] = value
                if len(values) == 2:
                    row["effect_status"] = "complete"
                    row["unavailable_reason"] = None
                    row["delta_length_minus_count_only"] = (
                        values["shapemix_length"] - values["shapemix_count_only"]
                    )
                    row["n_inner_available"] = 1
                else:
                    row["effect_status"] = "unavailable"
                    row["unavailable_reason"] = "; ".join(reasons)
                    row["delta_length_minus_count_only"] = np.nan
                    row["n_inner_available"] = 0
                pair_rows.append(row)

    pair_frame = pd.DataFrame(pair_rows)
    outer_rows: list[dict[str, Any]] = []
    group_columns = ["condition", "outer_split_seed", *dimension_columns, "metric"]
    for group_key, group in pair_frame.groupby(group_columns, sort=True, dropna=False):
        key_values = group_key if isinstance(group_key, tuple) else (group_key,)
        fields = dict(zip(group_columns, key_values))
        complete = group[group["effect_status"] == "complete"]
        valid = (
            set(group["inner_mixture_seed"]) == set(PRIMARY_INNER_MIXTURE_SEEDS)
            and len(group) == len(PRIMARY_INNER_MIXTURE_SEEDS)
            and len(complete) == len(PRIMARY_INNER_MIXTURE_SEEDS)
        )
        unavailable = group[group["effect_status"] != "complete"]
        outer_rows.append(
            {
                "effect_level": "outer_mean",
                **fields,
                "inner_mixture_seed": np.nan,
                "shapemix_length_status": "aggregated" if valid else "unavailable",
                "shapemix_length_run_id": None,
                "shapemix_length_value": (
                    float(complete["shapemix_length_value"].mean()) if valid else np.nan
                ),
                "shapemix_count_only_status": "aggregated" if valid else "unavailable",
                "shapemix_count_only_run_id": None,
                "shapemix_count_only_value": (
                    float(complete["shapemix_count_only_value"].mean()) if valid else np.nan
                ),
                "effect_status": "complete" if valid else "unavailable",
                "delta_length_minus_count_only": (
                    float(complete["delta_length_minus_count_only"].mean())
                    if valid
                    else np.nan
                ),
                "n_inner_expected": len(PRIMARY_INNER_MIXTURE_SEEDS),
                "n_inner_available": len(complete),
                "unavailable_reason": (
                    None
                    if valid
                    else "; ".join(
                        unavailable["unavailable_reason"].dropna().astype(str).tolist()
                    )
                    or "not all frozen inner-mixture effects are complete"
                ),
            }
        )
    return pd.concat([pair_frame, pd.DataFrame(outer_rows)], ignore_index=True)


def compute_cell_type_paired_effects(
    records: pd.DataFrame,
    cell_type_metrics: pd.DataFrame,
) -> pd.DataFrame:
    dimensions = [(cell_type,) for cell_type in FROZEN_CELL_TYPES]
    return build_secondary_nested_effects(
        records,
        cell_type_metrics,
        dimension_columns=("cell_type",),
        dimension_values=dimensions,
        metric_ids=("mae_v1", "rmse_v1", "pearson_r_v1", "spearman_rho_v1"),
    )


def compute_rare_cell_paired_effects(
    records: pd.DataFrame,
    rare_cell_metrics: pd.DataFrame,
) -> pd.DataFrame:
    tidy_rows: list[dict[str, Any]] = []
    for row in rare_cell_metrics.to_dict("records"):
        for metric in ("precision", "recall", "f1", "auprc"):
            tidy_rows.append(
                {
                    "run_id": row["run_id"],
                    "aggregation": row["aggregation"],
                    "cell_type": row["cell_type"],
                    "metric": metric,
                    "value": row[metric],
                    "defined": row[f"{metric}_defined"],
                    "undefined_reason": row[f"{metric}_undefined_reason"],
                }
            )
    tidy = pd.DataFrame(
        tidy_rows,
        columns=[
            "run_id",
            "aggregation",
            "cell_type",
            "metric",
            "value",
            "defined",
            "undefined_reason",
        ],
    )
    dimensions = [
        *(('per_type', cell_type) for cell_type in RARE_REFERENCE_TYPES),
        ("pooled_micro", None),
        ("macro_evaluable_types", None),
    ]
    return build_secondary_nested_effects(
        records,
        tidy,
        dimension_columns=("aggregation", "cell_type"),
        dimension_values=dimensions,
        metric_ids=("precision", "recall", "f1", "auprc"),
    )


def compute_performance(records: pd.DataFrame) -> pd.DataFrame:
    """Report convergence, runtime, and peak memory when it was recorded."""

    rows: list[dict[str, Any]] = []
    for record in records.to_dict("records"):
        base = {field: record[field] for field in RUN_ID_COLUMNS}
        row: dict[str, Any] = {
            **base,
            "status": record["status"],
            "error": _safe_error(record.get("error")),
            "runtime_seconds": float(record["wall_runtime_seconds"]),
            "runtime_source": "run.yaml:performance.wall_runtime_seconds",
            "peak_memory_bytes": int(record["peak_rss_bytes"]),
            "peak_memory_mb": float(record["peak_rss_mb"]),
            "peak_memory_source": "run.yaml:performance.peak_rss_bytes",
            "peak_memory_available": True,
            "rss_measurement": record["rss_measurement"],
            "execution_action": record["execution_action"],
            "converged": np.nan,
            "selected_restart": np.nan,
            "count_log_likelihood": np.nan,
            "shape_log_likelihood": np.nan,
            "abundance_log_prior": np.nan,
            "total_log_objective": np.nan,
        }
        if record["status"] == "success":
            diagnostics_path = Path(str(record["run_dir"])) / "results" / "diagnostics.json"
            if not diagnostics_path.exists():
                raise FileNotFoundError(f"Successful run is missing {diagnostics_path}.")
            with diagnostics_path.open() as handle:
                diagnostics = json.load(handle)
            fit = diagnostics.get("fit")
            if not isinstance(fit, Mapping):
                raise ValueError(f"{diagnostics_path} lacks fit diagnostics.")
            row["converged"] = bool(fit.get("success", False))
            row["selected_restart"] = fit.get("selected_restart")
            objective = fit.get("objective")
            if isinstance(objective, Mapping):
                for field in (
                    "count_log_likelihood",
                    "shape_log_likelihood",
                    "abundance_log_prior",
                    "total_log_objective",
                ):
                    row[field] = objective.get(field, np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def collect_reconstruction(records: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for record in records.to_dict("records"):
        if record["status"] != "success":
            continue
        path = (
            Path(str(record["run_dir"]))
            / "results"
            / "raw_method_output"
            / "reconstruction_summary.csv"
        )
        if not path.exists():
            raise FileNotFoundError(f"Successful ShapeMix run is missing {path}.")
        frame = pd.read_csv(path)
        for field in reversed(RUN_ID_COLUMNS):
            frame.insert(0, field, record[field])
        rows.append(frame)
    if rows:
        return pd.concat(rows, ignore_index=True)
    return pd.DataFrame(
        columns=list(RUN_ID_COLUMNS)
        + [
            "component",
            "bin_name",
            "layer_name",
            "entries",
            "observed_nonzero_entries",
            "expected_positive_entries",
            "observed_total",
            "expected_total",
            "residual_total",
            "absolute_error_sum",
            "squared_error_sum",
            "mean_absolute_error",
            "root_mean_squared_error",
        ]
    )


def collect_failures(records: pd.DataFrame) -> pd.DataFrame:
    columns = list(RUN_ID_COLUMNS) + ["status", "error"]
    failures = records[records["status"] != "success"].copy()
    return failures.loc[:, columns].reset_index(drop=True)


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False)


def _summary_manifest(
    run_groups: Sequence[Path],
    records: pd.DataFrame,
    provenance: pd.DataFrame,
    paired: pd.DataFrame,
    outer: pd.DataFrame,
    primary_summary: pd.DataFrame,
) -> dict[str, Any]:
    first_provenance = provenance.iloc[0]
    support = {
        condition: bool(
            primary_summary.loc[
                primary_summary["condition"] == condition, "directional_support_rule_met"
            ].all()
        )
        for condition in PRIMARY_CONDITIONS
    }
    return {
        "schema_version": 1,
        "benchmark_protocol_version": PROTOCOL_VERSION,
        "analysis": "paired ShapeMix length-minus-count-only ablation",
        "reporting_scope": REPORTING_SCOPE,
        "biological_replication": False,
        "donors": 1,
        "primary_metrics": list(PRIMARY_METRICS),
        "conditions": list(PRIMARY_CONDITIONS),
        "outer_split_seeds": list(PRIMARY_OUTER_SPLIT_SEEDS),
        "inner_mixture_seeds": list(PRIMARY_INNER_MIXTURE_SEEDS),
        "pairing_key": list(PAIR_KEY),
        "effect_definition": "shapemix_length minus shapemix_count_only; negative favors ShapeMix",
        "inner_aggregation": "mean of exactly two paired inner-mixture effects within outer split",
        "resampling_unit": "outer_split_seed",
        "bootstrap": {
            "method": "percentile bootstrap of mean outer effect",
            "seed_tuple": list(BOOTSTRAP_SEED_TUPLE),
            "replicates": BOOTSTRAP_REPLICATES,
            "interval": 0.95,
        },
        "sign_flip_test": "exact two-sided paired sign-flip across five outer effects; exploratory",
        "rare_reference_types": list(RARE_REFERENCE_TYPES),
        "rare_type_definition": "pre-split prevalence below 2% of the 9500-cell universe",
        "true_presence_rule": "truth proportion > 0",
        "predicted_presence_threshold": PREDICTED_PRESENCE_THRESHOLD,
        "run_groups": [str(Path(value).resolve()) for value in run_groups],
        "execution_provenance": {
            "experiment_config": {
                "path": first_provenance["experiment_config_path"],
                "source_sha256": first_provenance["experiment_config_source_sha256"],
                "resolved_sha256": first_provenance["experiment_config_resolved_sha256"],
            },
            "benchmark_protocol": {
                "path": first_provenance["benchmark_protocol_path"],
                "source_sha256": first_provenance["benchmark_protocol_sha256"],
            },
            "registry": {
                "path": first_provenance["registry_path"],
                "source_sha256": first_provenance["registry_source_sha256"],
            },
            "dataset_config_sha256_by_dataset": {
                str(dataset_id): str(group["dataset_config_sha256"].iloc[0])
                for dataset_id, group in provenance.groupby("dataset_id", sort=True)
            },
            "code": {
                "git_commit": first_provenance["git_commit"],
                "worktree_dirty": bool(first_provenance["git_worktree_dirty"]),
                "manifest_sha256": first_provenance["code_manifest_sha256"],
            },
            "shards": {
                "count": int(first_provenance["shard_count"]),
                "indices": sorted(provenance["shard_index"].astype(int).unique().tolist()),
            },
            "all_run_output_manifests_rehashed": True,
        },
        "counts": {
            "runs": len(records),
            "failed_runs": int((records["status"] != "success").sum()),
            "nested_pairs": int(paired[list(PAIR_KEY)].drop_duplicates().shape[0]),
            "unavailable_pair_metrics": int((paired["pair_status"] != "complete").sum()),
            "unavailable_outer_metrics": int((outer["outer_status"] != "complete").sum()),
        },
        "directional_support_rule": (
            "for each condition, both co-primary mean effects are below zero and at least "
            "four of five outer effects favor ShapeMix for each endpoint"
        ),
        "directional_support_rule_met": support,
        "scientific_interpretation_limit": (
            "Intervals and exploratory tests quantify conditional resampling variability "
            "within one donor and cannot establish donor-level or population-level generalization."
        ),
        "outputs": {
            "run_metrics": "run_metrics.csv",
            "paired_effects": "paired_effects.csv",
            "outer_effects": "outer_effects.csv",
            "primary_summary": "primary_summary.csv",
            "cell_type_metrics": "cell_type_metrics.csv",
            "cell_type_paired_effects": "cell_type_paired_effects.csv",
            "rare_cell_metrics": "rare_cell_metrics.csv",
            "rare_cell_paired_effects": "rare_cell_paired_effects.csv",
            "performance": "performance.csv",
            "reconstruction": "reconstruction.csv",
            "failures": "failures.csv",
            "provenance": "provenance.csv",
        },
    }


def summarize_benchmark(
    run_groups: Sequence[Path],
    output_dir: Path,
    *,
    registry: Path = ROOT / "data" / "registry" / "datasets.yaml",
    project_root: Path = ROOT,
) -> dict[str, Any]:
    """Validate, combine, summarize, and write the frozen primary benchmark."""

    groups = [Path(value).resolve() for value in run_groups]
    records = load_run_groups(
        groups,
        registry=Path(registry).resolve(),
        project_root=Path(project_root).resolve(),
    )
    validate_primary_design(records)
    validate_paired_run_contract(records)
    provenance = collect_and_validate_provenance(
        records,
        project_root=Path(project_root).resolve(),
    )
    comparison = _load_and_validate_comparisons(groups, records)
    run_metrics = compute_primary_metrics(records)
    verify_comparison_values(comparison, run_metrics)
    paired = build_paired_effects(records, run_metrics)
    outer = average_inner_effects(paired)
    primary_summary = summarize_outer_effects(outer)
    cell_type_metrics = compute_cell_type_metrics(records)
    cell_type_paired_effects = compute_cell_type_paired_effects(
        records, cell_type_metrics
    )
    rare_cell_metrics = compute_rare_cell_metrics(records)
    rare_cell_paired_effects = compute_rare_cell_paired_effects(
        records, rare_cell_metrics
    )
    performance = compute_performance(records)
    reconstruction = collect_reconstruction(records)
    failures = collect_failures(records)

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(run_metrics, output / "run_metrics.csv")
    _write_csv(paired, output / "paired_effects.csv")
    _write_csv(outer, output / "outer_effects.csv")
    _write_csv(primary_summary, output / "primary_summary.csv")
    _write_csv(cell_type_metrics, output / "cell_type_metrics.csv")
    _write_csv(cell_type_paired_effects, output / "cell_type_paired_effects.csv")
    _write_csv(rare_cell_metrics, output / "rare_cell_metrics.csv")
    _write_csv(rare_cell_paired_effects, output / "rare_cell_paired_effects.csv")
    _write_csv(performance, output / "performance.csv")
    _write_csv(reconstruction, output / "reconstruction.csv")
    _write_csv(failures, output / "failures.csv")
    public_provenance = provenance.drop(columns=["_code_files"])
    _write_csv(public_provenance, output / "provenance.csv")

    manifest = _summary_manifest(
        groups,
        records,
        provenance,
        paired,
        outer,
        primary_summary,
    )
    with (output / "summary.yaml").open("w") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize the frozen paired ShapeMix primary benchmark across one or more "
            "run-group directories."
        )
    )
    parser.add_argument("--run-groups", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--registry",
        type=Path,
        default=ROOT / "data" / "registry" / "datasets.yaml",
    )
    args = parser.parse_args()
    manifest = summarize_benchmark(
        args.run_groups,
        args.output_dir,
        registry=args.registry,
    )
    print(args.output_dir / "summary.yaml")
    if manifest["counts"]["failed_runs"]:
        print(
            "WARNING: failures are retained; formal summaries are unavailable where pairing is incomplete.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
