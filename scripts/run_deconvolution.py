#!/usr/bin/env python
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import platform
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Union

import pandas as pd
import psutil
import yaml

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data import load_deconvolution_input
from deconvatac.data.registry import get_dataset_config
from deconvatac.metrics import (
    PROPORTION_CONTRACT_VERSION,
    PROPORTION_ROW_SUM_ATOL,
    available_proportion_metrics,
    declared_cell_types_from_metadata,
    evaluate_proportion_metric,
    resolve_proportion_metric,
)
from deconvatac.methods import get_method


DEFAULT_METHOD_CONFIG_DIR = ROOT / "configs" / "methods"
BENCHMARK_PROTOCOL_PATH = ROOT / "docs" / "ShapeMix" / "benchmark_protocol.md"
OUTPUT_MANIFEST_FILENAME = "output_sha256.yaml"
OUTPUT_MANIFEST_SCHEMA_VERSION = 1
PROVENANCE_SCHEMA_VERSION = 1
RUN_REQUIRED_OUTPUTS = {
    "environment.txt",
    "inputs.yaml",
    "results/diagnostics.json",
    "results/proportions.csv",
}
EVALUATION_MODES = {"exact_truth", "prediction_only"}
CODE_SOURCE_PATTERNS = (
    "scripts/*.py",
    "src/deconvatac/**/*.py",
)
THREAD_ENVIRONMENT_VARIABLES = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
RESOURCE_GUARD_ENV = "DECONVATAC_RESOURCE_GUARD"
RESOURCE_MAX_ONE_MINUTE_LOAD = 6.0
RESOURCE_MIN_AVAILABLE_MEMORY_BYTES = 4 * 1024**3
RESOURCE_MAX_GPU_MEMORY_USED_MIB = 2048
RESOURCE_MAX_GPU_TEMPERATURE_C = 79
RESOURCE_MAX_DISPLAY_PROCESS_MEMORY_MIB = 512
RESOURCE_RECHECK_SECONDS = 30


def read_yaml(path: Optional[Union[str, Path]]) -> dict[str, Any]:
    if path is None:
        return {}
    with Path(path).open() as handle:
        return yaml.safe_load(handle) or {}


def build_run_id(
    dataset: str,
    modality: str,
    feature_set: str,
    method: str,
    method_run_id: Optional[str] = None,
) -> str:
    date = datetime.now().strftime("%Y-%m-%d")
    variant = method if method_run_id is None else method_run_id
    return f"{date}_{dataset}_{modality}_{feature_set}_{variant}"


def format_run_id(
    template: str,
    dataset: str,
    modality: str,
    feature_set: str,
    method: str,
    method_run_id: Optional[str] = None,
) -> str:
    variant = method if method_run_id is None else method_run_id
    return template.format(
        dataset=dataset,
        modality=modality,
        feature_set=feature_set,
        method=method,
        method_run_id=variant,
    )


def collect_environment(method: Optional[str] = None) -> dict[str, Any]:
    packages = {
        "deconvatac": "deconvatac",
        "anndata": "anndata",
        "muon": "muon",
        "numpy": "numpy",
        "pandas": "pandas",
        "scipy": "scipy",
        "yaml": "yaml",
    }
    # ShapeMix owns these optional dependencies. Do not import them merely to
    # inspect the environment for an unrelated method.
    if method is not None and method.lower() == "shapemix":
        packages.update({"torch": "torch", "pysam": "pysam"})
    versions: dict[str, str] = {}
    for package, module_name in packages.items():
        try:
            module = __import__(module_name)
            version = str(getattr(module, "__version__", "unknown"))
        except Exception as exc:
            version = f"unavailable: {exc}"
        versions[package] = version
    return {
        "python": sys.version.replace(chr(10), " "),
        "platform": platform.platform(),
        "packages": versions,
        "compute": collect_compute_environment(method),
    }


def collect_compute_environment(method: Optional[str] = None) -> dict[str, Any]:
    compute: dict[str, Any] = {
        "thread_environment": {
            name: os.environ.get(name) for name in THREAD_ENVIRONMENT_VARIABLES
        },
        "logical_cpu_count": psutil.cpu_count(logical=True),
        "physical_cpu_count": psutil.cpu_count(logical=False),
    }
    try:
        compute["cpu_affinity"] = psutil.Process().cpu_affinity()
    except (AttributeError, NotImplementedError, psutil.Error, OSError):
        compute["cpu_affinity"] = None
    if method is not None and method.lower() == "shapemix":
        try:
            import torch

            compute["torch_num_threads"] = int(torch.get_num_threads())
            compute["torch_num_interop_threads"] = int(torch.get_num_interop_threads())
        except Exception as exc:
            compute["torch_thread_error"] = str(exc)
    return compute


def write_environment(
    output_dir: Path,
    method: Optional[str] = None,
    metadata: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    environment = dict(metadata) if metadata is not None else collect_environment(method)
    lines = [
        f"python={environment['python']}",
        f"platform={environment['platform']}",
    ]
    lines.extend(f"{package}={version}" for package, version in environment["packages"].items())
    lines.append(f"compute={_canonical_config_json(environment.get('compute', {}))}")
    (output_dir / "environment.txt").write_text("\n".join(lines) + "\n")
    return environment


def validate_method_config(method: str, method_config: dict[str, Any]) -> None:
    configured_method = method_config.get("method")
    if configured_method is not None and configured_method != method:
        raise ValueError(f"Config method '{configured_method}' does not match method '{method}'.")


def _canonical_config_json(method_config: Mapping[str, Any]) -> str:
    """Return the stable representation used in hashes and CSV manifests."""
    return json.dumps(
        method_config,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def method_config_sha256(method_config: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_config_json(method_config).encode("utf-8")).hexdigest()


def _recorded_path(path: Optional[Union[str, Path]]) -> Optional[str]:
    if path is None:
        return None
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(ROOT.resolve()))
    except ValueError:
        return str(resolved)


def _source_sha256(path: Optional[Path]) -> Optional[str]:
    if path is None:
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_config_json(value).encode("utf-8")).hexdigest()


def _git_output(*arguments: str) -> Optional[str]:
    try:
        completed = subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()


def _campaign_code_provenance() -> dict[str, Any]:
    """Hash executable repository Python, including dirty/untracked files."""
    paths: set[Path] = set()
    for pattern in CODE_SOURCE_PATTERNS:
        paths.update(path for path in ROOT.glob(pattern) if path.is_file())
    files = {
        path.relative_to(ROOT).as_posix(): _source_sha256(path)
        for path in sorted(paths, key=lambda item: item.relative_to(ROOT).as_posix())
    }
    # Scope dirty-state reporting to executable/protocol/config inputs. New
    # result directories must not make concurrently launched shards disagree.
    status = _git_output(
        "status",
        "--short",
        "--untracked-files=normal",
        "--",
        "scripts",
        "src/deconvatac",
        "configs",
        "docs/ShapeMix",
        "requirements.txt",
        "pyproject.toml",
    )
    return {
        "git_commit": _git_output("rev-parse", "HEAD"),
        "git_worktree_dirty": None if status is None else bool(status),
        "git_status": None if status is None else status.splitlines(),
        "manifest_sha256": _canonical_sha256(files),
        "manifest_hash_encoding": "canonical_json_sorted_keys_compact_utf8",
        "selection": list(CODE_SOURCE_PATTERNS),
        "files": files,
    }


def _resolved_campaign_payload(
    experiment_config: Mapping[str, Any],
    jobs: Iterable[Mapping[str, Any]],
    metrics: Iterable[str],
    registry: Union[str, Path],
) -> dict[str, Any]:
    resolved_jobs = []
    for job in jobs:
        resolved_jobs.append(
            {
                "dataset": job["dataset"],
                "modality": job["modality"],
                "feature_set": job["feature_set"],
                "method": job["method"],
                "method_run_id": job["method_run_id"],
                "run_id": job["run_id"],
                "method_config": job["method_config"],
                "method_config_source": job["method_config_source"],
                "method_config_source_path": job["method_config_source_path"],
                "method_config_sha256": job["method_config_sha256"],
                "method_config_source_sha256": job["method_config_source_sha256"],
                "dataset_config_sha256": job["dataset_config_sha256"],
            }
        )
    registry_path = Path(registry).resolve()
    return {
        "experiment_config": copy.deepcopy(dict(experiment_config)),
        "resolved_jobs": resolved_jobs,
        "resolved_metrics": list(metrics),
        "registry": {
            "path": _recorded_path(registry_path),
            "source_sha256": _source_sha256(registry_path),
        },
    }


def build_execution_provenance(
    experiment_config_path: Path,
    experiment_config: Mapping[str, Any],
    jobs: list[Mapping[str, Any]],
    metrics: list[str],
    registry: Union[str, Path],
    shard: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    """Build one immutable provenance record shared by a campaign's runs."""
    resolved_payload = _resolved_campaign_payload(
        experiment_config,
        jobs,
        metrics,
        registry,
    )
    methods = {str(job["method"]).lower() for job in jobs}
    compute_method = "shapemix" if "shapemix" in methods else None
    registry_path = Path(registry).resolve()
    provenance = {
        "schema_version": PROVENANCE_SCHEMA_VERSION,
        "experiment_config": {
            "path": _recorded_path(experiment_config_path),
            "source_sha256": _source_sha256(experiment_config_path),
            "resolved_sha256": _canonical_sha256(resolved_payload),
            "resolved_hash_encoding": "canonical_json_sorted_keys_compact_utf8",
        },
        "benchmark_protocol": {
            "path": _recorded_path(BENCHMARK_PROTOCOL_PATH),
            "source_sha256": _source_sha256(BENCHMARK_PROTOCOL_PATH),
        },
        "registry": {
            "path": _recorded_path(registry_path),
            "source_sha256": _source_sha256(registry_path),
        },
        "code": _campaign_code_provenance(),
        "compute_environment": collect_compute_environment(compute_method),
    }
    if shard is not None:
        provenance["shard"] = copy.deepcopy(dict(shard))
    return provenance


def build_single_run_provenance(
    invocation: Mapping[str, Any],
    registry: Union[str, Path],
    method: str,
) -> dict[str, Any]:
    """Use the campaign schema for an auditable single-run invocation."""
    registry_path = Path(registry).resolve()
    return {
        "schema_version": PROVENANCE_SCHEMA_VERSION,
        "experiment_config": {
            "path": None,
            "source_sha256": None,
            "resolved_sha256": _canonical_sha256(dict(invocation)),
            "resolved_hash_encoding": "canonical_json_sorted_keys_compact_utf8",
            "mode": "single_run_cli",
        },
        "benchmark_protocol": {
            "path": _recorded_path(BENCHMARK_PROTOCOL_PATH),
            "source_sha256": _source_sha256(BENCHMARK_PROTOCOL_PATH),
        },
        "registry": {
            "path": _recorded_path(registry_path),
            "source_sha256": _source_sha256(registry_path),
        },
        "code": _campaign_code_provenance(),
        "compute_environment": collect_compute_environment(method),
    }


def _ru_maxrss_measurement() -> dict[str, Any]:
    try:
        import resource

        raw_value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    except (ImportError, OSError, ValueError):
        return {
            "raw_value": None,
            "units": None,
            "bytes": None,
            "semantics": "unavailable",
        }
    if sys.platform == "darwin":
        units = "bytes"
        byte_value = raw_value
    else:
        # Linux reports KiB. Other supported Unix Python builds generally use
        # KiB as well; retain the platform and raw units so this is auditable.
        units = "kibibytes"
        byte_value = raw_value * 1024
    return {
        "raw_value": raw_value,
        "units": units,
        "bytes": byte_value,
        "platform": sys.platform,
        "semantics": "process_lifetime_high_water_mark",
    }


class _PeakRSSMonitor:
    """Sample driver-plus-child RSS during one run without resetting state."""

    def __init__(self, interval_seconds: float = 0.02):
        self.interval_seconds = interval_seconds
        self._process = psutil.Process()
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._peak_bytes = 0
        self._cpu_seconds_start = 0.0

    def _cpu_seconds(self) -> float:
        processes = [self._process]
        try:
            processes.extend(self._process.children(recursive=True))
        except (psutil.Error, OSError):
            pass
        total = 0.0
        for process in processes:
            try:
                times = process.cpu_times()
                total += float(times.user) + float(times.system)
            except (psutil.Error, OSError):
                continue
        return total

    def _sample(self) -> None:
        processes = [self._process]
        try:
            processes.extend(self._process.children(recursive=True))
        except (psutil.Error, OSError):
            pass
        total = 0
        for process in processes:
            try:
                total += int(process.memory_info().rss)
            except (psutil.Error, OSError):
                continue
        self._peak_bytes = max(self._peak_bytes, total)

    def _run(self) -> None:
        while not self._stop_event.wait(self.interval_seconds):
            self._sample()

    def start(self) -> float:
        self._sample()
        self._cpu_seconds_start = self._cpu_seconds()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return time.perf_counter()

    def stop(self, started_at: float, scope: str) -> dict[str, Any]:
        elapsed = time.perf_counter() - started_at
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=max(1.0, self.interval_seconds * 4))
        self._sample()
        process_cpu_seconds = max(0.0, self._cpu_seconds() - self._cpu_seconds_start)
        ru_maxrss = _ru_maxrss_measurement()
        return {
            "wall_runtime_seconds": elapsed,
            "process_cpu_seconds": process_cpu_seconds,
            "average_process_cpu_cores": process_cpu_seconds / elapsed,
            "average_process_cpu_percent_of_one_core": 100.0
            * process_cpu_seconds
            / elapsed,
            "cpu_measurement": {
                "source": "psutil_process_tree_cpu_times",
                "semantics": "driver-plus-live-child user and system CPU time",
                "average_core_semantics": "CPU seconds divided by wall seconds",
            },
            "peak_rss_bytes": self._peak_bytes,
            "peak_rss_mb": self._peak_bytes / (1024.0 * 1024.0),
            "scope": scope,
            "clock": "time.perf_counter",
            "rss_measurement": {
                "source": "psutil_sampled_process_tree",
                "sample_interval_seconds": self.interval_seconds,
                "semantics": "maximum sampled sum of driver and recursive child RSS",
                "platform": sys.platform,
                "ru_maxrss": ru_maxrss,
            },
        }


def _query_gpu_resource_state() -> dict[str, Any]:
    applications = subprocess.run(
        [
            "nvidia-smi",
            "--query-compute-apps=pid,process_name,used_memory",
            "--format=csv,noheader,nounits",
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=15,
    ).stdout
    gpu_lines = subprocess.run(
        [
            "nvidia-smi",
            "--query-gpu=memory.used,memory.total,utilization.gpu,temperature.gpu",
            "--format=csv,noheader,nounits",
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=15,
    ).stdout.splitlines()
    if len(gpu_lines) != 1:
        raise RuntimeError(f"Expected one GPU from nvidia-smi, found {len(gpu_lines)}.")
    gpu_fields = [field.strip() for field in gpu_lines[0].split(",")]
    if len(gpu_fields) != 4:
        raise RuntimeError("Unexpected nvidia-smi GPU state response.")
    unrelated = []
    own_gpu_memory_mib = 0
    for line in applications.splitlines():
        if not line.strip():
            continue
        fields = [field.strip() for field in line.split(",", 2)]
        if len(fields) != 3:
            raise RuntimeError("Unexpected nvidia-smi application response.")
        pid, process_name, used_memory = int(fields[0]), fields[1], int(fields[2])
        if pid == os.getpid():
            own_gpu_memory_mib += used_memory
            continue
        if (
            Path(process_name).name == "gnome-remote-desktop-daemon"
            and used_memory <= RESOURCE_MAX_DISPLAY_PROCESS_MEMORY_MIB
        ):
            continue
        unrelated.append(
            {
                "pid": pid,
                "process_name": process_name,
                "used_memory_mib": used_memory,
            }
        )
    gpu_used_mib = int(gpu_fields[0])
    return {
        "gpu_memory_used_mib": gpu_used_mib,
        "gpu_external_memory_used_mib": max(0, gpu_used_mib - own_gpu_memory_mib),
        "own_gpu_memory_used_mib": own_gpu_memory_mib,
        "gpu_memory_total_mib": int(gpu_fields[1]),
        "gpu_utilization_percent": int(gpu_fields[2]),
        "gpu_temperature_c": int(gpu_fields[3]),
        "unrelated_gpu_processes": unrelated,
    }


def _resource_gate_state() -> dict[str, Any]:
    if os.environ.get(RESOURCE_GUARD_ENV) != "1":
        return {"enabled": False, "passed": True}
    state: dict[str, Any] = {
        "enabled": True,
        "captured_at": datetime.now().astimezone().isoformat(),
        "one_minute_load": float(os.getloadavg()[0]),
        "available_memory_bytes": int(psutil.virtual_memory().available),
        "process_nice": int(psutil.Process().nice()),
        "process_ionice": str(psutil.Process().ionice()),
        "thresholds": {
            "one_minute_load_max_exclusive": RESOURCE_MAX_ONE_MINUTE_LOAD,
            "available_memory_bytes_min": RESOURCE_MIN_AVAILABLE_MEMORY_BYTES,
            "gpu_external_memory_used_mib_max": RESOURCE_MAX_GPU_MEMORY_USED_MIB,
            "gpu_temperature_c_max": RESOURCE_MAX_GPU_TEMPERATURE_C,
        },
    }
    try:
        state.update(_query_gpu_resource_state())
        state["passed"] = bool(
            state["one_minute_load"] < RESOURCE_MAX_ONE_MINUTE_LOAD
            and state["available_memory_bytes"] >= RESOURCE_MIN_AVAILABLE_MEMORY_BYTES
            and state["gpu_external_memory_used_mib"] <= RESOURCE_MAX_GPU_MEMORY_USED_MIB
            and state["gpu_temperature_c"] <= RESOURCE_MAX_GPU_TEMPERATURE_C
            and not state["unrelated_gpu_processes"]
        )
    except Exception as exc:
        state["passed"] = False
        state["gpu_query_error"] = str(exc)
    return state


def _wait_for_resource_gate() -> dict[str, Any]:
    while True:
        state = _resource_gate_state()
        if state["passed"]:
            print(f"resource_job_preflight {json.dumps(state, sort_keys=True)}", flush=True)
            return state
        print(f"resource_job_wait {json.dumps(state, sort_keys=True)}", flush=True)
        time.sleep(RESOURCE_RECHECK_SECONDS)


def _execution_metadata(method_config: Mapping[str, Any]) -> dict[str, Any]:
    params = method_config.get("params", {})
    if not isinstance(params, Mapping):
        params = {}
    seed = params.get("seed")
    device = params.get("device")
    dtype = params.get("dtype")
    return {
        "seed": seed,
        "device": device,
        "dtype": dtype,
        "determinism": {
            "seed": seed,
            "device": device,
            "dtype": dtype,
            "source": "resolved_method_config",
        },
    }


def _configure_shapemix_torch_threads(method: str) -> None:
    """Enforce the co-tenant one-thread policy before ShapeMix imports work."""
    if str(method).lower() != "shapemix":
        return

    import torch

    torch.set_num_threads(1)
    if torch.get_num_interop_threads() == 1:
        return
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError as exc:
        if torch.get_num_interop_threads() != 1:
            raise RuntimeError(
                "ShapeMix requires torch inter-op threads=1 under the co-tenant policy."
            ) from exc

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
    method_run_id: Optional[str] = None,
    method_config_source_path: Optional[Union[str, Path]] = None,
    method_config_source: str = "provided",
    resolved_method_config_sha256: Optional[str] = None,
    method_config_source_sha256: Optional[str] = None,
    execution_provenance: Optional[Mapping[str, Any]] = None,
    resolved_dataset_config_sha256: Optional[str] = None,
) -> Path:
    method_config = copy.deepcopy(method_config or {})
    validate_method_config(method, method_config)
    resolved_method_run_id = method if method_run_id is None else method_run_id
    if not isinstance(resolved_method_run_id, str) or not resolved_method_run_id.strip():
        raise ValueError("method_run_id must be a nonempty string.")
    resolved_method_run_id = resolved_method_run_id.strip()
    config_sha256 = resolved_method_config_sha256 or method_config_sha256(method_config)
    source_path = _recorded_path(method_config_source_path)
    if method_config_source_sha256 is None and method_config_source_path is not None:
        method_config_source_sha256 = _source_sha256(Path(method_config_source_path).resolve())

    resolved_run_id = run_id or build_run_id(
        dataset,
        modality,
        feature_set,
        method,
        method_run_id=resolved_method_run_id,
    )
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

    _configure_shapemix_torch_threads(method)
    method_cls = get_method(method)
    method_instance = method_cls(**method_config.get("params", {}))
    result = method_instance.run(data)
    environment = collect_environment(method)

    extra_metadata = {
        "status": "pending_evaluation" if execution_provenance is not None else "success",
        "run_id": resolved_run_id,
        "method": method,
        "method_run_id": resolved_method_run_id,
        "method_config": method_config,
        "method_config_source": method_config_source,
        "method_config_source_path": source_path,
        "method_config_sha256": config_sha256,
        "method_config_source_sha256": method_config_source_sha256,
        "dataset_config_sha256": resolved_dataset_config_sha256
        or _canonical_sha256(data.metadata.get("dataset_config", {})),
        "environment": environment,
        "software_versions": environment["packages"],
        "proportion_evaluation": {
            "contract_version": PROPORTION_CONTRACT_VERSION,
            "row_sum_atol": PROPORTION_ROW_SUM_ATOL,
            "cell_types": data.cell_types,
            "evidence": "exact_truth" if data.truth is not None else "prediction_only",
        },
        "inputs": {
            "dataset_id": dataset,
            "modality": modality,
            "feature_set": feature_set,
            "registry": str(registry),
        },
        **_execution_metadata(method_config),
    }
    if execution_provenance is not None:
        extra_metadata["execution_provenance"] = copy.deepcopy(dict(execution_provenance))
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
    write_environment(output_dir, method=method, metadata=environment)

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


def _require_config_mapping(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{context} must resolve to a mapping.")
    return copy.deepcopy(dict(value))


def _read_config_mapping(path: Path, context: str) -> dict[str, Any]:
    return _require_config_mapping(read_yaml(path), context)


def _resolve_config_spec(
    method: str,
    spec: Any,
    *,
    source_label: str,
) -> dict[str, Any]:
    source_path: Optional[Path]
    source: str

    if isinstance(spec, Mapping):
        config = _require_config_mapping(spec, f"Config for method '{method}'")
        source_path = None
        source = f"inline_{source_label}"
    elif isinstance(spec, (str, Path)) and str(spec) != "default":
        source_path = resolve_config_path(spec).resolve()
        config = _read_config_mapping(source_path, f"Config file for method '{method}'")
        source = source_label
    elif spec is None or spec == "default":
        candidate = DEFAULT_METHOD_CONFIG_DIR / f"{method}.yaml"
        if candidate.exists():
            source_path = candidate.resolve()
            config = _read_config_mapping(source_path, f"Default config for method '{method}'")
            source = "default"
        else:
            source_path = None
            config = {}
            source = "implicit_empty_default"
    else:
        raise TypeError(
            f"Config for method '{method}' must be a path, an inline mapping, 'default', or null."
        )

    validate_method_config(method, config)
    return {
        "method_config": config,
        "method_config_source": source,
        "method_config_source_path": _recorded_path(source_path),
        "method_config_sha256": method_config_sha256(config),
        "method_config_source_sha256": _source_sha256(source_path),
        "_method_config_source_path_resolved": str(source_path) if source_path is not None else None,
    }


def resolve_method_config(method: str, experiment_config: dict[str, Any]) -> dict[str, Any]:
    """Resolve a legacy ``methods`` config, preserving the original return API."""
    method_configs = experiment_config.get("method_configs", {})
    if not isinstance(method_configs, Mapping):
        raise TypeError("method_configs must be a mapping from method names to configs.")
    method_spec = method_configs.get(method)
    return _resolve_config_spec(
        method,
        method_spec,
        source_label="legacy_method_configs",
    )["method_config"]


def resolve_config_path(path: Union[str, Path]) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return ROOT / path


def _nonempty_identifier(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a nonempty string.")
    return value.strip()


def resolve_method_runs(experiment_config: dict[str, Any]) -> list[dict[str, Any]]:
    """Normalize legacy and named method variants into fully resolved records."""
    has_methods = "methods" in experiment_config
    has_method_runs = "method_runs" in experiment_config
    if has_methods == has_method_runs:
        raise ValueError("Experiment config must define exactly one of 'methods' or 'method_runs'.")
    if has_method_runs and "method_configs" in experiment_config:
        raise ValueError(
            "method_configs is only valid with the legacy 'methods' schema; "
            "put each config on its method_runs entry."
        )

    resolved: list[dict[str, Any]] = []
    if has_methods:
        methods = _as_list(experiment_config.get("methods"))
        if not methods:
            raise ValueError("Experiment config must define at least one method.")
        method_configs = experiment_config.get("method_configs", {})
        if not isinstance(method_configs, Mapping):
            raise TypeError("method_configs must be a mapping from method names to configs.")
        for index, raw_method in enumerate(methods):
            method = _nonempty_identifier(raw_method, f"methods[{index}]").lower()
            config_record = _resolve_config_spec(
                method,
                method_configs.get(method),
                source_label="legacy_method_configs",
            )
            resolved.append(
                {
                    "method": method,
                    "method_run_id": method,
                    **config_record,
                }
            )
    else:
        raw_method_runs = experiment_config.get("method_runs")
        if not isinstance(raw_method_runs, list) or not raw_method_runs:
            raise ValueError("Experiment config must define at least one method_runs entry.")
        allowed_keys = {"id", "method", "config"}
        for index, raw_run in enumerate(raw_method_runs):
            if not isinstance(raw_run, Mapping):
                raise TypeError(f"method_runs[{index}] must be a mapping.")
            unknown_keys = sorted(set(raw_run).difference(allowed_keys))
            if unknown_keys:
                raise ValueError(
                    f"method_runs[{index}] has unknown key(s): {', '.join(map(str, unknown_keys))}."
                )
            method_run_id = _nonempty_identifier(raw_run.get("id"), f"method_runs[{index}].id")
            method = _nonempty_identifier(raw_run.get("method"), f"method_runs[{index}].method").lower()
            config_record = _resolve_config_spec(
                method,
                raw_run.get("config"),
                source_label="method_runs_config",
            )
            resolved.append(
                {
                    "method": method,
                    "method_run_id": method_run_id,
                    **config_record,
                }
            )

    method_run_ids = [record["method_run_id"] for record in resolved]
    duplicate_ids = sorted({value for value in method_run_ids if method_run_ids.count(value) > 1})
    if duplicate_ids:
        raise ValueError(f"Duplicate method run id(s): {', '.join(duplicate_ids)}.")

    # Validate registration while still in preflight. get_method imports only
    # the selected adapter; unrelated methods do not traverse ShapeMix.
    for method in dict.fromkeys(record["method"] for record in resolved):
        get_method(method)
    return resolved


def iter_experiment_jobs(
    experiment_config: dict[str, Any],
    registry: Union[str, Path],
) -> Iterable[dict[str, Any]]:
    datasets = _as_list(experiment_config.get("datasets"))
    method_runs = resolve_method_runs(experiment_config)
    configured_modalities = experiment_config.get("modalities")
    configured_feature_sets = experiment_config.get("feature_sets")
    skip_missing = experiment_config.get("skip_missing", True)

    if not datasets:
        raise ValueError("Experiment config must define at least one dataset.")
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

                for method_run in method_runs:
                    yield {
                        "dataset": str(dataset),
                        "modality": str(modality),
                        "feature_set": str(feature_set),
                        "dataset_config_sha256": _canonical_sha256(dataset_config),
                        **copy.deepcopy(method_run),
                    }


def read_run_metadata(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "run.yaml"
    if not path.exists():
        return {}
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def read_inputs_metadata(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "inputs.yaml"
    if not path.exists():
        return {}
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def _write_yaml(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w") as handle:
        yaml.safe_dump(dict(value), handle, sort_keys=False)
    temporary.replace(path)


def _artifact_paths(run_dir: Path) -> dict[str, Path]:
    excluded = {"run.yaml", OUTPUT_MANIFEST_FILENAME}
    artifacts: dict[str, Path] = {}
    for path in sorted(run_dir.rglob("*"), key=lambda item: str(item.relative_to(run_dir))):
        relative = path.relative_to(run_dir).as_posix()
        if relative in excluded or not path.is_file():
            continue
        if path.is_symlink():
            raise ValueError(f"Run output contains unsupported symbolic link: {path}")
        artifacts[relative] = path
    return artifacts


def _write_output_manifest(
    run_dir: Path,
    run_metadata_without_manifest: Mapping[str, Any],
) -> dict[str, Any]:
    files = {
        relative: _source_sha256(path)
        for relative, path in _artifact_paths(run_dir).items()
    }
    payload = {
        "schema_version": OUTPUT_MANIFEST_SCHEMA_VERSION,
        "algorithm": "sha256",
        "exclusions": ["run.yaml", OUTPUT_MANIFEST_FILENAME],
        "run_metadata_sha256": _canonical_sha256(run_metadata_without_manifest),
        "run_metadata_hash_encoding": "canonical_json_sorted_keys_compact_utf8",
        "files": files,
    }
    manifest_path = run_dir / OUTPUT_MANIFEST_FILENAME
    _write_yaml(manifest_path, payload)
    return {
        "path": OUTPUT_MANIFEST_FILENAME,
        "schema_version": OUTPUT_MANIFEST_SCHEMA_VERSION,
        "sha256": _source_sha256(manifest_path),
    }


def _finalize_run_directory(
    run_dir: Path,
    *,
    status: str,
    performance: Mapping[str, Any],
    execution_provenance: Mapping[str, Any],
    execution_action: str,
    fallback_metadata: Optional[Mapping[str, Any]] = None,
    error: Optional[str] = None,
) -> dict[str, Any]:
    """Finalize success/failure metadata, then hash every non-metadata artifact."""
    run_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_run_metadata(run_dir)
    if not metadata and fallback_metadata is not None:
        metadata = copy.deepcopy(dict(fallback_metadata))
    metadata.update(
        {
            "status": status,
            "execution_action": execution_action,
            "execution_provenance": copy.deepcopy(dict(execution_provenance)),
            "performance": copy.deepcopy(dict(performance)),
        }
    )
    if error is None:
        metadata.pop("error", None)
    else:
        metadata["error"] = error
    metadata.pop("output_manifest", None)
    _write_yaml(run_dir / "run.yaml", metadata)
    metadata["output_manifest"] = _write_output_manifest(run_dir, metadata)
    _write_yaml(run_dir / "run.yaml", metadata)
    return metadata


def _manifest_payload(run_dir: Path) -> dict[str, Any]:
    path = run_dir / OUTPUT_MANIFEST_FILENAME
    if not path.exists():
        raise ValueError(f"completed run is missing {path}")
    if path.is_symlink():
        raise ValueError(f"output manifest cannot be a symbolic link: {path}")
    with path.open() as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, Mapping):
        raise ValueError(f"output manifest must be a mapping: {path}")
    return dict(payload)


def _validate_output_manifest(run_dir: Path, metadata: Mapping[str, Any]) -> None:
    recorded = metadata.get("output_manifest")
    if not isinstance(recorded, Mapping):
        raise ValueError("run.yaml has no output_manifest mapping")
    if recorded.get("path") != OUTPUT_MANIFEST_FILENAME:
        raise ValueError("run.yaml records an unsupported output manifest path")
    if recorded.get("schema_version") != OUTPUT_MANIFEST_SCHEMA_VERSION:
        raise ValueError("run.yaml records an unsupported output manifest schema")
    manifest_path = run_dir / OUTPUT_MANIFEST_FILENAME
    if recorded.get("sha256") != _source_sha256(manifest_path):
        raise ValueError("output manifest SHA256 does not match run.yaml")

    payload = _manifest_payload(run_dir)
    expected_header = {
        "schema_version": OUTPUT_MANIFEST_SCHEMA_VERSION,
        "algorithm": "sha256",
        "exclusions": ["run.yaml", OUTPUT_MANIFEST_FILENAME],
    }
    for key, expected in expected_header.items():
        if payload.get(key) != expected:
            raise ValueError(f"output manifest has invalid {key}")
    metadata_without_manifest = copy.deepcopy(dict(metadata))
    metadata_without_manifest.pop("output_manifest", None)
    if payload.get("run_metadata_hash_encoding") != "canonical_json_sorted_keys_compact_utf8":
        raise ValueError("output manifest has invalid run metadata hash encoding")
    if payload.get("run_metadata_sha256") != _canonical_sha256(metadata_without_manifest):
        raise ValueError("run.yaml canonical metadata SHA256 mismatch")
    recorded_files = payload.get("files")
    if not isinstance(recorded_files, Mapping):
        raise ValueError("output manifest files must be a mapping")
    actual_paths = _artifact_paths(run_dir)
    if set(recorded_files) != set(actual_paths):
        missing = sorted(set(recorded_files).difference(actual_paths))
        extra = sorted(set(actual_paths).difference(recorded_files))
        raise ValueError(
            f"output manifest file set mismatch; missing={missing!r}, extra={extra!r}"
        )
    for relative, path in actual_paths.items():
        if recorded_files[relative] != _source_sha256(path):
            raise ValueError(f"output SHA256 mismatch for {relative}")


def _validate_resume_run(
    run_dir: Path,
    job: Mapping[str, Any],
    execution_provenance: Mapping[str, Any],
) -> dict[str, Any]:
    """Fail closed unless an existing run is exactly reusable."""
    if (run_dir / "run.yaml").is_symlink():
        raise ValueError("run.yaml cannot be a symbolic link")
    metadata = read_run_metadata(run_dir)
    if not metadata:
        raise ValueError("missing or empty run.yaml")
    expected_identity = {
        "run_id": job["run_id"],
        "dataset_id": job["dataset"],
        "modality": job["modality"],
        "feature_set": job["feature_set"],
        "method": job["method"],
        "method_run_id": job["method_run_id"],
        "method_config_sha256": job["method_config_sha256"],
        "method_config_source_sha256": job["method_config_source_sha256"],
        "dataset_config_sha256": job["dataset_config_sha256"],
    }
    for key, expected in expected_identity.items():
        if metadata.get(key) != expected:
            raise ValueError(
                f"run.yaml {key} mismatch: expected {expected!r}, found {metadata.get(key)!r}"
            )
    if metadata.get("status") != "success":
        raise ValueError(f"run status is {metadata.get('status')!r}, not 'success'")
    if _canonical_config_json(metadata.get("method_config", {})) != _canonical_config_json(
        job["method_config"]
    ):
        raise ValueError("resolved method config does not match")
    if metadata.get("execution_provenance") != execution_provenance:
        raise ValueError("execution provenance does not match the current campaign")

    missing_required = sorted(
        relative for relative in RUN_REQUIRED_OUTPUTS if not (run_dir / relative).is_file()
    )
    if missing_required:
        raise ValueError(f"required output(s) missing: {missing_required!r}")
    inputs_metadata = read_inputs_metadata(run_dir)
    truth_is_declared = inputs_metadata.get("truth") == "results/truth.csv"
    truth_path = run_dir / "results" / "truth.csv"
    if truth_is_declared != truth_path.is_file():
        raise ValueError(
            "inputs.yaml truth declaration and results/truth.csv presence do not agree"
        )
    performance = metadata.get("performance")
    if not isinstance(performance, Mapping):
        raise ValueError("run.yaml has no performance mapping")
    for field in ("wall_runtime_seconds", "peak_rss_bytes", "peak_rss_mb"):
        value = performance.get(field)
        if not isinstance(value, (int, float)) or isinstance(value, bool) or value < 0:
            raise ValueError(f"run.yaml performance.{field} is invalid")
    _validate_output_manifest(run_dir, metadata)
    return metadata


def _flatten_execution_metadata(metadata: Mapping[str, Any]) -> dict[str, Any]:
    provenance = metadata.get("execution_provenance")
    if not isinstance(provenance, Mapping):
        provenance = {}
    experiment = provenance.get("experiment_config")
    protocol = provenance.get("benchmark_protocol")
    registry = provenance.get("registry")
    code = provenance.get("code")
    shard = provenance.get("shard")
    performance = metadata.get("performance")
    output_manifest = metadata.get("output_manifest")
    if not isinstance(experiment, Mapping):
        experiment = {}
    if not isinstance(protocol, Mapping):
        protocol = {}
    if not isinstance(registry, Mapping):
        registry = {}
    if not isinstance(code, Mapping):
        code = {}
    if not isinstance(shard, Mapping):
        shard = {}
    if not isinstance(performance, Mapping):
        performance = {}
    if not isinstance(output_manifest, Mapping):
        output_manifest = {}
    rss_measurement = performance.get("rss_measurement", {})
    return {
        "experiment_config_path": experiment.get("path"),
        "experiment_config_source_sha256": experiment.get("source_sha256"),
        "experiment_config_resolved_sha256": experiment.get("resolved_sha256"),
        "benchmark_protocol_path": protocol.get("path"),
        "benchmark_protocol_sha256": protocol.get("source_sha256"),
        "registry_path": registry.get("path"),
        "registry_source_sha256": registry.get("source_sha256"),
        "git_commit": code.get("git_commit"),
        "git_worktree_dirty": code.get("git_worktree_dirty"),
        "code_manifest_sha256": code.get("manifest_sha256"),
        "shard_index": shard.get("index"),
        "shard_count": shard.get("count"),
        "shard_assignment": shard.get("assignment"),
        "wall_runtime_seconds": performance.get("wall_runtime_seconds"),
        "peak_rss_bytes": performance.get("peak_rss_bytes"),
        "peak_rss_mb": performance.get("peak_rss_mb"),
        "performance_scope": performance.get("scope"),
        "rss_measurement": _canonical_config_json(rss_measurement),
        "compute_environment": _canonical_config_json(
            performance.get("compute_environment", provenance.get("compute_environment", {}))
        ),
        "output_manifest_path": output_manifest.get("path"),
        "output_manifest_sha256": output_manifest.get("sha256"),
        "execution_action": metadata.get("execution_action"),
    }


def evaluate_run(
    run_dir: Path,
    metrics: list[str],
    evaluation_mode: str = "exact_truth",
) -> list[dict[str, Any]]:
    if evaluation_mode not in EVALUATION_MODES:
        raise ValueError(f"Unknown evaluation_mode {evaluation_mode!r}.")
    metadata = read_run_metadata(run_dir)
    inputs_metadata = read_inputs_metadata(run_dir)
    predicted = pd.read_csv(run_dir / "results" / "proportions.csv", index_col=0)
    cell_types = declared_cell_types_from_metadata(metadata, inputs_metadata)
    universe_json = _canonical_config_json(list(cell_types))

    expected_columns = list(cell_types)
    if list(predicted.columns) != expected_columns:
        raise ValueError(
            "Prediction columns must exactly match the declared cell_types order: "
            f"expected {expected_columns!r}, found {list(predicted.columns)!r}."
        )
    if predicted.empty:
        raise ValueError("Prediction output must contain at least one spot.")
    if not pd.Index(predicted.index).is_unique:
        raise ValueError("Prediction spot identifiers must be unique.")

    truth_is_declared = inputs_metadata.get("truth") == "results/truth.csv"
    truth_path = run_dir / "results" / "truth.csv"
    if truth_is_declared != truth_path.is_file():
        raise ValueError(
            "inputs.yaml truth declaration and results/truth.csv presence do not agree"
        )

    if evaluation_mode == "prediction_only":
        if metrics:
            raise ValueError("prediction_only evaluation cannot request exact-truth metrics.")
        if truth_is_declared:
            raise ValueError(
                "prediction_only evaluation requires a dataset descriptor without exact truth."
            )
        return []

    if not metrics:
        raise ValueError("exact_truth evaluation requires at least one metric.")
    if not truth_is_declared:
        raise ValueError(
            "exact_truth evaluation requires inputs.yaml to declare results/truth.csv."
        )
    truth = pd.read_csv(truth_path, index_col=0)

    rows = []
    for metric_name in metrics:
        evaluation = evaluate_proportion_metric(metric_name, truth, predicted, cell_types)
        rows.append(
            {
                "run_id": metadata.get("run_id", run_dir.name),
                "dataset_id": metadata.get("dataset_id"),
                "modality": metadata.get("modality"),
                "feature_set": metadata.get("feature_set"),
                "method": metadata.get("method"),
                "method_run_id": metadata.get("method_run_id", metadata.get("method")),
                "method_config": _canonical_config_json(metadata.get("method_config", {})),
                "method_config_source": metadata.get("method_config_source"),
                "method_config_source_path": metadata.get("method_config_source_path"),
                "method_config_sha256": metadata.get("method_config_sha256"),
                "method_config_source_sha256": metadata.get("method_config_source_sha256"),
                "dataset_config_sha256": metadata.get("dataset_config_sha256"),
                "seed": metadata.get("seed"),
                "device": metadata.get("device"),
                "dtype": metadata.get("dtype"),
                "determinism": _canonical_config_json(metadata.get("determinism", {})),
                **_flatten_execution_metadata(metadata),
                "status": "success",
                "metric": evaluation.metric_id,
                "metric_name": evaluation.metric_name,
                "metric_version": evaluation.metric_version,
                "evaluation_contract_version": evaluation.contract_version,
                "row_sum_atol": evaluation.row_sum_atol,
                "cell_type_universe": universe_json,
                "n_spots": evaluation.n_spots,
                "n_cell_types": evaluation.n_cell_types,
                "value": evaluation.value,
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


def _manifest_config_fields(job: Mapping[str, Any]) -> dict[str, Any]:
    execution = _execution_metadata(job["method_config"])
    return {
        "method_config": _canonical_config_json(job["method_config"]),
        "method_config_source": job["method_config_source"],
        "method_config_source_path": job["method_config_source_path"],
        "method_config_sha256": job["method_config_sha256"],
        "method_config_source_sha256": job["method_config_source_sha256"],
        "dataset_config_sha256": job["dataset_config_sha256"],
        "seed": execution["seed"],
        "device": execution["device"],
        "dtype": execution["dtype"],
        "determinism": _canonical_config_json(execution["determinism"]),
    }


def _manifest_row(
    job: Mapping[str, Any],
    batch_dir: Path,
    status: str,
    metadata: Optional[Mapping[str, Any]] = None,
    execution_action: Optional[str] = None,
) -> dict[str, Any]:
    row = {
        "dataset": job["dataset"],
        "dataset_id": job["dataset"],
        "modality": job["modality"],
        "feature_set": job["feature_set"],
        "method": job["method"],
        "method_run_id": job["method_run_id"],
        "run_id": job["run_id"],
        "run_dir": str(batch_dir / job["run_id"]),
        "status": status,
        **_manifest_config_fields(job),
    }
    if metadata is not None:
        row.update(_flatten_execution_metadata(metadata))
    if execution_action is not None:
        row["execution_action"] = execution_action
    return row


def resolve_experiment_jobs(
    experiment_config: dict[str, Any],
    registry: Union[str, Path],
    run_id_template: str,
) -> list[dict[str, Any]]:
    """Resolve the full campaign grid and its stable run IDs without writing."""
    if not isinstance(run_id_template, str) or not run_id_template:
        raise ValueError("run_id_template must be a nonempty string.")
    jobs = list(iter_experiment_jobs(experiment_config, registry=registry))
    if not jobs:
        raise ValueError("Experiment config did not produce any runnable jobs.")
    for job in jobs:
        run_id = format_run_id(
            run_id_template,
            dataset=job["dataset"],
            modality=job["modality"],
            feature_set=job["feature_set"],
            method=job["method"],
            method_run_id=job["method_run_id"],
        )
        if not run_id or not run_id.strip():
            raise ValueError("run_id_template resolved to an empty run ID.")
        run_id_path = Path(run_id)
        if run_id_path.is_absolute() or run_id_path.name != run_id or run_id in {".", ".."}:
            raise ValueError(
                f"Resolved run ID '{run_id}' must be a single relative directory name."
            )
        job["run_id"] = run_id
    run_ids = [job["run_id"] for job in jobs]
    duplicate_run_ids = sorted({value for value in run_ids if run_ids.count(value) > 1})
    if duplicate_run_ids:
        raise ValueError(f"Duplicate resolved run ID(s): {', '.join(duplicate_run_ids)}.")
    return jobs


def select_experiment_shard(
    jobs: list[dict[str, Any]],
    *,
    shard_index: Optional[int],
    shard_count: Optional[int],
) -> tuple[list[dict[str, Any]], Optional[dict[str, Any]]]:
    """Assign complete dataset/modality/feature units by stable sorted modulo."""
    if (shard_index is None) != (shard_count is None):
        raise ValueError("shard_index and shard_count must be provided together.")
    if shard_index is None:
        return jobs, None
    if isinstance(shard_index, bool) or not isinstance(shard_index, int):
        raise ValueError("shard_index must be an integer.")
    if isinstance(shard_count, bool) or not isinstance(shard_count, int) or shard_count < 1:
        raise ValueError("shard_count must be a positive integer.")
    if shard_index < 0 or shard_index >= shard_count:
        raise ValueError("shard_index must satisfy 0 <= shard_index < shard_count.")

    key_fields = ("dataset", "modality", "feature_set")
    units = sorted({tuple(str(job[field]) for field in key_fields) for job in jobs})
    selected_units = {
        unit for position, unit in enumerate(units) if position % shard_count == shard_index
    }
    selected = [
        job
        for job in jobs
        if tuple(str(job[field]) for field in key_fields) in selected_units
    ]
    if not selected:
        raise ValueError(
            f"Shard {shard_index} of {shard_count} contains no experiment jobs."
        )
    width = max(2, len(str(shard_count - 1)))
    suffix = f"__shard_{shard_index:0{width}d}_of_{shard_count:0{width}d}"
    return selected, {
        "index": shard_index,
        "count": shard_count,
        "index_base": 0,
        "assignment": "sorted_dataset_modality_feature_set_modulo",
        "unit_key_fields": list(key_fields),
        "selected_units": [list(unit) for unit in sorted(selected_units)],
        "run_group_suffix": suffix,
    }


def preflight_experiment_jobs(
    experiment_config: dict[str, Any],
    registry: Union[str, Path],
    run_id_template: str,
    batch_dir: Path,
    overwrite: bool,
    resume: bool = False,
    resolved_jobs: Optional[list[dict[str, Any]]] = None,
) -> list[dict[str, Any]]:
    """Resolve every job and reject run-directory collisions without writing."""
    jobs = (
        resolve_experiment_jobs(experiment_config, registry, run_id_template)
        if resolved_jobs is None
        else copy.deepcopy(resolved_jobs)
    )

    resolved_batch_dir = batch_dir.resolve()
    resolved_paths: dict[Path, str] = {}
    for job in jobs:
        run_path = (batch_dir / job["run_id"]).resolve()
        if run_path.parent != resolved_batch_dir:
            raise ValueError(f"Resolved run ID '{job['run_id']}' escapes the batch directory.")
        if run_path in resolved_paths:
            raise ValueError(
                "Resolved run-directory collision between "
                f"'{resolved_paths[run_path]}' and '{job['run_id']}'."
            )
        resolved_paths[run_path] = job["run_id"]
        job["_resume_existing"] = False
        if run_path.exists():
            if overwrite:
                continue
            if resume:
                job["_resume_existing"] = True
                continue
            raise FileExistsError(
                f"{run_path} already exists. Use --resume to validate/reuse it or --overwrite."
            )
    return jobs


def _write_batch_manifest(
    batch_dir: Path,
    *,
    run_group: str,
    status: str,
    metrics: list[str],
    execution_provenance: Mapping[str, Any],
    manifest_rows: list[Mapping[str, Any]],
) -> None:
    jobs = [
        {
            "run_id": row.get("run_id"),
            "method": row.get("method"),
            "method_run_id": row.get("method_run_id"),
            "status": row.get("status"),
            "execution_action": row.get("execution_action"),
            "run_dir": row.get("run_dir"),
            "output_manifest_sha256": row.get("output_manifest_sha256"),
            "wall_runtime_seconds": row.get("wall_runtime_seconds"),
            "peak_rss_bytes": row.get("peak_rss_bytes"),
        }
        for row in manifest_rows
    ]
    _write_yaml(
        batch_dir / "batch_manifest.yaml",
        {
            "schema_version": 1,
            "run_group": run_group,
            "status": status,
            "metrics": metrics,
            "execution_provenance": copy.deepcopy(dict(execution_provenance)),
            "runs": jobs,
        },
    )


def _validate_resume_batch_manifest(
    batch_dir: Path,
    execution_provenance: Mapping[str, Any],
) -> None:
    path = batch_dir / "batch_manifest.yaml"
    if not path.exists():
        return
    with path.open() as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, Mapping):
        raise ValueError(f"Cannot resume: batch manifest is invalid: {path}")
    if payload.get("execution_provenance") != execution_provenance:
        raise ValueError(
            "Cannot resume: batch manifest provenance does not match the current campaign."
        )


def _failure_fallback_metadata(
    job: Mapping[str, Any],
    registry: Union[str, Path],
    execution_provenance: Mapping[str, Any],
) -> dict[str, Any]:
    environment = collect_environment(str(job["method"]))
    return {
        "run_id": job["run_id"],
        "method": job["method"],
        "method_run_id": job["method_run_id"],
        "dataset_id": job["dataset"],
        "modality": job["modality"],
        "feature_set": job["feature_set"],
        "method_config": copy.deepcopy(job["method_config"]),
        "method_config_source": job["method_config_source"],
        "method_config_source_path": job["method_config_source_path"],
        "method_config_sha256": job["method_config_sha256"],
        "method_config_source_sha256": job["method_config_source_sha256"],
        "dataset_config_sha256": job["dataset_config_sha256"],
        "environment": environment,
        "software_versions": environment["packages"],
        "inputs": {
            "dataset_id": job["dataset"],
            "modality": job["modality"],
            "feature_set": job["feature_set"],
            "registry": str(registry),
        },
        "execution_provenance": copy.deepcopy(dict(execution_provenance)),
        **_execution_metadata(job["method_config"]),
    }


def _attach_compute_metadata(
    performance: dict[str, Any],
    execution_provenance: Mapping[str, Any],
    resource_preflight: Mapping[str, Any],
) -> dict[str, Any]:
    performance["compute_environment"] = copy.deepcopy(
        execution_provenance.get("compute_environment", {})
    )
    performance["resource_preflight"] = copy.deepcopy(dict(resource_preflight))
    return performance


def run_experiment(
    experiment_config_path: Union[str, Path],
    registry: Union[str, Path],
    output_root_override: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
    resume: bool = False,
    shard_index: Optional[int] = None,
    shard_count: Optional[int] = None,
) -> Path:
    experiment_config_path = resolve_config_path(experiment_config_path).resolve()
    experiment_config = read_yaml(experiment_config_path)
    effective_overwrite = overwrite or experiment_config.get("overwrite", False)
    if resume and effective_overwrite:
        raise ValueError("--resume and overwrite are mutually exclusive.")

    evaluation_mode = str(experiment_config.get("evaluation_mode", "exact_truth"))
    if evaluation_mode not in EVALUATION_MODES:
        raise ValueError(
            f"evaluation_mode must be one of {sorted(EVALUATION_MODES)!r}."
        )
    if evaluation_mode == "prediction_only":
        metrics = _as_list(experiment_config.get("metrics"))
        if metrics:
            raise ValueError("prediction_only evaluation requires metrics: [].")
    else:
        metrics = (
            _as_list(experiment_config["metrics"])
            if "metrics" in experiment_config
            else ["rmse_v1", "jsd_v2"]
        )
        if not metrics:
            raise ValueError("exact_truth evaluation requires at least one metric.")
    resolved_metric_ids: list[str] = []
    for metric_name in metrics:
        try:
            resolved_metric_ids.append(resolve_proportion_metric(metric_name).metric_id)
        except KeyError as exc:
            available = ", ".join(available_proportion_metrics())
            raise KeyError(f"Unknown metric '{metric_name}'. Available: {available}") from exc
    duplicate_metrics = sorted(
        {metric_id for metric_id in resolved_metric_ids if resolved_metric_ids.count(metric_id) > 1}
    )
    if duplicate_metrics:
        raise ValueError(
            "Metrics resolve to duplicate canonical endpoint(s): " + ", ".join(duplicate_metrics)
        )

    run_id_template = experiment_config.get(
        "run_id_template",
        "{dataset}__{modality}__{feature_set}__{method_run_id}",
    )
    all_jobs = resolve_experiment_jobs(experiment_config, registry, run_id_template)
    selected_jobs, shard_metadata = select_experiment_shard(
        all_jobs,
        shard_index=shard_index,
        shard_count=shard_count,
    )

    run_group = experiment_config.get("run_group")
    if run_group is None:
        if resume:
            raise ValueError("--resume requires experiment config field 'run_group'.")
        timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        run_group = f"{timestamp}_{experiment_config_path.stem}"
    run_group = str(run_group)
    if shard_metadata is not None:
        run_group += str(shard_metadata["run_group_suffix"])

    output_root = Path(output_root_override or experiment_config.get("output_root", ROOT / "results"))
    if not output_root.is_absolute():
        output_root = ROOT / output_root
    batch_dir = output_root / run_group

    if shard_metadata is not None and experiment_config.get("comparison_output") is not None:
        raise ValueError(
            "Sharded experiments require the default per-shard comparison output to avoid races."
        )
    continue_on_error = experiment_config.get("continue_on_error", False)
    jobs = preflight_experiment_jobs(
        experiment_config,
        registry=registry,
        run_id_template=run_id_template,
        batch_dir=batch_dir,
        overwrite=effective_overwrite,
        resume=resume,
        resolved_jobs=selected_jobs,
    )
    execution_provenance = build_execution_provenance(
        experiment_config_path,
        experiment_config,
        all_jobs,
        resolved_metric_ids,
        registry,
        shard=shard_metadata,
    )

    resumed_metadata: dict[str, dict[str, Any]] = {}
    if resume and batch_dir.exists():
        _validate_resume_batch_manifest(batch_dir, execution_provenance)
        for job in jobs:
            if not job["_resume_existing"]:
                continue
            run_dir = batch_dir / job["run_id"]
            try:
                metadata = _validate_resume_run(run_dir, job, execution_provenance)
                # Scientific evaluation is part of completion; artifact hashes
                # alone do not make a contract-invalid output reusable.
                evaluate_run(run_dir, metrics, evaluation_mode=evaluation_mode)
            except Exception as exc:
                raise ValueError(
                    f"Cannot resume run '{job['run_id']}': {exc}. "
                    "Failed or partial directories are preserved; inspect them and rerun "
                    "with an explicit overwrite only if replacement is intended."
                ) from exc
            resumed_metadata[job["run_id"]] = metadata

    # All schema, config, registry, and run-ID validation happens above this
    # first mutation of the output tree.
    batch_dir.mkdir(parents=True, exist_ok=True)

    successful_runs: list[Path] = []
    manifest_rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    _write_batch_manifest(
        batch_dir,
        run_group=run_group,
        status="running",
        metrics=resolved_metric_ids,
        execution_provenance=execution_provenance,
        manifest_rows=manifest_rows,
    )

    for job in jobs:
        run_id = job["run_id"]
        run_dir = batch_dir / run_id
        if job["_resume_existing"]:
            metadata = resumed_metadata[run_id]
            successful_runs.append(run_dir)
            manifest_rows.append(
                _manifest_row(
                    job,
                    batch_dir,
                    "success",
                    metadata=metadata,
                    execution_action="resumed",
                )
            )
            pd.DataFrame(manifest_rows).to_csv(batch_dir / "runs.csv", index=False)
            _write_batch_manifest(
                batch_dir,
                run_group=run_group,
                status="running",
                metrics=resolved_metric_ids,
                execution_provenance=execution_provenance,
                manifest_rows=manifest_rows,
            )
            continue

        resource_preflight = _wait_for_resource_gate()
        monitor = _PeakRSSMonitor()
        started_at = monitor.start()
        try:
            run_dir = run_one(
                dataset=job["dataset"],
                modality=job["modality"],
                feature_set=job["feature_set"],
                method=job["method"],
                method_config=job["method_config"],
                registry=registry,
                output_root=batch_dir,
                run_id=run_id,
                overwrite=effective_overwrite,
                method_run_id=job["method_run_id"],
                method_config_source_path=job["_method_config_source_path_resolved"],
                method_config_source=job["method_config_source"],
                resolved_method_config_sha256=job["method_config_sha256"],
                method_config_source_sha256=job["method_config_source_sha256"],
                execution_provenance=execution_provenance,
                resolved_dataset_config_sha256=job["dataset_config_sha256"],
            )
            # Evaluation is part of run validity. Contract violations and
            # non-finite endpoints are recorded as failures rather than
            # silently dropping spots, types, or metric rows.
            evaluate_run(run_dir, metrics, evaluation_mode=evaluation_mode)
            performance = _attach_compute_metadata(
                monitor.stop(
                    started_at,
                    "load_signature_fit_diagnostics_standard_write_and_evaluation",
                ),
                execution_provenance,
                resource_preflight,
            )
            metadata = _finalize_run_directory(
                run_dir,
                status="success",
                performance=performance,
                execution_provenance=execution_provenance,
                execution_action="executed",
            )
            successful_runs.append(run_dir)
            manifest_rows.append(
                _manifest_row(
                    job,
                    batch_dir,
                    "success",
                    metadata=metadata,
                    execution_action="executed",
                )
            )
        except Exception as exc:
            failure_diagnostics = getattr(exc, "diagnostics", None)
            if failure_diagnostics is not None and callable(
                getattr(failure_diagnostics, "to_dict", None)
            ):
                with (run_dir / "failure_diagnostics.json").open("w") as handle:
                    json.dump(
                        failure_diagnostics.to_dict(),
                        handle,
                        indent=2,
                        sort_keys=True,
                    )
                    handle.write("\n")
            performance = _attach_compute_metadata(
                monitor.stop(
                    started_at,
                    "load_signature_fit_diagnostics_standard_write_and_evaluation_until_failure",
                ),
                execution_provenance,
                resource_preflight,
            )
            metadata = _finalize_run_directory(
                run_dir,
                status="failed",
                performance=performance,
                execution_provenance=execution_provenance,
                execution_action="executed",
                fallback_metadata=_failure_fallback_metadata(
                    job,
                    registry,
                    execution_provenance,
                ),
                error=str(exc),
            )
            manifest_row = _manifest_row(
                job,
                batch_dir,
                "failed",
                metadata=metadata,
                execution_action="executed",
            )
            failure = {
                **manifest_row,
                "metric": None,
                "value": None,
                "error": str(exc),
            }
            failures.append(failure)
            manifest_rows.append(manifest_row)
            if not continue_on_error:
                pd.DataFrame(manifest_rows).to_csv(batch_dir / "runs.csv", index=False)
                pd.DataFrame(failures).to_csv(batch_dir / "failures.csv", index=False)
                _write_batch_manifest(
                    batch_dir,
                    run_group=run_group,
                    status="failed",
                    metrics=resolved_metric_ids,
                    execution_provenance=execution_provenance,
                    manifest_rows=manifest_rows,
                )
                raise

        pd.DataFrame(manifest_rows).to_csv(batch_dir / "runs.csv", index=False)
        if failures:
            pd.DataFrame(failures).to_csv(batch_dir / "failures.csv", index=False)
        _write_batch_manifest(
            batch_dir,
            run_group=run_group,
            status="running",
            metrics=resolved_metric_ids,
            execution_provenance=execution_provenance,
            manifest_rows=manifest_rows,
        )

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

    if evaluation_mode == "prediction_only":
        prediction_rows = [
            {
                **dict(row),
                "evaluation_mode": evaluation_mode,
                "metric": None,
                "value": None,
                "error": None,
            }
            for row in manifest_rows
            if row.get("status") == "success"
        ]
        prediction_rows.extend(failures)
        comparison_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(prediction_rows).to_csv(comparison_path, index=False)
    else:
        write_comparison(
            comparison_path,
            successful_runs=successful_runs,
            failures=failures,
            metrics=metrics,
        )
    _write_batch_manifest(
        batch_dir,
        run_group=run_group,
        status="completed_with_failures" if failures else "completed",
        metrics=resolved_metric_ids,
        execution_provenance=execution_provenance,
        manifest_rows=manifest_rows,
    )
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
    parser.add_argument("--method-run-id", help="Named method variant; defaults to --method.")
    parser.add_argument("--config", help="Method config for single-run mode.")
    parser.add_argument("--registry", default=str(ROOT / "data" / "registry" / "datasets.yaml"))
    parser.add_argument("--output-root")
    parser.add_argument("--run-id")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Validate and reuse completed experiment runs; reject failed/partial runs.",
    )
    parser.add_argument(
        "--shard-index",
        type=int,
        help="Zero-based deterministic experiment shard index (requires --shard-count).",
    )
    parser.add_argument(
        "--shard-count",
        type=int,
        help="Number of deterministic experiment shards (requires --shard-index).",
    )
    args = parser.parse_args()

    if args.experiment_config is not None:
        run_experiment(
            experiment_config_path=args.experiment_config,
            registry=args.registry,
            output_root_override=args.output_root,
            overwrite=args.overwrite,
            resume=args.resume,
            shard_index=args.shard_index,
            shard_count=args.shard_count,
        )
        return

    if args.resume or args.shard_index is not None or args.shard_count is not None:
        parser.error("--resume and --shard-index/--shard-count are experiment-mode options.")

    missing = [name for name in ("dataset", "modality", "method") if getattr(args, name) is None]
    if missing:
        parser.error(
            "single-run mode requires --dataset, --modality, and --method; "
            "batch mode requires --experiment-config."
        )

    config_path = resolve_config_path(args.config).resolve() if args.config is not None else None
    method_config = _read_config_mapping(config_path, "Single-run method config") if config_path else {}
    method = args.method.lower()
    method_run_id = args.method_run_id or method
    run_id = args.run_id or build_run_id(
        args.dataset,
        args.modality,
        args.feature_set,
        method,
        method_run_id=method_run_id,
    )
    output_root = Path(args.output_root or ROOT / "results")
    run_dir = output_root / run_id
    output_existed_before = run_dir.exists()
    dataset_config = get_dataset_config(
        args.dataset,
        registry_path=args.registry,
        project_root=ROOT,
    )
    dataset_config_sha256 = _canonical_sha256(dataset_config)
    config_sha256 = method_config_sha256(method_config)
    config_source_sha256 = _source_sha256(config_path)
    invocation = {
        "dataset": args.dataset,
        "modality": args.modality,
        "feature_set": args.feature_set,
        "method": method,
        "method_run_id": method_run_id,
        "run_id": run_id,
        "method_config": method_config,
        "method_config_source_path": _recorded_path(config_path),
        "method_config_source_sha256": config_source_sha256,
        "method_config_sha256": config_sha256,
        "dataset_config_sha256": dataset_config_sha256,
    }
    execution_provenance = build_single_run_provenance(invocation, args.registry, method)
    job = {
        **invocation,
        "method_config_source": (
            "single_run_config" if config_path else "implicit_empty_default"
        ),
    }
    resource_preflight = _wait_for_resource_gate()
    monitor = _PeakRSSMonitor()
    started_at = monitor.start()
    try:
        completed_run_dir = run_one(
            dataset=args.dataset,
            modality=args.modality,
            feature_set=args.feature_set,
            method=method,
            method_config=method_config,
            registry=args.registry,
            output_root=output_root,
            run_id=run_id,
            overwrite=args.overwrite,
            method_run_id=method_run_id,
            method_config_source_path=config_path,
            method_config_source=job["method_config_source"],
            resolved_method_config_sha256=config_sha256,
            method_config_source_sha256=config_source_sha256,
            execution_provenance=execution_provenance,
            resolved_dataset_config_sha256=dataset_config_sha256,
        )
        performance = _attach_compute_metadata(
            monitor.stop(started_at, "load_signature_fit_diagnostics_and_standard_write"),
            execution_provenance,
            resource_preflight,
        )
        _finalize_run_directory(
            completed_run_dir,
            status="success",
            performance=performance,
            execution_provenance=execution_provenance,
            execution_action="executed",
        )
    except Exception as exc:
        failure_diagnostics = getattr(exc, "diagnostics", None)
        if (
            run_dir.exists()
            and failure_diagnostics is not None
            and callable(getattr(failure_diagnostics, "to_dict", None))
        ):
            with (run_dir / "failure_diagnostics.json").open("w") as handle:
                json.dump(
                    failure_diagnostics.to_dict(),
                    handle,
                    indent=2,
                    sort_keys=True,
                )
                handle.write("\n")
        performance = _attach_compute_metadata(
            monitor.stop(started_at, "single_run_until_failure"),
            execution_provenance,
            resource_preflight,
        )
        if run_dir.exists() and (args.overwrite or not output_existed_before):
            _finalize_run_directory(
                run_dir,
                status="failed",
                performance=performance,
                execution_provenance=execution_provenance,
                execution_action="executed",
                fallback_metadata=_failure_fallback_metadata(
                    job,
                    args.registry,
                    execution_provenance,
                ),
                error=str(exc),
            )
        raise


if __name__ == "__main__":
    main()
