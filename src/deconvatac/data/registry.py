from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REGISTRY_PATH = PROJECT_ROOT / "data" / "registry" / "datasets.yaml"


def resolve_project_path(path: Union[str, Path], project_root: Optional[Union[str, Path]] = None) -> Path:
    """Resolve a path relative to the project root unless it is already absolute."""
    path = Path(path)
    if path.is_absolute():
        return path
    root = Path(project_root) if project_root is not None else PROJECT_ROOT
    return root / path


def read_yaml(path: Union[str, Path]) -> dict[str, Any]:
    """Read a YAML file and return an empty dict for empty files."""
    with Path(path).open() as handle:
        return yaml.safe_load(handle) or {}


def load_dataset_registry(registry_path: Optional[Union[str, Path]] = None) -> dict[str, Any]:
    """Load the dataset registry."""
    path = Path(registry_path) if registry_path is not None else DEFAULT_REGISTRY_PATH
    return read_yaml(path)


def get_dataset_config(
    dataset_id: str,
    registry_path: Optional[Union[str, Path]] = None,
    project_root: Optional[Union[str, Path]] = None,
) -> dict[str, Any]:
    """Resolve and load a dataset config from the registry."""
    registry = load_dataset_registry(registry_path)
    if dataset_id not in registry:
        available = ", ".join(sorted(registry))
        raise KeyError(f"Unknown dataset_id '{dataset_id}'. Available datasets: {available}")

    entry = registry[dataset_id]
    config_path = entry if isinstance(entry, str) else entry.get("config") or entry.get("path")
    if config_path is None:
        raise ValueError(f"Registry entry for '{dataset_id}' must define 'config' or 'path'.")

    resolved = resolve_project_path(config_path, project_root=project_root)
    config = read_yaml(resolved)
    config.setdefault("dataset_id", dataset_id)
    config.setdefault("config_path", str(resolved))
    return config
