"""Shared provenance helpers for generated ShapeMix data manifests."""

from __future__ import annotations

import hashlib
import importlib.metadata
import platform
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from deconvatac.data import FragmentShapeSpec


ROOT = Path(__file__).resolve().parents[1]
CODE_AND_PROTOCOL_PATHS = (
    "configs/data_sources/pbmc_granulocyte_sorted_10k_cellranger_arc_2.0.0.yaml",
    "docs/ShapeMix/model_specification.md",
    "docs/ShapeMix/benchmark_protocol.md",
    "scripts/shapemix_provenance.py",
    "scripts/audit_shapemix_signal.py",
    "scripts/prepare_shapemix_pbmc.py",
    "scripts/regenerate_shapemix_pbmc_simulations.py",
    "src/deconvatac/data/schemas.py",
    "src/deconvatac/data/validators.py",
    "src/deconvatac/data/loaders.py",
    "src/deconvatac/pp/fragment_shapes.py",
    "src/deconvatac/pp/feature_selection.py",
)


def _plain_data(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _plain_data(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return [_plain_data(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_plain_data(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_output(arguments: Sequence[str]) -> str | None:
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


def code_provenance() -> dict[str, Any]:
    """Record the base commit and content hashes of executable protocol inputs."""
    status = _git_output(("status", "--short", "--untracked-files=normal"))
    files = {
        relative: sha256_file(ROOT / relative)
        for relative in CODE_AND_PROTOCOL_PATHS
    }
    return {
        "git_commit": _git_output(("rev-parse", "HEAD")),
        "git_worktree_dirty": None if status is None else bool(status),
        "content_sha256": files,
    }


def software_versions() -> dict[str, str]:
    distributions = {
        "deconvatac": "deconvATAC",
        "anndata": "anndata",
        "numpy": "numpy",
        "pandas": "pandas",
        "scipy": "scipy",
        "pysam": "pysam",
        "pyyaml": "PyYAML",
    }
    versions = {"python": platform.python_version()}
    for key, distribution in distributions.items():
        try:
            versions[key] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[key] = "not-installed"
    return versions


def fragment_shape_declaration(metadata: Mapping[str, Any]) -> dict[str, Any]:
    spec = FragmentShapeSpec.from_mapping(metadata)
    return _plain_data({
        "schema_version": spec.schema_version,
        "axis": spec.axis,
        "count_unit": spec.count_unit,
        "read_support_policy": spec.read_support_policy,
        "peak_assignment": spec.peak_assignment,
        "left_cut_offset": spec.left_cut_offset,
        "right_cut_offset": spec.right_cut_offset,
        "bins": [bin_spec.to_dict(omit_none=True) for bin_spec in spec.bins],
    })


def matrix_summary(adata: Any, layer_names: Sequence[str]) -> dict[str, Any]:
    def nonzero_count(matrix: Any) -> int:
        return int(matrix.nnz) if hasattr(matrix, "nnz") else int(np.count_nonzero(matrix))

    return {
        "observations": int(adata.n_obs),
        "features": int(adata.n_vars),
        "x_nonzero": nonzero_count(adata.X),
        "x_total": int(adata.X.sum()),
        "layers": {
            str(layer_name): {
                "nonzero": nonzero_count(adata.layers[layer_name]),
                "total": int(adata.layers[layer_name].sum()),
            }
            for layer_name in layer_names
        },
    }
