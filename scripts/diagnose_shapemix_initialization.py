#!/usr/bin/env python
"""Record ShapeMix NNLS and deterministic-restart scales without MAP fitting."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Any

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from deconvatac.data import load_deconvolution_input
from deconvatac.shapemix.config import ShapeMixConfig
from deconvatac.shapemix.map import (
    ABUNDANCE_EPSILON,
    SEED_NAMESPACE,
    _CountSource,
    _nnls_initialization,
    _seeded_restart,
)
from deconvatac.shapemix.signatures import estimate_reference_signatures


DEFAULT_DATASET = (
    "gse129785_shapemix_physical_dilution_0p1_99p9_cd4mem_cd8naive"
)
DEFAULT_OUTPUT = (
    ROOT
    / "results/development/shapemix_gse129785_convergence_v1"
    / "initialization_diagnostics/scale.yaml"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite_summary(values: np.ndarray) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    if not np.isfinite(array).all() or (array <= 0).any():
        raise ValueError("Initialization diagnostic received invalid abundance values.")
    return {
        "minimum": float(array.min()),
        "median": float(np.median(array)),
        "maximum": float(array.max()),
        "sum": float(array.sum()),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default=DEFAULT_DATASET)
    parser.add_argument("--feature-set", default="selected_reference_peaks")
    parser.add_argument(
        "--method-config",
        type=Path,
        default=ROOT / "configs/methods/shapemix.yaml",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    method_config_path = args.method_config.resolve()
    with method_config_path.open() as handle:
        method_config: dict[str, Any] = yaml.safe_load(handle)
    config = ShapeMixConfig.from_mapping(method_config["params"])
    data = load_deconvolution_input(
        dataset_id=args.dataset,
        modality="atac",
        feature_set=args.feature_set,
        registry_path=ROOT / "data/registry/datasets.yaml",
        project_root=ROOT,
    )
    outer_seed = int(data.metadata["dataset_config"]["shapemix_seeds"]["outer_split_seed"])
    inner_seed = int(data.metadata["dataset_config"]["shapemix_seeds"]["inner_mixture_seed"])
    layer_names = data.fragment_shape.layer_names
    bin_names = tuple(item.name for item in data.fragment_shape.bins)
    signatures = estimate_reference_signatures(
        data.reference,
        data.labels_key,
        data.cell_types,
        outer_seed,
        config=config,
        layer_names=layer_names,
        bin_names=bin_names,
    )
    counts = _CountSource(
        tuple(data.spatial.layers[layer_name] for layer_name in layer_names)
    )
    initial, fallback_spots = _nnls_initialization(
        counts,
        signatures.A,
        tuple(data.spatial.obs_names.astype(str)),
    )
    restart_records = []
    for restart_index in range(config.restarts):
        seed_tuple = (
            SEED_NAMESPACE,
            outer_seed,
            inner_seed,
            config.seed,
            restart_index,
        )
        raw = _seeded_restart(initial, seed_tuple)
        abundance = np.exp(raw) + ABUNDANCE_EPSILON
        relative = abundance / initial
        restart_records.append(
            {
                "restart_index": restart_index,
                "seed_tuple": list(seed_tuple),
                "abundance": finite_summary(abundance),
                "relative_to_nnls": finite_summary(relative),
                "maximum_absolute_displacement": float(
                    np.max(np.abs(abundance - initial))
                ),
            }
        )
    observed_total = 0.0
    for spot_index in range(counts.shape[0]):
        observed_total += float(counts.collapsed_row(spot_index).sum())
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "diagnostic": "shapemix_initialization_scale_without_map_optimization",
        "dataset_id": args.dataset,
        "feature_set": args.feature_set,
        "dimensions": {
            "spots": counts.shape[0],
            "peaks": counts.shape[1],
            "bins": counts.shape[2],
            "cell_types": signatures.A.shape[0],
        },
        "observed_total_cut_sites": observed_total,
        "nnls_abundance": finite_summary(initial),
        "nnls_fallback_spots": list(fallback_spots),
        "deterministic_restarts": restart_records,
        "method_config": str(method_config_path.relative_to(ROOT)),
        "method_config_sha256": sha256(method_config_path),
        "map_optimization_performed": False,
    }
    with output.open("w") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)
    print(output)


if __name__ == "__main__":
    main()
