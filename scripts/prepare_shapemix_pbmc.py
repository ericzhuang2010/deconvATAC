#!/usr/bin/env python
"""Prepare cell-disjoint PBMC inputs for the ShapeMix benchmark.

The expensive fragments-file pass is represented by a reusable, full-universe
shape cache.  Every outer split then performs reference-only peak selection and
slices that cache by rows and ranked features.  This keeps the raw source-run
QC immutable while giving each derived matrix its own exact counters and hash.
"""

from __future__ import annotations

import argparse
import copy
import gc
import hashlib
import importlib.metadata
import os
import shutil
import sys
import tempfile
from collections.abc import Callable, Iterable, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import scipy
import yaml
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    # Repository scripts must prefer the working tree over a stale non-editable
    # installation in the active environment.
    sys.path.insert(0, str(SRC))

from deconvatac.data import ordered_feature_sha256  # noqa: E402
from deconvatac.pp import (  # noqa: E402
    FRAGMENT_SHAPE_LAYER_NAMES,
    FragmentRecord,
    PeakInterval,
    build_fragment_shape_anndata,
    count_fragment_shapes,
    count_fragment_shapes_from_records,
    read_peaks_bed,
)
from scripts.audit_shapemix_signal import (  # noqa: E402
    audit_shape_signal,
    write_signal_audit,
)
from scripts.shapemix_provenance import (  # noqa: E402
    code_provenance,
    fragment_shape_declaration,
    matrix_summary,
    software_versions,
)

PROTOCOL_VERSION = 1
RNG_NAMESPACE = 20260822
DEFAULT_OUTER_SPLIT_SEED = 0
PRIMARY_OUTER_SPLIT_SEEDS = (1_103, 2_203, 3_301, 4_409, 5_501)
CACHE_OUTER_SPLIT_SEEDS = (DEFAULT_OUTER_SPLIT_SEED, *PRIMARY_OUTER_SPLIT_SEEDS)
DEFAULT_N_PEAKS = 5_000
DEFAULT_MIN_REFERENCE_CELLS = 10
REFERENCE_FRACTION = 0.70

SOURCE_ROOT = (
    ROOT
    / "data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0"
)
DEFAULT_INPUT_PATHS = {
    "labels": ROOT / "data/raw/sources/snapatac2/pbmc10k_multiome/cell_type_mapping.csv",
    "label_summary": ROOT / "data/raw/sources/snapatac2/pbmc10k_multiome/cell_type_summary.csv",
    "reference": ROOT / "data/raw/references/pbmc_granulocyte_sorted_10k_multiome/atac/reference.h5ad",
    "peaks": SOURCE_ROOT / "pbmc_granulocyte_sorted_10k_atac_peaks.bed",
    "filtered_feature_matrix": SOURCE_ROOT / "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
    "fragments": SOURCE_ROOT / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
    "tabix_index": SOURCE_ROOT / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
}
PINNED_INPUT_SHA256 = {
    "labels": "3aa94d5f636c01c0159324984930cb14e775e49611102cb4dd02a491b63bf298",
    "label_summary": "809a3b285d996024daed95e0c2a8a17299ba93d1e5f0e884988dac9362150f4d",
    "reference": "cdaefffbfd5dd3cb36318f68158d4bec0df1d2b2bf63562b14065f1055cf6ee6",
    "peaks": "3975a4057f9caa3fb69ddaecc6ae9e530e77551717a1464c2d93ac9d73cb60ab",
    "filtered_feature_matrix": "f6824171378787baab244f559b8b438f79db2eb39f78d17b2196f7ecd2c03549",
    "fragments": "5075e32a0e9c6dded35b060bf90d6144375b150e131ffb0be121a93e3b5e1e38",
    "tabix_index": "3a516291d0e6e5ddf9f651f470b6312b83eeb26f153ea12cfa9d082760a5e7f5",
}

FROZEN_CELL_TYPE_COUNTS = {
    "CD14 Mono": 2_551,
    "CD4 Naive": 1_382,
    "CD8 Naive": 1_353,
    "CD4 TCM": 1_113,
    "CD16 Mono": 442,
    "NK": 403,
    "CD8 TEM_1": 322,
    "CD8 TEM_2": 315,
    "Intermediate B": 300,
    "Memory B": 298,
    "CD4 TEM": 286,
    "cDC": 180,
    "Treg": 157,
    "gdT": 143,
    "MAIT": 130,
    "Naive B": 125,
}
FROZEN_CELL_TYPES = tuple(FROZEN_CELL_TYPE_COUNTS)
SMOKE_CELL_TYPES = FROZEN_CELL_TYPES[:3]
SMOKE_N_PEAKS = 200
EXPECTED_REFERENCE_SHAPE = (9_627, 143_887)
EXPECTED_EXECUTABLE_CELLS = 9_500

DEFAULT_OUTPUT_ROOT = ROOT / "data/processed/shapemix/pbmc_granulocyte_sorted_10k"
DEFAULT_SHAPE_CACHE = DEFAULT_OUTPUT_ROOT / "full_universe_shape_cache.h5ad"


@dataclass
class PreparedSplit:
    """In-memory products and audit data for one outer split."""

    membership: pd.DataFrame
    split_sha256: str
    peak_selection: Any
    reference_cells: ad.AnnData
    heldout_test_cells: ad.AnnData


def sha256_file(path: Union[str, Path], chunk_bytes: int = 8 * 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest without loading the file into memory."""
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_bytes), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_pinned_files(
    paths: Mapping[str, Union[str, Path]],
    expected_sha256: Mapping[str, str],
) -> dict[str, str]:
    """Fail closed unless every named input exists and matches its pinned hash."""
    missing_specs = sorted(set(expected_sha256).difference(paths))
    unexpected_specs = sorted(set(paths).difference(expected_sha256))
    if missing_specs or unexpected_specs:
        raise ValueError(
            "Pinned input names differ: "
            f"missing={missing_specs}, unexpected={unexpected_specs}."
        )

    observed: dict[str, str] = {}
    for name in expected_sha256:
        path = Path(paths[name])
        if not path.is_file():
            raise FileNotFoundError(f"Pinned input '{name}' is missing: {path}")
        digest = sha256_file(path)
        expected = expected_sha256[name].lower()
        if digest != expected:
            raise ValueError(
                f"Pinned input '{name}' has SHA-256 {digest}, expected {expected}: {path}"
            )
        observed[name] = digest
    return observed


def _canonicalize_universe(
    universe: pd.DataFrame,
    cell_types: Sequence[str],
) -> pd.DataFrame:
    required = {"barcode", "cell_type"}
    if not required.issubset(universe.columns):
        raise ValueError("Cell universe must contain barcode and cell_type columns.")
    canonical = universe.loc[:, ["barcode", "cell_type"]].copy()
    if canonical.isna().any().any():
        raise ValueError("Cell universe cannot contain missing barcodes or cell types.")
    canonical["barcode"] = canonical["barcode"].astype(str)
    canonical["cell_type"] = canonical["cell_type"].astype(str)
    if canonical["barcode"].str.len().eq(0).any():
        raise ValueError("Cell universe cannot contain empty barcodes.")
    if canonical["barcode"].duplicated().any():
        raise ValueError("Cell universe barcodes must be unique.")

    ordered_types = tuple(str(value) for value in cell_types)
    if not ordered_types or len(set(ordered_types)) != len(ordered_types):
        raise ValueError("cell_types must be a non-empty ordered unique sequence.")
    observed_types = set(canonical["cell_type"])
    if observed_types != set(ordered_types):
        raise ValueError("Cell universe must contain exactly the declared cell types.")

    type_order = {cell_type: index for index, cell_type in enumerate(ordered_types)}
    canonical["_type_order"] = canonical["cell_type"].map(type_order)
    canonical = canonical.sort_values(
        ["_type_order", "barcode"], kind="mergesort"
    ).drop(columns="_type_order")
    return canonical.reset_index(drop=True)


def build_cell_universe(
    reference_barcodes: Sequence[str],
    label_mapping: pd.DataFrame,
    *,
    cell_types: Sequence[str],
    expected_counts: Optional[Mapping[str, int]] = None,
    reference_labels: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Intersect reference cells with labels and return canonical eligible cells."""
    reference_index = pd.Index([str(value) for value in reference_barcodes], name="barcode")
    if reference_index.has_duplicates:
        raise ValueError("Reference barcodes must be unique.")
    if "barcode" not in label_mapping or "cell_type" not in label_mapping:
        raise ValueError("Label mapping must contain barcode and cell_type columns.")

    labels = label_mapping.loc[:, ["barcode", "cell_type"]].copy()
    if labels[["barcode", "cell_type"]].isna().any().any():
        raise ValueError("Label mapping contains missing barcode or cell_type values.")
    labels["barcode"] = labels["barcode"].astype(str)
    labels["cell_type"] = labels["cell_type"].astype(str)
    if labels["barcode"].duplicated().any():
        raise ValueError("Label mapping barcodes must be unique.")
    labels = labels[labels["barcode"].isin(reference_index)]

    if reference_labels is not None:
        if len(reference_labels) != len(reference_index):
            raise ValueError("reference_labels must align one-to-one with reference_barcodes.")
        expected_by_barcode = pd.Series(
            [str(value) for value in reference_labels], index=reference_index
        )
        observed_by_barcode = labels.set_index("barcode")["cell_type"]
        compared = expected_by_barcode.loc[observed_by_barcode.index]
        mismatched = compared.to_numpy() != observed_by_barcode.to_numpy()
        if np.any(mismatched):
            barcode = str(observed_by_barcode.index[np.flatnonzero(mismatched)[0]])
            raise ValueError(f"Pinned mapping disagrees with reference label for {barcode}.")

    declared = tuple(str(value) for value in cell_types)
    universe = labels[labels["cell_type"].isin(declared)]
    universe = _canonicalize_universe(universe, declared)

    if expected_counts is not None:
        if tuple(expected_counts) != declared:
            raise ValueError("expected_counts keys must exactly follow cell_types order.")
        observed_counts = universe["cell_type"].value_counts()
        for cell_type, expected in expected_counts.items():
            observed = int(observed_counts.get(cell_type, 0))
            if observed != int(expected):
                raise ValueError(
                    f"Post-intersection count for '{cell_type}' is {observed}, expected {expected}."
                )
    return universe


def build_frozen_pbmc_universe(
    reference: ad.AnnData,
    label_mapping: pd.DataFrame,
) -> pd.DataFrame:
    """Validate the pinned reference and construct the exact 16-type universe."""
    if reference.shape != EXPECTED_REFERENCE_SHAPE:
        raise ValueError(
            f"Prepared ATAC reference has shape {reference.shape}, expected {EXPECTED_REFERENCE_SHAPE}."
        )
    if "cell_type" not in reference.obs:
        raise KeyError("Prepared ATAC reference is missing obs['cell_type'].")
    universe = build_cell_universe(
        reference.obs_names,
        label_mapping,
        cell_types=FROZEN_CELL_TYPES,
        expected_counts=FROZEN_CELL_TYPE_COUNTS,
        reference_labels=reference.obs["cell_type"].astype(str).tolist(),
    )
    if len(universe) != EXPECTED_EXECUTABLE_CELLS:
        raise ValueError(
            f"Executable PBMC universe has {len(universe)} cells, expected {EXPECTED_EXECUTABLE_CELLS}."
        )
    return universe


def subset_cell_universe(
    universe: pd.DataFrame,
    cell_types: Sequence[str],
) -> pd.DataFrame:
    """Create a smoke-specific universe while preserving canonical type order."""
    requested = tuple(str(value) for value in cell_types)
    unknown = sorted(set(requested).difference(FROZEN_CELL_TYPES))
    if unknown:
        raise ValueError(f"Unknown frozen PBMC cell types: {unknown}")
    canonical_requested = tuple(value for value in FROZEN_CELL_TYPES if value in requested)
    if requested != canonical_requested:
        raise ValueError("Smoke cell types must follow the frozen canonical order.")
    return _canonicalize_universe(
        universe[universe["cell_type"].isin(requested)], requested
    )


def canonical_universe_csv(universe: pd.DataFrame, cell_types: Sequence[str]) -> bytes:
    """Serialize the ordered cell universe using the documented CSV hash form."""
    canonical = _canonicalize_universe(universe, cell_types)
    return canonical.to_csv(index=False, lineterminator="\n").encode("utf-8")


def universe_sha256(universe: pd.DataFrame, cell_types: Sequence[str]) -> str:
    """Hash the canonical barcode-to-cell-type universe."""
    return hashlib.sha256(canonical_universe_csv(universe, cell_types)).hexdigest()


def canonicalize_split_membership(
    membership: pd.DataFrame,
    cell_types: Sequence[str],
) -> pd.DataFrame:
    """Validate basic split schema and restore protocol output order."""
    required = {"barcode", "cell_type", "pool"}
    if not required.issubset(membership.columns):
        raise ValueError("Split membership must contain barcode, cell_type, and pool columns.")
    canonical = membership.loc[:, ["barcode", "cell_type", "pool"]].copy()
    if canonical.isna().any().any():
        raise ValueError("Split membership cannot contain missing values.")
    for column in canonical.columns:
        canonical[column] = canonical[column].astype(str)
    if canonical["barcode"].duplicated().any():
        raise ValueError("A source barcode cannot occur in both split pools.")
    if set(canonical["pool"]) != {"reference", "heldout"}:
        raise ValueError("Split membership must contain exactly reference and heldout pools.")

    ordered_types = tuple(str(value) for value in cell_types)
    if set(canonical["cell_type"]) != set(ordered_types):
        raise ValueError("Split membership does not match the declared cell-type universe.")
    type_order = {cell_type: index for index, cell_type in enumerate(ordered_types)}
    canonical["_type_order"] = canonical["cell_type"].map(type_order)
    canonical = canonical.sort_values(
        ["_type_order", "barcode"], kind="mergesort"
    ).drop(columns="_type_order")
    return canonical.reset_index(drop=True)


def validate_split_membership(
    membership: pd.DataFrame,
    universe: pd.DataFrame,
    cell_types: Sequence[str],
    reference_fraction: float = REFERENCE_FRACTION,
) -> pd.DataFrame:
    """Validate disjointness, complete coverage, labels, and per-type 70/30 sizes."""
    if not 0 < reference_fraction < 1:
        raise ValueError("reference_fraction must be strictly between zero and one.")
    canonical = canonicalize_split_membership(membership, cell_types)
    canonical_universe = _canonicalize_universe(universe, cell_types)
    observed_lookup = canonical.set_index("barcode")["cell_type"]
    expected_lookup = canonical_universe.set_index("barcode")["cell_type"]
    if set(observed_lookup.index) != set(expected_lookup.index):
        raise ValueError("Split membership must cover the executable universe exactly once.")
    observed_lookup = observed_lookup.loc[expected_lookup.index]
    if not np.array_equal(observed_lookup.to_numpy(), expected_lookup.to_numpy()):
        raise ValueError("Split membership changes one or more frozen cell labels.")

    for cell_type in cell_types:
        type_rows = canonical[canonical["cell_type"] == cell_type]
        n_cells = len(type_rows)
        expected_reference = int(np.floor(reference_fraction * n_cells))
        observed_reference = int((type_rows["pool"] == "reference").sum())
        if observed_reference != expected_reference:
            raise ValueError(
                f"Reference pool for '{cell_type}' has {observed_reference} cells, "
                f"expected {expected_reference}."
            )
        if not {"reference", "heldout"}.issubset(set(type_rows["pool"])):
            raise ValueError(f"Both split pools must retain cell type '{cell_type}'.")
    return canonical


def stratified_reference_heldout_split(
    universe: pd.DataFrame,
    outer_split_seed: int,
    *,
    cell_types: Sequence[str],
    namespace: int = RNG_NAMESPACE,
    reference_fraction: float = REFERENCE_FRACTION,
) -> pd.DataFrame:
    """Create the frozen per-type PCG64/SeedSequence 70/30 split."""
    if isinstance(outer_split_seed, (bool, np.bool_)) or not isinstance(
        outer_split_seed, (int, np.integer)
    ):
        raise TypeError("outer_split_seed must be an integer.")
    if outer_split_seed < 0:
        raise ValueError("outer_split_seed must be nonnegative.")
    canonical = _canonicalize_universe(universe, cell_types)
    rows: list[dict[str, str]] = []
    for cell_type_order, cell_type in enumerate(cell_types, start=1):
        barcodes = canonical.loc[
            canonical["cell_type"] == cell_type, "barcode"
        ].sort_values(kind="mergesort").to_numpy(dtype=object)
        rng = np.random.Generator(
            np.random.PCG64(
                np.random.SeedSequence(
                    [int(namespace), int(outer_split_seed), cell_type_order]
                )
            )
        )
        permuted = rng.permutation(barcodes)
        n_reference = int(np.floor(reference_fraction * len(barcodes)))
        reference_barcodes = set(str(value) for value in permuted[:n_reference])
        rows.extend(
            {
                "barcode": str(barcode),
                "cell_type": str(cell_type),
                "pool": "reference" if str(barcode) in reference_barcodes else "heldout",
            }
            for barcode in barcodes
        )
    return validate_split_membership(
        pd.DataFrame(rows), canonical, cell_types, reference_fraction
    )


def canonical_split_csv(
    membership: pd.DataFrame,
    cell_types: Sequence[str],
) -> bytes:
    """Return the exact UTF-8 bytes written to split.csv and hashed."""
    canonical = canonicalize_split_membership(membership, cell_types)
    return canonical.to_csv(index=False, lineterminator="\n").encode("utf-8")


def split_membership_sha256(
    membership: pd.DataFrame,
    cell_types: Sequence[str],
) -> str:
    """Hash canonical split.csv bytes so the recorded hash is independently checkable."""
    return hashlib.sha256(canonical_split_csv(membership, cell_types)).hexdigest()


def select_reference_only_peaks(
    reference: ad.AnnData,
    membership: pd.DataFrame,
    *,
    cell_types: Sequence[str],
    cell_type_key: str = "cell_type",
    n_top_peaks: int = DEFAULT_N_PEAKS,
    min_reference_cells: int = DEFAULT_MIN_REFERENCE_CELLS,
    selector: Optional[Callable[..., Any]] = None,
) -> Any:
    """Call the frozen peak selector with reference-pool rows only."""
    canonical = canonicalize_split_membership(membership, cell_types)
    reference_barcodes = canonical.loc[
        canonical["pool"] == "reference", "barcode"
    ].tolist()
    missing = sorted(set(reference_barcodes).difference(reference.obs_names))
    if missing:
        raise ValueError(f"Reference AnnData is missing split barcodes, beginning with {missing[0]!r}.")
    if cell_type_key not in reference.obs:
        raise KeyError(f"Reference AnnData is missing obs[{cell_type_key!r}].")
    training = reference[reference_barcodes, :].copy()
    observed_labels = training.obs[cell_type_key].astype(str)
    expected_labels = canonical.set_index("barcode").loc[reference_barcodes, "cell_type"]
    if not np.array_equal(observed_labels.to_numpy(), expected_labels.to_numpy()):
        raise ValueError("Reference AnnData labels disagree with frozen split membership.")

    if selector is None:
        # Lazy import keeps this preparation module usable before optional peak
        # selection dependencies or future implementations are imported.
        from deconvatac.pp import select_reference_peaks

        selector = select_reference_peaks
    result = selector(
        training,
        cell_type_key,
        cell_types=list(cell_types),
        n_top_peaks=n_top_peaks,
        min_reference_cells=min_reference_cells,
        scale=1.0e4,
    )
    required = (
        "indices",
        "peak_ids",
        "candidate_feature_sha256",
        "selected_feature_sha256",
        "to_frame",
    )
    missing_result = [name for name in required if not hasattr(result, name)]
    if missing_result:
        raise TypeError(f"Peak selector result is missing fields: {missing_result}")
    peak_ids = tuple(str(value) for value in result.peak_ids)
    if len(peak_ids) != n_top_peaks or len(set(peak_ids)) != len(peak_ids):
        raise ValueError("Peak selector did not return the requested number of unique ranked peaks.")
    return result


def select_peaks_for_outer_splits(
    reference: ad.AnnData,
    universe: pd.DataFrame,
    outer_split_seeds: Sequence[int],
    *,
    cell_types: Sequence[str],
    n_top_peaks: int = DEFAULT_N_PEAKS,
    min_reference_cells: int = DEFAULT_MIN_REFERENCE_CELLS,
    selector: Optional[Callable[..., Any]] = None,
) -> dict[int, tuple[pd.DataFrame, Any]]:
    """Create split memberships and reference-only rankings for several seeds."""
    if len(outer_split_seeds) == 0 or len(set(outer_split_seeds)) != len(
        outer_split_seeds
    ):
        raise ValueError("outer_split_seeds must be a non-empty unique sequence.")
    selected: dict[int, tuple[pd.DataFrame, Any]] = {}
    for seed in outer_split_seeds:
        membership = stratified_reference_heldout_split(
            universe, seed, cell_types=cell_types
        )
        selection = select_reference_only_peaks(
            reference,
            membership,
            cell_types=cell_types,
            n_top_peaks=n_top_peaks,
            min_reference_cells=min_reference_cells,
            selector=selector,
        )
        selected[int(seed)] = (membership, selection)
    return selected


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def shape_provenance(
    source_sha256: Mapping[str, str],
    membership_sha256: str,
) -> dict[str, Any]:
    """Return resolved, validation-ready fragment-shape provenance."""
    if not {"fragments", "tabix_index"}.issubset(source_sha256):
        raise ValueError("source_sha256 must include fragments and tabix_index digests.")
    return {
        "source_sha256": dict(source_sha256),
        "split_sha256": str(membership_sha256),
        "coordinate_validation": {
            "selected_right_cut_offset": 0,
            "matrix_match": "exact",
            "mismatched_entries": 0,
            "absolute_error": 0,
        },
        "software_versions": {
            "deconvatac": _package_version("deconvATAC"),
            "anndata": ad.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "pysam": _package_version("pysam"),
        },
    }


def _shape_cache_from_result(
    result: Any,
    universe: pd.DataFrame,
    source_sha256: Mapping[str, str],
    cell_types: Sequence[str],
) -> ad.AnnData:
    canonical = _canonicalize_universe(universe, cell_types)
    universe_hash = universe_sha256(canonical, cell_types)
    obs = canonical.set_index("barcode")
    obs["cell_type"] = pd.Categorical(
        obs["cell_type"], categories=list(cell_types), ordered=True
    )
    cache = build_fragment_shape_anndata(
        result,
        obs=obs,
        provenance=shape_provenance(source_sha256, universe_hash),
    )
    cache.uns["shapemix_shape_cache"] = {
        "schema_version": PROTOCOL_VERSION,
        "scope": "provided_cells_and_ordered_peak_axis",
        "universe_sha256": universe_hash,
        "cell_types": np.asarray(list(cell_types), dtype=str),
    }
    return cache


def official_order_peak_union(
    selections: Sequence[Any],
    official_peaks: Sequence[PeakInterval],
) -> tuple[PeakInterval, ...]:
    """Return the selected-peak union in official BED order, never rank order."""
    if not selections:
        raise ValueError("At least one peak selection is required to construct a cache.")
    selected_ids: set[str] = set()
    for selection in selections:
        if not hasattr(selection, "peak_ids"):
            raise TypeError("Every peak selection must expose ranked peak_ids.")
        selected_ids.update(str(value) for value in selection.peak_ids)
    official_ids = tuple(peak.name for peak in official_peaks)
    missing = sorted(selected_ids.difference(official_ids))
    if missing:
        raise ValueError(f"Selected peak is absent from the official BED: {missing[0]!r}.")
    union = tuple(peak for peak in official_peaks if peak.name in selected_ids)
    if len(union) != len(selected_ids):
        raise AssertionError("Official-order peak union lost selected peaks.")
    return union


def build_shape_cache_from_records(
    records: Iterable[Union[str, bytes, FragmentRecord]],
    universe: pd.DataFrame,
    peaks: Sequence[Union[PeakInterval, str, Sequence[Any], Mapping[str, Any]]],
    *,
    source_sha256: Mapping[str, str],
    cell_types: Sequence[str],
    chunk_size: int = 1_000_000,
) -> ad.AnnData:
    """Build a full shape cache from an iterable (used by tests and small inputs)."""
    canonical = _canonicalize_universe(universe, cell_types)
    result = count_fragment_shapes_from_records(
        records,
        canonical["barcode"].tolist(),
        peaks,
        right_cut_offset=0,
        chunk_size=chunk_size,
    )
    return _shape_cache_from_result(result, canonical, source_sha256, cell_types)


def build_shape_cache_from_fragments(
    fragments_path: Union[str, Path],
    universe: pd.DataFrame,
    peaks: Sequence[Union[PeakInterval, str, Sequence[Any], Mapping[str, Any]]],
    *,
    source_sha256: Mapping[str, str],
    cell_types: Sequence[str],
    chunk_size: int = 1_000_000,
) -> ad.AnnData:
    """Stream the tabix-indexed fragments file once into the reusable cache."""
    canonical = _canonicalize_universe(universe, cell_types)
    result = count_fragment_shapes(
        fragments_path,
        canonical["barcode"].tolist(),
        peaks,
        right_cut_offset=0,
        chunk_size=chunk_size,
    )
    return _shape_cache_from_result(result, canonical, source_sha256, cell_types)


def _matrix_totals(adata: ad.AnnData) -> tuple[sparse.csr_matrix, dict[str, int]]:
    total = sparse.csr_matrix(adata.shape, dtype=np.int64)
    layer_totals: dict[str, int] = {}
    for layer_name in FRAGMENT_SHAPE_LAYER_NAMES:
        if layer_name not in adata.layers:
            raise ValueError(f"Shape cache is missing layer {layer_name!r}.")
        layer = sparse.csr_matrix(adata.layers[layer_name], dtype=np.int64)
        layer.sum_duplicates()
        layer.eliminate_zeros()
        layer.sort_indices()
        adata.layers[layer_name] = layer
        total = (total + layer).tocsr()
        layer_totals[layer_name] = int(layer.sum())
    total.sum_duplicates()
    total.eliminate_zeros()
    total.sort_indices()
    return total, layer_totals


def validate_shape_cache(
    cache: ad.AnnData,
    universe: pd.DataFrame,
    peak_ids: Sequence[str],
    *,
    source_sha256: Mapping[str, str],
    cell_types: Sequence[str],
) -> None:
    """Validate cache axes, conservation, provenance, and source identity."""
    canonical = _canonicalize_universe(universe, cell_types)
    if tuple(cache.obs_names.astype(str)) != tuple(canonical["barcode"]):
        raise ValueError("Shape-cache cells do not match the canonical executable universe.")
    if tuple(cache.var_names.astype(str)) != tuple(str(value) for value in peak_ids):
        raise ValueError("Shape-cache features do not match the official peak order.")
    if "cell_type" not in cache.obs or not np.array_equal(
        cache.obs["cell_type"].astype(str).to_numpy(),
        canonical["cell_type"].to_numpy(),
    ):
        raise ValueError("Shape-cache cell labels do not match the executable universe.")
    metadata = cache.uns.get("fragment_shape")
    if not isinstance(metadata, Mapping):
        raise ValueError("Shape cache is missing fragment_shape metadata.")
    if metadata.get("right_cut_offset") != 0:
        raise ValueError("Shape cache must use the validated right_cut_offset=0 convention.")
    coordinate_validation = metadata.get("coordinate_validation")
    if not isinstance(coordinate_validation, Mapping) or any(
        coordinate_validation.get(name) != expected
        for name, expected in {
            "selected_right_cut_offset": 0,
            "matrix_match": "exact",
            "mismatched_entries": 0,
            "absolute_error": 0,
        }.items()
    ):
        raise ValueError("Shape cache is missing the exact coordinate-validation result.")
    if not isinstance(metadata.get("preprocessing_counters"), Mapping):
        raise ValueError("Shape cache is missing immutable source-run preprocessing counters.")
    if not isinstance(metadata.get("software_versions"), Mapping):
        raise ValueError("Shape cache is missing software-version provenance.")
    observed_sources = metadata.get("source_sha256")
    if not isinstance(observed_sources, Mapping):
        raise ValueError("Shape cache is missing source hashes.")
    for name, expected in source_sha256.items():
        if observed_sources.get(name) != expected:
            raise ValueError(f"Shape-cache source hash differs for {name!r}.")
    expected_feature_hash = ordered_feature_sha256(str(value) for value in peak_ids)
    if metadata.get("feature_sha256") != expected_feature_hash:
        raise ValueError("Shape-cache feature hash does not match its ordered feature axis.")

    expected_total, layer_totals = _matrix_totals(cache)
    if not matrices_equal(cache.X, expected_total):
        raise ValueError("Shape-cache X does not equal the exact sum of its layers.")
    matrix_counters = metadata.get("matrix_counters")
    if not isinstance(matrix_counters, Mapping):
        raise ValueError("Shape cache is missing matrix counters.")
    if matrix_counters.get("assigned_cut_sites") != int(expected_total.sum()):
        raise ValueError("Shape-cache assigned_cut_sites counter is stale.")
    for layer_name, count in layer_totals.items():
        if matrix_counters.get(f"cut_sites_per_bin.{layer_name}") != count:
            raise ValueError(f"Shape-cache matrix counter is stale for {layer_name!r}.")
    cache_metadata = cache.uns.get("shapemix_shape_cache")
    expected_universe_hash = universe_sha256(canonical, cell_types)
    if not isinstance(cache_metadata, Mapping) or cache_metadata.get(
        "universe_sha256"
    ) != expected_universe_hash:
        raise ValueError("Shape-cache universe hash is missing or stale.")
    if metadata.get("split_sha256") != expected_universe_hash:
        raise ValueError("Shape-cache membership hash is missing or stale.")
    cached_cell_types = tuple(str(value) for value in cache_metadata.get("cell_types", ()))
    if cached_cell_types != tuple(cell_types):
        raise ValueError("Shape-cache cell-type order is missing or stale.")


def validate_shape_cache_against_official(
    cache: ad.AnnData,
    reference: ad.AnnData,
    universe: pd.DataFrame,
) -> None:
    """Require exact cache-X reconstruction on its complete cached peak union."""
    barcodes = universe["barcode"].astype(str).tolist()
    peak_ids = cache.var_names.astype(str).tolist()
    if not matrices_equal(cache.X, reference[barcodes, peak_ids].X):
        raise ValueError(
            "Shape-cache X does not reconstruct the official reference on its peak union."
        )


def derive_shape_pool(
    cache: ad.AnnData,
    barcodes: Sequence[str],
    ranked_peak_ids: Sequence[str],
    *,
    split_sha256: str,
    pool: str,
) -> ad.AnnData:
    """Slice a cache and refresh only feature/split/matrix-specific metadata."""
    if pool not in {"reference", "heldout"}:
        raise ValueError("pool must be 'reference' or 'heldout'.")
    ordered_barcodes = tuple(str(value) for value in barcodes)
    ordered_peaks = tuple(str(value) for value in ranked_peak_ids)
    if not ordered_barcodes or len(set(ordered_barcodes)) != len(ordered_barcodes):
        raise ValueError("Derived pool barcodes must be non-empty and unique.")
    if not ordered_peaks or len(set(ordered_peaks)) != len(ordered_peaks):
        raise ValueError("Ranked peaks must be non-empty and unique.")
    missing_barcodes = sorted(set(ordered_barcodes).difference(cache.obs_names))
    missing_peaks = sorted(set(ordered_peaks).difference(cache.var_names))
    if missing_barcodes:
        raise ValueError(f"Shape cache is missing barcode {missing_barcodes[0]!r}.")
    if missing_peaks:
        raise ValueError(f"Shape cache is missing peak {missing_peaks[0]!r}.")

    derived = cache[list(ordered_barcodes), list(ordered_peaks)].copy()
    metadata = copy.deepcopy(dict(derived.uns.get("fragment_shape", {})))
    preprocessing_counters = copy.deepcopy(metadata.get("preprocessing_counters"))
    total, layer_totals = _matrix_totals(derived)
    derived.X = total
    metadata["feature_sha256"] = ordered_feature_sha256(ordered_peaks)
    metadata["split_sha256"] = str(split_sha256)
    metadata["matrix_counters"] = {
        "assigned_cut_sites": int(total.sum()),
        **{
            f"cut_sites_per_bin.{layer_name}": count
            for layer_name, count in layer_totals.items()
        },
    }
    # The full source scan's counters are deliberately not rewritten after a
    # row/feature slice.  matrix_counters describe this derived object.
    metadata["preprocessing_counters"] = preprocessing_counters
    derived.uns["fragment_shape"] = metadata
    derived.uns["shapemix_preparation"] = {
        "protocol_version": PROTOCOL_VERSION,
        "pool": pool,
        "split_sha256": str(split_sha256),
        "source_shape_cache_universe_sha256": str(
            cache.uns["shapemix_shape_cache"]["universe_sha256"]
        ),
    }
    derived.obs["split_pool"] = pool
    return derived


def matrices_equal(left: Any, right: Any) -> bool:
    """Return exact equality for sparse or dense two-dimensional matrices."""
    left_csr = sparse.csr_matrix(left)
    right_csr = sparse.csr_matrix(right)
    if left_csr.shape != right_csr.shape:
        return False
    difference = left_csr - right_csr
    difference.eliminate_zeros()
    return difference.nnz == 0


def prepare_split_from_objects(
    reference: ad.AnnData,
    shape_cache: ad.AnnData,
    universe: pd.DataFrame,
    outer_split_seed: int,
    *,
    cell_types: Sequence[str],
    n_top_peaks: int = DEFAULT_N_PEAKS,
    min_reference_cells: int = DEFAULT_MIN_REFERENCE_CELLS,
    selector: Optional[Callable[..., Any]] = None,
    membership: Optional[pd.DataFrame] = None,
    peak_selection: Optional[Any] = None,
) -> PreparedSplit:
    """Build and validate one split using already loaded reference/cache objects."""
    if membership is None:
        membership = stratified_reference_heldout_split(
            universe, outer_split_seed, cell_types=cell_types
        )
    else:
        membership = validate_split_membership(membership, universe, cell_types)
    split_hash = split_membership_sha256(membership, cell_types)
    selection = peak_selection
    if selection is None:
        selection = select_reference_only_peaks(
            reference,
            membership,
            cell_types=cell_types,
            n_top_peaks=n_top_peaks,
            min_reference_cells=min_reference_cells,
            selector=selector,
        )
    elif len(tuple(selection.peak_ids)) != n_top_peaks:
        raise ValueError("Precomputed peak selection does not match n_top_peaks.")
    ranked_peaks = tuple(str(value) for value in selection.peak_ids)

    reference_barcodes = membership.loc[
        membership["pool"] == "reference", "barcode"
    ].tolist()
    heldout_barcodes = membership.loc[
        membership["pool"] == "heldout", "barcode"
    ].tolist()
    if set(reference_barcodes).intersection(heldout_barcodes):
        raise AssertionError("Reference and held-out barcode sets overlap.")

    reference_cells = derive_shape_pool(
        shape_cache,
        reference_barcodes,
        ranked_peaks,
        split_sha256=split_hash,
        pool="reference",
    )
    heldout_cells = derive_shape_pool(
        shape_cache,
        heldout_barcodes,
        ranked_peaks,
        split_sha256=split_hash,
        pool="heldout",
    )
    for derived in (reference_cells, heldout_cells):
        derived.uns["shapemix_preparation"]["outer_split_seed"] = int(
            outer_split_seed
        )

    official_reference = reference[reference_barcodes, list(ranked_peaks)].X
    official_heldout = reference[heldout_barcodes, list(ranked_peaks)].X
    if not matrices_equal(reference_cells.X, official_reference):
        raise ValueError("Reference shape layers do not reconstruct official selected-peak counts.")
    if not matrices_equal(heldout_cells.X, official_heldout):
        raise ValueError("Held-out shape layers do not reconstruct official selected-peak counts.")
    return PreparedSplit(
        membership=membership,
        split_sha256=split_hash,
        peak_selection=selection,
        reference_cells=reference_cells,
        heldout_test_cells=heldout_cells,
    )


def atomic_write_bytes(path: Union[str, Path], content: bytes) -> None:
    """Atomically replace one file using a temporary sibling."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_write_yaml(path: Union[str, Path], value: Mapping[str, Any]) -> None:
    """Write deterministic YAML through an atomic file replacement."""
    content = yaml.safe_dump(dict(value), sort_keys=False).encode("utf-8")
    atomic_write_bytes(path, content)


def atomic_write_h5ad(adata: ad.AnnData, path: Union[str, Path]) -> None:
    """Write an H5AD to a temporary sibling and atomically replace the target."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.stem}.", suffix=".h5ad"
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        adata.write_h5ad(temporary, compression="gzip")
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


@contextmanager
def atomic_output_directory(path: Union[str, Path]) -> Iterator[Path]:
    """Stage a new output tree and publish it with one atomic rename."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Output already exists and will not be overwritten: {path}")
    staging = Path(tempfile.mkdtemp(dir=path.parent, prefix=f".{path.name}.staging."))
    try:
        yield staging
        os.replace(staging, path)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def _output_record(path: Path, root: Path, role: str) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(root)),
        "role": role,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def display_path(path: Union[str, Path]) -> str:
    """Prefer a portable project-relative path when the target is in the repo."""
    path = Path(path)
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path)


def write_split_outputs(
    output_dir: Union[str, Path],
    prepared: PreparedSplit,
    *,
    outer_split_seed: int,
    cell_types: Sequence[str],
    input_paths: Mapping[str, Union[str, Path]],
    source_sha256: Mapping[str, str],
    n_top_peaks: int,
    min_reference_cells: int,
    shape_cache_path: Union[str, Path],
) -> Path:
    """Atomically write H5ADs, split/peak audits, and their manifest."""
    output_dir = Path(output_dir)
    canonical_membership = validate_split_membership(
        prepared.membership,
        prepared.membership.loc[:, ["barcode", "cell_type"]],
        cell_types,
    )
    if split_membership_sha256(canonical_membership, cell_types) != prepared.split_sha256:
        raise ValueError("Prepared split hash is stale.")
    selection_frame = prepared.peak_selection.to_frame().copy()
    ranked_peaks = tuple(str(value) for value in prepared.peak_selection.peak_ids)

    with atomic_output_directory(output_dir) as staging:
        reference_path = staging / "reference_cells.h5ad"
        heldout_path = staging / "heldout_test_cells.h5ad"
        split_path = staging / "split.csv"
        selected_path = staging / "selected_peaks.txt"
        selection_audit_path = staging / "peak_selection.csv"

        prepared.reference_cells.write_h5ad(reference_path, compression="gzip")
        prepared.heldout_test_cells.write_h5ad(heldout_path, compression="gzip")
        atomic_write_bytes(split_path, canonical_split_csv(canonical_membership, cell_types))
        atomic_write_bytes(
            selected_path, ("\n".join(ranked_peaks) + "\n").encode("utf-8")
        )
        atomic_write_bytes(
            selection_audit_path,
            selection_frame.to_csv(index=False, lineterminator="\n").encode("utf-8"),
        )
        signal_table, signal_summary = audit_shape_signal(
            prepared.reference_cells,
            labels_key="cell_type",
            cell_types=cell_types,
            split_seed=outer_split_seed,
        )
        signal_audit_path, signal_summary_path = write_signal_audit(
            signal_table,
            signal_summary,
            staging,
        )

        outputs = [
            _output_record(reference_path, staging, "training_reference_fragment_shapes"),
            _output_record(heldout_path, staging, "heldout_source_fragment_shapes"),
            _output_record(split_path, staging, "canonical_split_membership"),
            _output_record(selected_path, staging, "ranked_selected_peak_ids"),
            _output_record(selection_audit_path, staging, "ranked_peak_selection_audit"),
            _output_record(signal_audit_path, staging, "training_reference_shape_signal_audit"),
            _output_record(signal_summary_path, staging, "training_reference_shape_signal_summary"),
        ]
        counts = canonical_membership.groupby(["cell_type", "pool"], sort=False).size()
        fragment_metadata = prepared.reference_cells.uns["fragment_shape"]
        manifest = {
            "schema_version": PROTOCOL_VERSION,
            "dataset": "pbmc_granulocyte_sorted_10k",
            "benchmark_scope": (
                "primary_compatible"
                if tuple(cell_types) == FROZEN_CELL_TYPES and n_top_peaks == DEFAULT_N_PEAKS
                else "development_smoke"
            ),
            "limitation": "Conditional resampling from one donor; not donor-level generalization.",
            "code_provenance": code_provenance(),
            "software_versions": software_versions(),
            "inputs": [
                {
                    "name": name,
                    "path": display_path(input_paths[name]),
                    "sha256": source_sha256[name],
                }
                for name in source_sha256
            ],
            "cell_types": list(cell_types),
            "rng": {
                "numpy_version": np.__version__,
                "bit_generator": "PCG64",
                "seed_sequence": "SeedSequence([20260822, outer_split_seed, one_based_cell_type_order])",
                "namespace": RNG_NAMESPACE,
                "outer_split_seed": int(outer_split_seed),
            },
            "split": {
                "reference_fraction": REFERENCE_FRACTION,
                "sha256": prepared.split_sha256,
                "hash_encoding": "UTF-8 canonical split.csv bytes with LF line endings",
                "reference_cells": int((canonical_membership["pool"] == "reference").sum()),
                "heldout_cells": int((canonical_membership["pool"] == "heldout").sum()),
                "counts_by_cell_type": {
                    cell_type: {
                        pool: int(counts.get((cell_type, pool), 0))
                        for pool in ("reference", "heldout")
                    }
                    for cell_type in cell_types
                },
            },
            "peak_selection": {
                "training_pool_only": True,
                "n_top_peaks": int(n_top_peaks),
                "min_reference_cells": int(min_reference_cells),
                "scale": 1.0e4,
                "candidate_feature_sha256": str(
                    prepared.peak_selection.candidate_feature_sha256
                ),
                "selected_feature_sha256": str(
                    prepared.peak_selection.selected_feature_sha256
                ),
            },
            "fragment_shape": {
                **fragment_shape_declaration(fragment_metadata),
                "feature_sha256": str(fragment_metadata["feature_sha256"]),
                "split_sha256": str(fragment_metadata["split_sha256"]),
            },
            "matrices": {
                "reference_cells": matrix_summary(
                    prepared.reference_cells,
                    FRAGMENT_SHAPE_LAYER_NAMES,
                ),
                "heldout_test_cells": matrix_summary(
                    prepared.heldout_test_cells,
                    FRAGMENT_SHAPE_LAYER_NAMES,
                ),
            },
            "shape_cache": {
                "path": display_path(shape_cache_path),
                "raw_fragments_streamed_once_for_reuse": True,
                "preprocessing_counters_scope": "immutable full-universe source run",
                "derived_matrix_counters_scope": "this split and selected peak axis",
            },
            "signal_audit": {
                "scope": "training_reference_only",
                "heldout_counts_or_labels_read": False,
                "review_required_before_step_4": True,
                "summary_path": "signal_audit_summary.yaml",
            },
            "outputs": outputs,
        }
        atomic_write_yaml(staging / "manifest.yaml", manifest)
    return output_dir


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outer-split-seed",
        type=int,
        action="append",
        dest="outer_split_seeds",
        help="Repeat for multiple splits; defaults to development seed 0.",
    )
    parser.add_argument(
        "--n-peaks",
        type=int,
        help="Defaults to 5000 for the full benchmark or 200 for the frozen smoke universe.",
    )
    parser.add_argument(
        "--min-reference-cells", type=int, default=DEFAULT_MIN_REFERENCE_CELLS
    )
    parser.add_argument(
        "--cell-types",
        nargs="+",
        help="Use the frozen first-three-type development smoke universe; default is all 16 types.",
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--shape-cache", type=Path, default=DEFAULT_SHAPE_CACHE)
    parser.add_argument(
        "--build-shape-cache",
        action="store_true",
        help="Build the reusable full-universe cache if it does not already exist.",
    )
    parser.add_argument("--fragment-chunk-size", type=int, default=1_000_000)
    for name, default in DEFAULT_INPUT_PATHS.items():
        parser.add_argument(f"--{name.replace('_', '-')}", type=Path, default=default)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    input_paths = {
        name: Path(getattr(args, name)) for name in DEFAULT_INPUT_PATHS
    }
    source_hashes = validate_pinned_files(input_paths, PINNED_INPUT_SHA256)
    labels = pd.read_csv(input_paths["labels"])
    peaks = read_peaks_bed(input_paths["peaks"])
    peak_ids = tuple(peak.name for peak in peaks)

    reference = ad.read_h5ad(input_paths["reference"])
    full_universe = build_frozen_pbmc_universe(reference, labels)
    if tuple(reference.var_names.astype(str)) != peak_ids:
        raise ValueError("Prepared reference feature axis does not match the pinned peak BED.")

    if not args.shape_cache.exists():
        if not args.build_shape_cache:
            raise FileNotFoundError(
                f"Shape cache is missing: {args.shape_cache}. Pass --build-shape-cache "
                "once; later splits will reuse it without rescanning fragments."
            )
        cache_selections = select_peaks_for_outer_splits(
            reference,
            full_universe,
            CACHE_OUTER_SPLIT_SEEDS,
            cell_types=FROZEN_CELL_TYPES,
            n_top_peaks=DEFAULT_N_PEAKS,
            min_reference_cells=DEFAULT_MIN_REFERENCE_CELLS,
        )
        smoke_universe = subset_cell_universe(full_universe, SMOKE_CELL_TYPES)
        smoke_selections = select_peaks_for_outer_splits(
            reference,
            smoke_universe,
            (DEFAULT_OUTER_SPLIT_SEED,),
            cell_types=SMOKE_CELL_TYPES,
            n_top_peaks=SMOKE_N_PEAKS,
            min_reference_cells=DEFAULT_MIN_REFERENCE_CELLS,
        )
        cache_peaks = official_order_peak_union(
            [
                *(value[1] for value in cache_selections.values()),
                smoke_selections[DEFAULT_OUTER_SPLIT_SEED][1],
            ],
            peaks,
        )
        del reference
        gc.collect()
        cache = build_shape_cache_from_fragments(
            input_paths["fragments"],
            full_universe,
            cache_peaks,
            source_sha256=source_hashes,
            cell_types=FROZEN_CELL_TYPES,
            chunk_size=args.fragment_chunk_size,
        )
        cache.uns["shapemix_shape_cache"].update(
            {
                "scope": "official_order_union_of_reference_only_selected_peaks",
                "outer_split_seeds": np.asarray(CACHE_OUTER_SPLIT_SEEDS, dtype=np.int64),
                "n_top_peaks_per_split": DEFAULT_N_PEAKS,
                "min_reference_cells": DEFAULT_MIN_REFERENCE_CELLS,
                "smoke_outer_split_seed": DEFAULT_OUTER_SPLIT_SEED,
                "smoke_cell_types": np.asarray(SMOKE_CELL_TYPES, dtype=str),
                "smoke_n_top_peaks": SMOKE_N_PEAKS,
            }
        )
        reference = ad.read_h5ad(input_paths["reference"])
        validate_shape_cache_against_official(cache, reference, full_universe)
        atomic_write_h5ad(cache, args.shape_cache)
        del cache
        gc.collect()

    shape_cache = ad.read_h5ad(args.shape_cache)
    cached_peak_ids = tuple(shape_cache.var_names.astype(str))
    official_positions = {peak_id: index for index, peak_id in enumerate(peak_ids)}
    if any(peak_id not in official_positions for peak_id in cached_peak_ids) or list(
        cached_peak_ids
    ) != sorted(cached_peak_ids, key=official_positions.get):
        raise ValueError("Shape-cache features are not an official-order peak subset.")
    validate_shape_cache(
        shape_cache,
        full_universe,
        cached_peak_ids,
        source_sha256=source_hashes,
        cell_types=FROZEN_CELL_TYPES,
    )
    validate_shape_cache_against_official(shape_cache, reference, full_universe)

    cell_types = tuple(args.cell_types or FROZEN_CELL_TYPES)
    if args.cell_types and cell_types != SMOKE_CELL_TYPES:
        raise ValueError(
            "The version-1 smoke universe is frozen to the first three canonical types: "
            f"{list(SMOKE_CELL_TYPES)}."
        )
    n_top_peaks = args.n_peaks
    if n_top_peaks is None:
        n_top_peaks = SMOKE_N_PEAKS if args.cell_types else DEFAULT_N_PEAKS
    if args.cell_types and (
        n_top_peaks != SMOKE_N_PEAKS
        or args.min_reference_cells != DEFAULT_MIN_REFERENCE_CELLS
    ):
        raise ValueError(
            "The version-1 smoke selector is frozen to 200 peaks and "
            "min_reference_cells=10."
        )
    universe = (
        full_universe
        if cell_types == FROZEN_CELL_TYPES
        else subset_cell_universe(full_universe, cell_types)
    )
    outer_seeds = args.outer_split_seeds or [DEFAULT_OUTER_SPLIT_SEED]
    if len(set(outer_seeds)) != len(outer_seeds):
        raise ValueError("Repeated --outer-split-seed values are not allowed.")
    if cell_types == SMOKE_CELL_TYPES and outer_seeds != [DEFAULT_OUTER_SPLIT_SEED]:
        raise ValueError("The version-1 smoke dataset uses development outer split seed 0 only.")
    requested_selections = select_peaks_for_outer_splits(
        reference,
        universe,
        outer_seeds,
        cell_types=cell_types,
        n_top_peaks=n_top_peaks,
        min_reference_cells=args.min_reference_cells,
    )
    for outer_seed in outer_seeds:
        membership, peak_selection = requested_selections[outer_seed]
        if not set(peak_selection.peak_ids).issubset(shape_cache.var_names):
            raise ValueError(
                "Requested peak ranking is not covered by the reusable cache; "
                "use the frozen cache seed/selector contract or build a separately versioned cache."
            )
        prepared = prepare_split_from_objects(
            reference,
            shape_cache,
            universe,
            outer_seed,
            cell_types=cell_types,
            n_top_peaks=n_top_peaks,
            min_reference_cells=args.min_reference_cells,
            membership=membership,
            peak_selection=peak_selection,
        )
        output_dir = args.output_root / f"split_{outer_seed:03d}"
        write_split_outputs(
            output_dir,
            prepared,
            outer_split_seed=outer_seed,
            cell_types=cell_types,
            input_paths=input_paths,
            source_sha256=source_hashes,
            n_top_peaks=n_top_peaks,
            min_reference_cells=args.min_reference_cells,
            shape_cache_path=args.shape_cache,
        )
        print(f"Wrote {output_dir}", flush=True)


if __name__ == "__main__":
    main()
