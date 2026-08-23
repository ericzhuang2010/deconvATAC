#!/usr/bin/env python
"""Read-only validation of the canonical ShapeMix Step 3 data campaign.

The validator discovers the frozen campaign from the dataset registry, loads
every dataset through the maintained data loader, verifies all recorded file
and code hashes, and independently reconstructs every pseudo-spot from its
held-out source-cell provenance.  It never writes to the repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from scipy import sparse  # noqa: E402

from deconvatac.data import (  # noqa: E402
    FragmentShapeSpec,
    load_deconvolution_input,
    ordered_feature_sha256,
)
from deconvatac.data.registry import (  # noqa: E402
    DEFAULT_REGISTRY_PATH,
    load_dataset_registry,
    resolve_project_path,
)
from deconvatac.data.validators import (  # noqa: E402
    validate_fragment_shape_feature_axis,
)
from scripts.shapemix_provenance import (  # noqa: E402
    code_provenance,
    fragment_shape_declaration,
    matrix_summary,
)


DATASET_PREFIX = "pbmc_granulocyte_sorted_10k_shapemix_"
SEED_NAMESPACE = 20260822
PRIMARY_OUTER_SPLIT_SEEDS = (1103, 2203, 3301, 4409, 5501)
PRIMARY_INNER_MIXTURE_SEEDS = (101, 211)
CONDITION_INDICES = {"equal_celltype": 0, "observed_abundance": 1}
SHAPE_LAYERS = (
    "fragment_length_lt_100",
    "fragment_length_100_249",
    "fragment_length_ge_250",
)
FROZEN_CELL_TYPE_COUNTS: tuple[tuple[str, int], ...] = (
    ("CD14 Mono", 2551),
    ("CD4 Naive", 1382),
    ("CD8 Naive", 1353),
    ("CD4 TCM", 1113),
    ("CD16 Mono", 442),
    ("NK", 403),
    ("CD8 TEM_1", 322),
    ("CD8 TEM_2", 315),
    ("Intermediate B", 300),
    ("Memory B", 298),
    ("CD4 TEM", 286),
    ("cDC", 180),
    ("Treg", 157),
    ("gdT", 143),
    ("MAIT", 130),
    ("Naive B", 125),
)
FROZEN_CELL_TYPES = tuple(cell_type for cell_type, _ in FROZEN_CELL_TYPE_COUNTS)
SMOKE_CELL_TYPES = FROZEN_CELL_TYPES[:3]
PROCESSED_SPLIT_ROOT = (
    ROOT / "data/processed/shapemix/pbmc_granulocyte_sorted_10k"
)


class Step3ValidationError(RuntimeError):
    """Raised when a generated Step 3 artifact violates its frozen contract."""


@dataclass(frozen=True)
class DatasetExpectation:
    dataset_id: str
    condition: str
    outer_split_seed: int
    inner_mixture_seed: int
    smoke: bool

    @property
    def cell_types(self) -> tuple[str, ...]:
        return SMOKE_CELL_TYPES if self.smoke else FROZEN_CELL_TYPES

    @property
    def num_spots(self) -> int:
        return 32 if self.smoke else 1024

    @property
    def num_features(self) -> int:
        return 200 if self.smoke else 5000

    @property
    def grid_shape(self) -> tuple[int, int]:
        return (4, 8) if self.smoke else (32, 32)

    @property
    def dataset_scope(self) -> str:
        return "development_smoke" if self.smoke else "primary"

    @property
    def split_scope(self) -> str:
        return "development_smoke" if self.smoke else "primary_compatible"

    @property
    def split_dir(self) -> Path:
        if self.smoke:
            return PROCESSED_SPLIT_ROOT / "smoke/split_000"
        return PROCESSED_SPLIT_ROOT / f"split_{self.outer_split_seed:03d}"


@dataclass
class SplitBundle:
    expectation: DatasetExpectation
    manifest: dict[str, Any]
    reference: ad.AnnData
    heldout: ad.AnnData
    membership: pd.DataFrame
    selected_peaks: tuple[str, ...]


class HashCache:
    """Hash each unique path at most once during a campaign validation."""

    def __init__(self) -> None:
        self._digests: dict[Path, str] = {}
        self.bytes_read = 0

    def sha256(self, path: Path) -> str:
        resolved = Path(path).resolve()
        cached = self._digests.get(resolved)
        if cached is not None:
            return cached
        digest = hashlib.sha256()
        with resolved.open("rb") as handle:
            for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                self.bytes_read += len(chunk)
                digest.update(chunk)
        value = digest.hexdigest()
        self._digests[resolved] = value
        return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise Step3ValidationError(message)


def _read_yaml(path: Path) -> dict[str, Any]:
    with Path(path).open() as handle:
        value = yaml.safe_load(handle) or {}
    _require(isinstance(value, dict), f"YAML root must be a mapping: {path}")
    return value


def _resolve(path: str | Path, project_root: Path) -> Path:
    return resolve_project_path(path, project_root=project_root).resolve()


def _expected_campaign() -> tuple[DatasetExpectation, ...]:
    smoke = DatasetExpectation(
        dataset_id=(
            "pbmc_granulocyte_sorted_10k_shapemix_"
            "equal_celltype_split_000_mix_000_smoke"
        ),
        condition="equal_celltype",
        outer_split_seed=0,
        inner_mixture_seed=0,
        smoke=True,
    )
    primary = tuple(
        DatasetExpectation(
            dataset_id=(
                f"{DATASET_PREFIX}{condition}_split_{outer_seed:03d}_"
                f"mix_{inner_seed:03d}"
            ),
            condition=condition,
            outer_split_seed=outer_seed,
            inner_mixture_seed=inner_seed,
            smoke=False,
        )
        for outer_seed in PRIMARY_OUTER_SPLIT_SEEDS
        for inner_seed in PRIMARY_INNER_MIXTURE_SEEDS
        for condition in CONDITION_INDICES
    )
    return (smoke, *primary)


def discover_campaign_ids(
    registry: Mapping[str, Any],
) -> tuple[DatasetExpectation, ...]:
    """Require exactly the frozen one-smoke plus twenty-primary registry IDs."""
    expected = _expected_campaign()
    expected_ids = {item.dataset_id for item in expected}
    observed_ids = {
        str(dataset_id)
        for dataset_id in registry
        if str(dataset_id).startswith(DATASET_PREFIX)
    }
    missing = sorted(expected_ids - observed_ids)
    unexpected = sorted(observed_ids - expected_ids)
    _require(
        not missing and not unexpected,
        f"ShapeMix registry campaign differs: missing={missing}, unexpected={unexpected}.",
    )
    _require(len(observed_ids) == 21, f"Expected 21 ShapeMix registry IDs, found {len(observed_ids)}.")
    return expected


def _validate_file_record(
    record: Mapping[str, Any],
    *,
    base_dir: Path,
    hash_cache: HashCache,
    context: str,
    expected_path: Path | None = None,
) -> Path:
    _require(isinstance(record, Mapping), f"{context} file record must be a mapping.")
    raw_path = record.get("path")
    _require(isinstance(raw_path, str) and raw_path, f"{context} has no path.")
    path = _resolve(raw_path, base_dir)
    if expected_path is not None:
        _require(path == expected_path.resolve(), f"{context} path is {path}, expected {expected_path}.")
    _require(path.is_file(), f"{context} file is missing: {path}")
    observed_bytes = path.stat().st_size
    expected_bytes = record.get("bytes")
    if expected_bytes is not None:
        _require(
            isinstance(expected_bytes, int)
            and not isinstance(expected_bytes, bool)
            and expected_bytes == observed_bytes,
            f"{context} byte count is {observed_bytes}, expected {expected_bytes}.",
        )
    expected_hash = record.get("sha256")
    _require(isinstance(expected_hash, str) and len(expected_hash) == 64, f"{context} has no SHA-256.")
    observed_hash = hash_cache.sha256(path)
    _require(observed_hash == expected_hash, f"{context} SHA-256 mismatch for {path}.")
    return path


def _validate_code_provenance(
    manifest: Mapping[str, Any],
    current_code: Mapping[str, Any],
    context: str,
) -> None:
    recorded = manifest.get("code_provenance")
    _require(isinstance(recorded, Mapping), f"{context} is missing code_provenance.")
    recorded_hashes = recorded.get("content_sha256")
    current_hashes = current_code.get("content_sha256")
    _require(
        recorded_hashes == current_hashes,
        f"{context} code/protocol content hashes differ from the current working tree.",
    )
    _require(
        recorded.get("git_commit") == current_code.get("git_commit"),
        f"{context} base Git commit differs from the current checkout.",
    )


def _canonical_csr(matrix: Any) -> sparse.csr_matrix:
    result = matrix.tocsr(copy=True) if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    result.sum_duplicates()
    result.eliminate_zeros()
    result.sort_indices()
    return result


def _assert_sparse_equal(observed: Any, expected: Any, context: str) -> None:
    left = _canonical_csr(observed)
    right = _canonical_csr(expected)
    _require(left.shape == right.shape, f"{context} shape differs: {left.shape} versus {right.shape}.")
    difference = left - right
    difference.eliminate_zeros()
    if difference.nnz:
        absolute_error = int(np.abs(difference.data).sum())
        raise Step3ValidationError(
            f"{context} differs in {difference.nnz} entries (absolute error {absolute_error})."
        )


def _validate_shape_matrix(adata: ad.AnnData, context: str) -> tuple[str, ...]:
    validate_fragment_shape_feature_axis(adata, context)
    spec = FragmentShapeSpec.from_mapping(adata.uns["fragment_shape"])
    _require(spec.layer_names == SHAPE_LAYERS, f"{context} has the wrong shape-layer order.")
    total = sparse.csr_matrix(adata.shape, dtype=np.int64)
    for layer_name in spec.layer_names:
        _require(layer_name in adata.layers, f"{context} is missing layer {layer_name}.")
        total = _canonical_csr(total + _canonical_csr(adata.layers[layer_name]))
    _assert_sparse_equal(adata.X, total, f"{context}.X versus layer sum")
    return spec.layer_names


def _validate_split_manifest_files(
    manifest: Mapping[str, Any],
    split_dir: Path,
    hash_cache: HashCache,
) -> dict[str, Path]:
    inputs = manifest.get("inputs")
    _require(isinstance(inputs, list) and len(inputs) == 7, f"Split {split_dir} must record seven inputs.")
    input_names: set[str] = set()
    for record in inputs:
        _require(isinstance(record, Mapping), f"Split input record is malformed in {split_dir}.")
        name = record.get("name")
        _require(isinstance(name, str) and name not in input_names, f"Duplicate split input {name!r}.")
        input_names.add(name)
        _validate_file_record(
            record,
            base_dir=ROOT,
            hash_cache=hash_cache,
            context=f"split {split_dir.name} input {name}",
        )
    _require(
        input_names
        == {"labels", "label_summary", "reference", "peaks", "filtered_feature_matrix", "fragments", "tabix_index"},
        f"Split {split_dir} input names differ from the frozen source set.",
    )

    outputs = manifest.get("outputs")
    _require(isinstance(outputs, list) and len(outputs) == 7, f"Split {split_dir} must record seven outputs.")
    by_role: dict[str, Path] = {}
    for record in outputs:
        _require(isinstance(record, Mapping), f"Split output record is malformed in {split_dir}.")
        role = record.get("role")
        _require(isinstance(role, str) and role not in by_role, f"Duplicate split output role {role!r}.")
        by_role[role] = _validate_file_record(
            record,
            base_dir=split_dir,
            hash_cache=hash_cache,
            context=f"split {split_dir.name} output {role}",
        )
    expected_roles = {
        "training_reference_fragment_shapes",
        "heldout_source_fragment_shapes",
        "canonical_split_membership",
        "ranked_selected_peak_ids",
        "ranked_peak_selection_audit",
        "training_reference_shape_signal_audit",
        "training_reference_shape_signal_summary",
    }
    _require(set(by_role) == expected_roles, f"Split {split_dir} output roles differ from the contract.")
    return by_role


def _validate_signal_outputs(
    *,
    manifest: Mapping[str, Any],
    output_paths: Mapping[str, Path],
    expectation: DatasetExpectation,
    reference: ad.AnnData,
    selected_peaks: Sequence[str],
    hash_cache: HashCache,
) -> None:
    summary_path = output_paths["training_reference_shape_signal_summary"]
    audit_path = output_paths["training_reference_shape_signal_audit"]
    summary = _read_yaml(summary_path)
    _require(summary.get("scope") == "training_reference_only", f"Signal scope is wrong in {summary_path}.")
    _require(summary.get("split_seed") == expectation.outer_split_seed, f"Signal split seed is wrong in {summary_path}.")
    _require(summary.get("cell_types") == list(expectation.cell_types), f"Signal cell types are wrong in {summary_path}.")
    _require(summary.get("cells") == reference.n_obs, f"Signal cell count is wrong in {summary_path}.")
    _require(summary.get("peaks") == reference.n_vars, f"Signal peak count is wrong in {summary_path}.")
    _require(summary.get("positive_peaks") == reference.n_vars, f"Signal audit found a zero selected peak in {summary_path}.")
    _require(summary.get("layer_names") == list(SHAPE_LAYERS), f"Signal layers are wrong in {summary_path}.")
    _require(
        summary.get("feature_sha256") == manifest["fragment_shape"]["feature_sha256"],
        f"Signal feature hash is wrong in {summary_path}.",
    )
    _require(
        summary.get("split_sha256") == manifest["split"]["sha256"],
        f"Signal split hash is wrong in {summary_path}.",
    )
    expected_streams = {
        cell_type: [SEED_NAMESPACE, expectation.outer_split_seed, 23, type_index]
        for type_index, cell_type in enumerate(expectation.cell_types, start=1)
    }
    rng = summary.get("random_number_generator") or {}
    _require(rng.get("bit_generator") == "PCG64", f"Signal bit generator is wrong in {summary_path}.")
    _require(rng.get("seed_streams_by_cell_type") == expected_streams, f"Signal RNG streams are wrong in {summary_path}.")

    global_counts = summary.get("global_bin_counts")
    _require(isinstance(global_counts, Mapping), f"Signal global counts are missing in {summary_path}.")
    for layer_name in SHAPE_LAYERS:
        _require(
            global_counts.get(layer_name) == int(reference.layers[layer_name].sum()),
            f"Signal global count for {layer_name} is stale in {summary_path}.",
        )
    embedded_files = summary.get("files") or {}
    embedded_audit = embedded_files.get("signal_audit.csv")
    _validate_file_record(
        {"path": str(audit_path), **dict(embedded_audit or {})},
        base_dir=ROOT,
        hash_cache=hash_cache,
        context=f"signal audit embedded record for split {expectation.outer_split_seed}",
        expected_path=audit_path,
    )

    audit = pd.read_csv(audit_path)
    _require(audit["peak_id"].astype(str).tolist() == list(selected_peaks), f"Signal peak order is wrong in {audit_path}.")
    _require(
        (audit["nonzero_reference_cells"] >= 10).all(),
        f"Signal audit contains a selected peak below minimum coverage in {audit_path}.",
    )
    numeric = audit.select_dtypes(include=[np.number]).to_numpy(dtype=np.float64)
    _require(np.isfinite(numeric).all(), f"Signal audit contains non-finite values in {audit_path}.")


def validate_split_bundle(
    expectation: DatasetExpectation,
    *,
    current_code: Mapping[str, Any],
    hash_cache: HashCache,
) -> SplitBundle:
    split_dir = expectation.split_dir
    manifest_path = split_dir / "manifest.yaml"
    _require(manifest_path.is_file(), f"Split manifest is missing: {manifest_path}")
    manifest = _read_yaml(manifest_path)
    context = f"split {expectation.outer_split_seed}"
    _validate_code_provenance(manifest, current_code, context)
    _require(manifest.get("benchmark_scope") == expectation.split_scope, f"{context} scope is wrong.")
    _require(manifest.get("cell_types") == list(expectation.cell_types), f"{context} cell types are wrong.")
    rng = manifest.get("rng") or {}
    _require(rng.get("namespace") == SEED_NAMESPACE, f"{context} RNG namespace is wrong.")
    _require(rng.get("bit_generator") == "PCG64", f"{context} bit generator is wrong.")
    _require(rng.get("outer_split_seed") == expectation.outer_split_seed, f"{context} outer seed is wrong.")

    paths = _validate_split_manifest_files(manifest, split_dir, hash_cache)
    membership_path = paths["canonical_split_membership"]
    membership = pd.read_csv(membership_path, dtype=str)
    _require(list(membership.columns) == ["barcode", "cell_type", "pool"], f"{context} split columns are wrong.")
    _require(not membership.isna().any().any(), f"{context} split contains missing values.")
    _require(not membership["barcode"].duplicated().any(), f"{context} split contains duplicate barcodes.")
    _require(set(membership["pool"]) == {"reference", "heldout"}, f"{context} pools are wrong.")
    expected_order = {cell_type: index for index, cell_type in enumerate(expectation.cell_types)}
    ordered = membership.assign(_order=membership["cell_type"].map(expected_order)).sort_values(
        ["_order", "barcode"], kind="stable"
    )
    _require(
        membership.reset_index(drop=True).equals(ordered.drop(columns="_order").reset_index(drop=True)),
        f"{context} split rows are not in canonical type/barcode order.",
    )

    count_lookup = dict(FROZEN_CELL_TYPE_COUNTS)
    expected_counts_by_type = {
        cell_type: {
            "reference": int(math.floor(0.70 * count_lookup[cell_type])),
            "heldout": count_lookup[cell_type] - int(math.floor(0.70 * count_lookup[cell_type])),
        }
        for cell_type in expectation.cell_types
    }
    _require(
        manifest.get("split", {}).get("counts_by_cell_type") == expected_counts_by_type,
        f"{context} per-type split counts differ from the frozen 70/30 counts.",
    )
    _require(
        hash_cache.sha256(membership_path) == manifest["split"]["sha256"],
        f"{context} canonical split hash is stale.",
    )

    selected_peaks = tuple(
        line.strip()
        for line in paths["ranked_selected_peak_ids"].read_text().splitlines()
        if line.strip()
    )
    _require(len(selected_peaks) == expectation.num_features, f"{context} selected peak count is wrong.")
    _require(len(set(selected_peaks)) == len(selected_peaks), f"{context} selected peaks are duplicated.")
    expected_feature_hash = ordered_feature_sha256(selected_peaks)
    _require(
        expected_feature_hash == manifest["peak_selection"]["selected_feature_sha256"],
        f"{context} selected-peak hash is stale.",
    )
    _require(manifest["peak_selection"]["n_top_peaks"] == expectation.num_features, f"{context} n_top_peaks is wrong.")
    _require(manifest["peak_selection"]["min_reference_cells"] == 10, f"{context} coverage threshold is wrong.")

    peak_audit = pd.read_csv(paths["ranked_peak_selection_audit"])
    _require(peak_audit["rank"].tolist() == list(range(1, expectation.num_features + 1)), f"{context} peak ranks are wrong.")
    _require(peak_audit["peak_id"].astype(str).tolist() == list(selected_peaks), f"{context} peak audit order is wrong.")
    _require((peak_audit["nonzero_reference_cells"] >= 10).all(), f"{context} selected peak is under-covered.")

    reference = ad.read_h5ad(paths["training_reference_fragment_shapes"])
    heldout = ad.read_h5ad(paths["heldout_source_fragment_shapes"])
    for pool, adata in (("reference", reference), ("heldout", heldout)):
        expected_rows = membership.loc[membership["pool"] == pool]
        expected_names = expected_rows["barcode"].tolist()
        _require(adata.obs_names.astype(str).tolist() == expected_names, f"{context} {pool} observation order is wrong.")
        _require(adata.var_names.astype(str).tolist() == list(selected_peaks), f"{context} {pool} feature order is wrong.")
        _require(
            adata.obs["cell_type"].astype(str).tolist() == expected_rows["cell_type"].tolist(),
            f"{context} {pool} labels differ from split.csv.",
        )
        _require(set(adata.obs["split_pool"].astype(str)) == {pool}, f"{context} {pool} split_pool is wrong.")
        preparation = adata.uns.get("shapemix_preparation") or {}
        _require(preparation.get("pool") == pool, f"{context} {pool} preparation scope is wrong.")
        _require(preparation.get("outer_split_seed") == expectation.outer_split_seed, f"{context} {pool} outer seed is wrong.")
        _require(preparation.get("split_sha256") == manifest["split"]["sha256"], f"{context} {pool} split hash is wrong.")
        _validate_shape_matrix(adata, f"{context} {pool}")
        _require(
            matrix_summary(adata, SHAPE_LAYERS) == manifest["matrices"][f"{pool}_cells" if pool == "reference" else "heldout_test_cells"],
            f"{context} {pool} matrix summary is stale.",
        )

    expected_reference = sum(values["reference"] for values in expected_counts_by_type.values())
    expected_heldout = sum(values["heldout"] for values in expected_counts_by_type.values())
    _require(reference.shape == (expected_reference, expectation.num_features), f"{context} reference shape is wrong.")
    _require(heldout.shape == (expected_heldout, expectation.num_features), f"{context} heldout shape is wrong.")
    _require(
        manifest["fragment_shape"]["feature_sha256"] == expected_feature_hash,
        f"{context} fragment-shape feature hash is stale.",
    )
    _require(
        fragment_shape_declaration(reference.uns["fragment_shape"])
        == {key: manifest["fragment_shape"][key] for key in fragment_shape_declaration(reference.uns["fragment_shape"])},
        f"{context} fragment-shape declaration is wrong.",
    )
    _validate_signal_outputs(
        manifest=manifest,
        output_paths=paths,
        expectation=expectation,
        reference=reference,
        selected_peaks=selected_peaks,
        hash_cache=hash_cache,
    )
    return SplitBundle(
        expectation=expectation,
        manifest=manifest,
        reference=reference,
        heldout=heldout,
        membership=membership,
        selected_peaks=selected_peaks,
    )


def build_provenance_assignment(
    provenance_path: Path,
    *,
    spot_names: Sequence[str],
    heldout: ad.AnnData,
    cell_types: Sequence[str],
) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    """Build spot-by-heldout assignment and exact type counts from JSONL."""
    spot_names = tuple(str(value) for value in spot_names)
    cell_types = tuple(str(value) for value in cell_types)
    heldout_names = heldout.obs_names.astype(str)
    heldout_labels = heldout.obs["cell_type"].astype(str).to_numpy()
    type_index = {cell_type: index for index, cell_type in enumerate(cell_types)}
    pool_sizes = {
        cell_type: int(np.count_nonzero(heldout_labels == cell_type))
        for cell_type in cell_types
    }
    row_indices: list[int] = []
    column_indices: list[int] = []
    truth_counts = np.zeros((len(spot_names), len(cell_types)), dtype=np.int64)
    cell_counts = np.zeros(len(spot_names), dtype=np.int64)
    seen_spots: set[str] = set()

    with Path(provenance_path).open() as handle:
        rows = [json.loads(line) for line in handle if line.strip()]
    _require(len(rows) == len(spot_names), f"Provenance row count is wrong in {provenance_path}.")
    for spot_index, (expected_spot, row) in enumerate(zip(spot_names, rows)):
        _require(isinstance(row, Mapping), f"Malformed provenance row {spot_index} in {provenance_path}.")
        _require(set(row) == {"spot_id", "cell_count", "cell_id", "cell_type"}, f"Provenance keys are wrong for {expected_spot}.")
        spot_id = row["spot_id"]
        _require(spot_id == expected_spot and spot_id not in seen_spots, f"Provenance spot order is wrong at {expected_spot}.")
        seen_spots.add(str(spot_id))
        source_ids = row["cell_id"]
        source_types = row["cell_type"]
        declared_count = row["cell_count"]
        _require(
            isinstance(declared_count, int)
            and not isinstance(declared_count, bool)
            and declared_count > 0,
            f"Invalid cell count for {expected_spot}.",
        )
        _require(
            isinstance(source_ids, list)
            and isinstance(source_types, list)
            and len(source_ids) == len(source_types) == declared_count,
            f"Source provenance lengths disagree for {expected_spot}.",
        )
        positions = heldout_names.get_indexer(pd.Index([str(value) for value in source_ids]))
        _require((positions >= 0).all(), f"Non-heldout source barcode found for {expected_spot}.")
        observed_types = heldout_labels[positions]
        _require(
            observed_types.tolist() == [str(value) for value in source_types],
            f"Source barcode/type pairing is wrong for {expected_spot}.",
        )
        for cell_type in cell_types:
            selected_ids = [
                str(source_id)
                for source_id, source_type in zip(source_ids, source_types)
                if str(source_type) == cell_type
            ]
            requested = len(selected_ids)
            if requested <= pool_sizes[cell_type]:
                _require(
                    len(set(selected_ids)) == requested,
                    f"Source barcode was reused unnecessarily within {expected_spot} for {cell_type}.",
                )
        row_indices.extend([spot_index] * declared_count)
        column_indices.extend(positions.tolist())
        for cell_type in source_types:
            key = str(cell_type)
            _require(key in type_index, f"Unknown source cell type {key!r} for {expected_spot}.")
            truth_counts[spot_index, type_index[key]] += 1
        cell_counts[spot_index] = declared_count

    assignment = sparse.coo_matrix(
        (
            np.ones(len(row_indices), dtype=np.int64),
            (np.asarray(row_indices), np.asarray(column_indices)),
        ),
        shape=(len(spot_names), heldout.n_obs),
        dtype=np.int64,
    ).tocsr()
    assignment.sum_duplicates()
    _require(
        np.array_equal(np.asarray(assignment.sum(axis=1)).ravel(), cell_counts),
        f"Assignment matrix cell counts are wrong for {provenance_path}.",
    )
    return assignment, truth_counts, cell_counts


def _condition_probabilities(expectation: DatasetExpectation) -> np.ndarray:
    if expectation.condition == "equal_celltype":
        return np.full(len(expectation.cell_types), 1.0 / len(expectation.cell_types))
    counts = dict(FROZEN_CELL_TYPE_COUNTS)
    values = np.asarray([counts[cell_type] for cell_type in expectation.cell_types], dtype=np.float64)
    return values / values.sum()


def _validate_simulation_file_records(
    *,
    manifest: Mapping[str, Any],
    expectation: DatasetExpectation,
    dataset_dir: Path,
    split_bundle: SplitBundle,
    hash_cache: HashCache,
) -> dict[str, Path]:
    inputs = manifest.get("inputs")
    _require(isinstance(inputs, Mapping), f"Simulation inputs are missing for {expectation.dataset_id}.")
    expected_inputs = {
        "reference_h5ad": expectation.split_dir / "reference_cells.h5ad",
        "heldout_h5ad": expectation.split_dir / "heldout_test_cells.h5ad",
        "split_manifest": expectation.split_dir / "manifest.yaml",
    }
    _require(set(inputs) == set(expected_inputs), f"Simulation input roles are wrong for {expectation.dataset_id}.")
    for role, expected_path in expected_inputs.items():
        _validate_file_record(
            inputs[role],
            base_dir=ROOT,
            hash_cache=hash_cache,
            context=f"{expectation.dataset_id} input {role}",
            expected_path=expected_path,
        )

    outputs = manifest.get("outputs")
    _require(isinstance(outputs, Mapping), f"Simulation outputs are missing for {expectation.dataset_id}.")
    expected_output_paths = {
        "spatial_h5ad": dataset_dir / "atac/spatial.h5ad",
        "selected_peaks": dataset_dir / "atac/features/highly_variable.txt",
        "truth_proportions": dataset_dir / "truth/proportions.csv",
        "source_cells_by_spot": dataset_dir / "simulation/source_cells_by_spot.jsonl",
        "dataset_yaml": dataset_dir / "dataset.yaml",
    }
    if expectation.smoke:
        expected_output_paths["reference_h5ad"] = dataset_dir / "atac/reference_cells.h5ad"
    _require(set(outputs) == set(expected_output_paths), f"Simulation output roles are wrong for {expectation.dataset_id}.")
    validated: dict[str, Path] = {}
    for role, expected_path in expected_output_paths.items():
        validated[role] = _validate_file_record(
            outputs[role],
            base_dir=ROOT,
            hash_cache=hash_cache,
            context=f"{expectation.dataset_id} output {role}",
            expected_path=expected_path,
        )
    return validated


def validate_dataset(
    expectation: DatasetExpectation,
    *,
    registry_path: Path,
    registry: Mapping[str, Any],
    project_root: Path,
    split_bundle: SplitBundle,
    current_code: Mapping[str, Any],
    hash_cache: HashCache,
) -> dict[str, Any]:
    started = time.perf_counter()
    entry = registry[expectation.dataset_id]
    config_value = entry if isinstance(entry, str) else entry.get("config") or entry.get("path")
    _require(isinstance(config_value, str), f"Registry entry is malformed for {expectation.dataset_id}.")
    config_path = _resolve(config_value, project_root)
    dataset_dir = config_path.parent
    expected_dataset_dir = ROOT / "data/processed/datasets" / expectation.dataset_id
    _require(dataset_dir == expected_dataset_dir.resolve(), f"Registry path is wrong for {expectation.dataset_id}.")
    config = _read_yaml(config_path)
    _require(config.get("dataset_id") == expectation.dataset_id, f"Dataset ID is wrong in {config_path}.")
    _require(config.get("benchmark_scope") == expectation.dataset_scope, f"Dataset scope is wrong for {expectation.dataset_id}.")
    simulation_config = config.get("simulation") or {}
    _require(simulation_config.get("condition") == expectation.condition, f"Config condition is wrong for {expectation.dataset_id}.")
    _require(simulation_config.get("outer_split_seed") == expectation.outer_split_seed, f"Config outer seed is wrong for {expectation.dataset_id}.")
    _require(simulation_config.get("inner_mixture_seed") == expectation.inner_mixture_seed, f"Config inner seed is wrong for {expectation.dataset_id}.")
    _require(simulation_config.get("primary_dataset") is (not expectation.smoke), f"Config primary flag is wrong for {expectation.dataset_id}.")
    _require(simulation_config.get("depth_retain_probability") is None, f"Primary/smoke depth thinning is enabled for {expectation.dataset_id}.")

    manifest_path = _resolve(simulation_config.get("manifest", ""), project_root)
    _require(manifest_path == dataset_dir / "simulation/manifest.yaml", f"Simulation manifest path is wrong for {expectation.dataset_id}.")
    manifest = _read_yaml(manifest_path)
    _validate_code_provenance(manifest, current_code, expectation.dataset_id)
    _require(manifest.get("dataset_id") == expectation.dataset_id, f"Manifest ID is wrong for {expectation.dataset_id}.")
    _require(manifest.get("status") == "complete", f"Manifest status is not complete for {expectation.dataset_id}.")
    _require(manifest.get("benchmark_scope") == expectation.dataset_scope, f"Manifest scope is wrong for {expectation.dataset_id}.")
    output_paths = _validate_simulation_file_records(
        manifest=manifest,
        expectation=expectation,
        dataset_dir=dataset_dir,
        split_bundle=split_bundle,
        hash_cache=hash_cache,
    )

    simulation = manifest.get("simulation") or {}
    _require(simulation.get("condition") == expectation.condition, f"Manifest condition is wrong for {expectation.dataset_id}.")
    _require(simulation.get("num_spots") == expectation.num_spots, f"Spot count is wrong for {expectation.dataset_id}.")
    _require(simulation.get("num_features") == expectation.num_features, f"Feature count is wrong for {expectation.dataset_id}.")
    _require(simulation.get("grid_shape") == list(expectation.grid_shape), f"Grid is wrong for {expectation.dataset_id}.")
    _require(simulation.get("mean_cells_per_spot") == 10.0, f"Mean cells/spot is wrong for {expectation.dataset_id}.")
    _require(simulation.get("cell_types") == list(expectation.cell_types), f"Manifest cell types are wrong for {expectation.dataset_id}.")
    probabilities = np.asarray(
        [simulation["sampling_probabilities"][cell_type] for cell_type in expectation.cell_types],
        dtype=np.float64,
    )
    _require(
        np.allclose(probabilities, _condition_probabilities(expectation), rtol=0.0, atol=1.0e-15),
        f"Sampling probabilities are wrong for {expectation.dataset_id}.",
    )
    rng = simulation.get("random_number_generator") or {}
    condition_index = CONDITION_INDICES[expectation.condition]
    expected_rng_streams = {
        "cell_counts_and_types": [
            SEED_NAMESPACE,
            expectation.outer_split_seed,
            expectation.inner_mixture_seed,
            condition_index,
        ],
        "source_cells": [
            SEED_NAMESPACE,
            expectation.outer_split_seed,
            expectation.inner_mixture_seed,
            condition_index,
            1,
        ],
    }
    _require(rng.get("bit_generator") == "PCG64", f"Bit generator is wrong for {expectation.dataset_id}.")
    _require(rng.get("seed_namespace") == SEED_NAMESPACE, f"RNG namespace is wrong for {expectation.dataset_id}.")
    _require(rng.get("outer_split_seed") == expectation.outer_split_seed, f"Manifest outer seed is wrong for {expectation.dataset_id}.")
    _require(rng.get("inner_mixture_seed") == expectation.inner_mixture_seed, f"Manifest inner seed is wrong for {expectation.dataset_id}.")
    _require(rng.get("seed_streams") == expected_rng_streams, f"RNG streams are wrong for {expectation.dataset_id}.")
    thinning = simulation.get("depth_thinning") or {}
    _require(thinning == {"enabled": False, "retain_probability": None, "primary_dataset": not expectation.smoke}, f"Depth-thinning contract is wrong for {expectation.dataset_id}.")

    data = load_deconvolution_input(
        expectation.dataset_id,
        "atac",
        feature_set="all",
        registry_path=registry_path,
        project_root=project_root,
    )
    _require(data.spatial.shape == (expectation.num_spots, expectation.num_features), f"Loaded spatial shape is wrong for {expectation.dataset_id}.")
    _require(data.reference.n_vars == expectation.num_features, f"Loaded reference feature count is wrong for {expectation.dataset_id}.")
    _require(data.cell_types == list(expectation.cell_types), f"Loaded cell-type order is wrong for {expectation.dataset_id}.")
    _require(data.truth is not None, f"Truth is missing for {expectation.dataset_id}.")
    expected_spots = [f"spot_{index:04d}" for index in range(expectation.num_spots)]
    _require(data.spatial.obs_names.astype(str).tolist() == expected_spots, f"Spot order is wrong for {expectation.dataset_id}.")
    _require(
        np.array_equal(
            np.asarray(data.spatial.obsm["spatial"]),
            np.asarray(
                [
                    (index % expectation.grid_shape[1], index // expectation.grid_shape[1])
                    for index in range(expectation.num_spots)
                ],
                dtype=np.int64,
            ),
        ),
        f"Spatial row-major coordinates are wrong for {expectation.dataset_id}.",
    )
    _require(data.spatial.var_names.equals(split_bundle.heldout.var_names), f"Heldout/spatial feature axes differ for {expectation.dataset_id}.")

    provenance_path = output_paths["source_cells_by_spot"]
    assignment, truth_counts, cell_counts = build_provenance_assignment(
        provenance_path,
        spot_names=expected_spots,
        heldout=split_bundle.heldout,
        cell_types=expectation.cell_types,
    )
    _require(
        np.array_equal(data.spatial.obs["cell_count"].to_numpy(dtype=np.int64), cell_counts),
        f"Spatial cell counts differ from provenance for {expectation.dataset_id}.",
    )
    expected_truth = truth_counts / cell_counts[:, np.newaxis]
    loaded_truth = data.truth.to_numpy(dtype=np.float64)
    _require(
        np.allclose(loaded_truth, expected_truth, rtol=0.0, atol=1.0e-15)
        and np.array_equal(
            np.rint(loaded_truth * cell_counts[:, np.newaxis]).astype(np.int64),
            truth_counts,
        ),
        (
            "Truth CSV does not serialize the exact source-cell fractions/counts for "
            f"{expectation.dataset_id}."
        ),
    )
    _require(
        np.array_equal(np.asarray(data.spatial.obsm["proportions"]), expected_truth),
        f"Spatial embedded truth differs from provenance for {expectation.dataset_id}.",
    )

    expected_x = sparse.csr_matrix(data.spatial.shape, dtype=np.int64)
    for layer_name in SHAPE_LAYERS:
        expected_layer = _canonical_csr(
            assignment @ _canonical_csr(split_bundle.heldout.layers[layer_name])
        )
        _assert_sparse_equal(
            data.spatial.layers[layer_name],
            expected_layer,
            f"{expectation.dataset_id} layer {layer_name}",
        )
        expected_x = _canonical_csr(expected_x + expected_layer)
    _assert_sparse_equal(data.spatial.X, expected_x, f"{expectation.dataset_id}.X reconstructed from heldout cells")
    _require(
        matrix_summary(data.spatial, SHAPE_LAYERS) == simulation["matrix"],
        f"Simulation matrix summary is stale for {expectation.dataset_id}.",
    )
    _require(
        manifest["fragment_shape"]["feature_sha256"] == ordered_feature_sha256(data.spatial.var_names),
        f"Simulation feature hash is stale for {expectation.dataset_id}.",
    )
    _require(
        manifest["fragment_shape"]["split_sha256"] == split_bundle.manifest["split"]["sha256"],
        f"Simulation split hash is stale for {expectation.dataset_id}.",
    )
    return {
        "dataset_id": expectation.dataset_id,
        "scope": expectation.dataset_scope,
        "spots": expectation.num_spots,
        "features": expectation.num_features,
        "source_assignments": int(cell_counts.sum()),
        "seconds": time.perf_counter() - started,
    }


def validate_campaign(
    *,
    registry_path: Path = DEFAULT_REGISTRY_PATH,
    project_root: Path = ROOT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate all canonical Step 3 splits and datasets without writing files."""
    started = time.perf_counter()
    registry_path = Path(registry_path).resolve()
    project_root = Path(project_root).resolve()
    registry = load_dataset_registry(registry_path)
    campaign = discover_campaign_ids(registry)
    current_code = code_provenance()
    hash_cache = HashCache()
    results: list[dict[str, Any]] = []

    groups: dict[tuple[bool, int], list[DatasetExpectation]] = defaultdict(list)
    for expectation in campaign:
        groups[(expectation.smoke, expectation.outer_split_seed)].append(expectation)

    for group_key in sorted(groups, key=lambda value: (not value[0], value[1])):
        expectations = groups[group_key]
        split_started = time.perf_counter()
        split_bundle = validate_split_bundle(
            expectations[0],
            current_code=current_code,
            hash_cache=hash_cache,
        )
        if verbose:
            print(
                f"PASS split {expectations[0].outer_split_seed:04d} "
                f"({split_bundle.reference.n_obs} reference, "
                f"{split_bundle.heldout.n_obs} heldout, "
                f"{split_bundle.heldout.n_vars} peaks) "
                f"[{time.perf_counter() - split_started:.2f}s]",
                flush=True,
            )
        for expectation in expectations:
            result = validate_dataset(
                expectation,
                registry_path=registry_path,
                registry=registry,
                project_root=project_root,
                split_bundle=split_bundle,
                current_code=current_code,
                hash_cache=hash_cache,
            )
            results.append(result)
            if verbose:
                print(
                    f"PASS {expectation.dataset_id} "
                    f"({result['source_assignments']} source assignments) "
                    f"[{result['seconds']:.2f}s]",
                    flush=True,
                )
        del split_bundle

    elapsed = time.perf_counter() - started
    report = {
        "status": "pass",
        "datasets": len(results),
        "primary_datasets": sum(result["scope"] == "primary" for result in results),
        "smoke_datasets": sum(result["scope"] == "development_smoke" for result in results),
        "splits": len(groups),
        "source_assignments": sum(result["source_assignments"] for result in results),
        "unique_hashed_bytes": hash_cache.bytes_read,
        "seconds": elapsed,
        "dataset_results": results,
    }
    if verbose:
        print(
            "PASS ShapeMix Step 3 campaign: "
            f"{report['smoke_datasets']} smoke + {report['primary_datasets']} primary, "
            f"{report['splits']} splits, {report['source_assignments']} source assignments, "
            f"{report['unique_hashed_bytes'] / (1024 ** 3):.2f} GiB hashed "
            f"[{elapsed:.2f}s]",
            flush=True,
        )
    return report


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY_PATH)
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--json", action="store_true", help="Print the final report as JSON.")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    report = validate_campaign(
        registry_path=args.registry,
        project_root=args.project_root,
        verbose=not args.quiet,
    )
    if args.json:
        print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
