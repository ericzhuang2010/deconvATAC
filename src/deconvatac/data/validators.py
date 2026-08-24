from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable, Mapping
from numbers import Integral
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from .schemas import DeconvolutionInput, FragmentShapeBin, FragmentShapeSpec


_CANONICAL_FRAGMENT_SHAPE_BINS = (
    FragmentShapeBin(
        name="short",
        min_inclusive=0,
        max_exclusive=100,
        layer="fragment_length_lt_100",
    ),
    FragmentShapeBin(
        name="mono",
        min_inclusive=100,
        max_exclusive=250,
        layer="fragment_length_100_249",
    ),
    FragmentShapeBin(
        name="long",
        min_inclusive=250,
        max_exclusive=None,
        layer="fragment_length_ge_250",
    ),
)
_SHA256_PATTERN = re.compile(r"^[0-9a-fA-F]{64}$")
_FRAGMENT_SHAPE_COUNTER_FIELDS = {
    "total_rows",
    "header_rows",
    "invalid_schema_rows",
    "invalid_coordinate_rows",
    "unknown_barcodes",
    "filtered_contigs",
    "valid_rows",
    "retained_fragments",
    "fragments_with_assigned_cut_sites",
    "cut_sites_outside_peaks",
    "assigned_cut_sites",
    "read_support_total",
}


def ordered_feature_sha256(feature_names: Iterable[str]) -> str:
    """Hash an ordered feature axis using UTF-8 names with 8-byte length prefixes.

    Length prefixes make boundaries unambiguous: axes ``["ab", "c"]`` and
    ``["a", "bc"]`` cannot produce the same byte stream before hashing.
    """
    digest = hashlib.sha256()
    for feature_name in feature_names:
        if not isinstance(feature_name, str):
            raise TypeError("Feature names must be strings when computing feature_sha256.")
        encoded = feature_name.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, byteorder="big", signed=False))
        digest.update(encoded)
    return digest.hexdigest()


def normalize_proportions(values: pd.DataFrame, zero_policy: str = "zeros") -> pd.DataFrame:
    """Row-normalize non-negative abundance-like values into proportions."""
    if zero_policy not in {"zeros", "uniform"}:
        raise ValueError("zero_policy must be one of {'zeros', 'uniform'}.")

    numeric = values.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    numeric = numeric.clip(lower=0)
    row_sums = numeric.sum(axis=1)

    proportions = numeric.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)
    if zero_policy == "uniform":
        zero_rows = row_sums == 0
        if zero_rows.any() and proportions.shape[1] > 0:
            proportions.loc[zero_rows, :] = 1.0 / proportions.shape[1]

    return proportions


def _is_integer(value: Any) -> bool:
    return isinstance(value, Integral) and not isinstance(value, (bool, np.bool_))


def validate_fragment_shape_spec(spec: FragmentShapeSpec) -> None:
    """Validate the versioned parent-fragment-length declaration."""
    if not isinstance(spec, FragmentShapeSpec):
        raise TypeError("fragment_shape must be a FragmentShapeSpec.")
    if not _is_integer(spec.schema_version) or spec.schema_version != 1:
        raise ValueError("fragment_shape.schema_version must be integer 1.")
    if spec.axis != "parent_fragment_length_bp":
        raise ValueError("fragment_shape.axis must be 'parent_fragment_length_bp' for schema version 1.")
    if spec.count_unit != "deduplicated_cut_sites":
        raise ValueError("fragment_shape.count_unit must be 'deduplicated_cut_sites'.")
    if spec.read_support_policy != "ignore":
        raise ValueError("fragment_shape.read_support_policy must be 'ignore'.")
    if spec.peak_assignment != "containing_nonoverlapping_peak":
        raise ValueError(
            "fragment_shape.peak_assignment must be 'containing_nonoverlapping_peak'."
        )
    if spec.bins != _CANONICAL_FRAGMENT_SHAPE_BINS:
        raise ValueError(
            "fragment_shape.bins must be the ordered schema-version-1 bins "
            "[0, 100), [100, 250), and [250, infinity) with canonical layer names."
        )
    if spec.left_cut_offset is not None and (
        not _is_integer(spec.left_cut_offset) or spec.left_cut_offset != 0
    ):
        raise ValueError("fragment_shape.left_cut_offset must be 0 when declared.")
    if spec.right_cut_offset is not None and (
        not _is_integer(spec.right_cut_offset) or spec.right_cut_offset not in (-1, 0)
    ):
        raise ValueError("fragment_shape.right_cut_offset must be the validated value 0 or -1 when declared.")


def _matrix_values(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.data)
    return np.asarray(matrix)


def _validate_count_matrix(matrix: Any, name: str) -> None:
    values = _matrix_values(matrix)
    try:
        finite = np.isfinite(values)
    except TypeError as error:
        raise ValueError(f"{name} must contain numeric counts.") from error
    if not finite.all():
        raise ValueError(f"{name} contains non-finite counts.")
    if (values < 0).any():
        raise ValueError(f"{name} contains negative counts.")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError(f"{name} contains non-integer counts.")


def _matrix_equals_sparse_sum(matrix: Any, expected: sparse.csr_matrix) -> bool:
    expected = expected.copy()
    expected.sum_duplicates()
    expected.eliminate_zeros()

    if sparse.issparse(matrix):
        observed = matrix.tocsr(copy=True)
        observed.sum_duplicates()
        observed.eliminate_zeros()
        difference = observed - expected
        difference.eliminate_zeros()
        return difference.nnz == 0

    observed = np.asarray(matrix)
    if observed.shape != expected.shape or np.count_nonzero(observed) != expected.nnz:
        return False
    rows, columns = expected.nonzero()
    return np.array_equal(observed[rows, columns], expected.data)


def _validate_sha256(value: Any, field_name: str) -> None:
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"fragment_shape.{field_name} must be a 64-character SHA-256 hex digest.")


def _validate_source_sha256(value: Any) -> None:
    if not isinstance(value, Mapping) or not value:
        raise ValueError("fragment_shape.source_sha256 must be a non-empty digest mapping.")
    missing = {"fragments", "tabix_index"}.difference(value)
    if missing:
        raise ValueError(
            "fragment_shape.source_sha256 is missing required sources: "
            f"{', '.join(sorted(missing))}."
        )
    for source_name, digest in value.items():
        if not isinstance(source_name, str) or not source_name:
            raise ValueError("fragment_shape.source_sha256 keys must be non-empty strings.")
        _validate_sha256(digest, f"source_sha256.{source_name}")


def _bin_counter_values(
    counters: Mapping[str, Any],
    layer_names: tuple[str, ...],
    counter_name: str,
) -> dict[str, Any]:
    nested = counters.get("cut_sites_per_bin")
    if nested is not None and not isinstance(nested, Mapping):
        raise ValueError(f"fragment_shape.{counter_name}.cut_sites_per_bin must be a mapping.")
    bin_counts = dict(nested or {})

    prefix = "cut_sites_per_bin."
    flattened = {
        key[len(prefix) :]: value
        for key, value in counters.items()
        if isinstance(key, str) and key.startswith(prefix)
    }
    for layer_name, count in flattened.items():
        if layer_name in bin_counts and bin_counts[layer_name] != count:
            raise ValueError(
                f"Nested and flattened {counter_name} bin counters disagree for "
                f"'{layer_name}'."
            )
        bin_counts[layer_name] = count

    if set(bin_counts) != set(layer_names):
        raise ValueError(
            f"fragment_shape.{counter_name} must provide exactly one cut-site count "
            "for every declared layer."
        )
    return bin_counts


def _validate_preprocessing_counters(
    counters: Any,
    layer_names: tuple[str, ...],
) -> dict[str, int]:
    if not isinstance(counters, Mapping) or not counters:
        raise ValueError("fragment_shape.preprocessing_counters must be a non-empty mapping.")

    missing = _FRAGMENT_SHAPE_COUNTER_FIELDS.difference(counters)
    if missing:
        raise ValueError(
            "fragment_shape.preprocessing_counters is missing required fields: "
            f"{', '.join(sorted(missing))}."
        )
    for field_name in _FRAGMENT_SHAPE_COUNTER_FIELDS:
        count = counters[field_name]
        if not _is_integer(count) or count < 0:
            raise ValueError(
                f"fragment_shape.preprocessing_counters.{field_name} must be a nonnegative integer."
            )

    bin_counts = _bin_counter_values(counters, layer_names, "preprocessing_counters")
    if any(not _is_integer(count) or count < 0 for count in bin_counts.values()):
        raise ValueError("Fragment-shape per-layer counters must be nonnegative integers.")

    total_rows = counters["total_rows"]
    valid_rows = counters["valid_rows"]
    invalid_schema_rows = counters["invalid_schema_rows"]
    invalid_coordinate_rows = counters["invalid_coordinate_rows"]
    retained_fragments = counters["retained_fragments"]
    assigned_cut_sites = counters["assigned_cut_sites"]
    outside_cut_sites = counters["cut_sites_outside_peaks"]
    assigned_fragments = counters["fragments_with_assigned_cut_sites"]

    if valid_rows + invalid_schema_rows + invalid_coordinate_rows != total_rows:
        raise ValueError(
            "Fragment-shape row counters must satisfy valid_rows + invalid_schema_rows + "
            "invalid_coordinate_rows == total_rows."
        )
    if assigned_cut_sites + outside_cut_sites != 2 * retained_fragments:
        raise ValueError(
            "Fragment-shape cut-site counters must satisfy assigned_cut_sites + "
            "cut_sites_outside_peaks == 2 * retained_fragments."
        )
    if sum(bin_counts.values()) != assigned_cut_sites:
        raise ValueError("Fragment-shape per-layer counters must sum to assigned_cut_sites.")
    if retained_fragments > valid_rows:
        raise ValueError("retained_fragments cannot exceed valid_rows.")
    if counters["unknown_barcodes"] > valid_rows or counters["filtered_contigs"] > valid_rows:
        raise ValueError("Unknown-barcode and filtered-contig counters cannot exceed valid_rows.")
    if counters["read_support_total"] < valid_rows:
        raise ValueError("read_support_total must be at least valid_rows.")
    if assigned_fragments > retained_fragments:
        raise ValueError("fragments_with_assigned_cut_sites cannot exceed retained_fragments.")
    if assigned_fragments > assigned_cut_sites or assigned_cut_sites > 2 * assigned_fragments:
        raise ValueError(
            "assigned_cut_sites must be between one and two per fragment with assigned cut sites."
        )

    return {layer_name: int(bin_counts[layer_name]) for layer_name in layer_names}


def _validate_matrix_counters(
    counters: Any,
    layer_names: tuple[str, ...],
    layer_totals: Mapping[str, int],
    x_total: int,
    role: str,
) -> None:
    if not isinstance(counters, Mapping) or not counters:
        raise ValueError(f"{role}.uns['fragment_shape'].matrix_counters must be a non-empty mapping.")
    assigned_cut_sites = counters.get("assigned_cut_sites")
    if not _is_integer(assigned_cut_sites) or assigned_cut_sites <= 0:
        raise ValueError("fragment_shape.matrix_counters.assigned_cut_sites must be a positive integer.")

    bin_counts = _bin_counter_values(counters, layer_names, "matrix_counters")
    if any(not _is_integer(count) or count < 0 for count in bin_counts.values()):
        raise ValueError("Fragment-shape matrix per-layer counters must be nonnegative integers.")
    if sum(bin_counts.values()) != assigned_cut_sites:
        raise ValueError("Fragment-shape matrix per-layer counters must sum to assigned_cut_sites.")
    for layer_name in layer_names:
        if bin_counts[layer_name] != layer_totals[layer_name]:
            raise ValueError(
                "fragment_shape.matrix_counters does not match stored layer "
                f"'{layer_name}'."
            )
    if assigned_cut_sites != x_total:
        raise ValueError("fragment_shape.matrix_counters.assigned_cut_sites does not match stored X.")


def _validate_coordinate_provenance(
    value: Any, role: str, right_cut_offset: int
) -> None:
    if not isinstance(value, Mapping) or not value:
        raise ValueError(
            f"{role}.uns['fragment_shape'].coordinate_validation must be a non-empty mapping."
        )
    selected_offset = value.get("selected_right_cut_offset")
    if not _is_integer(selected_offset) or int(selected_offset) not in (-1, 0):
        raise ValueError(
            "fragment_shape.coordinate_validation.selected_right_cut_offset must be 0 or -1."
        )
    if int(selected_offset) != right_cut_offset:
        raise ValueError(
            "fragment_shape.coordinate_validation.selected_right_cut_offset must "
            "match fragment_shape.right_cut_offset."
        )
    expected_numeric = {
        "mismatched_entries": 0,
        "absolute_error": 0,
    }
    required_fields = set(expected_numeric) | {"matrix_match"}
    missing = required_fields - value.keys()
    if missing:
        raise ValueError(
            "fragment_shape.coordinate_validation is missing required fields: "
            f"{', '.join(sorted(missing))}."
        )
    for field_name, expected_value in expected_numeric.items():
        observed = value[field_name]
        if not _is_integer(observed) or observed != expected_value:
            raise ValueError(
                f"fragment_shape.coordinate_validation.{field_name} must be {expected_value}."
            )
    if value["matrix_match"] != "exact":
        raise ValueError("fragment_shape.coordinate_validation.matrix_match must be 'exact'.")


def validate_fragment_shape_feature_axis(adata: Any, role: str = "AnnData") -> str:
    """Validate and return the digest of an object's current ordered feature axis."""
    metadata = adata.uns.get("fragment_shape")
    if not isinstance(metadata, Mapping):
        raise ValueError(f"{role}.uns is missing valid 'fragment_shape' metadata.")
    stored_hash = metadata.get("feature_sha256")
    _validate_sha256(stored_hash, "feature_sha256")
    observed_hash = ordered_feature_sha256(adata.var_names)
    if stored_hash.lower() != observed_hash:
        raise ValueError(
            f"{role}.uns['fragment_shape'].feature_sha256 does not match the current ordered feature axis."
        )
    return observed_hash


def _validate_provenance(spec: FragmentShapeSpec, role: str) -> None:
    if spec.left_cut_offset is None or spec.right_cut_offset is None:
        raise ValueError(f"{role}.uns['fragment_shape'] must contain resolved cut-site offsets.")
    _validate_source_sha256(spec.source_sha256)
    _validate_sha256(spec.feature_sha256, "feature_sha256")
    _validate_sha256(spec.split_sha256, "split_sha256")

    _validate_coordinate_provenance(
        spec.coordinate_validation, role, int(spec.right_cut_offset)
    )
    if not isinstance(spec.software_versions, Mapping) or not spec.software_versions:
        raise ValueError(f"{role}.uns['fragment_shape'].software_versions must be a non-empty mapping.")
    if any(
        not isinstance(name, str) or not name or not isinstance(version, str) or not version
        for name, version in spec.software_versions.items()
    ):
        raise ValueError("fragment_shape.software_versions must map non-empty names to versions.")

    _validate_preprocessing_counters(spec.preprocessing_counters, spec.layer_names)


def _declared_fields_match(declared: FragmentShapeSpec, stored: FragmentShapeSpec, role: str) -> None:
    required_fields = (
        "schema_version",
        "axis",
        "count_unit",
        "read_support_policy",
        "peak_assignment",
        "bins",
    )
    optional_fields = (
        "left_cut_offset",
        "right_cut_offset",
        "source_sha256",
        "feature_sha256",
        "split_sha256",
        "coordinate_validation",
        "software_versions",
        "preprocessing_counters",
    )
    for field_name in required_fields:
        if getattr(declared, field_name) != getattr(stored, field_name):
            raise ValueError(
                f"{role}.uns['fragment_shape'].{field_name} does not match the dataset declaration."
            )
    for field_name in optional_fields:
        declared_value = getattr(declared, field_name)
        if declared_value is not None and declared_value != getattr(stored, field_name):
            raise ValueError(
                f"{role}.uns['fragment_shape'].{field_name} does not match the dataset declaration."
            )


def _validate_shape_object(
    adata,
    declared: FragmentShapeSpec,
    role: str,
) -> FragmentShapeSpec:
    if "fragment_shape" not in adata.uns:
        raise ValueError(f"{role}.uns is missing 'fragment_shape' metadata.")
    try:
        stored = FragmentShapeSpec.from_mapping(adata.uns["fragment_shape"])
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid {role}.uns['fragment_shape'] metadata: {error}") from error

    validate_fragment_shape_spec(stored)
    _declared_fields_match(declared, stored, role)
    _validate_provenance(stored, role)
    validate_fragment_shape_feature_axis(adata, role)

    if adata.X is None:
        raise ValueError(f"{role}.X is missing.")
    _validate_count_matrix(adata.X, f"{role}.X")

    layer_sum = sparse.csr_matrix(adata.shape, dtype=np.int64)
    layer_totals: dict[str, int] = {}
    for layer_name in declared.layer_names:
        if layer_name not in adata.layers:
            raise ValueError(f"{role}.layers is missing declared layer '{layer_name}'.")
        layer = adata.layers[layer_name]
        if not sparse.isspmatrix_csr(layer):
            raise ValueError(f"{role}.layers['{layer_name}'] must be a scipy CSR matrix.")
        if layer.shape != adata.shape:
            raise ValueError(f"{role}.layers['{layer_name}'] does not match the AnnData axes.")
        _validate_count_matrix(layer, f"{role}.layers['{layer_name}']")
        layer_sum = layer_sum + layer
        layer_totals[layer_name] = int(layer.sum())

    if not _matrix_equals_sparse_sum(adata.X, layer_sum):
        raise ValueError(f"{role}.X must equal the exact elementwise sum of fragment-shape layers.")
    x_total = int(adata.X.sum()) if sparse.issparse(adata.X) else int(np.asarray(adata.X).sum())
    _validate_matrix_counters(
        stored.matrix_counters,
        stored.layer_names,
        layer_totals,
        x_total,
        role,
    )
    return stored


def _validate_declared_cell_types(data: DeconvolutionInput) -> None:
    if data.cell_types is None:
        return
    if not isinstance(data.cell_types, list):
        raise TypeError("cell_types must be an ordered list when declared.")
    if not data.cell_types or any(not isinstance(value, str) or not value for value in data.cell_types):
        raise ValueError("cell_types must contain non-empty strings.")
    if len(set(data.cell_types)) != len(data.cell_types):
        raise ValueError("cell_types must not contain duplicates.")

    labels = data.reference.obs[data.labels_key]
    if labels.isna().any() or labels.astype(str).str.len().eq(0).any():
        raise ValueError(f"reference.obs['{data.labels_key}'] contains missing or empty labels.")
    observed_cell_types = set(pd.Index(labels).drop_duplicates())
    if observed_cell_types != set(data.cell_types):
        raise ValueError("cell_types must exactly match the reference label universe.")

    if data.truth is not None:
        if list(data.truth.columns) != data.cell_types:
            raise ValueError("truth columns must exactly match the declared cell_types order.")
        if data.truth.index.has_duplicates:
            raise ValueError("truth observation names must be unique.")
        if not data.truth.index.equals(data.spatial.obs_names):
            raise ValueError("truth rows must be aligned to spatial observations in the same order.")
        if data.truth.isna().any().any():
            raise ValueError("truth contains missing values.")


def validate_fragment_shape_input(data: DeconvolutionInput) -> None:
    """Validate opt-in ShapeMix layers, axes, metadata, and provenance."""
    if data.fragment_shape is None:
        return
    if data.modality != "atac":
        raise ValueError("fragment_shape metadata is supported only for ATAC inputs.")
    if data.cell_types is None:
        raise ValueError("cell_types must be declared for fragment-shape inputs.")

    validate_fragment_shape_spec(data.fragment_shape)
    reference_spec = _validate_shape_object(data.reference, data.fragment_shape, "reference")
    spatial_spec = _validate_shape_object(data.spatial, data.fragment_shape, "spatial")

    aligned_fields = (
        "schema_version",
        "axis",
        "count_unit",
        "read_support_policy",
        "peak_assignment",
        "bins",
        "left_cut_offset",
        "right_cut_offset",
        "feature_sha256",
        "coordinate_validation",
    )
    for field_name in aligned_fields:
        if getattr(reference_spec, field_name) != getattr(spatial_spec, field_name):
            raise ValueError(
                "reference and spatial fragment_shape metadata disagree for "
                f"'{field_name}'."
            )

    labels = data.reference.obs[data.labels_key]
    if labels.isna().any() or labels.astype(str).str.len().eq(0).any():
        raise ValueError(f"reference.obs['{data.labels_key}'] contains missing or empty labels.")


def validate_deconvolution_input(data: DeconvolutionInput) -> None:
    """Validate the shared input contract for deconvolution methods."""
    if data.modality not in {"atac", "rna"}:
        raise ValueError("modality must be 'atac' or 'rna'.")
    if data.reference.n_vars == 0:
        raise ValueError("reference contains no features.")
    if data.spatial.n_vars == 0:
        raise ValueError("spatial contains no features.")
    if list(data.reference.var_names) != list(data.spatial.var_names):
        raise ValueError("reference and spatial features must be aligned in the same order.")
    if data.labels_key not in data.reference.obs:
        raise ValueError(f"labels_key '{data.labels_key}' is missing from reference.obs.")
    if data.spatial_key not in data.spatial.obsm:
        raise ValueError(f"spatial coordinates '{data.spatial_key}' are missing from spatial.obsm.")

    if data.truth is not None:
        missing = data.spatial.obs_names.difference(data.truth.index)
        if len(missing) > 0:
            raise ValueError("truth is missing rows for one or more spatial observations.")

    _validate_declared_cell_types(data)
    if data.fragment_shape is not None:
        validate_fragment_shape_input(data)
