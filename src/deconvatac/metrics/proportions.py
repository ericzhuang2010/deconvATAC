"""Strict, versioned evaluation of cell-type proportion estimates.

The ShapeMix benchmark does not use the intersection-and-renormalize behavior
of the historical metric helpers. The declared cell-type universe and complete
truth spot set are part of the endpoint definition, so invalid or incomplete
output fails evaluation instead of disappearing from the denominator.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Callable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype
from scipy.spatial.distance import jensenshannon


PROPORTION_CONTRACT_VERSION = "v2"
PROPORTION_ROW_SUM_ATOL = 1.0e-6


@dataclass(frozen=True)
class ProportionMetric:
    """A registered scalar endpoint on already aligned proportion arrays."""

    metric_id: str
    name: str
    version: str
    evaluator: Callable[[np.ndarray, np.ndarray], float]
    lower_is_better: bool = True


@dataclass(frozen=True)
class MetricEvaluation:
    """A scalar value together with the contract information used to score it."""

    metric_id: str
    metric_name: str
    metric_version: str
    value: float
    contract_version: str
    row_sum_atol: float
    cell_types: tuple[str, ...]
    n_spots: int
    n_cell_types: int


def _rmse(true: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(true - predicted))))


def _jsd_v2(true: np.ndarray, predicted: np.ndarray) -> float:
    distances = jensenshannon(true, predicted, axis=1, base=2)
    return float(np.mean(np.square(distances)))


def _js_distance_v1(true: np.ndarray, predicted: np.ndarray) -> float:
    distances = jensenshannon(true, predicted, axis=1, base=2)
    return float(np.mean(distances))


_PROPORTION_METRICS = {
    "rmse_v1": ProportionMetric("rmse_v1", "rmse", "v1", _rmse),
    "jsd_v2": ProportionMetric("jsd_v2", "jsd", "v2", _jsd_v2),
    "js_distance_v1": ProportionMetric(
        "js_distance_v1",
        "js_distance",
        "v1",
        _js_distance_v1,
    ),
}

# These names occur in pre-Step-6 experiment YAML. They are accepted only at
# the registry boundary; output always records the canonical versioned ID.
_METRIC_ALIASES = {
    "rmse": "rmse_v1",
    "jsd": "jsd_v2",
}

PROPORTION_METRICS: Mapping[str, ProportionMetric] = MappingProxyType(_PROPORTION_METRICS)
PROPORTION_METRIC_ALIASES: Mapping[str, str] = MappingProxyType(_METRIC_ALIASES)


def available_proportion_metrics(*, include_aliases: bool = True) -> tuple[str, ...]:
    """Return accepted metric selectors in deterministic order."""

    names = set(PROPORTION_METRICS)
    if include_aliases:
        names.update(PROPORTION_METRIC_ALIASES)
    return tuple(sorted(names))


def resolve_proportion_metric(name: str) -> ProportionMetric:
    """Resolve an endpoint selector to its canonical registered definition."""

    if not isinstance(name, str) or not name:
        raise KeyError("Metric name must be a nonempty string.")
    metric_id = PROPORTION_METRIC_ALIASES.get(name, name)
    try:
        return PROPORTION_METRICS[metric_id]
    except KeyError as exc:
        available = ", ".join(available_proportion_metrics())
        raise KeyError(f"Unknown metric '{name}'. Available: {available}") from exc


def validate_cell_type_universe(cell_types: Sequence[str]) -> tuple[str, ...]:
    """Validate and freeze the dataset-declared ordered cell-type universe."""

    if isinstance(cell_types, (str, bytes)) or not isinstance(cell_types, Sequence):
        raise TypeError("cell_types must be a declared ordered sequence of strings.")
    universe = tuple(cell_types)
    if not universe:
        raise ValueError("cell_types must contain at least one declared cell type.")
    if any(not isinstance(cell_type, str) or not cell_type for cell_type in universe):
        raise ValueError("cell_types must contain nonempty strings.")
    if len(set(universe)) != len(universe):
        raise ValueError("cell_types must not contain duplicates.")
    return universe


def _validate_axis(frame: pd.DataFrame, label: str) -> None:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"{label} must be a pandas DataFrame.")
    if frame.index.has_duplicates:
        raise ValueError(f"{label} spot names must be unique.")
    if frame.columns.has_duplicates:
        raise ValueError(f"{label} cell-type columns must be unique.")
    if frame.index.hasnans:
        raise ValueError(f"{label} spot names must not contain missing values.")
    if frame.columns.hasnans:
        raise ValueError(f"{label} cell-type columns must not contain missing values.")


def _format_labels(values: pd.Index) -> str:
    return ", ".join(repr(value) for value in values.tolist()) or "none"


def _numeric_values(
    frame: pd.DataFrame,
    label: str,
    *,
    row_sum_atol: float,
) -> np.ndarray:
    nonnumeric = [
        str(column)
        for column, dtype in frame.dtypes.items()
        if not is_numeric_dtype(dtype) or is_bool_dtype(dtype)
    ]
    if nonnumeric:
        raise TypeError(f"{label} must be numeric; nonnumeric columns: {', '.join(nonnumeric)}.")

    values = frame.to_numpy(dtype=np.float64, copy=True)
    if values.size == 0 or values.shape[0] == 0:
        raise ValueError(f"{label} must contain at least one spot.")
    if not np.isfinite(values).all():
        raise ValueError(f"{label} values must all be finite.")
    if np.any(values < 0):
        raise ValueError(f"{label} values must all be nonnegative.")

    row_sums = values.sum(axis=1)
    zero_rows = np.flatnonzero(row_sums == 0)
    if zero_rows.size:
        spots = frame.index[zero_rows]
        raise ValueError(f"{label} contains all-zero rows: {_format_labels(spots)}.")
    invalid_sums = ~np.isclose(row_sums, 1.0, rtol=0.0, atol=row_sum_atol)
    if np.any(invalid_sums):
        spots = frame.index[np.flatnonzero(invalid_sums)]
        raise ValueError(
            f"{label} rows must sum to one within atol={row_sum_atol:g}; "
            f"invalid spots: {_format_labels(spots)}."
        )
    return values


def align_proportions(
    true: pd.DataFrame,
    predicted: pd.DataFrame,
    cell_types: Sequence[str],
    *,
    row_sum_atol: float = PROPORTION_ROW_SUM_ATOL,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate and align matrices under the frozen proportion contract.

    Truth columns must exactly equal ``cell_types`` in order. Prediction may
    omit a declared type (it is inserted as zero and therefore penalized), but
    it may not introduce an undeclared type. Prediction spots must be the exact
    truth spot set; only their order may differ.
    """

    if not isinstance(row_sum_atol, (int, float)) or not np.isfinite(row_sum_atol):
        raise TypeError("row_sum_atol must be a finite nonnegative number.")
    if row_sum_atol < 0:
        raise ValueError("row_sum_atol must be nonnegative.")

    universe = validate_cell_type_universe(cell_types)
    _validate_axis(true, "truth")
    _validate_axis(predicted, "prediction")

    if tuple(true.columns) != universe:
        raise ValueError(
            "Truth columns must exactly match the declared cell_types order: "
            f"expected {list(universe)!r}, observed {list(true.columns)!r}."
        )

    missing_spots = true.index.difference(predicted.index, sort=False)
    extra_spots = predicted.index.difference(true.index, sort=False)
    if len(missing_spots) or len(extra_spots):
        raise ValueError(
            "Prediction must contain exactly the truth spot set; "
            f"missing: {_format_labels(missing_spots)}; extra: {_format_labels(extra_spots)}."
        )

    unknown_types = predicted.columns.difference(pd.Index(universe), sort=False)
    if len(unknown_types):
        raise ValueError(
            "Prediction contains cell types outside the declared universe: "
            f"{_format_labels(unknown_types)}."
        )

    true_aligned = true.copy()
    predicted_aligned = predicted.reindex(index=true.index, columns=universe, fill_value=0.0)
    _numeric_values(true_aligned, "truth", row_sum_atol=float(row_sum_atol))
    _numeric_values(predicted_aligned, "prediction", row_sum_atol=float(row_sum_atol))
    return true_aligned, predicted_aligned


def evaluate_proportion_metric(
    metric_name: str,
    true: pd.DataFrame,
    predicted: pd.DataFrame,
    cell_types: Sequence[str],
    *,
    row_sum_atol: float = PROPORTION_ROW_SUM_ATOL,
) -> MetricEvaluation:
    """Evaluate one registered endpoint and return its full version metadata."""

    metric = resolve_proportion_metric(metric_name)
    true_aligned, predicted_aligned = align_proportions(
        true,
        predicted,
        cell_types,
        row_sum_atol=row_sum_atol,
    )
    true_values = true_aligned.to_numpy(dtype=np.float64, copy=False)
    predicted_values = predicted_aligned.to_numpy(dtype=np.float64, copy=False)
    value = float(metric.evaluator(true_values, predicted_values))
    if not np.isfinite(value):
        raise ValueError(f"Metric '{metric.metric_id}' produced a non-finite value.")
    if metric.metric_id == "jsd_v2" and not (0.0 <= value <= 1.0 + 1.0e-12):
        raise ValueError(f"Metric 'jsd_v2' produced an out-of-range value: {value}.")

    universe = validate_cell_type_universe(cell_types)
    return MetricEvaluation(
        metric_id=metric.metric_id,
        metric_name=metric.name,
        metric_version=metric.version,
        value=value,
        contract_version=PROPORTION_CONTRACT_VERSION,
        row_sum_atol=float(row_sum_atol),
        cell_types=universe,
        n_spots=len(true_aligned),
        n_cell_types=len(universe),
    )


def rmse_v1(
    true: pd.DataFrame,
    predicted: pd.DataFrame,
    cell_types: Sequence[str],
    *,
    row_sum_atol: float = PROPORTION_ROW_SUM_ATOL,
) -> float:
    """Return RMSE over every spot by declared-cell-type entry."""

    return evaluate_proportion_metric(
        "rmse_v1", true, predicted, cell_types, row_sum_atol=row_sum_atol
    ).value


def jsd_v2(
    true: pd.DataFrame,
    predicted: pd.DataFrame,
    cell_types: Sequence[str],
    *,
    row_sum_atol: float = PROPORTION_ROW_SUM_ATOL,
) -> float:
    """Return mean per-spot base-2 Jensen--Shannon divergence."""

    return evaluate_proportion_metric(
        "jsd_v2", true, predicted, cell_types, row_sum_atol=row_sum_atol
    ).value


def js_distance_v1(
    true: pd.DataFrame,
    predicted: pd.DataFrame,
    cell_types: Sequence[str],
    *,
    row_sum_atol: float = PROPORTION_ROW_SUM_ATOL,
) -> float:
    """Return the historical unsquared distance under the strict input contract."""

    return evaluate_proportion_metric(
        "js_distance_v1", true, predicted, cell_types, row_sum_atol=row_sum_atol
    ).value


# Source-level compatibility aliases. Unlike the pre-Step-6 implementations,
# these still require the declared universe and strict validation.
rmse = rmse_v1
jsd = jsd_v2


def declared_cell_types_from_metadata(
    run_metadata: Optional[Mapping[str, object]] = None,
    inputs_metadata: Optional[Mapping[str, object]] = None,
    override: Optional[Sequence[str]] = None,
) -> tuple[str, ...]:
    """Resolve, but never infer, the declared universe for a run.

    New runs record the universe directly in ``run.yaml``. The ``inputs.yaml``
    fallback makes pre-Step-6 ShapeMix runs evaluable when their captured
    dataset config already declared ``truth.cell_types``.
    """

    if override is not None:
        return validate_cell_type_universe(override)

    run_metadata = run_metadata or {}
    contract = run_metadata.get("proportion_evaluation")
    if isinstance(contract, Mapping) and contract.get("cell_types") is not None:
        return validate_cell_type_universe(contract["cell_types"])  # type: ignore[arg-type]
    if run_metadata.get("declared_cell_types") is not None:
        return validate_cell_type_universe(run_metadata["declared_cell_types"])  # type: ignore[arg-type]

    inputs_metadata = inputs_metadata or {}
    dataset_config = inputs_metadata.get("dataset_config")
    if isinstance(dataset_config, Mapping):
        modality = run_metadata.get("modality") or inputs_metadata.get("modality")
        modalities = dataset_config.get("modalities")
        modality_config = modalities.get(modality) if isinstance(modalities, Mapping) else None
        truth_spec = modality_config.get("truth") if isinstance(modality_config, Mapping) else None
        if not isinstance(truth_spec, Mapping):
            truth_spec = dataset_config.get("truth")
        if isinstance(truth_spec, Mapping) and truth_spec.get("cell_types") is not None:
            return validate_cell_type_universe(truth_spec["cell_types"])  # type: ignore[arg-type]

    raise ValueError(
        "No declared cell-type universe was found. Declare truth.cell_types in the dataset "
        "YAML, record it in run.yaml, or pass an explicit cell-type list."
    )
