"""Compact, JSON-safe diagnostics for ShapeMix MAP fits."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Optional, Sequence, Tuple

import torch
from torch import Tensor


def _json_float(value: Optional[float]) -> Optional[float]:
    """Return a finite built-in float or ``None`` for a non-finite value."""
    if value is None:
        return None
    normalized = float(value)
    return normalized if math.isfinite(normalized) else None


@dataclass(frozen=True)
class OptimizationStepRecord:
    """Immutable objective snapshot for one completed optimizer step."""

    step: int
    count_log_likelihood: Optional[float]
    shape_log_likelihood: Optional[float]
    abundance_log_prior: Optional[float]
    total_log_objective: Optional[float]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe training-history row."""
        return {
            "step": int(self.step),
            "count_log_likelihood": _json_float(self.count_log_likelihood),
            "shape_log_likelihood": _json_float(self.shape_log_likelihood),
            "abundance_log_prior": _json_float(self.abundance_log_prior),
            "total_log_objective": _json_float(self.total_log_objective),
        }


@dataclass(frozen=True)
class RestartRecord:
    """Immutable final record for one deterministic optimization restart."""

    restart_index: int
    seed_tuple: Tuple[int, ...]
    converged: bool
    stopping_reason: str
    steps: int
    count_log_likelihood: Optional[float]
    shape_log_likelihood: Optional[float]
    abundance_log_prior: Optional[float]
    total_log_objective: Optional[float]
    history: Tuple[OptimizationStepRecord, ...] = ()
    nonfinite_events: Tuple[str, ...] = ()
    runtime_seconds: Optional[float] = None
    best_step: Optional[int] = None
    final_gradient_norm: Optional[float] = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe mapping without tensors or non-finite numbers."""
        return {
            "restart_index": int(self.restart_index),
            "seed_tuple": [int(value) for value in self.seed_tuple],
            "converged": bool(self.converged),
            "stopping_reason": str(self.stopping_reason),
            "steps": int(self.steps),
            "count_log_likelihood": _json_float(self.count_log_likelihood),
            "shape_log_likelihood": _json_float(self.shape_log_likelihood),
            "abundance_log_prior": _json_float(self.abundance_log_prior),
            "total_log_objective": _json_float(self.total_log_objective),
            "history": [record.to_dict() for record in self.history],
            "nonfinite_events": [str(value) for value in self.nonfinite_events],
            "runtime_seconds": _json_float(self.runtime_seconds),
            "best_step": None if self.best_step is None else int(self.best_step),
            "final_gradient_norm": _json_float(self.final_gradient_norm),
        }


@dataclass(frozen=True)
class FitRecord:
    """Immutable summary of restart selection and the complete MAP fit."""

    use_shape: bool
    success: bool
    selected_restart: Optional[int]
    stopping_reason: str
    n_spots: int
    n_peaks: int
    n_bins: int
    n_cell_types: int
    dtype: str
    device: str
    restarts: Tuple[RestartRecord, ...]
    nonfinite_events: Tuple[str, ...] = ()
    runtime_seconds: Optional[float] = None
    count_cache_mode: str = "streamed_host_chunks"
    count_cache_bytes: int = 0
    count_cache_requested_bytes: int = 0
    device_total_memory_bytes: Optional[int] = None
    device_free_memory_bytes_before_cache: Optional[int] = None
    cuda_workspace_margin_bytes: Optional[int] = None
    device_index: Optional[int] = None
    device_name: Optional[str] = None
    device_compute_capability: Optional[str] = None
    torch_version: Optional[str] = None
    cuda_runtime_version: Optional[str] = None
    deterministic_algorithms: Optional[bool] = None
    float32_matmul_precision: Optional[str] = None
    torch_num_threads: Optional[int] = None
    torch_num_interop_threads: Optional[int] = None
    peak_device_memory_allocated_bytes: Optional[int] = None
    peak_device_memory_reserved_bytes: Optional[int] = None

    @property
    def selected_restart_record(self) -> Optional[RestartRecord]:
        """Return the selected restart, or ``None`` for a failed fit."""
        if self.selected_restart is None:
            return None
        matches = [
            record
            for record in self.restarts
            if record.restart_index == self.selected_restart
        ]
        if len(matches) != 1:
            raise ValueError(
                "selected_restart must identify exactly one restart record."
            )
        return matches[0]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe mapping suitable for ``diagnostics.json``."""
        return {
            "use_shape": bool(self.use_shape),
            "success": bool(self.success),
            "selected_restart": (
                None
                if self.selected_restart is None
                else int(self.selected_restart)
            ),
            "stopping_reason": str(self.stopping_reason),
            "dimensions": {
                "spots": int(self.n_spots),
                "peaks": int(self.n_peaks),
                "bins": int(self.n_bins),
                "cell_types": int(self.n_cell_types),
            },
            "dtype": str(self.dtype),
            "device": str(self.device),
            "restarts": [record.to_dict() for record in self.restarts],
            "nonfinite_events": [str(value) for value in self.nonfinite_events],
            "runtime_seconds": _json_float(self.runtime_seconds),
            "execution": {
                "count_cache_mode": str(self.count_cache_mode),
                "count_cache_bytes": int(self.count_cache_bytes),
                "count_cache_requested_bytes": int(
                    self.count_cache_requested_bytes
                ),
                "device_total_memory_bytes": self.device_total_memory_bytes,
                "device_free_memory_bytes_before_cache": (
                    self.device_free_memory_bytes_before_cache
                ),
                "cuda_workspace_margin_bytes": self.cuda_workspace_margin_bytes,
                "device_index": self.device_index,
                "device_name": self.device_name,
                "device_compute_capability": self.device_compute_capability,
                "torch_version": self.torch_version,
                "cuda_runtime_version": self.cuda_runtime_version,
                "deterministic_algorithms": self.deterministic_algorithms,
                "float32_matmul_precision": self.float32_matmul_precision,
                "torch_num_threads": self.torch_num_threads,
                "torch_num_interop_threads": self.torch_num_interop_threads,
                "peak_device_memory_allocated_bytes": (
                    self.peak_device_memory_allocated_bytes
                ),
                "peak_device_memory_reserved_bytes": (
                    self.peak_device_memory_reserved_bytes
                ),
            },
        }


# A descriptive alias for callers that prefer the diagnostics-oriented name.
FitDiagnostics = FitRecord


@dataclass(frozen=True)
class ReconstructionSummary:
    """Immutable count/bin reconstruction summary accumulated over chunks."""

    count_entries: int
    observed_nonzero_entries: int
    expected_positive_entries: int
    observed_total: float
    expected_total: float
    residual_total: float
    absolute_error_sum: float
    squared_error_sum: float
    observed_by_bin: Tuple[float, ...]
    expected_by_bin: Tuple[float, ...]
    residual_by_bin: Tuple[float, ...]
    absolute_error_by_bin: Tuple[float, ...]
    squared_error_by_bin: Tuple[float, ...]

    @property
    def mean_absolute_error(self) -> float:
        """Return mean absolute total-count reconstruction error."""
        if self.count_entries == 0:
            return 0.0
        return self.absolute_error_sum / self.count_entries

    @property
    def root_mean_squared_error(self) -> float:
        """Return root mean squared total-count reconstruction error."""
        if self.count_entries == 0:
            return 0.0
        return math.sqrt(self.squared_error_sum / self.count_entries)

    @property
    def observed_nonzero_fraction(self) -> float:
        """Return the fraction of observed spot/peak entries with nonzero counts."""
        if self.count_entries == 0:
            return 0.0
        return self.observed_nonzero_entries / self.count_entries

    @property
    def expected_positive_fraction(self) -> float:
        """Return the fraction of fitted spot/peak rates above zero."""
        if self.count_entries == 0:
            return 0.0
        return self.expected_positive_entries / self.count_entries

    def to_dict(self) -> dict[str, Any]:
        """Return only finite built-ins and lists for strict JSON writers."""
        count = {
            "entries": int(self.count_entries),
            "observed_nonzero_entries": int(self.observed_nonzero_entries),
            "observed_nonzero_fraction": float(self.observed_nonzero_fraction),
            "expected_positive_entries": int(self.expected_positive_entries),
            "expected_positive_fraction": float(self.expected_positive_fraction),
            "observed_total": float(self.observed_total),
            "expected_total": float(self.expected_total),
            "residual_total": float(self.residual_total),
            "absolute_error_sum": float(self.absolute_error_sum),
            "squared_error_sum": float(self.squared_error_sum),
            "mean_absolute_error": float(self.mean_absolute_error),
            "root_mean_squared_error": float(self.root_mean_squared_error),
        }
        shape = None
        if self.observed_by_bin:
            shape = {
                "observed_by_bin": [float(value) for value in self.observed_by_bin],
                "expected_by_bin": [float(value) for value in self.expected_by_bin],
                "residual_by_bin": [float(value) for value in self.residual_by_bin],
                "absolute_error_by_bin": [
                    float(value) for value in self.absolute_error_by_bin
                ],
                "squared_error_by_bin": [
                    float(value) for value in self.squared_error_by_bin
                ],
            }
        return {"count": count, "shape": shape}


class ReconstructionAccumulator:
    """Accumulate reconstruction diagnostics without retaining dense chunks.

    Call :meth:`update` once per non-overlapping spot/peak chunk.  Optional bin
    tensors are immediately reduced over every axis except the final bin axis,
    so memory usage remains ``O(B)`` after each call.
    """

    def __init__(self, n_bins: Optional[int] = None) -> None:
        if n_bins is not None:
            if isinstance(n_bins, bool) or not isinstance(n_bins, int):
                raise TypeError("n_bins must be an integer or None.")
            if n_bins < 1:
                raise ValueError("n_bins must be strictly positive.")
        self._n_bins = n_bins
        self._count_entries = 0
        self._observed_nonzero_entries = 0
        self._expected_positive_entries = 0
        self._observed_total = 0.0
        self._expected_total = 0.0
        self._absolute_error_sum = 0.0
        self._squared_error_sum = 0.0
        self._observed_by_bin = [0.0] * (n_bins or 0)
        self._expected_by_bin = [0.0] * (n_bins or 0)
        self._absolute_error_by_bin = [0.0] * (n_bins or 0)
        self._squared_error_by_bin = [0.0] * (n_bins or 0)

    @staticmethod
    def _as_finite_double(value: Tensor, name: str) -> Tensor:
        if not isinstance(value, Tensor):
            value = torch.as_tensor(value)
        detached = value.detach().to(device="cpu", dtype=torch.float64)
        if not bool(torch.isfinite(detached).all().item()):
            raise ValueError(f"{name} contains non-finite values.")
        return detached

    def update(
        self,
        observed_totals: Tensor,
        expected_totals: Tensor,
        *,
        observed_bins: Optional[Tensor] = None,
        expected_bins: Optional[Tensor] = None,
    ) -> None:
        """Reduce one chunk into the running count and optional bin summaries."""
        observed = self._as_finite_double(observed_totals, "observed_totals")
        expected = self._as_finite_double(expected_totals, "expected_totals")
        if observed.shape != expected.shape:
            raise ValueError("observed_totals and expected_totals must have equal shapes.")
        if observed.ndim != 2:
            raise ValueError("Total-count chunks must have spot-by-peak shape.")
        if bool((observed < 0).any().item()) or bool((expected < 0).any().item()):
            raise ValueError("Observed counts and expected rates must be nonnegative.")

        residual = observed - expected
        if (observed_bins is None) != (expected_bins is None):
            raise ValueError("observed_bins and expected_bins must be supplied together.")
        bins = None
        observed_sums = expected_sums = absolute_sums = squared_sums = None
        if observed_bins is not None:
            observed_shape = self._as_finite_double(observed_bins, "observed_bins")
            expected_shape = self._as_finite_double(expected_bins, "expected_bins")
            if observed_shape.shape != expected_shape.shape:
                raise ValueError("observed_bins and expected_bins must have equal shapes.")
            if observed_shape.ndim != observed.ndim + 1:
                raise ValueError("Bin tensors require exactly one additional final axis.")
            if observed_shape.shape[:-1] != observed.shape:
                raise ValueError("Bin tensor leading dimensions must match total tensors.")
            if bool((observed_shape < 0).any().item()) or bool(
                (expected_shape < 0).any().item()
            ):
                raise ValueError(
                    "Observed bin counts and expected bin rates must be nonnegative."
                )

            bins = observed_shape.shape[-1]
            if self._n_bins is not None and bins != self._n_bins:
                raise ValueError(
                    f"Expected {self._n_bins} bins, observed a chunk with {bins}."
                )
            if not torch.equal(observed_shape.sum(dim=-1), observed):
                raise ValueError("Observed bin counts must sum exactly to observed totals.")
            if not torch.allclose(
                expected_shape.sum(dim=-1), expected, rtol=1.0e-5, atol=1.0e-6
            ):
                raise ValueError("Expected bin rates must sum to expected total rates.")

            reduction_axes: Sequence[int] = tuple(range(observed_shape.ndim - 1))
            shape_residual = observed_shape - expected_shape
            observed_sums = observed_shape.sum(dim=reduction_axes).tolist()
            expected_sums = expected_shape.sum(dim=reduction_axes).tolist()
            absolute_sums = shape_residual.abs().sum(dim=reduction_axes).tolist()
            squared_sums = shape_residual.square().sum(dim=reduction_axes).tolist()

        # Mutate the accumulator only after the complete chunk is validated.
        self._count_entries += observed.numel()
        self._observed_nonzero_entries += int((observed > 0).sum().item())
        self._expected_positive_entries += int((expected > 0).sum().item())
        self._observed_total += float(observed.sum().item())
        self._expected_total += float(expected.sum().item())
        self._absolute_error_sum += float(residual.abs().sum().item())
        self._squared_error_sum += float(residual.square().sum().item())
        if bins is None:
            return
        if self._n_bins is None:
            self._n_bins = bins
            self._observed_by_bin = [0.0] * bins
            self._expected_by_bin = [0.0] * bins
            self._absolute_error_by_bin = [0.0] * bins
            self._squared_error_by_bin = [0.0] * bins
        assert observed_sums is not None
        assert expected_sums is not None
        assert absolute_sums is not None
        assert squared_sums is not None
        for bin_index in range(bins):
            self._observed_by_bin[bin_index] += float(observed_sums[bin_index])
            self._expected_by_bin[bin_index] += float(expected_sums[bin_index])
            self._absolute_error_by_bin[bin_index] += float(absolute_sums[bin_index])
            self._squared_error_by_bin[bin_index] += float(squared_sums[bin_index])

    def finalize(self) -> ReconstructionSummary:
        """Freeze the accumulated values into a JSON-safe summary record."""
        residual_by_bin = tuple(
            observed - expected
            for observed, expected in zip(
                self._observed_by_bin, self._expected_by_bin
            )
        )
        return ReconstructionSummary(
            count_entries=self._count_entries,
            observed_nonzero_entries=self._observed_nonzero_entries,
            expected_positive_entries=self._expected_positive_entries,
            observed_total=self._observed_total,
            expected_total=self._expected_total,
            residual_total=self._observed_total - self._expected_total,
            absolute_error_sum=self._absolute_error_sum,
            squared_error_sum=self._squared_error_sum,
            observed_by_bin=tuple(self._observed_by_bin),
            expected_by_bin=tuple(self._expected_by_bin),
            residual_by_bin=residual_by_bin,
            absolute_error_by_bin=tuple(self._absolute_error_by_bin),
            squared_error_by_bin=tuple(self._squared_error_by_bin),
        )
