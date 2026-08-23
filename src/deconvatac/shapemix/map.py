"""Deterministic, streaming MAP fitting for fixed ShapeMix signatures."""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Any, Mapping, Optional, Sequence, Union

import numpy as np
import torch
from scipy import sparse
from scipy.optimize import nnls

from .config import ShapeMixConfig
from .diagnostics import FitRecord, OptimizationStepRecord, RestartRecord
from .likelihood import gamma_log_prior, likelihood_components, positive_abundance


SEED_NAMESPACE = 20260822
ABUNDANCE_EPSILON = 1.0e-8
INITIAL_ABUNDANCE_FLOOR = 5.0e-2
RESTART_LOG_STANDARD_DEVIATION = 0.20


@dataclass(frozen=True)
class _ObjectiveValues:
    """Complete finite MAP objective split into its three terms."""

    count_log_likelihood: float
    shape_log_likelihood: float
    abundance_log_prior: float

    @property
    def total_log_objective(self) -> float:
        return (
            self.count_log_likelihood
            + self.shape_log_likelihood
            + self.abundance_log_prior
        )

    @property
    def finite(self) -> bool:
        return bool(
            np.isfinite(
                [
                    self.count_log_likelihood,
                    self.shape_log_likelihood,
                    self.abundance_log_prior,
                    self.total_log_objective,
                ]
            ).all()
        )


@dataclass(frozen=True)
class MAPFitResult:
    """Absolute effective abundance, proportions, and restart diagnostics."""

    abundance: np.ndarray
    proportions: np.ndarray
    selected_restart: int
    count_log_likelihood: float
    shape_log_likelihood: float
    abundance_log_prior: float
    diagnostics: FitRecord
    nnls_fallback_spots: tuple[str, ...]
    spot_names: tuple[str, ...]
    cell_types: tuple[str, ...]
    use_shape: bool

    def __post_init__(self) -> None:
        for name in ("abundance", "proportions"):
            array = np.ascontiguousarray(getattr(self, name), dtype="<f8")
            immutable = np.frombuffer(array.tobytes(order="C"), dtype="<f8").reshape(
                array.shape
            )
            object.__setattr__(self, name, immutable)

    @property
    def total_log_objective(self) -> float:
        """Return the selected restart's complete MAP log objective."""
        return (
            self.count_log_likelihood
            + self.shape_log_likelihood
            + self.abundance_log_prior
        )

    @property
    def restart_diagnostics(self) -> tuple[RestartRecord, ...]:
        """Expose the shared restart records for convenient callers."""
        return self.diagnostics.restarts

    def to_diagnostics_dict(self) -> dict[str, Any]:
        """Return compact run diagnostics without duplicating fitted matrices."""
        return {
            **self.diagnostics.to_dict(),
            "objective": {
                "count_log_likelihood": self.count_log_likelihood,
                "shape_log_likelihood": self.shape_log_likelihood,
                "abundance_log_prior": self.abundance_log_prior,
                "total_log_objective": self.total_log_objective,
            },
            "nnls_fallback_spots": list(self.nnls_fallback_spots),
        }


class MAPOptimizationError(RuntimeError):
    """Raised when no finite, converged restart can be selected."""

    def __init__(
        self,
        message: str,
        diagnostics: FitRecord,
    ) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics


class _CountSource:
    """Chunked access to dense counts or ordered sparse bin layers."""

    def __init__(self, values: Any) -> None:
        self._dense: Optional[np.ndarray] = None
        self._layers: Optional[tuple[Any, ...]] = None

        dense = self._as_dense_tensor(values)
        if dense is not None:
            self._validate_values(dense, "shape_counts")
            if dense.ndim != 3:
                raise ValueError("Dense shape_counts must have shape [spots, peaks, bins].")
            if min(dense.shape) <= 0:
                raise ValueError("shape_counts axes must all be non-empty.")
            self._dense = dense
            self.shape = tuple(int(value) for value in dense.shape)
            return

        if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
            raise TypeError(
                "shape_counts must be a dense [S,P,B] array or an ordered sequence "
                "of 2-D bin matrices."
            )
        layers = tuple(values)
        if not layers:
            raise ValueError("At least one ordered shape-count layer is required.")
        normalized = []
        matrix_shape = None
        for bin_index, layer in enumerate(layers):
            if sparse.issparse(layer):
                layer = layer.tocsr(copy=True)
                layer.sum_duplicates()
                layer.eliminate_zeros()
                layer.sort_indices()
                values_to_check = layer.data
            else:
                layer = np.asarray(layer)
                if layer.ndim != 2:
                    raise ValueError(f"Shape-count layer {bin_index} must be two-dimensional.")
                values_to_check = layer
            self._validate_values(values_to_check, f"shape_counts[{bin_index}]")
            if matrix_shape is None:
                matrix_shape = layer.shape
            elif layer.shape != matrix_shape:
                raise ValueError("All shape-count layers must have identical axes.")
            normalized.append(layer)
        if matrix_shape is None or min(matrix_shape) <= 0:
            raise ValueError("shape_counts axes must all be non-empty.")
        self._layers = tuple(normalized)
        self.shape = (int(matrix_shape[0]), int(matrix_shape[1]), len(normalized))

    @staticmethod
    def _as_dense_tensor(values: Any) -> Optional[np.ndarray]:
        if isinstance(values, torch.Tensor):
            if values.device.type != "cpu":
                raise ValueError("ShapeMix MAP fitting accepts CPU tensors only.")
            return values.detach().numpy()
        if isinstance(values, np.ndarray):
            return values
        return None

    @staticmethod
    def _validate_values(values: Any, name: str) -> None:
        array = np.asarray(values)
        if not np.issubdtype(array.dtype, np.number):
            raise ValueError(f"{name} must contain numeric counts.")
        if not np.isfinite(array).all():
            raise ValueError(f"{name} contains non-finite counts.")
        if (array < 0).any() or not np.equal(array, np.floor(array)).all():
            raise ValueError(f"{name} must contain nonnegative integer counts.")

    def reordered(self, spot_order: np.ndarray, peak_order: np.ndarray) -> "_CountSource":
        if self._dense is not None:
            return _CountSource(self._dense[spot_order][:, peak_order, :])
        assert self._layers is not None
        return _CountSource(
            [layer[spot_order][:, peak_order] for layer in self._layers]
        )

    def numpy_chunk(self, spot_slice: slice, peak_slice: slice) -> np.ndarray:
        if self._dense is not None:
            return np.asarray(self._dense[spot_slice, peak_slice, :], dtype=np.float32)
        assert self._layers is not None
        chunks = []
        for layer in self._layers:
            selected = layer[spot_slice, peak_slice]
            if sparse.issparse(selected):
                selected = selected.toarray()
            chunks.append(np.asarray(selected, dtype=np.float32))
        return np.stack(chunks, axis=-1)

    def collapsed_row(self, spot_index: int) -> np.ndarray:
        if self._dense is not None:
            return np.asarray(self._dense[spot_index].sum(axis=-1), dtype=np.float64)
        assert self._layers is not None
        result = np.zeros(self.shape[1], dtype=np.float64)
        for layer in self._layers:
            row = layer[spot_index]
            if sparse.issparse(row):
                row = row.toarray()
            result += np.asarray(row, dtype=np.float64).ravel()
        return result


def _canonical_order(
    names: Optional[Sequence[str]],
    length: int,
    prefix: str,
) -> tuple[tuple[str, ...], np.ndarray]:
    if names is None:
        normalized = tuple(f"{prefix}_{index:04d}" for index in range(length))
        return normalized, np.arange(length, dtype=np.int64)
    normalized = tuple(str(name) for name in names)
    if len(normalized) != length:
        raise ValueError(f"{prefix}_names must contain exactly {length} values.")
    if any(not name for name in normalized) or len(set(normalized)) != length:
        raise ValueError(f"{prefix}_names must be non-empty and unique.")
    order = np.asarray(
        sorted(range(length), key=lambda index: normalized[index].encode("utf-8")),
        dtype=np.int64,
    )
    return normalized, order


def _resolve_signatures(
    A: Any,
    omega: Any,
    phi_ref: Any,
    signatures: Any,
) -> tuple[np.ndarray, np.ndarray, float]:
    if signatures is not None:
        if any(value is not None for value in (A, omega, phi_ref)):
            raise ValueError("Pass either signatures or A/omega/phi_ref, not both.")
        try:
            A = signatures.A
            omega = signatures.omega
            phi_ref = signatures.phi_ref
        except AttributeError as error:
            raise TypeError("signatures must expose A, omega, and phi_ref.") from error
    if any(value is None for value in (A, omega, phi_ref)):
        raise ValueError("A, omega, and phi_ref are required fixed signatures.")

    if isinstance(A, torch.Tensor):
        A = A.detach().cpu().numpy()
    if isinstance(omega, torch.Tensor):
        omega = omega.detach().cpu().numpy()
    A_array = np.asarray(A, dtype=np.float64)
    omega_array = np.asarray(omega, dtype=np.float64)
    if A_array.ndim != 2:
        raise ValueError("A must have shape [cell_types, peaks].")
    if omega_array.ndim != 3:
        raise ValueError("omega must have shape [cell_types, peaks, bins].")
    if omega_array.shape[:2] != A_array.shape:
        raise ValueError("A and omega must have identical cell-type and peak axes.")
    if not np.isfinite(A_array).all() or (A_array <= 0).any():
        raise ValueError("Fixed model-version-1 A must contain finite strictly positive rates.")
    if not np.isfinite(omega_array).all() or (omega_array <= 0).any():
        raise ValueError("omega must contain finite strictly positive probabilities.")
    if not np.allclose(omega_array.sum(axis=-1), 1.0, rtol=0.0, atol=1.0e-6):
        raise ValueError("Every omega cell-type/peak row must sum to one.")
    if isinstance(phi_ref, (bool, np.bool_)):
        raise TypeError("phi_ref must be numeric.")
    phi = float(phi_ref)
    if not np.isfinite(phi) or phi <= 0:
        raise ValueError("phi_ref must be finite and positive.")
    return A_array, omega_array, phi


def _nnls_initialization(
    counts: _CountSource,
    A: np.ndarray,
    canonical_spot_names: tuple[str, ...],
) -> tuple[np.ndarray, tuple[str, ...]]:
    n_spots = counts.shape[0]
    n_cell_types = A.shape[0]
    initialization = np.empty((n_spots, n_cell_types), dtype=np.float64)
    fallback_spots = []
    design = np.asarray(A.T, dtype=np.float64)
    uniform_signature_total = float(A.sum())

    for spot_index in range(n_spots):
        totals = counts.collapsed_row(spot_index)
        use_fallback = totals.sum() <= 0
        estimate = None
        if not use_fallback:
            try:
                estimate, _ = nnls(design, totals)
                use_fallback = (
                    not np.isfinite(estimate).all() or float(estimate.sum()) <= 0
                )
            except (RuntimeError, ValueError, np.linalg.LinAlgError):
                use_fallback = True
        if use_fallback:
            if totals.sum() > 0 and uniform_signature_total > 0:
                uniform_value = float(totals.sum()) / uniform_signature_total
            else:
                # The Gamma(shape=2, rate=1) prior mode is the natural fallback
                # for an all-zero spot with no count-based initialization signal.
                uniform_value = 1.0
            estimate = np.full(n_cell_types, uniform_value, dtype=np.float64)
            fallback_spots.append(canonical_spot_names[spot_index])
        initialization[spot_index] = np.maximum(estimate, INITIAL_ABUNDANCE_FLOOR)
    return initialization, tuple(fallback_spots)


def _softplus_inverse(values: np.ndarray, epsilon: float) -> np.ndarray:
    adjusted = np.maximum(np.asarray(values, dtype=np.float64) - epsilon, epsilon)
    # x + log(1-exp(-x)) is stable for both small and large positive x.
    return adjusted + np.log(-np.expm1(-adjusted))


def _seeded_restart(
    initial_abundance: np.ndarray,
    seed_tuple: tuple[int, int, int, int, int],
) -> np.ndarray:
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed_tuple)))
    perturbation = rng.normal(
        loc=0.0,
        scale=RESTART_LOG_STANDARD_DEVIATION,
        size=initial_abundance.shape,
    )
    abundance = initial_abundance * np.exp(perturbation)
    return _softplus_inverse(np.maximum(abundance, INITIAL_ABUNDANCE_FLOOR), ABUNDANCE_EPSILON)


def _torch_signatures(
    A: np.ndarray,
    omega: np.ndarray,
) -> tuple[torch.Tensor, torch.Tensor]:
    return (
        torch.as_tensor(A, dtype=torch.float32, device="cpu"),
        torch.as_tensor(omega, dtype=torch.float32, device="cpu"),
    )


def _complete_objective(
    raw_z: torch.Tensor,
    counts: _CountSource,
    A: torch.Tensor,
    omega: torch.Tensor,
    phi_ref: float,
    config: ShapeMixConfig,
) -> _ObjectiveValues:
    count_value = 0.0
    shape_value = 0.0
    prior_value = 0.0
    with torch.no_grad():
        for spot_start in range(0, counts.shape[0], config.spot_batch_size):
            spot_stop = min(spot_start + config.spot_batch_size, counts.shape[0])
            spot_slice = slice(spot_start, spot_stop)
            z_batch = positive_abundance(raw_z[spot_slice], eps=ABUNDANCE_EPSILON)
            prior = gamma_log_prior(
                z_batch,
                shape=config.abundance_prior_shape,
                rate=config.abundance_prior_rate,
            )
            prior_value += float(prior.detach())
            for peak_start in range(0, counts.shape[1], config.peak_chunk_size):
                peak_stop = min(peak_start + config.peak_chunk_size, counts.shape[1])
                peak_slice = slice(peak_start, peak_stop)
                shape_counts = torch.as_tensor(
                    counts.numpy_chunk(spot_slice, peak_slice),
                    dtype=torch.float32,
                    device="cpu",
                )
                components = likelihood_components(
                    z_batch,
                    shape_counts.sum(dim=-1),
                    A[:, peak_slice],
                    phi_ref,
                    shape_counts=shape_counts,
                    omega=omega[:, peak_slice, :],
                    use_shape=config.use_shape,
                    total_likelihood=config.total_likelihood,
                    eps=ABUNDANCE_EPSILON,
                )
                count_value += float(components.count_log_likelihood.detach())
                shape_value += float(components.shape_log_likelihood.detach())
    return _ObjectiveValues(count_value, shape_value, prior_value)


def _streaming_backward(
    raw_z: torch.Tensor,
    counts: _CountSource,
    A: torch.Tensor,
    omega: torch.Tensor,
    phi_ref: float,
    config: ShapeMixConfig,
) -> Optional[str]:
    """Accumulate an exact full gradient while releasing every chunk graph."""
    for spot_start in range(0, counts.shape[0], config.spot_batch_size):
        spot_stop = min(spot_start + config.spot_batch_size, counts.shape[0])
        spot_slice = slice(spot_start, spot_stop)
        z_for_prior = positive_abundance(raw_z[spot_slice], eps=ABUNDANCE_EPSILON)
        prior = gamma_log_prior(
            z_for_prior,
            shape=config.abundance_prior_shape,
            rate=config.abundance_prior_rate,
        )
        if not torch.isfinite(prior):
            return f"nonfinite_prior_spots_{spot_start}_{spot_stop}"
        (-prior).backward()

        for peak_start in range(0, counts.shape[1], config.peak_chunk_size):
            peak_stop = min(peak_start + config.peak_chunk_size, counts.shape[1])
            peak_slice = slice(peak_start, peak_stop)
            # Recompute z so this small graph can be released immediately after
            # backward instead of retaining a graph for the full S x P x B fit.
            z_chunk = positive_abundance(raw_z[spot_slice], eps=ABUNDANCE_EPSILON)
            shape_counts = torch.as_tensor(
                counts.numpy_chunk(spot_slice, peak_slice),
                dtype=torch.float32,
                device="cpu",
            )
            components = likelihood_components(
                z_chunk,
                shape_counts.sum(dim=-1),
                A[:, peak_slice],
                phi_ref,
                shape_counts=shape_counts,
                omega=omega[:, peak_slice, :],
                use_shape=config.use_shape,
                total_likelihood=config.total_likelihood,
                eps=ABUNDANCE_EPSILON,
            )
            objective = components.total_log_objective
            if not torch.isfinite(objective):
                return (
                    f"nonfinite_likelihood_spots_{spot_start}_{spot_stop}_"
                    f"peaks_{peak_start}_{peak_stop}"
                )
            (-objective).backward()
    return None


def _history_record(step: int, summary: _ObjectiveValues) -> OptimizationStepRecord:
    return OptimizationStepRecord(
        step=step,
        count_log_likelihood=summary.count_log_likelihood,
        shape_log_likelihood=summary.shape_log_likelihood,
        abundance_log_prior=summary.abundance_log_prior,
        total_log_objective=summary.total_log_objective,
    )


def fit_shapemix_map(
    shape_counts: Any,
    A: Any = None,
    omega: Any = None,
    phi_ref: Any = None,
    *,
    signatures: Any = None,
    config: Optional[Union[ShapeMixConfig, Mapping[str, Any]]] = None,
    outer_split_seed: int,
    inner_mixture_seed: int,
    spot_names: Optional[Sequence[str]] = None,
    feature_names: Optional[Sequence[str]] = None,
    cell_types: Optional[Sequence[str]] = None,
) -> MAPFitResult:
    """Fit positive effective abundances with deterministic streamed Adam MAP.

    ``shape_counts`` may be a small dense ``[S,P,B]`` tensor/array or an
    ordered sequence of two-dimensional SciPy CSR bin layers.  Sparse inputs
    are densified only one configured spot/peak chunk at a time.
    """
    if config is None:
        config = ShapeMixConfig()
    elif isinstance(config, Mapping):
        config = ShapeMixConfig.from_mapping(config)
    elif not isinstance(config, ShapeMixConfig):
        raise TypeError("config must be a ShapeMixConfig or strict parameter mapping.")
    if isinstance(outer_split_seed, (bool, np.bool_)) or not isinstance(
        outer_split_seed, (int, np.integer)
    ):
        raise TypeError("outer_split_seed must be an integer.")
    if isinstance(inner_mixture_seed, (bool, np.bool_)) or not isinstance(
        inner_mixture_seed, (int, np.integer)
    ):
        raise TypeError("inner_mixture_seed must be an integer.")
    outer_split_seed = int(outer_split_seed)
    inner_mixture_seed = int(inner_mixture_seed)
    if outer_split_seed < 0 or inner_mixture_seed < 0:
        raise ValueError("Split and mixture seeds must be nonnegative.")

    if signatures is not None:
        try:
            signature_cell_types = tuple(signatures.cell_types)
            signature_peak_ids = tuple(signatures.peak_ids)
        except AttributeError as error:
            raise TypeError(
                "signatures must expose ordered cell_types and peak_ids axes."
            ) from error
        if cell_types is None:
            cell_types = signature_cell_types
        elif tuple(cell_types) != signature_cell_types:
            raise ValueError("Explicit cell_types do not match signatures.cell_types order.")
        if feature_names is None:
            feature_names = signature_peak_ids
        elif tuple(feature_names) != signature_peak_ids:
            raise ValueError("Explicit feature_names do not match signatures.peak_ids order.")

    counts = _CountSource(shape_counts)
    A_array, omega_array, phi = _resolve_signatures(A, omega, phi_ref, signatures)
    if A_array.shape[1] != counts.shape[1] or omega_array.shape[2] != counts.shape[2]:
        raise ValueError("Observed counts and fixed signatures have incompatible peak/bin axes.")
    n_cell_types = A_array.shape[0]
    if cell_types is None:
        ordered_cell_types = tuple(f"cell_type_{index:02d}" for index in range(n_cell_types))
    else:
        ordered_cell_types = tuple(str(value) for value in cell_types)
        if (
            len(ordered_cell_types) != n_cell_types
            or any(not value for value in ordered_cell_types)
            or len(set(ordered_cell_types)) != n_cell_types
        ):
            raise ValueError("cell_types must exactly name the unique A cell-type axis.")

    original_spot_names, spot_order = _canonical_order(spot_names, counts.shape[0], "spot")
    _, peak_order = _canonical_order(feature_names, counts.shape[1], "feature")
    canonical_spot_names = tuple(original_spot_names[index] for index in spot_order)
    counts = counts.reordered(spot_order, peak_order)
    A_array = A_array[:, peak_order]
    omega_array = omega_array[:, peak_order, :]
    initial_abundance, fallback_spots = _nnls_initialization(
        counts,
        A_array,
        canonical_spot_names,
    )
    A_tensor, omega_tensor = _torch_signatures(A_array, omega_array)

    fit_started = time.perf_counter()
    diagnostics: list[RestartRecord] = []
    converged_candidates: list[tuple[float, int, np.ndarray, _ObjectiveValues]] = []
    for restart_index in range(config.restarts):
        restart_started = time.perf_counter()
        seed_tuple = (
            SEED_NAMESPACE,
            outer_split_seed,
            inner_mixture_seed,
            config.seed,
            restart_index,
        )
        raw_initial = _seeded_restart(initial_abundance, seed_tuple)
        raw_z = torch.nn.Parameter(torch.as_tensor(raw_initial, dtype=torch.float32, device="cpu"))
        optimizer = torch.optim.Adam([raw_z], lr=config.learning_rate)

        initial_summary = _complete_objective(
            raw_z, counts, A_tensor, omega_tensor, phi, config
        )
        history = [_history_record(0, initial_summary)]
        nonfinite_events: list[str] = []
        best_summary = initial_summary
        best_raw = raw_z.detach().clone()
        best_step = 0
        # Center the conditional term at its restart-specific initial value.
        # Its absolute log probability contains data-only constants that must
        # not change the relative stopping scale between the nested arms.  Its
        # informative change during optimization is still included in full.
        initial_shape_objective = initial_summary.shape_log_likelihood
        anchor_objective = (
            initial_summary.count_log_likelihood
            + initial_summary.abundance_log_prior
        )
        stale_steps = 0
        steps_completed = 0
        converged = False
        stopping_reason = "max_steps"
        final_gradient_norm = math.nan

        if not initial_summary.finite:
            nonfinite_events.append("nonfinite_initial_complete_objective")
            stopping_reason = "nonfinite_initial_objective"
        else:
            for step in range(1, config.max_steps + 1):
                optimizer.zero_grad(set_to_none=True)
                nonfinite_chunk = _streaming_backward(
                    raw_z,
                    counts,
                    A_tensor,
                    omega_tensor,
                    phi,
                    config,
                )
                if nonfinite_chunk is not None:
                    nonfinite_events.append(f"step_{step}:{nonfinite_chunk}")
                    stopping_reason = "nonfinite_objective"
                    break
                if raw_z.grad is None:
                    nonfinite_events.append(f"step_{step}:missing_gradient")
                    stopping_reason = "missing_gradient"
                    break
                if not torch.isfinite(raw_z.grad).all():
                    nonfinite_events.append(f"step_{step}:nonfinite_gradient")
                    stopping_reason = "nonfinite_gradient"
                    break
                final_gradient_norm = float(
                    torch.linalg.vector_norm(raw_z.grad.detach()).cpu().item()
                )
                optimizer.step()
                steps_completed = step
                if not torch.isfinite(raw_z).all():
                    nonfinite_events.append(f"step_{step}:nonfinite_parameter")
                    stopping_reason = "nonfinite_parameter"
                    break

                summary = _complete_objective(
                    raw_z, counts, A_tensor, omega_tensor, phi, config
                )
                history.append(_history_record(step, summary))
                objective_value = summary.total_log_objective
                if not summary.finite:
                    nonfinite_events.append(f"step_{step}:nonfinite_complete_objective")
                    stopping_reason = "nonfinite_complete_objective"
                    break
                if objective_value > best_summary.total_log_objective:
                    best_summary = summary
                    best_raw = raw_z.detach().clone()
                    best_step = step

                stopping_objective = (
                    summary.count_log_likelihood
                    + summary.abundance_log_prior
                    + summary.shape_log_likelihood
                    - initial_shape_objective
                )
                required_improvement = config.tolerance * max(1.0, abs(anchor_objective))
                if stopping_objective > anchor_objective + required_improvement:
                    anchor_objective = stopping_objective
                    stale_steps = 0
                else:
                    stale_steps += 1
                if stale_steps >= config.patience:
                    converged = True
                    stopping_reason = "objective_patience"
                    break
                if final_gradient_norm <= config.tolerance:
                    converged = True
                    stopping_reason = "gradient_tolerance"
                    break

        finite = best_summary.finite and not nonfinite_events
        if finite and converged:
            with torch.no_grad():
                best_abundance = (
                    positive_abundance(best_raw, eps=ABUNDANCE_EPSILON)
                    .detach()
                    .cpu()
                    .numpy()
                    .astype(np.float64)
                )
            converged_candidates.append(
                (
                    best_summary.total_log_objective,
                    restart_index,
                    best_abundance,
                    best_summary,
                )
            )
        diagnostics.append(
            RestartRecord(
                restart_index=restart_index,
                seed_tuple=seed_tuple,
                converged=converged,
                stopping_reason=stopping_reason,
                steps=steps_completed,
                count_log_likelihood=best_summary.count_log_likelihood,
                shape_log_likelihood=best_summary.shape_log_likelihood,
                abundance_log_prior=best_summary.abundance_log_prior,
                total_log_objective=best_summary.total_log_objective,
                history=tuple(history),
                nonfinite_events=tuple(nonfinite_events),
                runtime_seconds=time.perf_counter() - restart_started,
                best_step=best_step,
                final_gradient_norm=final_gradient_norm,
            )
        )

    if not converged_candidates:
        reasons = ", ".join(
            f"restart {record.restart_index}: {record.stopping_reason}"
            for record in diagnostics
        )
        failure_diagnostics = FitRecord(
            use_shape=config.use_shape,
            success=False,
            selected_restart=None,
            stopping_reason="no_finite_converged_restart",
            n_spots=counts.shape[0],
            n_peaks=counts.shape[1],
            n_bins=counts.shape[2],
            n_cell_types=n_cell_types,
            dtype=config.dtype,
            device=config.device,
            restarts=tuple(diagnostics),
            nonfinite_events=tuple(
                f"restart_{record.restart_index}:{event}"
                for record in diagnostics
                for event in record.nonfinite_events
            ),
            runtime_seconds=time.perf_counter() - fit_started,
        )
        raise MAPOptimizationError(
            f"No finite converged ShapeMix MAP restart succeeded ({reasons}).",
            failure_diagnostics,
        )
    _, selected_restart, canonical_abundance, selected_summary = max(
        converged_candidates,
        key=lambda candidate: (candidate[0], -candidate[1]),
    )
    canonical_proportions = canonical_abundance / canonical_abundance.sum(axis=1, keepdims=True)
    if (
        not np.isfinite(canonical_abundance).all()
        or (canonical_abundance <= 0).any()
        or not np.isfinite(canonical_proportions).all()
        or (canonical_proportions < 0).any()
        or not np.allclose(canonical_proportions.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-6)
    ):
        invalid_diagnostics = FitRecord(
            use_shape=config.use_shape,
            success=False,
            selected_restart=None,
            stopping_reason="invalid_selected_output",
            n_spots=counts.shape[0],
            n_peaks=counts.shape[1],
            n_bins=counts.shape[2],
            n_cell_types=n_cell_types,
            dtype=config.dtype,
            device=config.device,
            restarts=tuple(diagnostics),
            runtime_seconds=time.perf_counter() - fit_started,
        )
        raise MAPOptimizationError(
            "Selected MAP restart produced invalid output proportions.",
            invalid_diagnostics,
        )

    abundance = np.empty_like(canonical_abundance)
    proportions = np.empty_like(canonical_proportions)
    abundance[spot_order] = canonical_abundance
    proportions[spot_order] = canonical_proportions
    fit_diagnostics = FitRecord(
        use_shape=config.use_shape,
        success=True,
        selected_restart=selected_restart,
        stopping_reason="selected_largest_finite_converged_objective",
        n_spots=counts.shape[0],
        n_peaks=counts.shape[1],
        n_bins=counts.shape[2],
        n_cell_types=n_cell_types,
        dtype=config.dtype,
        device=config.device,
        restarts=tuple(diagnostics),
        nonfinite_events=tuple(
            f"restart_{record.restart_index}:{event}"
            for record in diagnostics
            for event in record.nonfinite_events
        ),
        runtime_seconds=time.perf_counter() - fit_started,
    )
    return MAPFitResult(
        abundance=abundance,
        proportions=proportions,
        selected_restart=selected_restart,
        count_log_likelihood=selected_summary.count_log_likelihood,
        shape_log_likelihood=selected_summary.shape_log_likelihood,
        abundance_log_prior=selected_summary.abundance_log_prior,
        diagnostics=fit_diagnostics,
        nnls_fallback_spots=tuple(sorted(fallback_spots, key=lambda value: value.encode("utf-8"))),
        spot_names=original_spot_names,
        cell_types=ordered_cell_types,
        use_shape=config.use_shape,
    )
