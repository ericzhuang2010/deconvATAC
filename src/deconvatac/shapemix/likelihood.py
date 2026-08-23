"""Stable likelihood and prior primitives for ShapeMix model version 1.

All public objective values are *log* objectives to maximize.  The pure
``likelihood_components`` function operates on one spot/peak chunk and is the
primitive intended for memory-bounded gradient accumulation in the MAP
optimizer.  The convenience wrappers aggregate small problems without ever
constructing a dense expected tensor larger than the requested chunk.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Optional, Tuple

import torch
import torch.nn.functional as F
from torch import Tensor


DEFAULT_EPS = 1.0e-8
DEFAULT_PRIOR_SHAPE = 2.0
DEFAULT_PRIOR_RATE = 1.0

TotalLikelihood = Literal["negative_binomial", "poisson"]


@dataclass(frozen=True)
class ObjectiveComponents:
    """Scalar terms in the ShapeMix log objective.

    ``shape_log_likelihood`` is an exact scalar zero when shape is disabled.
    ``abundance_log_prior`` is an exact scalar zero for a chunk-level
    likelihood result; full-objective helpers add the prior once.
    """

    count_log_likelihood: Tensor
    shape_log_likelihood: Tensor
    abundance_log_prior: Tensor

    @property
    def total_log_objective(self) -> Tensor:
        """Return the complete scalar log objective to maximize."""
        return (
            self.count_log_likelihood
            + self.shape_log_likelihood
            + self.abundance_log_prior
        )

    def detached_values(self) -> dict[str, float]:
        """Return JSON-safe scalar values detached from autograd."""
        values = {
            "count_log_likelihood": self.count_log_likelihood,
            "shape_log_likelihood": self.shape_log_likelihood,
            "abundance_log_prior": self.abundance_log_prior,
            "total_log_objective": self.total_log_objective,
        }
        result: dict[str, float] = {}
        for name, value in values.items():
            if value.numel() != 1:
                raise ValueError(f"{name} must be scalar, observed shape {tuple(value.shape)}.")
            result[name] = float(value.detach().cpu().item())
        return result


def _validate_eps(eps: float) -> float:
    if isinstance(eps, bool) or not isinstance(eps, (int, float)):
        raise TypeError("eps must be a real number.")
    normalized = float(eps)
    if not math.isfinite(normalized) or normalized <= 0.0:
        raise ValueError("eps must be finite and strictly positive.")
    return normalized


def _require_floating_matrix(value: Tensor, name: str) -> None:
    if not isinstance(value, Tensor):
        raise TypeError(f"{name} must be a torch.Tensor.")
    if value.ndim != 2:
        raise ValueError(f"{name} must have two dimensions, observed {value.ndim}.")
    if not torch.is_floating_point(value):
        raise TypeError(f"{name} must have a floating-point dtype.")


def _check_same_device_dtype(left: Tensor, right: Tensor, left_name: str, right_name: str) -> None:
    if left.device != right.device:
        raise ValueError(f"{left_name} and {right_name} must be on the same device.")
    if left.dtype != right.dtype:
        raise TypeError(f"{left_name} and {right_name} must have the same dtype.")


def _validate_finite_nonnegative(value: Tensor, name: str) -> None:
    """Validate public probability-function inputs outside hot optimizer paths."""
    if not bool(torch.isfinite(value).all().item()):
        raise ValueError(f"{name} must contain only finite values.")
    if bool((value < 0).any().item()):
        raise ValueError(f"{name} must be nonnegative.")


def _validate_counts(value: Tensor, name: str) -> None:
    _validate_finite_nonnegative(value, name)
    if torch.is_floating_point(value) and bool(
        (value != torch.floor(value)).any().item()
    ):
        raise ValueError(f"{name} must contain integer-valued counts.")


def positive_abundance(raw_z: Tensor, eps: float = DEFAULT_EPS) -> Tensor:
    """Transform unconstrained spot abundances to strictly positive values."""
    _require_floating_matrix(raw_z, "raw_z")
    eps = _validate_eps(eps)
    return F.softplus(raw_z) + eps


def expected_total_rate(z: Tensor, accessibility: Tensor) -> Tensor:
    """Return ``mu[s,p] = sum_c z[s,c] * A[c,p]`` for one chunk."""
    _require_floating_matrix(z, "z")
    _require_floating_matrix(accessibility, "accessibility")
    _check_same_device_dtype(z, accessibility, "z", "accessibility")
    if z.shape[1] != accessibility.shape[0]:
        raise ValueError(
            "z cell-type dimension must equal the first accessibility dimension: "
            f"{z.shape[1]} != {accessibility.shape[0]}."
        )
    return z @ accessibility


def expected_bin_rates(z: Tensor, accessibility: Tensor, omega: Tensor) -> Tensor:
    """Return ``v[s,p,b]`` for one spot/peak chunk.

    Only the returned chunk has shape ``S_chunk x P_chunk x B``.  Callers must
    pass normalized ``omega[c,p,:]`` from the fixed-signature contract.
    """
    _require_floating_matrix(z, "z")
    _require_floating_matrix(accessibility, "accessibility")
    _check_same_device_dtype(z, accessibility, "z", "accessibility")
    if not isinstance(omega, Tensor):
        raise TypeError("omega must be a torch.Tensor.")
    if omega.ndim != 3:
        raise ValueError(f"omega must have three dimensions, observed {omega.ndim}.")
    _check_same_device_dtype(z, omega, "z", "omega")
    if omega.shape[:2] != accessibility.shape:
        raise ValueError(
            "omega cell-type/peak dimensions must equal accessibility shape: "
            f"{tuple(omega.shape[:2])} != {tuple(accessibility.shape)}."
        )
    if omega.shape[2] < 1:
        raise ValueError("omega must contain at least one bin.")
    if z.shape[1] != accessibility.shape[0]:
        raise ValueError(
            "z cell-type dimension must equal the first accessibility dimension: "
            f"{z.shape[1]} != {accessibility.shape[0]}."
        )
    weighted_shape = accessibility.unsqueeze(-1) * omega
    return torch.einsum("sc,cpb->spb", z, weighted_shape)


def negative_binomial_log_prob(
    counts: Tensor,
    mean: Tensor,
    inverse_dispersion: Tensor,
    eps: float = DEFAULT_EPS,
    *,
    validate_args: bool = True,
) -> Tensor:
    """Return elementwise NB log probabilities in mean/inverse-dispersion form.

    The parameterization has ``E[X] = mean`` and
    ``Var[X] = mean + mean**2 / inverse_dispersion``.  The calculation uses
    ``log1p(mean / inverse_dispersion)`` and evaluates the gamma-ratio portion
    in float64 when inputs are float32.  This avoids the severe cancellation
    that otherwise occurs at the protocol's near-Poisson dispersion floor;
    the returned tensor retains ``mean.dtype`` and gradients flow back to the
    original float32 parameters.
    """
    eps = _validate_eps(eps)
    if not isinstance(counts, Tensor) or not isinstance(mean, Tensor):
        raise TypeError("counts and mean must be torch.Tensor instances.")
    if not isinstance(inverse_dispersion, Tensor):
        inverse_dispersion = torch.as_tensor(
            inverse_dispersion, dtype=mean.dtype, device=mean.device
        )
    if not torch.is_floating_point(mean):
        raise TypeError("mean must have a floating-point dtype.")
    if inverse_dispersion.device != mean.device:
        raise ValueError("inverse_dispersion and mean must be on the same device.")

    if validate_args:
        _validate_counts(counts, "counts")
        _validate_finite_nonnegative(mean, "mean")
        if not bool(torch.isfinite(inverse_dispersion).all().item()) or bool(
            (inverse_dispersion <= 0).any().item()
        ):
            raise ValueError("inverse_dispersion must be finite and strictly positive.")

    output_dtype = mean.dtype
    work_dtype = (
        torch.float64
        if output_dtype in (torch.float16, torch.bfloat16, torch.float32)
        else output_dtype
    )
    n = counts.to(device=mean.device, dtype=work_dtype)
    mu = mean.to(dtype=work_dtype)
    dispersion = inverse_dispersion.to(dtype=work_dtype)
    n, mu, dispersion = torch.broadcast_tensors(n, mu, dispersion)

    safe_mu = mu.clamp_min(eps)
    safe_dispersion = dispersion.clamp_min(eps)
    log_one_plus_ratio = torch.log1p(mu / safe_dispersion)
    log_mu_term = torch.where(
        n > 0,
        n
        * (
            torch.log(safe_mu)
            - torch.log(safe_dispersion)
            - log_one_plus_ratio
        ),
        torch.zeros_like(n),
    )
    log_probability = (
        torch.lgamma(n + safe_dispersion)
        - torch.lgamma(safe_dispersion)
        - torch.lgamma(n + 1.0)
        - safe_dispersion * log_one_plus_ratio
        + log_mu_term
    )
    return log_probability.to(dtype=output_dtype)


def poisson_log_prob(
    counts: Tensor,
    rate: Tensor,
    eps: float = DEFAULT_EPS,
    *,
    validate_args: bool = True,
) -> Tensor:
    """Return elementwise Poisson log probabilities with exact zero handling."""
    eps = _validate_eps(eps)
    if not isinstance(counts, Tensor) or not isinstance(rate, Tensor):
        raise TypeError("counts and rate must be torch.Tensor instances.")
    if not torch.is_floating_point(rate):
        raise TypeError("rate must have a floating-point dtype.")
    if validate_args:
        _validate_counts(counts, "counts")
        _validate_finite_nonnegative(rate, "rate")
    n = counts.to(device=rate.device, dtype=rate.dtype)
    n, rate = torch.broadcast_tensors(n, rate)
    positive_rate = rate > 0
    safe_rate = torch.where(positive_rate, rate, torch.ones_like(rate))
    log_rate = torch.where(
        positive_rate,
        torch.log(safe_rate),
        torch.full_like(rate, -torch.inf),
    )
    selected_log_rate = torch.where(n > 0, log_rate, torch.zeros_like(log_rate))
    count_times_log_rate = n * selected_log_rate
    return count_times_log_rate - rate - torch.lgamma(n + 1.0)


def conditional_multinomial_log_prob(
    shape_counts: Tensor,
    bin_rates: Tensor,
    eps: float = DEFAULT_EPS,
    *,
    validate_args: bool = True,
) -> Tensor:
    """Return direct conditional-multinomial log probabilities per row.

    The final axis is the bin axis.  Probabilities are normalized with
    log-sum-exp.  Observed zero bins contribute an exact zero, and a row whose
    observed total is zero is explicitly replaced by an exact zero.
    """
    eps = _validate_eps(eps)
    if not isinstance(shape_counts, Tensor) or not isinstance(bin_rates, Tensor):
        raise TypeError("shape_counts and bin_rates must be torch.Tensor instances.")
    if shape_counts.shape != bin_rates.shape:
        raise ValueError(
            "shape_counts and bin_rates must have identical shapes: "
            f"{tuple(shape_counts.shape)} != {tuple(bin_rates.shape)}."
        )
    if shape_counts.ndim < 1 or shape_counts.shape[-1] < 1:
        raise ValueError("shape_counts and bin_rates require a nonempty final bin axis.")
    if not torch.is_floating_point(bin_rates):
        raise TypeError("bin_rates must have a floating-point dtype.")

    if validate_args:
        _validate_counts(shape_counts, "shape_counts")
        _validate_finite_nonnegative(bin_rates, "bin_rates")

    counts = shape_counts.to(device=bin_rates.device, dtype=bin_rates.dtype)
    positive_rates = bin_rates > 0
    # Taking log(1), followed by a constant -inf replacement, gives exact
    # support semantics at zero without a 0/0 gradient in LogBackward.
    safe_rates = torch.where(positive_rates, bin_rates, torch.ones_like(bin_rates))
    log_weights = torch.where(
        positive_rates,
        torch.log(safe_rates),
        torch.full_like(bin_rates, -torch.inf),
    )
    has_positive_rate = positive_rates.any(dim=-1, keepdim=True)
    normalization_weights = torch.where(
        has_positive_rate, log_weights, torch.zeros_like(log_weights)
    )
    log_probabilities = log_weights - torch.logsumexp(
        normalization_weights, dim=-1, keepdim=True
    )
    selected_log_probabilities = torch.where(
        counts > 0, log_probabilities, torch.zeros_like(log_probabilities)
    )
    observed_log_probability = (counts * selected_log_probabilities).sum(dim=-1)
    totals = counts.sum(dim=-1)
    combinatorial = torch.lgamma(totals + 1.0) - torch.lgamma(counts + 1.0).sum(
        dim=-1
    )
    result = combinatorial + observed_log_probability
    return torch.where(totals > 0, result, torch.zeros_like(result))


def gamma_log_prior(
    z: Tensor,
    shape: float = DEFAULT_PRIOR_SHAPE,
    rate: float = DEFAULT_PRIOR_RATE,
    eps: float = DEFAULT_EPS,
) -> Tensor:
    """Return the summed normalized Gamma(shape, rate) abundance log prior."""
    eps = _validate_eps(eps)
    if not isinstance(z, Tensor) or not torch.is_floating_point(z):
        raise TypeError("z must be a floating-point torch.Tensor.")
    for value, name in ((shape, "shape"), (rate, "rate")):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError(f"{name} must be a real number.")
        if not math.isfinite(float(value)) or float(value) <= 0.0:
            raise ValueError(f"{name} must be finite and strictly positive.")
    shape_tensor = z.new_tensor(float(shape))
    rate_tensor = z.new_tensor(float(rate))
    elementwise = (
        shape_tensor * torch.log(rate_tensor)
        - torch.lgamma(shape_tensor)
        + (shape_tensor - 1.0) * torch.log(z.clamp_min(eps))
        - rate_tensor * z
    )
    return elementwise.sum()


def _spot_inverse_dispersion(
    z: Tensor,
    phi_ref: Tensor,
    eps: float,
) -> Tensor:
    phi = torch.as_tensor(phi_ref, dtype=z.dtype, device=z.device)
    if phi.numel() != 1:
        raise ValueError("phi_ref must be scalar.")
    effective_abundance = z.sum(dim=-1, keepdim=True)
    return effective_abundance.clamp_min(eps) * phi.reshape(())


def likelihood_components(
    z: Tensor,
    total_counts: Tensor,
    accessibility: Tensor,
    phi_ref: Tensor,
    *,
    shape_counts: Optional[Tensor] = None,
    omega: Optional[Tensor] = None,
    use_shape: bool,
    total_likelihood: TotalLikelihood = "negative_binomial",
    eps: float = DEFAULT_EPS,
) -> ObjectiveComponents:
    """Evaluate one spot/peak chunk without adding the abundance prior.

    Parameters
    ----------
    z
        Positive abundance with shape ``S_chunk x C``.
    total_counts
        Observed collapsed counts with shape ``S_chunk x P_chunk``.
    accessibility
        Fixed total signature with shape ``C x P_chunk``.
    phi_ref
        Scalar reference inverse-dispersion.  It is ignored only by the
        Poisson golden-test likelihood.
    shape_counts, omega
        Required only when ``use_shape`` is true, with shapes
        ``S_chunk x P_chunk x B`` and ``C x P_chunk x B``.

    Notes
    -----
    The returned prior is zero so this function may be called once per peak
    chunk without multiplying the prior.  For streaming gradient accumulation,
    backpropagate each returned likelihood term, then backpropagate
    ``gamma_log_prior(z)`` exactly once before the optimizer step.
    """
    eps = _validate_eps(eps)
    if not isinstance(use_shape, bool):
        raise TypeError("use_shape must be a boolean.")
    _require_floating_matrix(z, "z")
    mean = expected_total_rate(z, accessibility)
    if not isinstance(total_counts, Tensor) or total_counts.ndim != 2:
        raise ValueError("total_counts must be a two-dimensional torch.Tensor.")
    if total_counts.shape != mean.shape:
        raise ValueError(
            "total_counts shape must match the expected total-rate shape: "
            f"{tuple(total_counts.shape)} != {tuple(mean.shape)}."
        )
    observed_totals = total_counts.to(device=mean.device, dtype=mean.dtype)

    if total_likelihood == "negative_binomial":
        inverse_dispersion = _spot_inverse_dispersion(z, phi_ref, eps)
        count_log_likelihood = negative_binomial_log_prob(
            observed_totals,
            mean,
            inverse_dispersion,
            eps=eps,
            validate_args=False,
        ).sum()
    elif total_likelihood == "poisson":
        count_log_likelihood = poisson_log_prob(
            observed_totals, mean, eps=eps, validate_args=False
        ).sum()
    else:
        raise ValueError(
            "total_likelihood must be 'negative_binomial' or the Poisson test oracle."
        )

    shape_log_likelihood = z.new_zeros(())
    if use_shape:
        if shape_counts is None or omega is None:
            raise ValueError("shape_counts and omega are required when use_shape=True.")
        # If every cell type has the exact same bin probabilities, the
        # conditional likelihood is mathematically independent of abundance.
        # Evaluate that case from the shared probabilities themselves.  A
        # zero-valued connection to ``z`` preserves autograd compatibility
        # while guaranteeing an exact zero gradient; routing this case through
        # expected_bin_rates would leave cancellation roundoff that adaptive
        # optimizers can amplify into a different fit.
        shared_omega = omega[0:1, :, :].expand_as(omega)
        if torch.equal(omega, shared_omega):
            zero_z_connection = z[:, :1].unsqueeze(-1) * 0.0
            bin_rates = omega[0, :, :].unsqueeze(0) + zero_z_connection
        else:
            bin_rates = expected_bin_rates(z, accessibility, omega)
        if shape_counts.shape != bin_rates.shape:
            raise ValueError(
                "shape_counts shape must equal the expected bin-rate shape: "
                f"{tuple(shape_counts.shape)} != {tuple(bin_rates.shape)}."
            )
        shape_log_likelihood = conditional_multinomial_log_prob(
            shape_counts, bin_rates, eps=eps, validate_args=False
        ).sum()

    return ObjectiveComponents(
        count_log_likelihood=count_log_likelihood,
        shape_log_likelihood=shape_log_likelihood,
        abundance_log_prior=z.new_zeros(()),
    )


def _observed_count_parts(observed_counts: Tensor) -> Tuple[Tensor, Optional[Tensor]]:
    if not isinstance(observed_counts, Tensor):
        raise TypeError("observed_counts must be a torch.Tensor.")
    if observed_counts.ndim == 2:
        return observed_counts, None
    if observed_counts.ndim == 3 and observed_counts.shape[-1] >= 1:
        return observed_counts.sum(dim=-1), observed_counts
    raise ValueError(
        "observed_counts must have shape S x P or S x P x B with at least one bin."
    )


def objective_components_from_abundance(
    z: Tensor,
    accessibility: Tensor,
    observed_counts: Tensor,
    phi_ref: Tensor,
    *,
    omega: Optional[Tensor] = None,
    use_shape: bool,
    total_likelihood: TotalLikelihood = "negative_binomial",
    prior_shape: float = DEFAULT_PRIOR_SHAPE,
    prior_rate: float = DEFAULT_PRIOR_RATE,
    eps: float = DEFAULT_EPS,
) -> ObjectiveComponents:
    """Evaluate the complete objective from already-positive abundance."""
    total_counts, shape_counts = _observed_count_parts(observed_counts)
    likelihood = likelihood_components(
        z,
        total_counts,
        accessibility,
        phi_ref,
        shape_counts=shape_counts,
        omega=omega,
        use_shape=use_shape,
        total_likelihood=total_likelihood,
        eps=eps,
    )
    prior = gamma_log_prior(
        z, shape=prior_shape, rate=prior_rate, eps=eps
    )
    return ObjectiveComponents(
        count_log_likelihood=likelihood.count_log_likelihood,
        shape_log_likelihood=likelihood.shape_log_likelihood,
        abundance_log_prior=prior,
    )


def objective_components_from_raw(
    raw_z: Tensor,
    accessibility: Tensor,
    observed_counts: Tensor,
    phi_ref: Tensor,
    *,
    omega: Optional[Tensor] = None,
    use_shape: bool,
    total_likelihood: TotalLikelihood = "negative_binomial",
    prior_shape: float = DEFAULT_PRIOR_SHAPE,
    prior_rate: float = DEFAULT_PRIOR_RATE,
    eps: float = DEFAULT_EPS,
) -> ObjectiveComponents:
    """Transform ``raw_z`` and evaluate the complete objective."""
    return objective_components_from_abundance(
        positive_abundance(raw_z, eps=eps),
        accessibility,
        observed_counts,
        phi_ref,
        omega=omega,
        use_shape=use_shape,
        total_likelihood=total_likelihood,
        prior_shape=prior_shape,
        prior_rate=prior_rate,
        eps=eps,
    )


def _normalize_chunk_size(value: Optional[int], total: int, name: str) -> int:
    if value is None:
        return total
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer or None.")
    if value <= 0:
        raise ValueError(f"{name} must be strictly positive.")
    return min(value, total)


def chunked_objective_components_from_abundance(
    z: Tensor,
    accessibility: Tensor,
    observed_counts: Tensor,
    phi_ref: Tensor,
    *,
    omega: Optional[Tensor] = None,
    use_shape: bool,
    spot_batch_size: Optional[int] = None,
    peak_chunk_size: Optional[int] = None,
    total_likelihood: TotalLikelihood = "negative_binomial",
    prior_shape: float = DEFAULT_PRIOR_SHAPE,
    prior_rate: float = DEFAULT_PRIOR_RATE,
    eps: float = DEFAULT_EPS,
) -> ObjectiveComponents:
    """Aggregate objective components over bounded spot and peak chunks.

    This convenience function is useful for evaluation and parity tests.  Its
    scalar result retains the autograd graph for all chunks.  Optimizers that
    need strict peak-bounded backward memory should call
    :func:`likelihood_components` on each chunk and backpropagate each scalar
    before adding the prior once.
    """
    _require_floating_matrix(z, "z")
    _require_floating_matrix(accessibility, "accessibility")
    total_counts, shape_counts = _observed_count_parts(observed_counts)
    n_spots, n_peaks = total_counts.shape
    if n_spots < 1 or n_peaks < 1:
        raise ValueError("observed_counts must contain at least one spot and one peak.")
    if z.shape[0] != n_spots:
        raise ValueError("z and observed_counts must have the same number of spots.")
    if accessibility.shape[1] != n_peaks:
        raise ValueError(
            "accessibility and observed_counts must have the same number of peaks."
        )
    if use_shape:
        if shape_counts is None or omega is None:
            raise ValueError(
                "Three-dimensional observed_counts and omega are required when use_shape=True."
            )
        if omega.shape != (
            accessibility.shape[0],
            n_peaks,
            shape_counts.shape[2],
        ):
            raise ValueError("omega dimensions do not match accessibility/observed_counts.")

    spot_batch_size = _normalize_chunk_size(
        spot_batch_size, n_spots, "spot_batch_size"
    )
    peak_chunk_size = _normalize_chunk_size(
        peak_chunk_size, n_peaks, "peak_chunk_size"
    )
    count_log_likelihood = z.new_zeros(())
    shape_log_likelihood = z.new_zeros(())
    for spot_start in range(0, n_spots, spot_batch_size):
        spot_stop = min(spot_start + spot_batch_size, n_spots)
        z_chunk = z[spot_start:spot_stop]
        for peak_start in range(0, n_peaks, peak_chunk_size):
            peak_stop = min(peak_start + peak_chunk_size, n_peaks)
            chunk = likelihood_components(
                z_chunk,
                total_counts[spot_start:spot_stop, peak_start:peak_stop],
                accessibility[:, peak_start:peak_stop],
                phi_ref,
                shape_counts=(
                    None
                    if shape_counts is None
                    else shape_counts[
                        spot_start:spot_stop, peak_start:peak_stop, :
                    ]
                ),
                omega=(
                    None
                    if omega is None
                    else omega[:, peak_start:peak_stop, :]
                ),
                use_shape=use_shape,
                total_likelihood=total_likelihood,
                eps=eps,
            )
            count_log_likelihood = (
                count_log_likelihood + chunk.count_log_likelihood
            )
            shape_log_likelihood = (
                shape_log_likelihood + chunk.shape_log_likelihood
            )
    return ObjectiveComponents(
        count_log_likelihood=count_log_likelihood,
        shape_log_likelihood=shape_log_likelihood,
        abundance_log_prior=gamma_log_prior(
            z, shape=prior_shape, rate=prior_rate, eps=eps
        ),
    )


def chunked_objective_components_from_raw(
    raw_z: Tensor,
    accessibility: Tensor,
    observed_counts: Tensor,
    phi_ref: Tensor,
    *,
    omega: Optional[Tensor] = None,
    use_shape: bool,
    spot_batch_size: Optional[int] = None,
    peak_chunk_size: Optional[int] = None,
    total_likelihood: TotalLikelihood = "negative_binomial",
    prior_shape: float = DEFAULT_PRIOR_SHAPE,
    prior_rate: float = DEFAULT_PRIOR_RATE,
    eps: float = DEFAULT_EPS,
) -> ObjectiveComponents:
    """Transform ``raw_z`` and aggregate the chunked complete objective."""
    return chunked_objective_components_from_abundance(
        positive_abundance(raw_z, eps=eps),
        accessibility,
        observed_counts,
        phi_ref,
        omega=omega,
        use_shape=use_shape,
        spot_batch_size=spot_batch_size,
        peak_chunk_size=peak_chunk_size,
        total_likelihood=total_likelihood,
        prior_shape=prior_shape,
        prior_rate=prior_rate,
        eps=eps,
    )


def factorized_poisson_log_likelihood(
    shape_counts: Tensor,
    bin_rates: Tensor,
    eps: float = DEFAULT_EPS,
) -> Tensor:
    """Return Poisson-total plus conditional-multinomial log likelihood."""
    if shape_counts.shape != bin_rates.shape:
        raise ValueError("shape_counts and bin_rates must have identical shapes.")
    total_counts = shape_counts.sum(dim=-1)
    total_rates = bin_rates.sum(dim=-1)
    return poisson_log_prob(total_counts, total_rates, eps=eps).sum() + (
        conditional_multinomial_log_prob(shape_counts, bin_rates, eps=eps).sum()
    )


def independent_poisson_bin_log_likelihood(
    shape_counts: Tensor,
    bin_rates: Tensor,
    eps: float = DEFAULT_EPS,
) -> Tensor:
    """Return the independent-bin Poisson oracle used only for golden tests."""
    if shape_counts.shape != bin_rates.shape:
        raise ValueError("shape_counts and bin_rates must have identical shapes.")
    return poisson_log_prob(shape_counts, bin_rates, eps=eps).sum()
