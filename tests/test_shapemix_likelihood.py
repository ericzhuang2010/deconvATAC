import json
import math

import pytest
import torch

from deconvatac.shapemix.diagnostics import (
    FitRecord,
    OptimizationStepRecord,
    ReconstructionAccumulator,
    RestartRecord,
)
from deconvatac.shapemix.likelihood import (
    chunked_objective_components_from_raw,
    conditional_multinomial_log_prob,
    expected_bin_rates,
    expected_total_rate,
    factorized_poisson_log_likelihood,
    gamma_log_prior,
    independent_poisson_bin_log_likelihood,
    likelihood_components,
    negative_binomial_log_prob,
    objective_components_from_raw,
    positive_abundance,
)


def _toy_inputs():
    raw_z = torch.tensor(
        [[0.2, -0.7], [1.1, 0.3], [-0.4, 0.9]], dtype=torch.float32
    )
    accessibility = torch.tensor(
        [[0.8, 1.4, 0.5, 2.0], [1.6, 0.4, 1.2, 0.7]], dtype=torch.float32
    )
    omega = torch.tensor(
        [
            [
                [0.70, 0.20, 0.10],
                [0.20, 0.30, 0.50],
                [0.55, 0.25, 0.20],
                [0.15, 0.70, 0.15],
            ],
            [
                [0.10, 0.35, 0.55],
                [0.65, 0.25, 0.10],
                [0.20, 0.30, 0.50],
                [0.50, 0.20, 0.30],
            ],
        ],
        dtype=torch.float32,
    )
    observed = torch.tensor(
        [
            [[2, 0, 1], [0, 0, 0], [1, 1, 0], [0, 3, 0]],
            [[0, 1, 2], [4, 0, 0], [0, 0, 0], [1, 0, 1]],
            [[1, 0, 0], [0, 2, 1], [2, 1, 2], [0, 0, 1]],
        ],
        dtype=torch.float32,
    )
    return raw_z, accessibility, omega, observed


def test_expected_rate_conservation():
    raw_z, accessibility, omega, observed = _toy_inputs()
    z = positive_abundance(raw_z)
    total_rates = expected_total_rate(z, accessibility)
    bin_rates = expected_bin_rates(z, accessibility, omega)

    assert total_rates.dtype == torch.float32
    assert bin_rates.shape == (*total_rates.shape, 3)
    torch.testing.assert_close(bin_rates.sum(dim=-1), total_rates)
    assert torch.equal(observed.sum(dim=-1), observed.sum(dim=-1))


def test_negative_binomial_parameterization_and_abundance_scaling():
    counts = torch.tensor([[0.0, 2.0], [4.0, 1.0]], dtype=torch.float32)
    mean = torch.tensor([[0.7, 2.3], [4.1, 1.2]], dtype=torch.float32)
    inverse_dispersion = torch.tensor([[3.0], [7.0]], dtype=torch.float32)

    actual = negative_binomial_log_prob(counts, mean, inverse_dispersion)
    distribution = torch.distributions.NegativeBinomial(
        total_count=inverse_dispersion,
        logits=torch.log(mean) - torch.log(inverse_dispersion),
    )
    torch.testing.assert_close(actual, distribution.log_prob(counts), atol=2e-6, rtol=2e-6)
    torch.testing.assert_close(distribution.mean, mean)
    torch.testing.assert_close(
        distribution.variance, mean + mean.square() / inverse_dispersion
    )

    z = torch.tensor([[1.0, 2.0], [3.0, 4.0]], dtype=torch.float32)
    accessibility = torch.tensor([[0.5, 0.8], [1.2, 0.4]], dtype=torch.float32)
    totals = torch.tensor([[1.0, 2.0], [3.0, 0.0]], dtype=torch.float32)
    phi_ref = torch.tensor(5.0)
    components = likelihood_components(
        z,
        totals,
        accessibility,
        phi_ref,
        use_shape=False,
    )
    expected = negative_binomial_log_prob(
        totals,
        expected_total_rate(z, accessibility),
        z.sum(dim=-1, keepdim=True) * phi_ref,
    ).sum()
    torch.testing.assert_close(components.count_log_likelihood, expected)


def test_gamma_shape_rate_prior_matches_torch():
    z = torch.tensor([[0.3, 1.4], [2.2, 0.8]], dtype=torch.float32)
    expected = torch.distributions.Gamma(
        torch.tensor(2.0), torch.tensor(1.0)
    ).log_prob(z).sum()
    torch.testing.assert_close(gamma_log_prior(z), expected)


def test_poisson_factorization_equals_independent_bins():
    shape_counts = torch.tensor(
        [[[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]], [[0.0, 0.0, 0.0], [2.0, 1.0, 1.0]]],
        dtype=torch.float32,
    )
    bin_rates = torch.tensor(
        [[[0.8, 0.0, 1.7], [0.0, 2.4, 0.0]], [[0.3, 0.2, 0.1], [1.5, 0.7, 0.9]]],
        dtype=torch.float32,
    )
    factorized = factorized_poisson_log_likelihood(shape_counts, bin_rates)
    independent = independent_poisson_bin_log_likelihood(shape_counts, bin_rates)
    torch.testing.assert_close(factorized, independent, atol=2e-6, rtol=2e-6)


def test_conditional_tiny_and_zero_probability_support_is_exact():
    tiny = torch.tensor([[1.0e-20, 1.0]], dtype=torch.float32)
    tiny_observed = torch.tensor([[1.0, 0.0]], dtype=torch.float32)
    tiny_result = conditional_multinomial_log_prob(tiny_observed, tiny)
    expected = math.log(1.0e-20 / (1.0 + 1.0e-20))
    assert tiny_result.item() == pytest.approx(expected, rel=2e-6)

    rates = torch.tensor(
        [[0.0, 2.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]],
        dtype=torch.float32,
    )
    counts = torch.tensor(
        [[0.0, 3.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 0.0]],
        dtype=torch.float32,
    )
    result = conditional_multinomial_log_prob(counts, rates)
    assert result[0].item() == 0.0
    assert torch.isneginf(result[1])
    assert result[2].item() == 0.0


def test_identical_omega_has_no_abundance_information():
    raw_z, accessibility, _, observed = _toy_inputs()
    shared = torch.tensor(
        [
            [0.2, 0.3, 0.5],
            [0.6, 0.1, 0.3],
            [0.1, 0.7, 0.2],
            [0.4, 0.4, 0.2],
        ],
        dtype=torch.float32,
    )
    omega = shared.unsqueeze(0).repeat(2, 1, 1)
    first_raw = raw_z.clone().requires_grad_(True)
    second_raw = (raw_z + 1.7).clone().requires_grad_(True)

    first_z = positive_abundance(first_raw)
    second_z = positive_abundance(second_raw)
    first_shape = likelihood_components(
        first_z,
        observed.sum(dim=-1),
        accessibility,
        torch.tensor(4.0),
        shape_counts=observed,
        omega=omega,
        use_shape=True,
    ).shape_log_likelihood
    second_shape = likelihood_components(
        second_z,
        observed.sum(dim=-1),
        accessibility,
        torch.tensor(4.0),
        shape_counts=observed,
        omega=omega,
        use_shape=True,
    ).shape_log_likelihood

    torch.testing.assert_close(first_shape, second_shape, atol=0, rtol=0)
    gradient = torch.autograd.grad(first_shape, first_raw)[0]
    assert torch.count_nonzero(gradient).item() == 0


def test_one_bin_conditional_is_exact_zero():
    raw_z = torch.tensor([[0.1, -0.4], [0.7, 0.3]], dtype=torch.float32)
    accessibility = torch.tensor([[1.0, 0.4], [0.7, 1.2]], dtype=torch.float32)
    observed = torch.tensor([[[2.0], [0.0]], [[1.0], [3.0]]], dtype=torch.float32)
    omega = torch.ones((2, 2, 1), dtype=torch.float32)

    shape = objective_components_from_raw(
        raw_z,
        accessibility,
        observed,
        torch.tensor(3.0),
        omega=omega,
        use_shape=True,
    )
    count = objective_components_from_raw(
        raw_z,
        accessibility,
        observed,
        torch.tensor(3.0),
        omega=omega,
        use_shape=False,
    )
    assert shape.shape_log_likelihood.item() == 0.0
    assert shape.total_log_objective.item() == count.total_log_objective.item()


def test_float32_objective_has_finite_gradients_at_large_dispersion():
    raw_z, accessibility, omega, observed = _toy_inputs()
    raw_z = raw_z.requires_grad_(True)
    components = objective_components_from_raw(
        raw_z,
        accessibility,
        observed,
        torch.tensor(1.0e8, dtype=torch.float32),
        omega=omega,
        use_shape=True,
    )
    assert components.total_log_objective.dtype == torch.float32
    assert torch.isfinite(components.total_log_objective)
    (-components.total_log_objective).backward()
    assert raw_z.grad is not None
    assert torch.isfinite(raw_z.grad).all()


@pytest.mark.parametrize("use_shape", [False, True])
def test_chunked_objective_and_gradient_parity(use_shape):
    raw_z, accessibility, omega, observed = _toy_inputs()
    full_raw = raw_z.clone().requires_grad_(True)
    chunked_raw = raw_z.clone().requires_grad_(True)

    full = objective_components_from_raw(
        full_raw,
        accessibility,
        observed,
        torch.tensor(6.5),
        omega=omega,
        use_shape=use_shape,
    )
    chunked = chunked_objective_components_from_raw(
        chunked_raw,
        accessibility,
        observed,
        torch.tensor(6.5),
        omega=omega,
        use_shape=use_shape,
        spot_batch_size=2,
        peak_chunk_size=2,
    )

    for field in (
        "count_log_likelihood",
        "shape_log_likelihood",
        "abundance_log_prior",
        "total_log_objective",
    ):
        torch.testing.assert_close(
            getattr(full, field), getattr(chunked, field), atol=1e-5, rtol=1e-6
        )
    full.total_log_objective.backward()
    chunked.total_log_objective.backward()
    torch.testing.assert_close(full_raw.grad, chunked_raw.grad, atol=2e-5, rtol=2e-5)


def test_reconstruction_accumulator_and_fit_records_are_json_safe():
    observed_bins = torch.tensor(
        [[[1.0, 0.0], [0.0, 2.0], [1.0, 1.0]], [[0.0, 0.0], [2.0, 1.0], [0.0, 3.0]]]
    )
    expected_bins = torch.tensor(
        [[[0.8, 0.2], [0.3, 1.5], [1.4, 0.7]], [[0.1, 0.2], [1.8, 1.3], [0.4, 2.5]]]
    )
    observed_totals = observed_bins.sum(dim=-1)
    expected_totals = expected_bins.sum(dim=-1)

    accumulator = ReconstructionAccumulator(n_bins=2)
    accumulator.update(
        observed_totals[:, :1],
        expected_totals[:, :1],
        observed_bins=observed_bins[:, :1],
        expected_bins=expected_bins[:, :1],
    )
    accumulator.update(
        observed_totals[:, 1:],
        expected_totals[:, 1:],
        observed_bins=observed_bins[:, 1:],
        expected_bins=expected_bins[:, 1:],
    )
    reconstruction = accumulator.finalize()
    assert reconstruction.count_entries == observed_totals.numel()
    assert reconstruction.observed_total == pytest.approx(observed_totals.sum().item())
    assert reconstruction.expected_total == pytest.approx(expected_totals.sum().item())
    assert reconstruction.observed_by_bin == pytest.approx(
        observed_bins.sum(dim=(0, 1)).tolist()
    )

    history = OptimizationStepRecord(
        step=1,
        count_log_likelihood=-11.0,
        shape_log_likelihood=-2.5,
        abundance_log_prior=-3.2,
        total_log_objective=-16.7,
    )
    restart = RestartRecord(
        restart_index=0,
        seed_tuple=(20260822, 1103, 101, 0, 0),
        converged=True,
        stopping_reason="tolerance",
        steps=37,
        count_log_likelihood=-10.0,
        shape_log_likelihood=-2.0,
        abundance_log_prior=-3.0,
        total_log_objective=-15.0,
        history=(history,),
        runtime_seconds=float("nan"),
    )
    fit = FitRecord(
        use_shape=True,
        success=True,
        selected_restart=0,
        stopping_reason="selected_finite_converged_restart",
        n_spots=2,
        n_peaks=3,
        n_bins=2,
        n_cell_types=2,
        dtype="float32",
        device="cpu",
        restarts=(restart,),
    )
    assert fit.selected_restart_record is restart
    payload = {"fit": fit.to_dict(), "reconstruction": reconstruction.to_dict()}
    assert payload["fit"]["restarts"][0]["runtime_seconds"] is None
    assert payload["fit"]["restarts"][0]["history"][0]["step"] == 1
    json.dumps(payload, allow_nan=False)
