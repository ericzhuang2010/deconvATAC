from __future__ import annotations

import numpy as np
import pytest
import torch

from deconvatac.shapemix import map as map_module
from deconvatac.shapemix.config import ShapeMixConfig
from deconvatac.shapemix.likelihood import likelihood_components, positive_abundance
from deconvatac.shapemix.map import fit_shapemix_map


pytestmark = pytest.mark.skipif(
    not torch.cuda.is_available(),
    reason="CUDA hardware is required for ShapeMix GPU qualification.",
)


def _toy() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    accessibility = np.asarray(
        [
            [30.0, 5.0, 25.0, 4.0, 20.0, 6.0],
            [5.0, 30.0, 4.0, 25.0, 6.0, 20.0],
        ],
        dtype=np.float64,
    )
    omega = np.empty((2, 6, 3), dtype=np.float64)
    omega[0, :, :] = (0.75, 0.15, 0.10)
    omega[1, :, :] = (0.10, 0.15, 0.75)
    abundance = np.asarray([[8.0, 2.0], [2.0, 8.0], [6.0, 4.0]])
    counts = np.rint(
        np.einsum("sc,cp,cpb->spb", abundance, accessibility, omega)
    ).astype(np.int64)
    return counts, accessibility, omega


def _config(device: str, cuda_count_cache: str = "auto") -> ShapeMixConfig:
    return ShapeMixConfig(
        use_shape=True,
        learning_rate=0.05,
        max_steps=500,
        patience=25,
        tolerance=5.0e-5,
        restarts=1,
        spot_batch_size=2,
        peak_chunk_size=3,
        device=device,
        cuda_count_cache=cuda_count_cache,
    )


def _protocol_config(device: str) -> ShapeMixConfig:
    return ShapeMixConfig(
        use_shape=True,
        restarts=1,
        spot_batch_size=2,
        peak_chunk_size=3,
        device=device,
    )


def test_cuda_objective_and_gradient_match_cpu() -> None:
    counts, accessibility, omega = _toy()
    raw_cpu = torch.tensor([[1.5, 0.5], [0.5, 1.5]], requires_grad=True)
    raw_cuda = raw_cpu.detach().to("cuda").requires_grad_(True)

    def evaluate(raw_z: torch.Tensor):
        device = raw_z.device
        z = positive_abundance(raw_z)
        count_tensor = torch.as_tensor(counts[:2], dtype=torch.float32, device=device)
        access_tensor = torch.as_tensor(accessibility, dtype=torch.float32, device=device)
        omega_tensor = torch.as_tensor(omega, dtype=torch.float32, device=device)
        components = likelihood_components(
            z,
            count_tensor.sum(dim=-1),
            access_tensor,
            100.0,
            shape_counts=count_tensor,
            omega=omega_tensor,
            use_shape=True,
        )
        loss = -components.total_log_objective
        loss.backward()
        return float(loss.detach().cpu()), raw_z.grad.detach().cpu().numpy()

    cpu_value, cpu_gradient = evaluate(raw_cpu)
    cuda_value, cuda_gradient = evaluate(raw_cuda)
    assert cuda_value == pytest.approx(cpu_value, rel=1.0e-5, abs=1.0e-4)
    np.testing.assert_allclose(cuda_gradient, cpu_gradient, rtol=1.0e-5, atol=1.0e-4)


def test_cuda_map_matches_cpu_and_repeats_deterministically() -> None:
    counts, accessibility, omega = _toy()
    common = {
        "outer_split_seed": 0,
        "inner_mixture_seed": 0,
        "spot_names": ("spot_b", "spot_a", "spot_c"),
        "feature_names": ("peak_3", "peak_1", "peak_5", "peak_0", "peak_4", "peak_2"),
        "cell_types": ("type_0", "type_1"),
    }
    cpu = fit_shapemix_map(
        counts, accessibility, omega, 100.0, config=_protocol_config("cpu"), **common
    )
    cuda = fit_shapemix_map(
        counts,
        accessibility,
        omega,
        100.0,
        config=_protocol_config("cuda:0"),
        **common,
    )
    repeat = fit_shapemix_map(
        counts,
        accessibility,
        omega,
        100.0,
        config=_protocol_config("cuda:0"),
        **common,
    )

    # CPU and CUDA can select adjacent checkpoints on the same flat,
    # deterministic convergence path because their transcendental kernels
    # round differently. Keep the toy-map numerical bound narrow while also
    # requiring the restart and convergence length to agree exactly.
    np.testing.assert_allclose(cuda.proportions, cpu.proportions, rtol=0.0, atol=2.0e-5)
    assert cuda.restart_diagnostics[0].steps == cpu.restart_diagnostics[0].steps
    np.testing.assert_allclose(repeat.proportions, cuda.proportions, rtol=0.0, atol=1.0e-7)
    assert repeat.selected_restart == cuda.selected_restart
    execution = cuda.to_diagnostics_dict()["execution"]
    assert execution["count_cache_mode"] == "full_cuda"
    assert execution["count_cache_bytes"] == counts.size * 4
    assert execution["peak_device_memory_allocated_bytes"] > 0
    assert execution["device_index"] == 0
    assert execution["device_name"]
    assert execution["device_compute_capability"] == "8.6"
    assert execution["cuda_runtime_version"]
    assert execution["deterministic_algorithms"] is True


def test_cuda_streamed_cache_matches_full_cache() -> None:
    counts, accessibility, omega = _toy()
    common = {
        "outer_split_seed": 0,
        "inner_mixture_seed": 0,
        "spot_names": ("spot_b", "spot_a", "spot_c"),
        "feature_names": ("peak_3", "peak_1", "peak_5", "peak_0", "peak_4", "peak_2"),
        "cell_types": ("type_0", "type_1"),
    }
    cached = fit_shapemix_map(
        counts, accessibility, omega, 100.0, config=_config("cuda:0"), **common
    )
    streamed = fit_shapemix_map(
        counts, accessibility, omega, 100.0, config=_config("cuda:0", "disabled"), **common
    )

    np.testing.assert_allclose(
        streamed.proportions, cached.proportions, rtol=0.0, atol=1.0e-7
    )
    assert streamed.to_diagnostics_dict()["execution"]["count_cache_mode"] == (
        "streamed_host_chunks"
    )


def test_cuda_request_fails_closed_when_cuda_is_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    with pytest.raises(RuntimeError, match="torch.cuda.is_available"):
        map_module._resolve_device(ShapeMixConfig(device="cuda:0"))


def test_cuda_cache_oom_falls_back_to_streaming(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    counts, _, _ = _toy()
    source = map_module._CountSource(counts)

    def fail_cache_allocation(*args, **kwargs):
        raise torch.OutOfMemoryError("qualification test")

    monkeypatch.setattr(torch, "empty", fail_cache_allocation)
    cache, metadata = map_module._prepare_count_cache(
        source,
        torch.device("cuda"),
        _config("cuda:0"),
    )

    assert cache is None
    assert metadata["count_cache_mode"] == "streamed_after_cache_oom"
    assert metadata["count_cache_bytes"] == 0
