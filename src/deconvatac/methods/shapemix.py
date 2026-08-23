"""Maintained-runner adapter for the fixed-signature ShapeMix MAP model."""

from __future__ import annotations

import json
import math
from importlib.metadata import PackageNotFoundError, version
from numbers import Integral
from pathlib import Path
from typing import Any, Mapping, Optional

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from deconvatac.data import (
    DeconvolutionInput,
    DeconvolutionResult,
    validate_deconvolution_input,
)
from deconvatac.shapemix.config import ShapeMixConfig

from .base import BaseDeconvolver


_MODEL_VERSION = 1
_NATIVE_FILENAMES = {
    "training_history": "training_history.csv",
    "restart_summary": "restart_summary.csv",
    "reconstruction_summary": "reconstruction_summary.csv",
    "residuals_by_peak_and_bin": "residuals_by_peak_and_bin.csv.gz",
    "signature_summary": "signature_summary.yaml",
}


def _require_seed_metadata(data: DeconvolutionInput) -> tuple[int, int]:
    metadata = data.metadata
    if not isinstance(metadata, Mapping):
        raise ValueError("ShapeMix requires mapping-valued dataset metadata.")
    dataset_config = metadata.get("dataset_config")
    if not isinstance(dataset_config, Mapping):
        raise ValueError(
            "ShapeMix requires data.metadata['dataset_config'] with simulation seeds."
        )
    configured_id = dataset_config.get("dataset_id")
    if configured_id is not None and configured_id != data.dataset_id:
        raise ValueError("dataset_config.dataset_id does not match data.dataset_id.")
    simulation = dataset_config.get("simulation")
    if not isinstance(simulation, Mapping):
        raise ValueError(
            "ShapeMix requires dataset_config.simulation with outer and inner seeds."
        )

    seeds = []
    for name in ("outer_split_seed", "inner_mixture_seed"):
        value = simulation.get(name)
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise ValueError(f"dataset_config.simulation.{name} must be an integer.")
        normalized = int(value)
        if normalized < 0:
            raise ValueError(
                f"dataset_config.simulation.{name} must be nonnegative."
            )
        seeds.append(normalized)
    return seeds[0], seeds[1]


def _dense_slice(matrix: Any, spot_slice: slice, peak_slice: slice) -> np.ndarray:
    selected = matrix[spot_slice, peak_slice]
    if sparse.issparse(selected):
        selected = selected.toarray()
    return np.asarray(selected, dtype=np.float64)


def _reconstruction_diagnostics(
    data: DeconvolutionInput,
    abundance: np.ndarray,
    signatures: Any,
    config: ShapeMixConfig,
    *,
    retain_peak_bin_residuals: bool,
) -> tuple[Any, Optional[pd.DataFrame]]:
    """Reduce fitted expectations chunkwise without retaining ``S x P x B``."""
    from deconvatac.shapemix.diagnostics import ReconstructionAccumulator

    if data.fragment_shape is None:
        raise ValueError("fragment_shape is required for reconstruction diagnostics.")
    layer_names = data.fragment_shape.layer_names
    bin_names = tuple(bin_spec.name for bin_spec in data.fragment_shape.bins)
    layers = tuple(data.spatial.layers[layer_name] for layer_name in layer_names)
    n_spots, n_peaks = data.spatial.shape
    n_bins = len(layer_names)
    accumulator = ReconstructionAccumulator(n_bins=n_bins)

    observed_peak_bin = expected_peak_bin = None
    absolute_error_peak_bin = squared_error_peak_bin = None
    observed_nonzero_peak_bin = None
    if retain_peak_bin_residuals:
        array_shape = (n_peaks, n_bins)
        observed_peak_bin = np.zeros(array_shape, dtype=np.float64)
        expected_peak_bin = np.zeros(array_shape, dtype=np.float64)
        absolute_error_peak_bin = np.zeros(array_shape, dtype=np.float64)
        squared_error_peak_bin = np.zeros(array_shape, dtype=np.float64)
        observed_nonzero_peak_bin = np.zeros(array_shape, dtype=np.int64)

    for spot_start in range(0, n_spots, config.spot_batch_size):
        spot_stop = min(spot_start + config.spot_batch_size, n_spots)
        spot_slice = slice(spot_start, spot_stop)
        z_chunk = abundance[spot_slice]
        for peak_start in range(0, n_peaks, config.peak_chunk_size):
            peak_stop = min(peak_start + config.peak_chunk_size, n_peaks)
            peak_slice = slice(peak_start, peak_stop)
            observed_bins = np.stack(
                [
                    _dense_slice(layer, spot_slice, peak_slice)
                    for layer in layers
                ],
                axis=-1,
            )
            accessibility = signatures.A[:, peak_slice]
            omega = signatures.omega[:, peak_slice, :]
            weighted_shape = accessibility[:, :, None] * omega
            expected_bins = np.einsum(
                "sc,cpb->spb", z_chunk, weighted_shape, optimize=True
            )
            observed_totals = observed_bins.sum(axis=-1, dtype=np.float64)
            expected_totals = z_chunk @ accessibility
            if not np.allclose(
                expected_bins.sum(axis=-1, dtype=np.float64),
                expected_totals,
                rtol=1.0e-10,
                atol=1.0e-10,
            ):
                raise RuntimeError(
                    "Fitted bin expectations do not conserve expected peak totals."
                )
            accumulator.update(
                observed_totals,
                expected_totals,
                observed_bins=observed_bins,
                expected_bins=expected_bins,
            )

            if observed_peak_bin is not None:
                residual = observed_bins - expected_bins
                observed_peak_bin[peak_slice] += observed_bins.sum(axis=0)
                expected_peak_bin[peak_slice] += expected_bins.sum(axis=0)
                absolute_error_peak_bin[peak_slice] += np.abs(residual).sum(axis=0)
                squared_error_peak_bin[peak_slice] += np.square(residual).sum(axis=0)
                observed_nonzero_peak_bin[peak_slice] += (observed_bins > 0).sum(
                    axis=0
                )

    summary = accumulator.finalize()
    if observed_peak_bin is None:
        return summary, None

    residual_frame = pd.DataFrame(
        {
            "peak_id": np.repeat(
                data.spatial.var_names.astype(str).to_numpy(), n_bins
            ),
            "bin_name": np.tile(np.asarray(bin_names, dtype=object), n_peaks),
            "layer_name": np.tile(np.asarray(layer_names, dtype=object), n_peaks),
            "observed": observed_peak_bin.reshape(-1),
            "expected": expected_peak_bin.reshape(-1),
            "residual": (observed_peak_bin - expected_peak_bin).reshape(-1),
            "absolute_error_sum": absolute_error_peak_bin.reshape(-1),
            "squared_error_sum": squared_error_peak_bin.reshape(-1),
            "observed_nonzero_spots": observed_nonzero_peak_bin.reshape(-1),
        }
    )
    return summary, residual_frame


def _training_history_frame(fit: Any) -> pd.DataFrame:
    rows = []
    for restart in fit.restart_diagnostics:
        for record in restart.history:
            rows.append(
                {
                    "restart_index": restart.restart_index,
                    "selected": restart.restart_index == fit.selected_restart,
                    **record.to_dict(),
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "restart_index",
            "selected",
            "step",
            "count_log_likelihood",
            "shape_log_likelihood",
            "abundance_log_prior",
            "total_log_objective",
        ],
    )


def _restart_summary_frame(fit: Any) -> pd.DataFrame:
    rows = []
    for restart in fit.restart_diagnostics:
        row = restart.to_dict()
        history = row.pop("history")
        row["history_rows"] = len(history)
        row["selected"] = restart.restart_index == fit.selected_restart
        row["seed_tuple"] = json.dumps(row["seed_tuple"], separators=(",", ":"))
        row["nonfinite_events"] = json.dumps(
            row["nonfinite_events"], separators=(",", ":")
        )
        rows.append(row)
    return pd.DataFrame(rows)


def _reconstruction_summary_frame(summary: Any, data: DeconvolutionInput) -> pd.DataFrame:
    if data.fragment_shape is None:
        raise ValueError("fragment_shape is required for reconstruction diagnostics.")
    rows = [
        {
            "component": "total_count",
            "bin_name": None,
            "layer_name": None,
            "entries": summary.count_entries,
            "observed_nonzero_entries": summary.observed_nonzero_entries,
            "expected_positive_entries": summary.expected_positive_entries,
            "observed_total": summary.observed_total,
            "expected_total": summary.expected_total,
            "residual_total": summary.residual_total,
            "absolute_error_sum": summary.absolute_error_sum,
            "squared_error_sum": summary.squared_error_sum,
            "mean_absolute_error": summary.mean_absolute_error,
            "root_mean_squared_error": summary.root_mean_squared_error,
        }
    ]
    for bin_index, bin_spec in enumerate(data.fragment_shape.bins):
        entries = summary.count_entries
        absolute_error = summary.absolute_error_by_bin[bin_index]
        squared_error = summary.squared_error_by_bin[bin_index]
        rows.append(
            {
                "component": "shape_bin",
                "bin_name": bin_spec.name,
                "layer_name": bin_spec.layer,
                "entries": entries,
                "observed_nonzero_entries": None,
                "expected_positive_entries": None,
                "observed_total": summary.observed_by_bin[bin_index],
                "expected_total": summary.expected_by_bin[bin_index],
                "residual_total": summary.residual_by_bin[bin_index],
                "absolute_error_sum": absolute_error,
                "squared_error_sum": squared_error,
                "mean_absolute_error": absolute_error / entries,
                "root_mean_squared_error": math.sqrt(squared_error / entries),
            }
        )
    return pd.DataFrame(rows)


def _compact_fit_diagnostics(fit: Any) -> dict[str, Any]:
    payload = fit.to_diagnostics_dict()
    for restart in payload["restarts"]:
        history = restart.pop("history")
        restart["history_rows"] = len(history)
    return payload


def _package_version(package: str) -> Optional[str]:
    try:
        return version(package)
    except PackageNotFoundError:
        return None


def _write_native_outputs(
    output_dir: Path,
    data: DeconvolutionInput,
    fit: Any,
    reconstruction: Any,
    residuals: pd.DataFrame,
    signature_summary: dict[str, Any],
) -> dict[str, str]:
    native_dir = output_dir / "results" / "raw_method_output"
    native_dir.mkdir(parents=True, exist_ok=True)
    _training_history_frame(fit).to_csv(
        native_dir / _NATIVE_FILENAMES["training_history"], index=False
    )
    _restart_summary_frame(fit).to_csv(
        native_dir / _NATIVE_FILENAMES["restart_summary"], index=False
    )
    _reconstruction_summary_frame(reconstruction, data).to_csv(
        native_dir / _NATIVE_FILENAMES["reconstruction_summary"], index=False
    )
    residuals.to_csv(
        native_dir / _NATIVE_FILENAMES["residuals_by_peak_and_bin"],
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    with (native_dir / _NATIVE_FILENAMES["signature_summary"]).open("w") as handle:
        yaml.safe_dump(signature_summary, handle, sort_keys=False)
    return {
        name: f"results/raw_method_output/{filename}"
        for name, filename in _NATIVE_FILENAMES.items()
    }


class ShapeMixDeconvolver(BaseDeconvolver):
    """ATAC-only adapter for both nested fixed-signature ShapeMix arms."""

    method_name = "shapemix"

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self.model_config = ShapeMixConfig.from_mapping(kwargs)

    def run(self, data: DeconvolutionInput) -> DeconvolutionResult:
        """Validate shape-aware ATAC input, fit MAP, and emit compact outputs."""
        if data.modality != "atac":
            raise ValueError("ShapeMix supports only ATAC inputs.")
        if data.fragment_shape is None:
            raise ValueError("ShapeMix requires declared fragment_shape metadata and layers.")
        if data.cell_types is None:
            raise ValueError("ShapeMix requires an ordered cell_types declaration.")
        validate_deconvolution_input(data)
        outer_split_seed, inner_mixture_seed = _require_seed_metadata(data)

        # Import the Torch-dependent core only after ShapeMix is selected.  The
        # registry can list or run every other method without importing Torch.
        from deconvatac.shapemix.map import fit_shapemix_map
        from deconvatac.shapemix.signatures import estimate_reference_signatures

        layer_names = data.fragment_shape.layer_names
        bin_names = tuple(bin_spec.name for bin_spec in data.fragment_shape.bins)
        signatures = estimate_reference_signatures(
            data.reference,
            data.labels_key,
            data.cell_types,
            outer_split_seed,
            config=self.model_config,
            layer_names=layer_names,
            bin_names=bin_names,
        )
        fit = fit_shapemix_map(
            tuple(data.spatial.layers[layer_name] for layer_name in layer_names),
            signatures=signatures,
            config=self.model_config,
            outer_split_seed=outer_split_seed,
            inner_mixture_seed=inner_mixture_seed,
            spot_names=data.spatial.obs_names.astype(str).tolist(),
            feature_names=data.spatial.var_names.astype(str).tolist(),
            cell_types=data.cell_types,
        )

        abundance = pd.DataFrame(
            fit.abundance,
            index=data.spatial.obs_names.copy(),
            columns=pd.Index(data.cell_types, dtype=object),
        )
        proportions = pd.DataFrame(
            fit.proportions,
            index=data.spatial.obs_names.copy(),
            columns=pd.Index(data.cell_types, dtype=object),
        )
        if (
            not np.isfinite(abundance.to_numpy()).all()
            or (abundance.to_numpy() <= 0).any()
            or not np.isfinite(proportions.to_numpy()).all()
            or (proportions.to_numpy() < 0).any()
            or not np.allclose(
                proportions.sum(axis=1).to_numpy(), 1.0, rtol=0.0, atol=1.0e-6
            )
        ):
            raise RuntimeError("ShapeMix returned invalid abundance or proportion rows.")

        should_write = data.output_dir is not None
        reconstruction, residuals = _reconstruction_diagnostics(
            data,
            fit.abundance,
            signatures,
            self.model_config,
            retain_peak_bin_residuals=should_write,
        )
        signature_summary = {
            "schema_version": 1,
            "method": self.method_name,
            "model_version": _MODEL_VERSION,
            "source": "training_reference_only",
            "outer_split_seed": outer_split_seed,
            "signature": signatures.to_metadata(),
            "fragment_shape": data.fragment_shape.to_dict(omit_none=True),
            "config_sha256": self.model_config.content_sha256,
            "config": self.model_config.to_dict(),
        }
        native_outputs = None
        if should_write:
            assert residuals is not None
            native_outputs = _write_native_outputs(
                Path(data.output_dir),
                data,
                fit,
                reconstruction,
                residuals,
                signature_summary,
            )

        diagnostics = {
            "method": self.method_name,
            "model_version": _MODEL_VERSION,
            "use_shape": self.model_config.use_shape,
            "config": self.model_config.to_dict(),
            "config_sha256": self.model_config.content_sha256,
            "seeds": {
                "outer_split_seed": outer_split_seed,
                "inner_mixture_seed": inner_mixture_seed,
                "method_seed": self.model_config.seed,
            },
            "fragment_shape": data.fragment_shape.to_dict(omit_none=True),
            "signature": signatures.to_metadata(),
            "fit": _compact_fit_diagnostics(fit),
            "reconstruction": reconstruction.to_dict(),
            "software_versions": {
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "torch": _package_version("torch"),
                "pysam": _package_version("pysam"),
            },
            "native_outputs": native_outputs,
        }
        return DeconvolutionResult(
            method=self.method_name,
            dataset_id=data.dataset_id,
            modality=data.modality,
            feature_set=data.feature_set,
            proportions=proportions,
            abundance=abundance,
            diagnostics=diagnostics,
            output_dir=Path(data.output_dir) if data.output_dir is not None else None,
        )
