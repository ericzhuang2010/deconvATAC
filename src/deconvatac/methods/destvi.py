from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from deconvatac.data import DeconvolutionInput, DeconvolutionResult, normalize_proportions

from .base import BaseDeconvolver


def read_destvi_proportions(path: Union[str, Path], zero_policy: str = "zeros") -> pd.DataFrame:
    """Read DestVI predictions in standardized proportion form."""
    values = pd.read_csv(path, index_col=0)
    values.index = values.index.astype(str)
    return normalize_proportions(values, zero_policy=zero_policy)


class DestVIDeconvolver(BaseDeconvolver):
    """Adapter for the existing DestVI wrapper."""

    method_name = "destvi"

    def run(self, data: DeconvolutionInput) -> DeconvolutionResult:
        from deconvatac.tl.destvi import destvi

        output_dir = Path(data.output_dir) if data.output_dir is not None else None
        temp_dir: Optional[tempfile.TemporaryDirectory[str]] = None
        if output_dir is None:
            temp_dir = tempfile.TemporaryDirectory()
            native_dir = Path(temp_dir.name)
        else:
            native_dir = output_dir / "results" / "raw_method_output"
            native_dir.mkdir(parents=True, exist_ok=True)

        try:
            params = {
                "labels_key": self.config.get("labels_key", data.labels_key),
                "layer_spatial": self.config.get("layer_spatial"),
                "layer_ref": self.config.get("layer_ref"),
                "use_gpu": self.config.get("use_gpu", True),
                "max_epochs_spatial": self.config.get("max_epochs_spatial", 2000),
                "max_epochs_ref": self.config.get("max_epochs_ref", 300),
                "return_adatas": False,
                "plots": self.config.get("plots", False),
                "results_path": str(native_dir),
                "model_ref_kwargs": self.config.get("model_ref_kwargs", {}),
                "train_ref_kwargs": self.config.get("train_ref_kwargs", {}),
                "model_spatial_kwargs": self.config.get("model_spatial_kwargs", {}),
                "train_spatial_kwargs": self.config.get("train_spatial_kwargs", {}),
            }
            destvi(data.spatial, data.reference, **params)

            result_path = native_dir / "predicted_proportions.csv"
            zero_policy = self.config.get("zero_policy", "zeros")
            proportions = read_destvi_proportions(result_path, zero_policy=zero_policy)
            diagnostics = {
                "method": self.method_name,
                "params": {key: value for key, value in params.items() if key != "results_path"},
                "native_output_dir": str(native_dir),
                "native_files": {"predicted_proportions": str(result_path)},
            }
            return DeconvolutionResult(
                method=self.method_name,
                dataset_id=data.dataset_id,
                modality=data.modality,
                feature_set=data.feature_set,
                proportions=proportions,
                diagnostics=diagnostics,
            )
        finally:
            if temp_dir is not None:
                temp_dir.cleanup()
