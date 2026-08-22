from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from deconvatac.data import DeconvolutionInput, DeconvolutionResult, normalize_proportions

from .base import BaseDeconvolver


def read_rctd_proportions(path: Union[str, Path]) -> pd.DataFrame:
    """Read RCTD's estimated proportions CSV in standardized form."""
    proportions = pd.read_csv(path, index_col=0)
    proportions.index = proportions.index.astype(str)
    return normalize_proportions(proportions)


class RCTDDeconvolver(BaseDeconvolver):
    """Adapter for the existing RCTD wrapper."""

    method_name = "rctd"

    def run(self, data: DeconvolutionInput) -> DeconvolutionResult:
        from deconvatac.tl.rctd import rctd

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
                "doublet_mode": self.config.get("doublet_mode", "full"),
                "r_lib_path": self.config.get("r_lib_path"),
                "results_path": str(native_dir),
                "create_rctd_kwargs": self.config.get("create_rctd_kwargs", {}),
            }
            rctd(data.spatial, data.reference, **params)

            result_path = native_dir / "estimated_proportions.csv"
            proportions = read_rctd_proportions(result_path)
            diagnostics = {
                "method": self.method_name,
                "params": {key: value for key, value in params.items() if key != "results_path"},
                "native_output_dir": str(native_dir),
                "native_files": {"estimated_proportions": str(result_path)},
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
