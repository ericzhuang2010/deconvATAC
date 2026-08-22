from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from deconvatac.data import DeconvolutionInput, DeconvolutionResult, normalize_proportions

from .base import BaseDeconvolver


def read_tangram_proportions(path: Union[str, Path], zero_policy: str = "zeros") -> pd.DataFrame:
    """Read Tangram cell-type predictions in standardized proportion form."""
    values = pd.read_csv(path, index_col=0)
    values.index = values.index.astype(str)
    return normalize_proportions(values, zero_policy=zero_policy)


class TangramDeconvolver(BaseDeconvolver):
    """Adapter for the existing Tangram wrapper."""

    method_name = "tangram"

    def run(self, data: DeconvolutionInput) -> DeconvolutionResult:
        from deconvatac.tl.tangram import tangram

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
                "run_rank_genes": self.config.get("run_rank_genes", False),
                "layer_rank_genes": self.config.get("layer_rank_genes"),
                "num_epochs": self.config.get("num_epochs", 1000),
                "device": self.config.get("device", "cpu"),
                "return_adatas": False,
                "result_path": str(native_dir),
            }
            params.update(self.config.get("map_kwargs", {}))
            tangram(data.spatial, data.reference, **params)

            result_path = native_dir / "tangram_ct_pred.csv"
            zero_policy = self.config.get("zero_policy", "zeros")
            proportions = read_tangram_proportions(result_path, zero_policy=zero_policy)
            diagnostics = {
                "method": self.method_name,
                "params": {key: value for key, value in params.items() if key != "result_path"},
                "native_output_dir": str(native_dir),
                "native_files": {"cell_type_predictions": str(result_path)},
                "proportion_source": "row-normalized tangram_ct_pred.csv",
            }
            return DeconvolutionResult(
                method=self.method_name,
                dataset_id=data.dataset_id,
                modality=data.modality,
                feature_set=data.feature_set,
                proportions=proportions,
                abundance=proportions,
                diagnostics=diagnostics,
            )
        finally:
            if temp_dir is not None:
                temp_dir.cleanup()
