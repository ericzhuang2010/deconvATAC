from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional

import pandas as pd

from deconvatac.data import DeconvolutionInput, DeconvolutionResult, normalize_proportions

from .base import BaseDeconvolver


def clean_cell2location_columns(columns: pd.Index) -> pd.Index:
    """Remove common Cell2location abundance prefixes from result columns."""
    prefixes = (
        "meanscell_abundance_w_sf_",
        "means_cell_abundance_w_sf_",
        "q05cell_abundance_w_sf_",
        "q05_cell_abundance_w_sf_",
    )
    cleaned = []
    for column in columns:
        name = str(column)
        for prefix in prefixes:
            if name.startswith(prefix):
                name = name[len(prefix) :]
                break
        cleaned.append(name)
    return pd.Index(cleaned)


class Cell2LocationDeconvolver(BaseDeconvolver):
    """Adapter for the existing Cell2location wrapper."""

    method_name = "cell2location"

    def run(self, data: DeconvolutionInput) -> DeconvolutionResult:
        from deconvatac.tl.cell2location import cell2location

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
                "N_cells_per_location": self.config.get("N_cells_per_location", 8),
                "detection_alpha": self.config.get("detection_alpha", 20),
                "labels_key": self.config.get("labels_key", data.labels_key),
                "layer_spatial": self.config.get("layer_spatial"),
                "layer_ref": self.config.get("layer_ref"),
                "use_gpu": self.config.get("use_gpu", True),
                "max_epochs_spatial": self.config.get("max_epochs_spatial", 30000),
                "max_epochs_ref": self.config.get("max_epochs_ref"),
                "return_adatas": False,
                "plots": self.config.get("plots", False),
                "results_path": str(native_dir),
                "setup_ref_kwargs": self.config.get("setup_ref_kwargs", {}),
                "train_ref_kwargs": self.config.get("train_ref_kwargs", {}),
                "setup_spatial_kwargs": self.config.get("setup_spatial_kwargs", {}),
                "train_spatial_kwargs": self.config.get("train_spatial_kwargs", {}),
            }
            cell2location(data.spatial, data.reference, **params)

            abundance_path = native_dir / "means_cell_abundance_w_sf.csv"
            q05_path = native_dir / "q05_cell_abundance_w_sf.csv"
            abundance = pd.read_csv(abundance_path, index_col=0)
            abundance.columns = clean_cell2location_columns(abundance.columns)
            proportions = normalize_proportions(
                abundance,
                zero_policy=self.config.get("zero_policy", "zeros"),
            )

            diagnostics = {
                "method": self.method_name,
                "params": {key: value for key, value in params.items() if key != "results_path"},
                "native_output_dir": str(native_dir),
                "native_files": {
                    "mean_abundance": str(abundance_path),
                    "q05_abundance": str(q05_path),
                },
                "proportion_source": "row-normalized means_cell_abundance_w_sf.csv",
            }
            return DeconvolutionResult(
                method=self.method_name,
                dataset_id=data.dataset_id,
                modality=data.modality,
                feature_set=data.feature_set,
                proportions=proportions,
                abundance=abundance,
                diagnostics=diagnostics,
            )
        finally:
            if temp_dir is not None:
                temp_dir.cleanup()
