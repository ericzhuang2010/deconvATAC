from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from deconvatac.data import DeconvolutionInput, DeconvolutionResult, normalize_proportions

from .base import BaseDeconvolver


_SPOT_COLUMNS = ("cell_ID", "spot", "spot_id", "barcode", "location", "location_id")
_CELL_TYPE_COLUMNS = ("cell_type", "celltype", "feat_ID", "feature", "annotation")
_VALUE_COLUMNS = ("proportion", "proportions", "enrichment", "DWLS", "estimate", "value")


def _drop_unnamed_columns(values: pd.DataFrame) -> pd.DataFrame:
    return values.loc[:, [column for column in values.columns if not str(column).startswith("Unnamed:")]]


def read_spatialdwls_proportions(path: Union[str, Path], zero_policy: str = "zeros") -> pd.DataFrame:
    """Read SpatialDWLS/Giotto output in standardized proportion form."""
    values = _drop_unnamed_columns(pd.read_csv(path))
    values.columns = [str(column) for column in values.columns]

    cell_type_column = next((column for column in _CELL_TYPE_COLUMNS if column in values.columns), None)
    value_column = next((column for column in _VALUE_COLUMNS if column in values.columns), None)
    spot_column = next((column for column in _SPOT_COLUMNS if column in values.columns), None)

    if cell_type_column is not None and value_column is not None and spot_column is not None:
        wide = values.pivot_table(
            index=spot_column,
            columns=cell_type_column,
            values=value_column,
            aggfunc="sum",
            fill_value=0.0,
        )
    else:
        if spot_column is not None:
            wide = values.set_index(spot_column)
        else:
            wide = values.copy()
            wide.index = wide.index.astype(str)
        metadata_columns = set(_SPOT_COLUMNS) | set(_CELL_TYPE_COLUMNS) | set(_VALUE_COLUMNS)
        wide = wide.drop(columns=[column for column in wide.columns if column in metadata_columns], errors="ignore")
        wide = wide.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all")

    wide.index = wide.index.astype(str)
    return normalize_proportions(wide, zero_policy=zero_policy)


class SpatialDWLSDeconvolver(BaseDeconvolver):
    """Adapter for the existing SpatialDWLS wrapper."""

    method_name = "spatialdwls"

    def run(self, data: DeconvolutionInput) -> DeconvolutionResult:
        from deconvatac.tl.spatialdwls import spatialdwls

        output_dir = Path(data.output_dir) if data.output_dir is not None else None
        temp_dir: Optional[tempfile.TemporaryDirectory[str]] = None
        if output_dir is None:
            temp_dir = tempfile.TemporaryDirectory()
            native_dir = Path(temp_dir.name)
        else:
            native_dir = output_dir / "results" / "raw_method_output"
            native_dir.mkdir(parents=True, exist_ok=True)

        try:
            default_cluster_key = "leiden_lsi" if data.modality == "atac" else "leiden_pca"
            cluster_key = self.config.get("cluster_key") or default_cluster_key
            configured_tfidf = self.config.get("tfidf")
            tfidf = data.modality == "atac" if configured_tfidf is None else configured_tfidf
            params = {
                "labels_key": self.config.get("labels_key", data.labels_key),
                "cluster_key": cluster_key,
                "n_cell": self.config.get("n_cell", 50),
                "tfidf": tfidf,
                "r_lib_path": self.config.get("r_lib_path"),
                "results_path": str(native_dir),
            }
            spatialdwls(data.spatial, data.reference, **params)

            result_path = native_dir / "proportions.csv"
            zero_policy = self.config.get("zero_policy", "zeros")
            proportions = read_spatialdwls_proportions(result_path, zero_policy=zero_policy)
            diagnostics = {
                "method": self.method_name,
                "params": {key: value for key, value in params.items() if key != "results_path"},
                "native_output_dir": str(native_dir),
                "native_files": {"proportions": str(result_path)},
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
