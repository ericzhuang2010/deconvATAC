import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

import anndata as ad
import pandas as pd
import yaml


@dataclass
class DeconvolutionInput:
    """Standard input passed to every deconvolution method."""

    dataset_id: str
    modality: str
    feature_set: str
    spatial: ad.AnnData
    reference: ad.AnnData
    labels_key: str
    spatial_key: str = "spatial"
    truth: Optional[pd.DataFrame] = None
    output_dir: Optional[Path] = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class DeconvolutionResult:
    """Standard output returned by every deconvolution method."""

    method: str
    dataset_id: str
    modality: str
    feature_set: str
    proportions: pd.DataFrame
    abundance: Optional[pd.DataFrame] = None
    diagnostics: dict[str, Any] = field(default_factory=dict)
    output_dir: Optional[Path] = None

    def write(self, output_dir: Union[str, Path], extra_metadata: Optional[dict[str, Any]] = None) -> Path:
        """Write standardized run outputs and return the output directory."""
        output_dir = Path(output_dir)
        results_dir = output_dir / "results"
        results_dir.mkdir(parents=True, exist_ok=True)

        self.proportions.to_csv(results_dir / "proportions.csv")
        if self.abundance is not None:
            self.abundance.to_csv(results_dir / "abundance.csv")

        with (results_dir / "diagnostics.json").open("w") as handle:
            json.dump(self.diagnostics, handle, indent=2, sort_keys=True, default=str)

        run_metadata: dict[str, Any] = {
            "method": self.method,
            "dataset_id": self.dataset_id,
            "modality": self.modality,
            "feature_set": self.feature_set,
            "results": {
                "proportions": "results/proportions.csv",
                "abundance": "results/abundance.csv" if self.abundance is not None else None,
                "diagnostics": "results/diagnostics.json",
            },
        }
        if extra_metadata:
            run_metadata.update(extra_metadata)

        with (output_dir / "run.yaml").open("w") as handle:
            yaml.safe_dump(run_metadata, handle, sort_keys=False)

        self.output_dir = output_dir
        return output_dir
