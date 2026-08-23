from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

import anndata as ad
import pandas as pd
import yaml


@dataclass(frozen=True)
class FragmentShapeBin:
    """One ordered bin in a fragment-shape axis."""

    name: str
    min_inclusive: int
    max_exclusive: Optional[int]
    layer: str

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> FragmentShapeBin:
        """Construct a bin from YAML or ``AnnData.uns`` metadata."""
        if not isinstance(value, Mapping):
            raise TypeError("Each fragment_shape bin must be a mapping.")

        required = {"name", "min_inclusive", "layer"}
        missing = sorted(required.difference(value))
        if missing:
            raise ValueError(f"fragment_shape bin is missing required fields: {', '.join(missing)}.")

        return cls(
            name=value["name"],
            min_inclusive=value["min_inclusive"],
            max_exclusive=value.get("max_exclusive"),
            layer=value["layer"],
        )

    def to_dict(self, omit_none: bool = False) -> dict[str, Any]:
        """Return a serializable representation in canonical field order."""
        value: dict[str, Any] = {
            "name": self.name,
            "min_inclusive": self.min_inclusive,
            "max_exclusive": self.max_exclusive,
            "layer": self.layer,
        }
        if omit_none:
            return {key: item for key, item in value.items() if item is not None}
        return value


def _fragment_shape_bins(value: Any) -> tuple[FragmentShapeBin, ...]:
    """Normalize YAML and H5AD-safe bin encodings while preserving order."""
    if isinstance(value, pd.DataFrame):
        records = value.to_dict(orient="records")
    elif isinstance(value, Mapping):
        column_names = {"name", "min_inclusive", "max_exclusive", "layer"}
        if column_names.issubset(value):
            columns = {key: list(value[key]) for key in column_names}
            lengths = {len(column) for column in columns.values()}
            if len(lengths) != 1:
                raise ValueError("Parallel fragment_shape bin fields must have equal lengths.")
            records = [
                {key: columns[key][index] for key in column_names}
                for index in range(next(iter(lengths), 0))
            ]
        else:
            # Numeric keys are the recommended H5AD-safe representation because
            # AnnData 0.9 cannot serialize a list of nested mappings.
            keys = list(value)
            if keys and all(str(key).isdigit() for key in keys):
                keys = sorted(keys, key=lambda key: int(str(key)))
            elif keys and all(isinstance(value[key], Mapping) and "order" in value[key] for key in keys):
                keys = sorted(keys, key=lambda key: value[key]["order"])
            records = []
            for key in keys:
                record = value[key]
                if not isinstance(record, Mapping):
                    raise TypeError("Each fragment_shape bin must be a mapping.")
                record = dict(record)
                record.setdefault("name", str(key))
                records.append(record)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        records = list(value)
    elif hasattr(value, "tolist"):
        records = value.tolist()
    else:
        raise TypeError("fragment_shape.bins must be an ordered sequence or mapping.")

    return tuple(FragmentShapeBin.from_mapping(record) for record in records)


@dataclass(frozen=True)
class FragmentShapeSpec:
    """Typed declaration of aligned fragment-shape layers and provenance."""

    schema_version: int
    axis: str
    count_unit: str
    read_support_policy: str
    peak_assignment: str
    bins: tuple[FragmentShapeBin, ...]
    left_cut_offset: Optional[int] = None
    right_cut_offset: Optional[int] = None
    source_sha256: Optional[dict[str, str]] = None
    feature_sha256: Optional[str] = None
    split_sha256: Optional[str] = None
    coordinate_validation: Optional[dict[str, Any]] = None
    software_versions: Optional[dict[str, str]] = None
    preprocessing_counters: Optional[dict[str, Any]] = None
    matrix_counters: Optional[dict[str, Any]] = None

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> FragmentShapeSpec:
        """Construct a typed specification from YAML or ``AnnData.uns``."""
        if not isinstance(value, Mapping):
            raise TypeError("fragment_shape metadata must be a mapping.")

        required = {
            "schema_version",
            "axis",
            "count_unit",
            "read_support_policy",
            "peak_assignment",
            "bins",
        }
        missing = sorted(required.difference(value))
        if missing:
            raise ValueError(f"fragment_shape is missing required fields: {', '.join(missing)}.")

        return cls(
            schema_version=value["schema_version"],
            axis=value["axis"],
            count_unit=value["count_unit"],
            read_support_policy=value["read_support_policy"],
            peak_assignment=value["peak_assignment"],
            bins=_fragment_shape_bins(value["bins"]),
            left_cut_offset=value.get("left_cut_offset"),
            right_cut_offset=value.get("right_cut_offset"),
            source_sha256=value.get("source_sha256"),
            feature_sha256=value.get("feature_sha256"),
            split_sha256=value.get("split_sha256"),
            coordinate_validation=value.get("coordinate_validation"),
            software_versions=value.get("software_versions"),
            preprocessing_counters=value.get("preprocessing_counters"),
            matrix_counters=value.get("matrix_counters"),
        )

    @property
    def layer_names(self) -> tuple[str, ...]:
        """Return layer names in declared bin order."""
        return tuple(bin_spec.layer for bin_spec in self.bins)

    def to_dict(self, omit_none: bool = False) -> dict[str, Any]:
        """Return the YAML-style representation of this specification."""
        value: dict[str, Any] = {
            "schema_version": self.schema_version,
            "axis": self.axis,
            "count_unit": self.count_unit,
            "read_support_policy": self.read_support_policy,
            "peak_assignment": self.peak_assignment,
            "left_cut_offset": self.left_cut_offset,
            "right_cut_offset": self.right_cut_offset,
            "bins": [bin_spec.to_dict(omit_none=omit_none) for bin_spec in self.bins],
            "source_sha256": self.source_sha256,
            "feature_sha256": self.feature_sha256,
            "split_sha256": self.split_sha256,
            "coordinate_validation": self.coordinate_validation,
            "software_versions": self.software_versions,
            "preprocessing_counters": self.preprocessing_counters,
            "matrix_counters": self.matrix_counters,
        }
        if omit_none:
            return {key: item for key, item in value.items() if item is not None}
        return value

    def to_uns(self) -> dict[str, Any]:
        """Return an H5AD-safe representation for ``uns['fragment_shape']``.

        AnnData 0.9 does not serialize lists of nested mappings. Numeric keys
        preserve the bin order without relying on HDF5 group iteration order.
        """
        value = self.to_dict(omit_none=True)
        value["bins"] = {
            str(index): bin_spec.to_dict(omit_none=True) for index, bin_spec in enumerate(self.bins)
        }
        return value


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
    fragment_shape: Optional[FragmentShapeSpec] = None
    cell_types: Optional[list[str]] = None


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
