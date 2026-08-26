#!/usr/bin/env python
"""Build the audited GSE244618 human hippocampus ShapeMix reference."""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import os
import platform
import re
import tempfile
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import scipy
import yaml

from deconvatac.data import FragmentShapeSpec, ordered_feature_sha256
from deconvatac.data.validators import (
    validate_fragment_shape_feature_axis,
    validate_fragment_shape_spec,
)
from deconvatac.pp.fragment_shapes import (
    FragmentRecord,
    build_fragment_shape_anndata,
    build_peak_index,
    count_fragment_shapes_from_records,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs/data_sources/shapemix_gse244618_reference.yaml"
LOCK_PATH = ROOT / "configs/data_sources/shapemix_gse244618_reference_lock.yaml"
CELL_TYPES = (
    "Inhibitory neurons",
    "Excitatory neurons",
    "Astrocytes",
    "Oligodendrocytes",
    "OPCs",
    "Microglia",
)
MIN_CELLS_PER_SAMPLE_TYPE = 20
MIN_REFERENCE_CELLS_PER_PEAK = 10
N_TOP_PEAKS = 5_000
COUNT_CHUNK_SIZE = 1_000_000
LABEL_VERSION = "broad6_v1"
SAMPLE_PATTERN = re.compile(r"_(MM_\d+)\.bedpe\.gz$")


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a mapping.")
    return value


CONFIG = _read_yaml(CONFIG_PATH)
FAMILY_ROOT = ROOT / str(CONFIG["processed_directory"])
LABEL_ROOT = FAMILY_ROOT / "labels" / LABEL_VERSION
AXIS_ROOT = FAMILY_ROOT / "feature_axes" / LABEL_VERSION
CACHE_ROOT = FAMILY_ROOT / "fragment_shape_cache" / LABEL_VERSION
REFERENCE_ROOT = ROOT / "data/processed/references" / str(CONFIG["standardized_reference_id"])
RAW_ROOT = ROOT / str(CONFIG["raw_directory"])
ANNOTATIONS_PATH = RAW_ROOT / "portal_metadata/Table_S3_nuclei.tsv.gz"
TAXONOMY_PATH = RAW_ROOT / "portal_metadata/Table_S4.xlsx"
CCRE_PATH = RAW_ROOT / "portal_metadata/Table_S6_ccres.bed.gz"


def _repository_path(path: Path) -> str:
    return path.absolute().relative_to(ROOT.absolute()).as_posix()


def _write_yaml(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
    )
    try:
        with os.fdopen(descriptor, "w") as handle:
            yaml.safe_dump(dict(value), handle, sort_keys=False)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _software_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }


def _lock_hashes() -> dict[str, str]:
    files = _read_yaml(LOCK_PATH).get("files")
    if not isinstance(files, list):
        raise ValueError("The GSE244618 source lock has no file records.")
    hashes = {
        str(record["name"]): str(record["sha256"])
        for record in files
        if isinstance(record, Mapping)
    }
    if len(hashes) != int(_read_yaml(LOCK_PATH)["resources"]):
        raise ValueError("The GSE244618 source lock is incomplete.")
    return hashes


def _sample_records() -> list[dict[str, Any]]:
    records = []
    for resource in CONFIG["resources"]:
        if resource.get("role") != "selected_hippocampal_parent_fragments":
            continue
        name = str(resource["name"])
        match = SAMPLE_PATTERN.search(name)
        if match is None:
            raise ValueError(f"Cannot derive author sample from {name}.")
        records.append(
            {
                "gsm": str(resource["gsm"]),
                "sample": match.group(1),
                "region": str(resource["region"]),
                "donor": int(resource["donor"]),
                "name": name,
                "path": RAW_ROOT / str(resource["destination"]).split(
                    "data/raw/sources/ncbi_geo/GSE244618/", 1
                )[-1],
            }
        )
    if len(records) != 9:
        raise ValueError("The frozen GSE244618 subset must contain nine BEDPE files.")
    if len({record["sample"] for record in records}) != 9:
        raise ValueError("The frozen GSE244618 author samples must be unique.")
    return records


def broad_label(cell_class: str, subclass: str) -> str | None:
    if cell_class == "GABA":
        return "Inhibitory neurons"
    if cell_class == "GLUT":
        return "Excitatory neurons"
    if cell_class != "NonN":
        return None
    if subclass in {"ACBGM", "ASCNT", "ASCT"}:
        return "Astrocytes"
    if subclass == "OGC":
        return "Oligodendrocytes"
    if subclass == "OPC":
        return "OPCs"
    if subclass == "MGC":
        return "Microglia"
    return None


def _source_labels() -> pd.DataFrame:
    samples = _sample_records()
    sample_metadata = pd.DataFrame.from_records(samples).set_index("sample")
    frame = pd.read_csv(
        ANNOTATIONS_PATH,
        sep="\t",
        usecols=["barcode", "sample", "cellclass", "subclass", "celltype"],
        dtype=str,
    )
    frame = frame[frame["sample"].isin(sample_metadata.index)].copy()
    if frame.empty:
        raise ValueError("No selected GSE244618 annotations were found.")
    frame["cell_type"] = [
        broad_label(cell_class, subclass)
        for cell_class, subclass in zip(frame["cellclass"], frame["subclass"])
    ]
    frame = frame[frame["cell_type"].notna()].copy()
    frame["cell_type"] = pd.Categorical(
        frame["cell_type"],
        categories=list(CELL_TYPES),
        ordered=True,
    )
    frame["gsm"] = frame["sample"].map(sample_metadata["gsm"])
    frame["region"] = frame["sample"].map(sample_metadata["region"])
    frame["donor"] = frame["sample"].map(sample_metadata["donor"]).astype(int)
    frame["fragment_barcode"] = frame["barcode"].astype(str)
    frame["cell_id"] = frame["sample"].astype(str) + ":" + frame["barcode"].astype(str)
    if frame["cell_id"].duplicated().any():
        raise ValueError("GSE244618 sample-plus-barcode cell IDs are not unique.")
    support = frame.groupby(
        ["sample", "cell_type"],
        observed=False,
    ).size().unstack(fill_value=0)
    support = support.reindex(
        index=[record["sample"] for record in samples],
        columns=list(CELL_TYPES),
        fill_value=0,
    )
    failures = np.argwhere(support.to_numpy() < MIN_CELLS_PER_SAMPLE_TYPE)
    if len(failures):
        details = [
            f"{support.index[row]}:{support.columns[column]}={support.iat[row, column]}"
            for row, column in failures
        ]
        raise ValueError(
            "GSE244618 broad-label support gate failed: " + ", ".join(details)
        )
    return frame.sort_values(["sample", "barcode"]).reset_index(drop=True)


def build_labels() -> pd.DataFrame:
    cells_path = LABEL_ROOT / "cells.tsv.gz"
    manifest_path = LABEL_ROOT / "manifest.yaml"
    if cells_path.is_file() and manifest_path.is_file():
        print("gse244618 labels status=reused", flush=True)
        return pd.read_csv(cells_path, sep="\t")
    if LABEL_ROOT.exists():
        raise FileExistsError(f"Partial immutable label directory: {LABEL_ROOT}")
    labels = _source_labels()
    LABEL_ROOT.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{LABEL_VERSION}.", dir=LABEL_ROOT.parent)
    )
    try:
        cells_output = temporary / "cells.tsv.gz"
        support_output = temporary / "support.tsv"
        labels.to_csv(cells_output, sep="\t", index=False, compression="gzip")
        labels.groupby(
            ["sample", "donor", "region", "cell_type"],
            observed=False,
        ).size().rename("cells").reset_index().to_csv(
            support_output,
            sep="\t",
            index=False,
        )
        hashes = _lock_hashes()
        _write_yaml(
            temporary / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "label_version": LABEL_VERSION,
                "cell_types": list(CELL_TYPES),
                "minimum_cells_per_sample_type": MIN_CELLS_PER_SAMPLE_TYPE,
                "retained_cells": len(labels),
                "samples": 9,
                "excluded_author_subclasses": ["EC", "PER", "SMC"],
                "inputs": {
                    "annotations": _repository_path(ANNOTATIONS_PATH),
                    "annotations_sha256": hashes[ANNOTATIONS_PATH.name],
                    "taxonomy": _repository_path(TAXONOMY_PATH),
                    "taxonomy_sha256": hashes[TAXONOMY_PATH.name],
                },
                "output_sha256": _sha256(cells_output),
            },
        )
        temporary.rename(LABEL_ROOT)
    except BaseException:
        import shutil

        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(f"gse244618 labels cells={len(labels)} status=completed", flush=True)
    return labels


def parse_bedpe_record(line: str | bytes) -> FragmentRecord:
    if isinstance(line, bytes):
        line = line.decode("utf-8")
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != 10 or any(field == "" for field in fields):
        raise ValueError("Expected ten non-empty BEDPE fields.")
    chrom1, start1_text, end1_text, chrom2, start2_text, end2_text = fields[:6]
    read_name, mapq_text, strand1, strand2 = fields[6:]
    if chrom1 != chrom2:
        raise ValueError("BEDPE mates must be on the same chromosome.")
    if {strand1, strand2} != {"+", "-"}:
        raise ValueError("BEDPE mates must have opposite strands.")
    try:
        start1, end1 = int(start1_text), int(end1_text)
        start2, end2 = int(start2_text), int(end2_text)
        mapq = int(mapq_text)
    except ValueError as exc:
        raise ValueError("BEDPE coordinate and MAPQ fields must be integers.") from exc
    if min(start1, start2, mapq) < 0 or end1 <= start1 or end2 <= start2:
        raise ValueError("BEDPE coordinates or MAPQ are invalid.")
    cut1 = start1 if strand1 == "+" else end1
    cut2 = start2 if strand2 == "+" else end2
    start, end = sorted((cut1, cut2))
    if end <= start:
        raise ValueError("BEDPE 5-prime cut sites do not define a positive fragment.")
    barcode = read_name.split(":", 1)[0]
    if not barcode:
        raise ValueError("BEDPE read name has no barcode prefix.")
    return FragmentRecord(
        chrom=chrom1,
        start=start,
        end=end,
        barcode=barcode,
        read_support=1,
    )


def iter_bedpe_fragments(path: Path) -> Iterator[FragmentRecord]:
    with gzip.open(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            try:
                yield parse_bedpe_record(line)
            except ValueError as exc:
                raise ValueError(f"{path}: invalid BEDPE row {line_number}: {exc}") from exc


def _candidate_ccres() -> pd.DataFrame:
    rows = []
    with gzip.open(CCRE_PATH, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"cCRE row {line_number} has fewer than four fields.")
            chrom, start_text, end_text, peak_id = fields[:4]
            start, end = int(start_text), int(end_text)
            rows.append((chrom, start, end, peak_id))
    frame = pd.DataFrame(rows, columns=["chrom", "start", "end", "peak_id"])
    if frame.empty or frame["peak_id"].duplicated().any():
        raise ValueError("The author cCRE axis must be non-empty and unique.")
    build_peak_index(
        list(frame[["chrom", "start", "end", "peak_id"]].itertuples(index=False, name=None))
    )
    return frame


def _aggregate_feature_statistics(
    labels: pd.DataFrame,
    candidates: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    peaks = build_peak_index(
        list(candidates[["chrom", "start", "end", "peak_id"]].itertuples(index=False, name=None))
    )
    counts = np.zeros((len(CELL_TYPES), len(candidates)), dtype=np.int64)
    coverage = np.zeros((len(CELL_TYPES), len(candidates)), dtype=np.int32)
    type_index = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    audit_records: list[dict[str, Any]] = []
    for source in _sample_records():
        sample_labels = labels[labels["sample"] == source["sample"]]
        barcode_types = {
            str(barcode): type_index[str(cell_type)]
            for barcode, cell_type in zip(
                sample_labels["fragment_barcode"],
                sample_labels["cell_type"],
            )
        }
        current_barcode: str | None = None
        current_group: int | None = None
        current_peaks: set[int] = set()
        closed_barcodes: set[str] = set()
        rows = 0
        assigned_cut_sites = 0

        def flush_cell() -> None:
            if current_group is not None and current_peaks:
                indices = np.fromiter(current_peaks, dtype=np.int64)
                coverage[current_group, indices] += 1

        for record in iter_bedpe_fragments(source["path"]):
            rows += 1
            if record.barcode != current_barcode:
                flush_cell()
                if current_barcode is not None:
                    closed_barcodes.add(current_barcode)
                if record.barcode in closed_barcodes:
                    raise ValueError(
                        f"{source['name']} is not grouped contiguously by barcode."
                    )
                current_barcode = record.barcode
                current_group = barcode_types.get(record.barcode)
                current_peaks = set()
            if current_group is not None:
                for coordinate in (record.start, record.end):
                    peak = peaks.assign(record.chrom, coordinate)
                    if peak is not None:
                        counts[current_group, peak] += 1
                        current_peaks.add(peak)
                        assigned_cut_sites += 1
            if rows % 5_000_000 == 0:
                print(
                    f"gse244618 feature_scan sample={source['sample']} rows={rows}",
                    flush=True,
                )
        flush_cell()
        audit_records.append(
            {
                "sample": source["sample"],
                "gsm": source["gsm"],
                "rows": rows,
                "selected_annotation_barcodes": len(barcode_types),
                "observed_barcode_blocks": len(closed_barcodes) + int(current_barcode is not None),
                "assigned_cut_sites": assigned_cut_sites,
            }
        )
        print(
            f"gse244618 feature_scan sample={source['sample']} rows={rows} status=completed",
            flush=True,
        )
    return counts, coverage, audit_records


def build_feature_axis() -> pd.DataFrame:
    selected_path = AXIS_ROOT / "selected_ccres.tsv.gz"
    manifest_path = AXIS_ROOT / "manifest.yaml"
    if selected_path.is_file() and manifest_path.is_file():
        print("gse244618 feature_axis status=reused", flush=True)
        return pd.read_csv(selected_path, sep="\t")
    if AXIS_ROOT.exists():
        raise FileExistsError(f"Partial immutable feature axis: {AXIS_ROOT}")
    labels = build_labels()
    candidates = _candidate_ccres()
    counts, coverage_by_type, audit_records = _aggregate_feature_statistics(
        labels,
        candidates,
    )
    type_totals = counts.sum(axis=1, dtype=np.float64)
    if (type_totals <= 0).any() or not np.isfinite(type_totals).all():
        raise ValueError("At least one GSE244618 broad type has no assigned cut sites.")
    normalized = np.log2(1.0 + 1.0e4 * counts.astype(np.float64) / type_totals[:, None])
    score = np.var(normalized, axis=0, ddof=0)
    coverage = coverage_by_type.sum(axis=0, dtype=np.int64)
    total_counts = counts.sum(axis=0, dtype=np.int64)
    ranked = candidates.copy()
    ranked["score"] = score
    ranked["nonzero_reference_cells"] = coverage
    ranked["total_reference_counts"] = total_counts
    ranked = ranked[ranked["nonzero_reference_cells"] >= MIN_REFERENCE_CELLS_PER_PEAK]
    ranked = ranked.sort_values(
        ["score", "nonzero_reference_cells", "total_reference_counts", "peak_id"],
        ascending=[False, False, False, True],
        kind="mergesort",
    )
    if len(ranked) < N_TOP_PEAKS:
        raise ValueError(f"Only {len(ranked)} author cCREs pass the support gate.")
    selected = ranked.head(N_TOP_PEAKS).copy()
    selected.insert(0, "rank", np.arange(1, N_TOP_PEAKS + 1))
    AXIS_ROOT.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{LABEL_VERSION}.", dir=AXIS_ROOT.parent)
    )
    try:
        output = temporary / "selected_ccres.tsv.gz"
        selected.to_csv(output, sep="\t", index=False, compression="gzip")
        (temporary / "selected_peaks.txt").write_text(
            "\n".join(selected["peak_id"].astype(str)) + "\n"
        )
        _write_yaml(
            temporary / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "axis_version": LABEL_VERSION,
                "label_ontology": list(CELL_TYPES),
                "candidate_ccres": len(candidates),
                "eligible_ccres": len(ranked),
                "selected_ccres": N_TOP_PEAKS,
                "feature_sha256": ordered_feature_sha256(selected["peak_id"].astype(str)),
                "selector": {
                    "minimum_nonzero_reference_cells": MIN_REFERENCE_CELLS_PER_PEAK,
                    "score": "population_variance_log2_1_plus_scaled_type_rate",
                    "scale": 1.0e4,
                    "tie_breaks": [
                        "score_desc",
                        "coverage_desc",
                        "total_count_desc",
                        "peak_id_asc",
                    ],
                },
                "source_scans": audit_records,
                "inputs": {
                    "labels": _repository_path(LABEL_ROOT / "cells.tsv.gz"),
                    "candidate_ccres": _repository_path(CCRE_PATH),
                },
                "software_versions": _software_versions(),
            },
        )
        temporary.rename(AXIS_ROOT)
    except BaseException:
        import shutil

        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print("gse244618 feature_axis peaks=5000 status=completed", flush=True)
    return selected


def _combined_metadata(
    values: Iterable[Mapping[str, Any]],
    layers: Mapping[str, sparse.csr_matrix],
    feature_names: Iterable[str],
    split_sha256: str,
    source_hashes: Mapping[str, str],
) -> dict[str, Any]:
    metadata_values = [copy.deepcopy(dict(value)) for value in values]
    if not metadata_values:
        raise ValueError("Cannot combine zero fragment-shape metadata records.")
    metadata = metadata_values[0]
    counters: dict[str, int] = {}
    for value in metadata_values:
        for key, count in dict(value.get("preprocessing_counters", {})).items():
            if isinstance(count, (int, np.integer)):
                counters[key] = counters.get(key, 0) + int(count)
    layer_totals = {name: int(matrix.sum()) for name, matrix in layers.items()}
    metadata["preprocessing_counters"] = counters
    metadata["matrix_counters"] = {
        "assigned_cut_sites": int(sum(layer_totals.values())),
        **{f"cut_sites_per_bin.{name}": count for name, count in layer_totals.items()},
    }
    metadata["feature_sha256"] = ordered_feature_sha256(feature_names)
    metadata["split_sha256"] = split_sha256
    metadata["source_sha256"] = dict(source_hashes)
    return metadata


def build_reference() -> Path:
    reference_path = REFERENCE_ROOT / "atac/reference.h5ad"
    reference_manifest = REFERENCE_ROOT / "reference.yaml"
    if reference_path.is_file() and reference_manifest.is_file():
        print("gse244618 reference status=reused", flush=True)
        return reference_path
    if REFERENCE_ROOT.exists():
        raise FileExistsError(f"Partial immutable reference directory: {REFERENCE_ROOT}")
    labels = build_labels()
    selected = build_feature_axis()
    peaks = list(
        selected[["chrom", "start", "end", "peak_id"]].itertuples(index=False, name=None)
    )
    var = selected.set_index("peak_id")[["chrom", "start", "end"]].copy()
    label_sha256 = _sha256(LABEL_ROOT / "cells.tsv.gz")
    locked_hashes = _lock_hashes()
    cache_objects: list[Path] = []
    for source in _sample_records():
        output_dir = CACHE_ROOT / source["sample"]
        output_path = output_dir / "cells.h5ad"
        if output_path.is_file():
            cache_objects.append(output_path)
            print(f"gse244618 fragment_cache sample={source['sample']} status=reused", flush=True)
            continue
        if output_dir.exists():
            raise FileExistsError(f"Partial GSE244618 cache directory: {output_dir}")
        sample_labels = labels[labels["sample"] == source["sample"]].copy()
        sample_labels = sample_labels.sort_values("fragment_barcode")
        barcodes = sample_labels["fragment_barcode"].astype(str).tolist()
        result = count_fragment_shapes_from_records(
            iter_bedpe_fragments(source["path"]),
            barcodes,
            peaks,
            right_cut_offset=0,
            chunk_size=COUNT_CHUNK_SIZE,
        )
        obs = sample_labels.set_index(sample_labels["fragment_barcode"].astype(str))
        obs.index.name = "fragment_barcode"
        source_hashes = {
            source["name"]: locked_hashes[source["name"]],
            ANNOTATIONS_PATH.name: locked_hashes[ANNOTATIONS_PATH.name],
            TAXONOMY_PATH.name: locked_hashes[TAXONOMY_PATH.name],
            CCRE_PATH.name: locked_hashes[CCRE_PATH.name],
        }
        shape = build_fragment_shape_anndata(
            result,
            obs=obs,
            var=var,
            provenance={
                "split_sha256": label_sha256,
                "source_sha256": source_hashes,
                "coordinate_validation": {
                    "selected_right_cut_offset": 0,
                    "matrix_match": "not_available",
                    "validation_method": "deposited_strand_aware_bedpe_5prime",
                    "semantic_match": "exact",
                },
                "software_versions": _software_versions(),
            },
        )
        shape.obs_names = pd.Index(sample_labels["cell_id"].astype(str), name="cell_id")
        output_dir.mkdir(parents=True)
        temporary = output_dir / ".cells.h5ad.tmp"
        shape.write_h5ad(temporary, compression="gzip")
        temporary.replace(output_path)
        _write_yaml(
            output_dir / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "sample": source["sample"],
                "gsm": source["gsm"],
                "cells": shape.n_obs,
                "peaks": shape.n_vars,
                "source_sha256": source_hashes,
                "preprocessing_counters": result.qc.to_dict(),
                "output": _repository_path(output_path),
            },
        )
        cache_objects.append(output_path)
        print(
            f"gse244618 fragment_cache sample={source['sample']} status=completed",
            flush=True,
        )

    obs_parts = []
    layer_parts: dict[str, list[sparse.csr_matrix]] = {}
    metadata_parts = []
    source_hashes: dict[str, str] = {
        ANNOTATIONS_PATH.name: locked_hashes[ANNOTATIONS_PATH.name],
        TAXONOMY_PATH.name: locked_hashes[TAXONOMY_PATH.name],
        CCRE_PATH.name: locked_hashes[CCRE_PATH.name],
    }
    for source, path in zip(_sample_records(), cache_objects):
        cache = ad.read_h5ad(path)
        obs_parts.append(cache.obs.copy())
        spec = FragmentShapeSpec.from_mapping(cache.uns["fragment_shape"])
        for layer_name in spec.layer_names:
            layer_parts.setdefault(layer_name, []).append(
                sparse.csr_matrix(cache.layers[layer_name], dtype=np.int64)
            )
        metadata_parts.append(cache.uns["fragment_shape"])
        source_hashes[source["name"]] = locked_hashes[source["name"]]
    obs = pd.concat(obs_parts, axis=0)
    if obs.index.duplicated().any() or set(obs["cell_type"]) != set(CELL_TYPES):
        raise ValueError("The assembled GSE244618 reference cell axis is invalid.")
    layers = {
        name: sparse.vstack(parts, format="csr", dtype=np.int64)
        for name, parts in layer_parts.items()
    }
    x = sparse.csr_matrix((len(obs), len(var)), dtype=np.int64)
    for matrix in layers.values():
        x = (x + matrix).tocsr()
    reference = ad.AnnData(X=x, obs=obs, var=var.copy())
    for name, matrix in layers.items():
        reference.layers[name] = matrix
    reference.uns["fragment_shape"] = _combined_metadata(
        metadata_parts,
        layers,
        reference.var_names.astype(str),
        label_sha256,
        source_hashes,
    )
    validate_fragment_shape_spec(
        FragmentShapeSpec.from_mapping(reference.uns["fragment_shape"])
    )
    validate_fragment_shape_feature_axis(reference, "GSE244618 reference")
    temporary_root = Path(
        tempfile.mkdtemp(
            prefix=f".{CONFIG['standardized_reference_id']}.",
            dir=REFERENCE_ROOT.parent,
        )
    )
    try:
        (temporary_root / "atac").mkdir()
        temporary_path = temporary_root / "atac/reference.h5ad"
        reference.write_h5ad(temporary_path, compression="gzip")
        _write_yaml(
            temporary_root / "reference.yaml",
            {
                "schema_version": 1,
                "reference_id": CONFIG["standardized_reference_id"],
                "source_dataset_id": CONFIG["source_dataset_id"],
                "description": (
                    "Nine-sample adult human hippocampus snATAC reference on a "
                    "reference-only 5,000-cCRE GRCh38 axis."
                ),
                "labels_key": "cell_type",
                "genome_build": "GRCh38",
                "cell_types": list(CELL_TYPES),
                "counts": {"cells": reference.n_obs, "peaks": reference.n_vars},
                "modalities": {
                    "atac": {
                        "path": _repository_path(
                            REFERENCE_ROOT / "atac/reference.h5ad"
                        ),
                        "feature_type": "cCREs",
                    }
                },
                "provenance": {
                    "source_lock": _repository_path(LOCK_PATH),
                    "label_manifest": _repository_path(LABEL_ROOT / "manifest.yaml"),
                    "feature_manifest": _repository_path(AXIS_ROOT / "manifest.yaml"),
                },
            },
        )
        temporary_root.rename(REFERENCE_ROOT)
    except BaseException:
        import shutil

        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    print(
        f"gse244618 reference cells={reference.n_obs} peaks={reference.n_vars} status=completed",
        flush=True,
    )
    return reference_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("labels", "feature-axis", "reference", "all"),
        default="all",
    )
    args = parser.parse_args()
    if args.stage in {"labels", "all"}:
        build_labels()
    if args.stage in {"feature-axis", "all"}:
        build_feature_axis()
    if args.stage in {"reference", "all"}:
        build_reference()


if __name__ == "__main__":
    main()
