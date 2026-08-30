#!/usr/bin/env python3
"""Materialize audited GSE205055/GSE263333 real-spatial ShapeMix datasets."""

from __future__ import annotations

import argparse
import copy
import hashlib
import os
import platform
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

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
from deconvatac.pp.fragment_shapes import build_fragment_shape_anndata, count_fragment_shapes


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TEMPLATES = (
    ROOT / "configs/datasets/shapemix_gse205055_real_spatial_v1.yaml",
    ROOT / "configs/datasets/shapemix_gse263333_real_spatial_v1.yaml",
)
REGISTRY_PATH = ROOT / "data/registry/datasets.yaml"
COUNT_CHUNK_SIZE = 1_000_000


def repository_path(path: Path) -> str:
    """Return a repository-relative path without resolving the data symlink."""
    return path.absolute().relative_to(ROOT.absolute()).as_posix()


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_yaml(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(descriptor, "w") as handle:
            yaml.safe_dump(dict(value), handle, sort_keys=False)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def software_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }


def require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {repository_path(path)}")
    return path


def require_under(path: Path, expected_parent: Path, label: str) -> None:
    path_text = path.absolute().as_posix()
    parent_text = expected_parent.absolute().as_posix().rstrip("/") + "/"
    if not path_text.startswith(parent_text):
        raise ValueError(
            f"{label} violates the canonical layout: {repository_path(path)}; "
            f"expected under {repository_path(expected_parent)}"
        )


def validate_audited_path(
    record: Mapping[str, Any],
    field: str,
    bytes_field: str,
    *,
    expected_parent: Path,
    label: str,
) -> Path:
    path = ROOT / str(record[field])
    require_under(path, expected_parent, label)
    require_file(path, label)
    expected_bytes = int(record[bytes_field])
    observed_bytes = path.stat().st_size
    if observed_bytes != expected_bytes:
        raise ValueError(
            f"{label} size changed: expected {expected_bytes}, observed {observed_bytes}: "
            f"{repository_path(path)}"
        )
    return path


def fragment_records(audit: Mapping[str, Any], gsm: str, assay: str | None = None) -> list[dict[str, Any]]:
    records = [
        dict(value)
        for value in audit.get("files", [])
        if str(value.get("gsm")) == gsm
        and (assay is None or str(value.get("assay")) == assay)
    ]
    if not records:
        suffix = f" assay={assay}" if assay is not None else ""
        raise KeyError(f"No audited fragment record for {gsm}{suffix}")
    return records


def read_barcodes(record: Mapping[str, Any]) -> list[str]:
    path = require_file(ROOT / str(record["barcode_path"]), "audited barcode list")
    if file_digest(path) != str(record["barcode_sha256"]):
        raise ValueError(f"Audited barcode-list digest changed: {repository_path(path)}")
    with path.open() as handle:
        values = [line.strip() for line in handle if line.strip()]
    if len(values) != int(record["unique_barcodes"]) or len(set(values)) != len(values):
        raise ValueError(f"Audited barcode-list cardinality changed: {repository_path(path)}")
    return values


def canonicalize_barcodes(values: Iterable[str], policy: Mapping[str, Any]) -> list[str]:
    values = [str(value) for value in values]
    suffix = policy.get("fragment_terminal_suffix_to_strip")
    if suffix is None:
        result = values
    else:
        suffix = str(suffix)
        require_suffix = bool(policy.get("require_suffix_on_every_fragment_barcode", False))
        if require_suffix and any(not value.endswith(suffix) for value in values):
            raise ValueError(f"A fragment barcode lacks required terminal suffix {suffix!r}")
        result = [value[: -len(suffix)] if value.endswith(suffix) else value for value in values]
    if len(set(result)) != len(result):
        raise ValueError("Barcode canonicalization created duplicates")
    return result


def alignment_manifest(template: Mapping[str, Any], section: Mapping[str, Any]) -> tuple[Path, dict[str, Any]]:
    path = ROOT / str(template["reusable_root"]) / "manifests/cross_modality_alignment" / f"{section['group']}.yaml"
    record = read_yaml(require_file(path, "cross-modality alignment manifest"))
    if str(record.get("group")) != str(section["group"]):
        raise ValueError(f"Alignment manifest group mismatch: {repository_path(path)}")
    return path, record


def load_coordinates(
    work_root: Path,
    gsm: str,
    canonical_barcodes: list[str],
) -> pd.DataFrame:
    path = require_file(
        work_root / "spatial_coordinates" / gsm / "coordinates.csv",
        "audited spatial coordinates",
    )
    coordinates = pd.read_csv(path, dtype={"barcode": str})
    required = {
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pixel_row_fullres",
        "pixel_col_fullres",
    }
    missing = sorted(required.difference(coordinates.columns))
    if missing:
        raise ValueError(f"Coordinate table lacks columns {missing}: {repository_path(path)}")
    coordinates["barcode"] = coordinates["barcode"].astype(str)
    if coordinates["barcode"].duplicated().any():
        raise ValueError(f"Coordinate barcodes are not unique: {repository_path(path)}")
    indexed = coordinates.set_index("barcode")
    absent = sorted(set(canonical_barcodes).difference(indexed.index))
    if absent:
        raise ValueError(
            f"{len(absent)} fragment barcodes lack coordinates for {gsm}; first={absent[0]}"
        )
    selected = indexed.loc[canonical_barcodes].copy()
    spatial = selected[["pixel_col_fullres", "pixel_row_fullres"]].to_numpy(dtype=np.float64)
    if not np.isfinite(spatial).all():
        raise ValueError(f"Non-finite spatial coordinates for {gsm}")
    selected.index = pd.Index(canonical_barcodes, name="spot")
    return selected


def reference_inputs(reference_id: str) -> tuple[Path, dict[str, Any], pd.DataFrame, FragmentShapeSpec]:
    manifest_path = ROOT / "data/processed/references" / reference_id / "reference.yaml"
    manifest = read_yaml(require_file(manifest_path, "standardized reference manifest"))
    if str(manifest.get("reference_id")) != reference_id:
        raise ValueError(f"Reference ID mismatch in {repository_path(manifest_path)}")
    reference_path = ROOT / str(manifest["modalities"]["atac"]["path"])
    require_under(reference_path, ROOT / "data/processed/references", "reference H5AD")
    require_file(reference_path, "standardized reference H5AD")
    reference = ad.read_h5ad(reference_path, backed="r")
    try:
        var = reference.var.copy()
        metadata = copy.deepcopy(dict(reference.uns["fragment_shape"]))
        feature_names = reference.var_names.astype(str).tolist()
    finally:
        reference.file.close()
    required = {"chrom", "start", "end"}
    missing = sorted(required.difference(var.columns))
    if missing:
        raise ValueError(f"Reference feature axis lacks columns {missing}: {reference_id}")
    var.index = pd.Index(feature_names, name=var.index.name or "peak")
    if len(var) != 5_000:
        raise ValueError(f"Frozen real-spatial protocol requires exactly 5,000 reference features: {reference_id}")
    spec = FragmentShapeSpec.from_mapping(metadata)
    validate_fragment_shape_spec(spec)
    if spec.right_cut_offset is None:
        raise ValueError(f"Reference lacks a selected right-cut offset: {reference_id}")
    return reference_path, manifest, var, spec


def peak_records(var: pd.DataFrame) -> list[tuple[str, int, int, str]]:
    return [
        (str(chrom), int(start), int(end), str(name))
        for name, chrom, start, end in zip(
            var.index, var["chrom"], var["start"], var["end"], strict=True
        )
    ]


def bin_records(spec: FragmentShapeSpec) -> list[dict[str, Any]]:
    return [value.to_dict(omit_none=False) for value in spec.bins]


def shape_provenance(
    record: Mapping[str, Any],
    audit_path: Path,
    alignment_path: Path,
    var: pd.DataFrame,
    spec: FragmentShapeSpec,
) -> dict[str, Any]:
    return {
        "source_sha256": {
            Path(str(record["normalized"])).name: str(record["normalized_sha256"]),
            Path(repository_path(audit_path)).name: file_digest(audit_path),
        },
        "feature_sha256": ordered_feature_sha256(var.index.astype(str)),
        # This field is part of the strictly aligned model-input contract and
        # must equal the reference metadata byte-for-value. Spatial alignment
        # provenance belongs in the dataset/section manifests instead.
        "coordinate_validation": copy.deepcopy(spec.coordinate_validation),
        "split_sha256": file_digest(alignment_path),
        "software_versions": software_versions(),
    }


def build_atac_spatial(
    record: Mapping[str, Any],
    fragments_path: Path,
    raw_barcodes: list[str],
    canonical_barcodes: list[str],
    coordinates: pd.DataFrame,
    var: pd.DataFrame,
    spec: FragmentShapeSpec,
    provenance: Mapping[str, Any],
) -> tuple[ad.AnnData, dict[str, Any]]:
    result = count_fragment_shapes(
        fragments_path,
        raw_barcodes,
        peak_records(var),
        right_cut_offset=int(spec.right_cut_offset),
        bins=bin_records(spec),
        chunk_size=COUNT_CHUNK_SIZE,
    )
    obs = coordinates.copy()
    obs.index = pd.Index(raw_barcodes, name="fragment_barcode")
    obs.insert(0, "fragment_barcode", raw_barcodes)
    spatial = build_fragment_shape_anndata(
        result,
        obs=obs,
        var=var,
        provenance=provenance,
    )
    spatial.obs_names = pd.Index(canonical_barcodes, name="spot")
    spatial.obsm["spatial"] = coordinates[
        ["pixel_col_fullres", "pixel_row_fullres"]
    ].to_numpy(dtype=np.float64)
    validate_fragment_shape_feature_axis(spatial, "real-spatial ATAC")
    observed = FragmentShapeSpec.from_mapping(spatial.uns["fragment_shape"])
    validate_fragment_shape_spec(observed)
    if observed.layer_names != spec.layer_names:
        raise ValueError("Spatial fragment-shape layers differ from the reference")
    return spatial, result.qc.to_dict()


def build_epigenome_validation(
    record: Mapping[str, Any],
    fragments_path: Path,
    raw_barcodes: list[str],
    canonical_barcodes: list[str],
    coordinates: pd.DataFrame,
    var: pd.DataFrame,
    spec: FragmentShapeSpec,
) -> tuple[ad.AnnData, dict[str, Any]]:
    result = count_fragment_shapes(
        fragments_path,
        raw_barcodes,
        peak_records(var),
        right_cut_offset=int(spec.right_cut_offset),
        bins=bin_records(spec),
        chunk_size=COUNT_CHUNK_SIZE,
    )
    obs = coordinates.loc[canonical_barcodes].copy()
    obs.insert(0, "fragment_barcode", raw_barcodes)
    validation = ad.AnnData(
        X=result.X.copy(),
        obs=obs,
        var=var.copy(),
    )
    validation.obsm["spatial"] = obs[
        ["pixel_col_fullres", "pixel_row_fullres"]
    ].to_numpy(dtype=np.float64)
    validation.uns["validation_modality"] = {
        "schema_version": 1,
        "assay": str(record["assay"]),
        "gsm": str(record["gsm"]),
        "interpretation": "orthogonal_epigenome_validation_only",
        "not_a_shapemix_input": True,
        "not_composition_truth": True,
        "source_sha256": str(record["normalized_sha256"]),
        "preprocessing_counters": result.qc.to_dict(),
        "feature_sha256": ordered_feature_sha256(var.index.astype(str)),
    }
    return validation, result.qc.to_dict()


def copy_validation_h5ad(
    source: Path,
    destination: Path,
    spatial_barcodes: pd.Index,
) -> dict[str, Any]:
    require_file(source, "preprocessed validation H5AD")
    value = ad.read_h5ad(source, backed="r")
    try:
        barcodes = value.obs_names.astype(str)
        observations = int(value.n_obs)
        features = int(value.n_vars)
    finally:
        value.file.close()
    if barcodes.duplicated().any():
        raise ValueError(f"Validation H5AD has duplicate barcodes: {repository_path(source)}")
    overlap = len(spatial_barcodes.intersection(barcodes))
    if overlap == 0:
        raise ValueError(f"Validation H5AD has no spatial-barcode overlap: {repository_path(source)}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return {
        "path": repository_path(destination),
        "sha256": file_digest(destination),
        "observations": observations,
        "features": features,
        "spatial_barcode_overlap": overlap,
    }


def write_features(path: Path, names: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name in names:
            handle.write(f"{name}\n")


def descriptor(
    template: Mapping[str, Any],
    section: Mapping[str, Any],
    reference_path: Path,
    reference_manifest: Mapping[str, Any],
    spatial_path: Path,
    feature_path: Path,
    fragment_shape: FragmentShapeSpec,
    validation: Mapping[str, Any],
) -> dict[str, Any]:
    dataset_id = str(section["dataset_id"])
    accession = "GSE205055" if "gse205055" in str(template["template_id"]) else "GSE263333"
    evaluation_design: dict[str, Any] = {
        "comparison_key": "section",
        "group": str(section["group"]),
        "atac_gsm": str(section["atac_gsm"]),
        "reference_id": str(section["reference_id"]),
        "truth_limitation": (
            "No exact per-spot cell-type composition truth is available; RNA, protein, "
            "histone, anatomy, and spatial-continuity endpoints are orthogonal validation only."
        ),
        "reference_caveat": section.get("reference_caveat"),
    }
    evaluation_design = {key: value for key, value in evaluation_design.items() if value is not None}
    declared_shape = fragment_shape.to_dict(omit_none=True)
    # Per-object source, split, software, and counter provenance differs
    # legitimately between the reference and spatial H5ADs. The dataset-level
    # declaration contains only fields that form their aligned model contract.
    for key in (
        "source_sha256",
        "split_sha256",
        "software_versions",
        "preprocessing_counters",
        "matrix_counters",
    ):
        declared_shape.pop(key, None)
    return {
        "dataset_id": dataset_id,
        "source": f"{accession} spatial ATAC section {section['atac_gsm']}",
        "description": (
            "Real spatial ATAC section projected onto a frozen, reference-only 5,000-feature "
            "axis for ShapeMix prediction and orthogonal cross-modality validation."
        ),
        "labels_key": "cell_type",
        "spatial_key": "spatial",
        "benchmark_scope": "real_spatial_orthogonal_validation",
        "evaluation_design": evaluation_design,
        "modalities": {
            "atac": {
                "reference": {"path": repository_path(reference_path)},
                "spatial": {"path": repository_path(spatial_path)},
                "labels_key": "cell_type",
                "spatial_key": "spatial",
                "cell_types": [str(value) for value in reference_manifest["cell_types"]],
                "fragment_shape": declared_shape,
                "feature_sets": {
                    "all": {"mode": "all"},
                    "selected_reference_features": {"path": repository_path(feature_path)},
                },
            }
        },
        "validation": copy.deepcopy(dict(validation)),
    }


def section_paths(template: Mapping[str, Any]) -> tuple[Path, Path, Path]:
    family_root = ROOT / str(template["reusable_root"])
    family_name = family_root.name
    work_root = ROOT / "data/work/preprocessing" / family_name / "source_preprocessing_v1"
    audit_path = family_root / "source_audit/fragments.yaml"
    return family_root, work_root, audit_path


def materialize_section(template: Mapping[str, Any], section: Mapping[str, Any]) -> Path:
    dataset_id = str(section["dataset_id"])
    final_root = ROOT / "data/processed/datasets" / dataset_id
    completed = (final_root / "manifest.yaml").is_file() and (final_root / "dataset.yaml").is_file()
    if final_root.exists():
        if not completed:
            raise FileExistsError(f"Partial immutable dataset directory: {repository_path(final_root)}")
        print(f"real_spatial dataset={dataset_id} status=reused", flush=True)
        return final_root

    family_root, work_root, audit_path = section_paths(template)
    audit = read_yaml(require_file(audit_path, "fragment source audit"))
    alignment_path, alignment = alignment_manifest(template, section)
    barcode_policy = dict(alignment.get("barcode_policy", {}))
    atac_gsm = str(section["atac_gsm"])
    atac_records = fragment_records(audit, atac_gsm, "atac")
    if len(atac_records) != 1:
        raise ValueError(f"Expected one ATAC fragment record for {atac_gsm}; found {len(atac_records)}")
    atac_record = atac_records[0]
    fragments_path = validate_audited_path(
        atac_record,
        "normalized",
        "normalized_bytes",
        expected_parent=family_root / "normalized_fragments",
        label="normalized ATAC fragments",
    )
    validate_audited_path(
        atac_record,
        "tabix",
        "tabix_bytes",
        expected_parent=family_root / "normalized_fragments",
        label="ATAC tabix index",
    )
    raw_barcodes = read_barcodes(atac_record)
    canonical_barcodes = canonicalize_barcodes(raw_barcodes, barcode_policy)
    coordinates = load_coordinates(work_root, atac_gsm, canonical_barcodes)
    reference_path, reference_manifest, var, reference_spec = reference_inputs(
        str(section["reference_id"])
    )
    provenance = shape_provenance(
        atac_record,
        audit_path,
        alignment_path,
        var,
        reference_spec,
    )
    spatial, atac_qc = build_atac_spatial(
        atac_record,
        fragments_path,
        raw_barcodes,
        canonical_barcodes,
        coordinates,
        var,
        reference_spec,
        provenance,
    )
    spatial_spec = FragmentShapeSpec.from_mapping(spatial.uns["fragment_shape"])

    final_root.parent.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(tempfile.mkdtemp(prefix=f".{dataset_id}.", dir=final_root.parent))
    try:
        spatial_path = temporary_root / "atac/spatial.h5ad"
        spatial_path.parent.mkdir(parents=True)
        spatial.write_h5ad(spatial_path, compression="gzip")
        feature_path = temporary_root / "atac/features/selected_reference_features.txt"
        write_features(feature_path, spatial.var_names.astype(str))

        validation: dict[str, Any] = {
            "truth_status": "orthogonal_validation_not_exact_composition_truth",
            "exact_truth": False,
            "alignment_manifest": {
                "path": repository_path(final_root / "validation/alignment.yaml"),
                "sha256": file_digest(alignment_path),
            },
            "epigenome": [],
        }
        alignment_destination = temporary_root / "validation/alignment.yaml"
        alignment_destination.parent.mkdir(parents=True)
        shutil.copy2(alignment_path, alignment_destination)

        rna_gsm = section.get("rna_gsm")
        if rna_gsm:
            rna_source = work_root / "validation_modalities/rna" / str(rna_gsm) / "rna.h5ad"
            validation["rna"] = {
                "gsm": str(rna_gsm),
                **copy_validation_h5ad(
                    rna_source,
                    temporary_root / "validation/rna/rna.h5ad",
                    spatial.obs_names,
                ),
            }
            validation["rna"]["path"] = repository_path(final_root / "validation/rna/rna.h5ad")

        protein_gsm = section.get("protein_gsm")
        if protein_gsm:
            protein_source = work_root / "validation_modalities/protein" / str(protein_gsm) / "protein.h5ad"
            validation["protein"] = {
                "gsm": str(protein_gsm),
                **copy_validation_h5ad(
                    protein_source,
                    temporary_root / "validation/protein/protein.h5ad",
                    spatial.obs_names,
                ),
            }
            validation["protein"]["path"] = repository_path(final_root / "validation/protein/protein.h5ad")

        epigenome_qc: dict[str, Any] = {}
        for gsm_value in section.get("validation_epigenome_gsms", []):
            gsm = str(gsm_value)
            for record in fragment_records(audit, gsm):
                if str(record.get("assay")) == "atac":
                    continue
                source = validate_audited_path(
                    record,
                    "normalized",
                    "normalized_bytes",
                    expected_parent=work_root / "validation_modalities/epigenome",
                    label="normalized epigenome fragments",
                )
                validate_audited_path(
                    record,
                    "tabix",
                    "tabix_bytes",
                    expected_parent=work_root / "validation_modalities/epigenome",
                    label="epigenome tabix index",
                )
                epigenome_raw = read_barcodes(record)
                epigenome_canonical = canonicalize_barcodes(epigenome_raw, barcode_policy)
                absent = sorted(set(epigenome_canonical).difference(coordinates.index))
                if absent:
                    raise ValueError(
                        f"Epigenome barcodes are not aligned to ATAC coordinates for {dataset_id}: {absent[0]}"
                    )
                validation_object, qc = build_epigenome_validation(
                    record,
                    source,
                    epigenome_raw,
                    epigenome_canonical,
                    coordinates,
                    var,
                    reference_spec,
                )
                assay = str(record["assay"])
                filename = f"{assay}__{gsm}.h5ad"
                destination = temporary_root / "validation/epigenome" / filename
                destination.parent.mkdir(parents=True, exist_ok=True)
                validation_object.write_h5ad(destination, compression="gzip")
                final_destination = final_root / "validation/epigenome" / filename
                validation["epigenome"].append(
                    {
                        "gsm": gsm,
                        "assay": assay,
                        "path": repository_path(final_destination),
                        "sha256": file_digest(destination),
                        "spots": validation_object.n_obs,
                        "features": validation_object.n_vars,
                        "interpretation": "orthogonal_validation_only",
                    }
                )
                epigenome_qc[f"{gsm}:{assay}"] = qc

        descriptor_value = descriptor(
            template,
            section,
            reference_path,
            reference_manifest,
            final_root / "atac/spatial.h5ad",
            final_root / "atac/features/selected_reference_features.txt",
            spatial_spec,
            validation,
        )
        atomic_yaml(temporary_root / "dataset.yaml", descriptor_value)
        atomic_yaml(
            temporary_root / "manifest.yaml",
            {
                "schema_version": 1,
                "status": "complete",
                "completed_at": datetime.now(timezone.utc).isoformat(),
                "dataset_id": dataset_id,
                "template": repository_path(
                    ROOT / "configs/datasets" / f"{template['template_id']}.yaml"
                ),
                "source_accession": str(audit["accession"]),
                "source_fragment_audit": repository_path(audit_path),
                "source_fragment_audit_sha256": file_digest(audit_path),
                "reference_id": str(section["reference_id"]),
                "reference_manifest_sha256": file_digest(reference_path.parents[1] / "reference.yaml"),
                "atac_gsm": atac_gsm,
                "atac_source_sha256": str(atac_record["normalized_sha256"]),
                "spots": spatial.n_obs,
                "features": spatial.n_vars,
                "atac_preprocessing_counters": atac_qc,
                "epigenome_preprocessing_counters": epigenome_qc,
                "outputs": {
                    "dataset": repository_path(final_root / "dataset.yaml"),
                    "spatial": repository_path(final_root / "atac/spatial.h5ad"),
                    "features": repository_path(
                        final_root / "atac/features/selected_reference_features.txt"
                    ),
                },
                "truth_used": False,
                "outcome_data_used_for_reference_or_feature_selection": False,
            },
        )
        temporary_root.rename(final_root)
    except BaseException:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    print(
        f"real_spatial dataset={dataset_id} spots={spatial.n_obs} features={spatial.n_vars} "
        "status=completed",
        flush=True,
    )
    return final_root


def register_datasets(roots: Iterable[Path]) -> None:
    registry = read_yaml(REGISTRY_PATH) if REGISTRY_PATH.is_file() else {}
    changed = False
    for root in roots:
        dataset_id = root.name
        expected = {"config": repository_path(root / "dataset.yaml")}
        current = registry.get(dataset_id)
        if current is not None and current != expected:
            raise ValueError(f"Registry entry conflicts for {dataset_id}: {current}")
        if current is None:
            registry[dataset_id] = expected
            changed = True
    if changed:
        atomic_yaml(REGISTRY_PATH, registry)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--template",
        action="append",
        type=Path,
        help="Template YAML; may be repeated. Defaults to both frozen spatial templates.",
    )
    parser.add_argument("--dataset-id", action="append", help="Materialize only selected dataset IDs.")
    return parser.parse_args()


def main() -> None:
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run this materializer through scripts/run_shapemix_low_impact.sh")
    args = parse_args()
    template_paths = tuple(args.template) if args.template else DEFAULT_TEMPLATES
    selected_ids = set(args.dataset_id or [])
    roots: list[Path] = []
    observed_ids: set[str] = set()
    for path in template_paths:
        path = path if path.is_absolute() else ROOT / path
        template = read_yaml(require_file(path, "real-spatial template"))
        for section in template.get("sections", []):
            dataset_id = str(section["dataset_id"])
            observed_ids.add(dataset_id)
            if selected_ids and dataset_id not in selected_ids:
                continue
            roots.append(materialize_section(template, section))
    missing = sorted(selected_ids.difference(observed_ids))
    if missing:
        raise ValueError(f"Requested dataset IDs are absent from the templates: {missing}")
    register_datasets(roots)
    print(f"real_spatial registered={len(roots)} status=complete", flush=True)


if __name__ == "__main__":
    main()
