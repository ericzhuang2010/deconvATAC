#!/usr/bin/env python3
"""Update spatial preprocessing audits after the canonical same-disk path migration."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import yaml


ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ("gse205055_spatial", "gse263333_spatial_mux")


def replacements(family: str) -> dict[str, str]:
    processed = f"data/processed/shapemix/{family}"
    work = f"data/work/preprocessing/{family}/source_preprocessing_v1"
    return {
        f"{processed}/extracted_payload/": f"{work}/extracted_payload/",
        f"{processed}/normalized_atac_fragments/": f"{processed}/normalized_fragments/",
        f"{processed}/validation_modalities/": f"{work}/validation_modalities/",
        f"{processed}/spatial_coordinates/": f"{work}/spatial_coordinates/",
        f"{processed}/cross_modality_alignment/": (
            f"{processed}/manifests/cross_modality_alignment/"
        ),
    }


def rewrite_value(value: Any, mapping: Mapping[str, str]) -> tuple[Any, list[tuple[str, str]]]:
    changes: list[tuple[str, str]] = []
    if isinstance(value, str):
        for old, new in mapping.items():
            if value.startswith(old):
                replacement = new + value[len(old) :]
                return replacement, [(value, replacement)]
        return value, changes
    if isinstance(value, list):
        output = []
        for item in value:
            updated, item_changes = rewrite_value(item, mapping)
            output.append(updated)
            changes.extend(item_changes)
        return output, changes
    if isinstance(value, dict):
        output = {}
        for key, item in value.items():
            updated, item_changes = rewrite_value(item, mapping)
            output[key] = updated
            changes.extend(item_changes)
        return output, changes
    return value, changes


def atomic_yaml(path: Path, value: Mapping[str, Any]) -> None:
    descriptor, temporary_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(descriptor, "w") as handle:
            yaml.safe_dump(dict(value), handle, sort_keys=False)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def migrate_yaml(path: Path, mapping: Mapping[str, str]) -> list[tuple[str, str]]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    updated, changes = rewrite_value(value, mapping)
    if changes:
        validate_targets(changes)
        atomic_yaml(path, updated)
    return changes


def migrate_inventory(path: Path, mapping: Mapping[str, str]) -> list[tuple[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames
    if not fieldnames or "path" not in fieldnames:
        raise ValueError(f"Payload inventory lacks a path column: {path}")
    changes: list[tuple[str, str]] = []
    for row in rows:
        updated, row_changes = rewrite_value(row["path"], mapping)
        row["path"] = updated
        changes.extend(row_changes)
    if changes:
        validate_targets(changes)
        descriptor, temporary_name = tempfile.mkstemp(
            dir=path.parent, prefix=f".{path.name}."
        )
        try:
            with os.fdopen(descriptor, "w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            os.replace(temporary_name, path)
        except BaseException:
            Path(temporary_name).unlink(missing_ok=True)
            raise
    return changes


def validate_targets(changes: list[tuple[str, str]]) -> None:
    missing = sorted({new for _, new in changes if not (ROOT / new).exists()})
    if missing:
        raise FileNotFoundError(f"Migrated manifest target(s) are absent: {missing[:10]}")


def migrate_family(family: str) -> dict[str, Any]:
    family_root = ROOT / "data" / "processed" / "shapemix" / family
    mapping = replacements(family)
    yaml_paths = sorted(
        path
        for root in (family_root / "source_audit", family_root / "manifests")
        if root.exists()
        for path in root.rglob("*.yaml")
        if path.name != "layout_migration_v1.yaml"
    )
    changes_by_file: dict[str, int] = {}
    all_changes: list[tuple[str, str]] = []
    for path in yaml_paths:
        changes = migrate_yaml(path, mapping)
        if changes:
            changes_by_file[str(path.relative_to(ROOT))] = len(changes)
            all_changes.extend(changes)
    inventory = family_root / "source_audit" / "payload_inventory.tsv"
    if inventory.is_file():
        changes = migrate_inventory(inventory, mapping)
        if changes:
            changes_by_file[str(inventory.relative_to(ROOT))] = len(changes)
            all_changes.extend(changes)
    validate_targets(all_changes)
    record = {
        "schema_version": 1,
        "status": "complete",
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "family": family,
        "mapping": mapping,
        "changed_files": changes_by_file,
        "changed_path_references": len(all_changes),
        "validation": "every rewritten target exists",
    }
    atomic_yaml(family_root / "manifests" / "layout_migration_v1.yaml", record)
    return record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", action="append", choices=FAMILIES)
    return parser.parse_args()


def main() -> None:
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run this migration through run_shapemix_low_impact.sh")
    args = parse_args()
    for family in args.family or FAMILIES:
        record = migrate_family(family)
        print(
            f"layout_migration family={family} changed_paths="
            f"{record['changed_path_references']} status=completed",
            flush=True,
        )


if __name__ == "__main__":
    main()
