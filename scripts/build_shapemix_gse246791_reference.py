#!/usr/bin/env python3
"""Build frozen broad labels for the GSE246791 adult mouse-brain reference."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import shutil
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import h5py
import numpy as np
import yaml

from scripts.prepare_shapemix_gse246791 import CONFIG, LOCK_PATH, selected_sources


ROOT = Path(__file__).resolve().parents[1]
LABEL_VERSION = "broad9_v1"
CELL_TYPES = (
    "Excitatory neurons",
    "Inhibitory neurons",
    "Other neurons",
    "Astroglia/ependymal",
    "Oligodendrocytes",
    "OPCs",
    "Microglia/macrophages",
    "Other immune",
    "Vascular/stromal",
)
MIN_CELLS_PER_TYPE = {
    cell_type: 20 if cell_type == "Other immune" else 100
    for cell_type in CELL_TYPES
}
METADATA_PATH = (
    ROOT
    / "data/raw/sources/ncbi_geo/GSE246791/portal_metadata/SI_Table_2_nuclei.txt"
)
EXPECTED_METADATA_COLUMNS = (
    "CellID",
    "Sample",
    "Barcode",
    "# of Fragments",
    "TSSe",
    "L1",
    "L2",
    "L3",
    "L4",
    "pL4",
    "NeuronTransmitter",
    "Subclass",
)


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


def locked_metadata_record() -> Mapping[str, Any]:
    lock = load_yaml(LOCK_PATH)
    matches = [
        record
        for record in lock.get("files", [])
        if isinstance(record, Mapping) and record.get("name") == METADATA_PATH.name
    ]
    if len(matches) != 1:
        raise ValueError("Expected one locked GSE246791 cell-metadata record")
    return matches[0]


def subclass_id(value: str) -> int:
    token = value.strip().split(" ", 1)[0]
    try:
        result = int(token)
    except ValueError as exc:
        raise ValueError(f"Subclass lacks a numeric author ID: {value!r}") from exc
    if result < 1 or result > 339:
        raise ValueError(f"Subclass author ID is out of range: {result}")
    return result


def broad_label(neurotransmitter: str, subclass: str) -> str | None:
    identifier = subclass_id(subclass)
    transmitter = neurotransmitter.strip()
    lowered = transmitter.lower()
    if transmitter not in {"NN", "LQ", ""}:
        if "glut" in lowered:
            return "Excitatory neurons"
        if "gaba" in lowered or "gly" in lowered:
            return "Inhibitory neurons"
        return "Other neurons"
    if 316 <= identifier <= 325:
        return "Astroglia/ependymal"
    if identifier == 326:
        return "OPCs"
    if identifier in {327, 328}:
        return "Oligodendrocytes"
    if 329 <= identifier <= 333:
        return "Vascular/stromal"
    if identifier in {334, 335}:
        return "Microglia/macrophages"
    if 336 <= identifier <= 338:
        return "Other immune"
    if identifier == 339:
        return None
    raise ValueError(
        f"Non-neuronal transmitter {transmitter!r} has unexpected subclass {subclass!r}"
    )


def h5ad_barcodes_and_counts(path: Path) -> tuple[list[str], np.ndarray]:
    with h5py.File(path, "r") as handle:
        barcodes = [
            value.decode("ascii") if isinstance(value, bytes) else str(value)
            for value in handle["obs/index"][:]
        ]
        counts = np.asarray(handle["obs/n_fragment"][:], dtype=np.int64)
    if len(barcodes) != len(counts) or len(barcodes) != len(set(barcodes)):
        raise ValueError(f"Invalid deposited H5AD observation axis: {path}")
    if np.any(counts <= 0):
        raise ValueError(f"Deposited H5AD contains nonpositive fragment counts: {path}")
    return barcodes, counts


def build_labels() -> Path:
    family_root = ROOT / str(CONFIG["processed_directory"])
    output_root = family_root / "labels" / LABEL_VERSION
    cells_path = output_root / "cells.tsv.gz"
    manifest_path = output_root / "manifest.yaml"
    if cells_path.is_file() and manifest_path.is_file():
        print("gse246791 labels status=reused", flush=True)
        return cells_path
    if output_root.exists():
        raise FileExistsError(f"Partial immutable label directory: {output_root}")

    sources = selected_sources()
    sample_metadata: dict[str, dict[str, Any]] = {}
    expected: dict[tuple[str, str], int] = {}
    for source in sources:
        if not source.output_path.is_file():
            raise FileNotFoundError(
                f"Run prepare_shapemix_gse246791.py first: {source.output_path}"
            )
        barcodes, counts = h5ad_barcodes_and_counts(source.output_path)
        sample_metadata[source.sample] = {
            "gsm": source.gsm,
            "major_region": source.major_region,
            "h5ad": source.output_path,
            "barcodes": set(barcodes),
            "cells": len(barcodes),
        }
        expected.update(
            {
                (source.sample, barcode): int(count)
                for barcode, count in zip(barcodes, counts, strict=True)
            }
        )
    if len(expected) != sum(value["cells"] for value in sample_metadata.values()):
        raise ValueError("Selected sample-plus-barcode keys are not unique")

    metadata_record = locked_metadata_record()
    if (
        not METADATA_PATH.is_file()
        or METADATA_PATH.stat().st_size != int(metadata_record["bytes"])
        or file_digest(METADATA_PATH) != str(metadata_record["sha256"])
    ):
        raise ValueError("Locked GSE246791 cell metadata changed or is absent")

    output_root.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{LABEL_VERSION}.", dir=output_root.parent))
    retained_keys: set[tuple[str, str]] = set()
    global_support: Counter[str] = Counter()
    sample_support: Counter[tuple[str, str]] = Counter()
    excluded_lq = 0
    fragment_relative_errors: list[float] = []
    try:
        with METADATA_PATH.open(newline="") as input_handle, gzip.open(
            temporary / "cells.tsv.gz", "wt", newline=""
        ) as output_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != EXPECTED_METADATA_COLUMNS:
                raise ValueError(f"Unexpected SI Table 2 columns: {reader.fieldnames}")
            fieldnames = (
                "cell_id",
                "sample",
                "gsm",
                "major_region",
                "barcode",
                "cell_type",
                "author_fragments",
                "h5ad_n_fragment",
                "tsse",
                "l1",
                "l2",
                "l3",
                "l4",
                "pl4",
                "neuron_transmitter",
                "author_subclass",
            )
            writer = csv.DictWriter(output_handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row_number, row in enumerate(reader, start=2):
                sample = str(row["Sample"])
                metadata = sample_metadata.get(sample)
                if metadata is None:
                    continue
                barcode = str(row["Barcode"])
                key = (sample, barcode)
                h5_count = expected.get(key)
                if h5_count is None:
                    continue
                if key in retained_keys:
                    raise ValueError(f"Duplicate selected metadata key at row {row_number}: {key}")
                retained_keys.add(key)
                label = broad_label(str(row["NeuronTransmitter"]), str(row["Subclass"]))
                if label is None:
                    excluded_lq += 1
                    continue
                author_count = int(str(row["# of Fragments"]))
                if author_count <= 0:
                    raise ValueError(f"Nonpositive author fragment count at row {row_number}")
                fragment_relative_errors.append(abs(author_count - h5_count) / h5_count)
                cell_id = str(row["CellID"])
                expected_cell_id = f"{sample}.{barcode}"
                if cell_id != expected_cell_id:
                    raise ValueError(
                        f"Author cell ID does not equal sample.barcode at row {row_number}"
                    )
                writer.writerow(
                    {
                        "cell_id": cell_id,
                        "sample": sample,
                        "gsm": metadata["gsm"],
                        "major_region": metadata["major_region"],
                        "barcode": barcode,
                        "cell_type": label,
                        "author_fragments": author_count,
                        "h5ad_n_fragment": h5_count,
                        "tsse": row["TSSe"],
                        "l1": row["L1"],
                        "l2": row["L2"],
                        "l3": row["L3"],
                        "l4": row["L4"],
                        "pl4": row["pL4"],
                        "neuron_transmitter": row["NeuronTransmitter"],
                        "author_subclass": row["Subclass"],
                    }
                )
                global_support[label] += 1
                sample_support[(sample, label)] += 1

        missing_keys = sorted(set(expected).difference(retained_keys))
        if missing_keys:
            raise ValueError(
                f"{len(missing_keys)} deposited H5AD cells lack exact SI Table 2 annotations"
            )
        failures = {
            cell_type: global_support[cell_type]
            for cell_type in CELL_TYPES
            if global_support[cell_type] < MIN_CELLS_PER_TYPE[cell_type]
        }
        if failures:
            raise ValueError(f"Adult-brain broad-label support gate failed: {failures}")
        with (temporary / "support.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("sample", "gsm", "major_region", "cell_type", "cells"))
            for source in sources:
                for cell_type in CELL_TYPES:
                    writer.writerow(
                        (
                            source.sample,
                            source.gsm,
                            source.major_region,
                            cell_type,
                            sample_support[(source.sample, cell_type)],
                        )
                    )
        errors = np.asarray(fragment_relative_errors, dtype=float)
        manifest = {
            "schema_version": 1,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "label_version": LABEL_VERSION,
            "cell_types": list(CELL_TYPES),
            "minimum_cells_per_type": dict(MIN_CELLS_PER_TYPE),
            "retained_cells": int(sum(global_support.values())),
            "excluded_lq_cells": excluded_lq,
            "samples": len(sources),
            "support": {cell_type: global_support[cell_type] for cell_type in CELL_TYPES},
            "mapping": {
                "neuronal": (
                    "author NeuronTransmitter: Glut -> excitatory; Gaba/Gly -> "
                    "inhibitory; remaining non-NN transmitters -> other neurons"
                ),
                "non_neuronal_author_subclass_ids": {
                    "316-325": "Astroglia/ependymal",
                    "326": "OPCs",
                    "327-328": "Oligodendrocytes",
                    "329-333": "Vascular/stromal",
                    "334-335": "Microglia/macrophages",
                    "336-338": "Other immune",
                    "339": "excluded low quality",
                },
                "outcome_data_used": False,
            },
            "fragment_count_crosscheck": {
                "median_absolute_relative_error": float(np.median(errors)),
                "p95_absolute_relative_error": float(np.quantile(errors, 0.95)),
            },
            "inputs": {
                "metadata": str(METADATA_PATH.relative_to(ROOT)),
                "metadata_sha256": str(metadata_record["sha256"]),
                "h5ad_source_objects": [
                    str(source.output_path.relative_to(ROOT)) for source in sources
                ],
            },
            "outputs": {"cells": "cells.tsv.gz", "support": "support.tsv"},
        }
        with (temporary / "manifest.yaml").open("w") as handle:
            yaml.safe_dump(manifest, handle, sort_keys=False)
        temporary.rename(output_root)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(
        f"gse246791 labels cells={sum(global_support.values())} "
        f"excluded_lq={excluded_lq} status=completed",
        flush=True,
    )
    return cells_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("labels",), default="labels")
    return parser.parse_args()


def main() -> None:
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run this builder through run_shapemix_low_impact.sh")
    parse_args()
    build_labels()


if __name__ == "__main__":
    main()
