#!/usr/bin/env python3
"""Stream the GSE216371 workbook and audit all author-annotated E13.5 cells."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import shutil
import tempfile
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator, Mapping

import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs/data_sources/shapemix_gse216371_reference.yaml"
LOCK_PATH = ROOT / "configs/data_sources/shapemix_gse216371_reference_lock.yaml"
LABEL_VERSION = "author_e13_5_v1"
CCRE_VERSION = "author_ccres_v1"
REQUIRED_COLUMNS = (
    "Cell_Barcodes",
    "Developmental_stages",
    "Embryo_ID",
    "Round4_Barcode",
    "FRIP",
    "Fragments",
    "Main_cluster_name",
    "Main_cluster_number",
    "Sub_cluster",
    "Sub_cluster_Annotation",
)
OUTPUT_COLUMNS = (
    "cell_id",
    "developmental_stage",
    "embryo_id",
    "round4_barcode",
    "frip",
    "fragments",
    "tss_enrichment",
    "author_main_cluster",
    "author_main_cluster_number",
    "author_subcluster",
    "author_subcluster_annotation",
)
CCRE_COLUMNS = (
    "cCREs_name",
    "cCREs_type",
    "nearest_TSS",
    "distance_To_TSS (bp)",
    "CTCF binding",
    "94_Modules_in_Fig2g",
    "merged_Module_for_GREAT_and_Motif analysis_in_Fig2ij",
)
EXPECTED_CCRE_COUNT = 830_873


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


def project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def excel_column_index(reference: str) -> int:
    letters = "".join(character for character in reference if character.isalpha())
    if not letters:
        raise ValueError(f"Cell reference has no column letters: {reference!r}")
    result = 0
    for character in letters.upper():
        if not ("A" <= character <= "Z"):
            raise ValueError(f"Invalid Excel column reference: {reference!r}")
        result = result * 26 + ord(character) - ord("A") + 1
    return result - 1


def load_shared_strings(workbook: zipfile.ZipFile) -> list[str]:
    values: list[str] = []
    with workbook.open("xl/sharedStrings.xml") as handle:
        for _, element in ET.iterparse(handle, events=("end",)):
            if local_name(element.tag) != "si":
                continue
            values.append(
                "".join(
                    node.text or ""
                    for node in element.iter()
                    if local_name(node.tag) == "t"
                )
            )
            element.clear()
    return values


def cell_text(cell: ET.Element, shared_strings: list[str]) -> str:
    value_node = next(
        (node for node in cell.iter() if local_name(node.tag) == "v"),
        None,
    )
    if cell.attrib.get("t") == "inlineStr":
        return "".join(
            node.text or "" for node in cell.iter() if local_name(node.tag) == "t"
        )
    if value_node is None or value_node.text is None:
        return ""
    value = value_node.text
    if cell.attrib.get("t") == "s":
        index = int(value)
        if index < 0 or index >= len(shared_strings):
            raise ValueError(f"Shared-string index is out of range: {index}")
        return shared_strings[index]
    return value


def iter_sheet_rows(
    workbook: zipfile.ZipFile,
    member: str,
    shared_strings: list[str],
) -> Iterator[tuple[int, dict[int, str]]]:
    with workbook.open(member) as handle:
        for _, element in ET.iterparse(handle, events=("end",)):
            if local_name(element.tag) != "row":
                continue
            row_number = int(element.attrib.get("r", "0"))
            values: dict[int, str] = {}
            for cell in element:
                if local_name(cell.tag) != "c":
                    continue
                reference = cell.attrib.get("r")
                if reference is None:
                    raise ValueError(f"Row {row_number} contains a cell without a reference")
                values[excel_column_index(reference)] = cell_text(cell, shared_strings)
            yield row_number, values
            element.clear()


def normalize_stage(value: str) -> str:
    normalized = value.strip().upper().replace(" ", "")
    if normalized in {"E13.5", "13.5", "E13.50", "13.50"}:
        return "E13.5"
    return value.strip()


def parse_ccre_name(value: str) -> tuple[str, int, int, str]:
    try:
        chrom, start_text, end_text = value.strip().rsplit("_", 2)
        start, end = int(start_text), int(end_text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid author cCRE identifier: {value!r}") from exc
    if not chrom.startswith("chr") or start < 0 or end <= start:
        raise ValueError(f"Invalid author cCRE coordinates: {value!r}")
    return chrom, start, end, f"{chrom}:{start}-{end}"


def workbook_source(config: Mapping[str, Any]) -> Path:
    return (
        project_path(str(config["raw_directory"]))
        / "processed_downloads"
        / "GSE216371_Cell_annotation_and_cCREs_in_MOPA.xlsx"
    )


def locked_source_record(path: Path) -> Mapping[str, Any]:
    lock = load_yaml(LOCK_PATH)
    matches = [
        record
        for record in lock.get("files", [])
        if isinstance(record, Mapping) and record.get("name") == path.name
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one source-lock record for {path.name}")
    return matches[0]


def audit_labels(config: Mapping[str, Any]) -> Path:
    family_root = project_path(str(config["processed_directory"]))
    output_root = family_root / "labels" / LABEL_VERSION
    cells_path = output_root / "cells.tsv.gz"
    manifest_path = output_root / "manifest.yaml"
    if cells_path.is_file() and manifest_path.is_file():
        print("gse216371 labels status=reused", flush=True)
        return cells_path
    if output_root.exists():
        raise FileExistsError(f"Partial immutable label audit: {output_root}")

    source = workbook_source(config)
    record = locked_source_record(source)
    if not source.is_file() or source.stat().st_size != int(record["bytes"]):
        raise ValueError(f"Locked workbook byte count changed or is absent: {source}")
    if file_digest(source) != str(record["sha256"]):
        raise ValueError(f"Locked workbook SHA-256 changed: {source}")

    output_root.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{LABEL_VERSION}.", dir=output_root.parent))
    seen: set[str] = set()
    main_support: Counter[str] = Counter()
    subcluster_support: Counter[tuple[str, str]] = Counter()
    stage_support: Counter[str] = Counter()
    retained = 0
    missing_annotation = 0
    try:
        with zipfile.ZipFile(source) as workbook:
            shared_strings = load_shared_strings(workbook)
            rows = iter_sheet_rows(workbook, "xl/worksheets/sheet1.xml", shared_strings)
            header_row, header_values = next(rows)
            headers = {value.strip(): index for index, value in header_values.items()}
            missing = [name for name in REQUIRED_COLUMNS if name not in headers]
            tss_columns = [name for name in headers if name.upper().startswith("TSS")]
            if header_row != 1 or missing or len(tss_columns) != 1:
                raise ValueError(
                    f"Unexpected GSE216371 metadata header; row={header_row} "
                    f"missing={missing} tss_columns={tss_columns}"
                )
            with gzip.open(temporary / "cells.tsv.gz", "wt", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
                writer.writeheader()
                for row_number, values in rows:
                    stage = normalize_stage(values.get(headers["Developmental_stages"], ""))
                    stage_support[stage or "<missing>"] += 1
                    if stage != "E13.5":
                        continue
                    cell_id = values.get(headers["Cell_Barcodes"], "").strip()
                    main_cluster = values.get(headers["Main_cluster_name"], "").strip()
                    annotation = values.get(headers["Sub_cluster_Annotation"], "").strip()
                    if not cell_id:
                        raise ValueError(f"E13.5 row {row_number} has no cell barcode")
                    if cell_id in seen:
                        raise ValueError(f"Duplicate E13.5 cell barcode: {cell_id}")
                    seen.add(cell_id)
                    if not main_cluster or not annotation:
                        missing_annotation += 1
                        continue
                    row = {
                        "cell_id": cell_id,
                        "developmental_stage": stage,
                        "embryo_id": values.get(headers["Embryo_ID"], "").strip(),
                        "round4_barcode": values.get(headers["Round4_Barcode"], "").strip(),
                        "frip": values.get(headers["FRIP"], "").strip(),
                        "fragments": values.get(headers["Fragments"], "").strip(),
                        "tss_enrichment": values.get(headers[tss_columns[0]], "").strip(),
                        "author_main_cluster": main_cluster,
                        "author_main_cluster_number": values.get(
                            headers["Main_cluster_number"], ""
                        ).strip(),
                        "author_subcluster": values.get(headers["Sub_cluster"], "").strip(),
                        "author_subcluster_annotation": annotation,
                    }
                    writer.writerow(row)
                    main_support[main_cluster] += 1
                    subcluster_support[(main_cluster, annotation)] += 1
                    retained += 1
                    if retained % 100_000 == 0:
                        print(f"gse216371 label_scan retained={retained}", flush=True)

        if retained == 0 or not main_support:
            raise ValueError("No annotated E13.5 cells were retained from the workbook")
        with (temporary / "main_cluster_support.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("author_main_cluster", "cells"))
            writer.writerows(sorted(main_support.items()))
        with (temporary / "subcluster_support.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                ("author_main_cluster", "author_subcluster_annotation", "cells")
            )
            for (main_cluster, annotation), count in sorted(subcluster_support.items()):
                writer.writerow((main_cluster, annotation, count))
        manifest = {
            "schema_version": 1,
            "status": "author_labels_audited_broad_ontology_pending",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "label_version": LABEL_VERSION,
            "retained_stage": "E13.5",
            "retained_annotated_cells": retained,
            "e13_5_cells_missing_author_annotation": missing_annotation,
            "main_clusters": len(main_support),
            "subcluster_annotations": len(subcluster_support),
            "stage_support": dict(sorted(stage_support.items())),
            "source": str(source.relative_to(ROOT)),
            "source_bytes": source.stat().st_size,
            "source_sha256": str(record["sha256"]),
            "parameters": {
                "sheet": "MOPA cell metadata",
                "sheet_member": "xl/worksheets/sheet1.xml",
                "streaming_xml": True,
                "outcome_data_used": False,
            },
            "outputs": {
                "cells": "cells.tsv.gz",
                "main_cluster_support": "main_cluster_support.tsv",
                "subcluster_support": "subcluster_support.tsv",
            },
        }
        with (temporary / "manifest.yaml").open("w") as handle:
            yaml.safe_dump(manifest, handle, sort_keys=False)
        temporary.rename(output_root)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(
        f"gse216371 labels retained={retained} main_clusters={len(main_support)} "
        "status=completed",
        flush=True,
    )
    return cells_path


def audit_ccres(config: Mapping[str, Any]) -> Path:
    family_root = project_path(str(config["processed_directory"]))
    output_root = family_root / "feature_axes" / CCRE_VERSION
    ccres_path = output_root / "candidate_ccres.tsv.gz"
    manifest_path = output_root / "manifest.yaml"
    if ccres_path.is_file() and manifest_path.is_file():
        print("gse216371 ccres status=reused", flush=True)
        return ccres_path
    if output_root.exists():
        raise FileExistsError(f"Partial immutable cCRE audit: {output_root}")

    source = workbook_source(config)
    record = locked_source_record(source)
    if not source.is_file() or source.stat().st_size != int(record["bytes"]):
        raise ValueError(f"Locked workbook byte count changed or is absent: {source}")
    if file_digest(source) != str(record["sha256"]):
        raise ValueError(f"Locked workbook SHA-256 changed: {source}")

    output_root.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{CCRE_VERSION}.", dir=output_root.parent))
    seen: set[str] = set()
    intervals: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    type_support: Counter[str] = Counter()
    ctcf_support: Counter[str] = Counter()
    rows_written = 0
    try:
        with zipfile.ZipFile(source) as workbook:
            shared_strings = load_shared_strings(workbook)
            rows = iter_sheet_rows(workbook, "xl/worksheets/sheet2.xml", shared_strings)
            header_row, header_values = next(rows)
            ordered_headers = tuple(
                value.strip() for _, value in sorted(header_values.items())
            )
            if header_row != 1 or ordered_headers != CCRE_COLUMNS:
                raise ValueError(
                    f"Unexpected GSE216371 cCRE header; row={header_row} "
                    f"columns={ordered_headers}"
                )
            headers = {value: index for index, value in header_values.items()}
            fieldnames = (
                "source_index",
                "peak_id",
                "author_ccre_name",
                "chrom",
                "start",
                "end",
                "ccre_type",
                "nearest_tss",
                "distance_to_tss_bp",
                "ctcf_binding",
                "module_94",
                "module_merged",
            )
            with gzip.open(temporary / "candidate_ccres.tsv.gz", "wt", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                for row_number, values in rows:
                    author_name = values.get(headers["cCREs_name"], "").strip()
                    chrom, start, end, peak_id = parse_ccre_name(author_name)
                    if end - start != 500:
                        raise ValueError(f"Author cCRE is not 500 bp at row {row_number}: {author_name}")
                    if peak_id in seen:
                        raise ValueError(f"Duplicate author cCRE at row {row_number}: {peak_id}")
                    seen.add(peak_id)
                    ccre_type = values.get(headers["cCREs_type"], "").strip()
                    nearest_tss = values.get(headers["nearest_TSS"], "").strip()
                    distance_text = values.get(headers["distance_To_TSS (bp)"], "").strip()
                    try:
                        distance = int(float(distance_text))
                    except ValueError as exc:
                        raise ValueError(
                            f"Invalid cCRE-to-TSS distance at row {row_number}: {distance_text!r}"
                        ) from exc
                    ctcf = values.get(headers["CTCF binding"], "").strip()
                    if not ccre_type or not nearest_tss or ctcf not in {"Yes", "No"}:
                        raise ValueError(f"Incomplete cCRE annotation at row {row_number}")
                    writer.writerow(
                        {
                            "source_index": row_number - 2,
                            "peak_id": peak_id,
                            "author_ccre_name": author_name,
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "ccre_type": ccre_type,
                            "nearest_tss": nearest_tss,
                            "distance_to_tss_bp": distance,
                            "ctcf_binding": ctcf,
                            "module_94": values.get(
                                headers["94_Modules_in_Fig2g"], ""
                            ).strip(),
                            "module_merged": values.get(
                                headers[
                                    "merged_Module_for_GREAT_and_Motif analysis_in_Fig2ij"
                                ],
                                "",
                            ).strip(),
                        }
                    )
                    intervals[chrom].append((start, end, peak_id))
                    type_support[ccre_type] += 1
                    ctcf_support[ctcf] += 1
                    rows_written += 1
                    if rows_written % 100_000 == 0:
                        print(f"gse216371 ccre_scan rows={rows_written}", flush=True)

        if rows_written != EXPECTED_CCRE_COUNT:
            raise ValueError(
                f"Expected {EXPECTED_CCRE_COUNT} author cCREs; found {rows_written}"
            )
        overlaps = []
        for chrom, chrom_intervals in intervals.items():
            chrom_intervals.sort(key=lambda value: (value[0], value[1], value[2]))
            previous_end = -1
            previous_name = ""
            for start, end, peak_id in chrom_intervals:
                if start < previous_end:
                    overlaps.append((chrom, previous_name, peak_id))
                    if len(overlaps) == 10:
                        break
                if end > previous_end:
                    previous_end = end
                    previous_name = peak_id
            if len(overlaps) == 10:
                break
        if overlaps:
            raise ValueError(f"Author cCREs overlap; first examples: {overlaps}")
        manifest = {
            "schema_version": 1,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "axis_version": CCRE_VERSION,
            "candidate_intervals": rows_written,
            "interval_width_bp": 500,
            "nonoverlapping": True,
            "contigs": {chrom: len(values) for chrom, values in sorted(intervals.items())},
            "ccre_type_support": dict(sorted(type_support.items())),
            "ctcf_support": dict(sorted(ctcf_support.items())),
            "source": str(source.relative_to(ROOT)),
            "source_bytes": source.stat().st_size,
            "source_sha256": str(record["sha256"]),
            "parameters": {
                "sheet": "cCREs in MOPA",
                "sheet_member": "xl/worksheets/sheet2.xml",
                "streaming_xml": True,
                "outcome_data_used": False,
            },
            "outputs": {"candidate_ccres": "candidate_ccres.tsv.gz"},
        }
        with (temporary / "manifest.yaml").open("w") as handle:
            yaml.safe_dump(manifest, handle, sort_keys=False)
        temporary.rename(output_root)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(f"gse216371 ccres candidates={rows_written} status=completed", flush=True)
    return ccres_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    return parser.parse_args()


def main() -> None:
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run this audit through run_shapemix_low_impact.sh")
    args = parse_args()
    config = load_yaml(args.config)
    audit_labels(config)
    audit_ccres(config)


if __name__ == "__main__":
    main()
