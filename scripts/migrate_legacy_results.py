#!/usr/bin/env python
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEGACY_DIRS = [
    ROOT / "cell2location_results",
    ROOT / "rctd_results",
    ROOT / "data" / "deconvolution_results",
]


def relpath(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Archive legacy top-level result folders.")
    parser.add_argument("--archive-root", default=str(ROOT / "data" / "archive" / "legacy_results"))
    parser.add_argument("--execute", action="store_true", help="Copy legacy folders. Default is dry-run only.")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing archived copies.")
    args = parser.parse_args()

    archive_root = Path(args.archive_root)
    planned: list[tuple[Path, Path, str]] = []
    for source in LEGACY_DIRS:
        if source.exists():
            destination = archive_root / source.name
            status = "planned"
            if destination.exists():
                status = "overwrite" if args.overwrite else "exists"
            planned.append((source, destination, status))

    archive_root.mkdir(parents=True, exist_ok=True)
    manifest_path = archive_root / "manifest.tsv"
    if not planned:
        print("No legacy result source folders found.")
        if manifest_path.exists():
            print(manifest_path)
            return

    with manifest_path.open("w") as manifest:
        manifest.write("source\tdestination\tstatus\n")
        for source, destination, status in planned:
            manifest.write(f"{relpath(source)}\t{relpath(destination)}\t{status}\n")

    for source, destination, status in planned:
        print(f"{source} -> {destination} [{status}]")
        if args.execute:
            if destination.exists():
                if not args.overwrite:
                    continue
                shutil.rmtree(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(source, destination, dirs_exist_ok=False)

    if not args.execute:
        print("Dry-run only. Re-run with --execute to copy these folders.")
    print(manifest_path)


if __name__ == "__main__":
    main()
