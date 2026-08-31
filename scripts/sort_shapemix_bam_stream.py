#!/usr/bin/env python3
"""Read SAM from stdin and write a one-thread name-sorted BAM with pinned pysam."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run the streaming sorter under run_shapemix_low_impact.sh")
    args = parse_args()
    output = args.output.absolute()
    if output.exists():
        raise FileExistsError(f"Streaming BAM output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    pysam.sort("-n", "-@", "0", "-O", "BAM", "-o", str(output), "-")


if __name__ == "__main__":
    main()
