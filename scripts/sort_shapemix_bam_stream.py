#!/usr/bin/env python3
"""Read SAM from stdin and write a resource-configurable name-sorted BAM."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=0)
    args = parser.parse_args()
    if args.threads < 0:
        parser.error("--threads must be nonnegative")
    return args


def main() -> None:
    if os.environ.get("DECONVATAC_RESOURCE_GUARD") != "1":
        raise RuntimeError("Run the streaming sorter under run_shapemix_low_impact.sh")
    args = parse_args()
    output = args.output.absolute()
    if output.exists():
        raise FileExistsError(f"Streaming BAM output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    pysam.sort("-n", "-@", str(args.threads), "-O", "BAM", "-o", str(output), "-")


if __name__ == "__main__":
    main()
