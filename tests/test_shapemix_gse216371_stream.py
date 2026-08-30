from __future__ import annotations

import gzip
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/shapemix_gse216371_stream.cpp"
EVENT_DTYPE = np.dtype(
    [("cell", "<u4"), ("feature", "<u4"), ("layer", "u1")], align=False
)


@pytest.fixture(scope="module")
def streamer(tmp_path_factory: pytest.TempPathFactory) -> Path:
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("g++ is unavailable")
    output = tmp_path_factory.mktemp("gse216371_stream") / "streamer"
    subprocess.run(
        [
            compiler,
            "-O2",
            "-std=c++20",
            "-Wall",
            "-Wextra",
            "-Werror",
            str(SOURCE),
            "-lz",
            "-o",
            str(output),
        ],
        check=True,
        cwd=ROOT,
    )
    return output


def write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def source_payload() -> bytes:
    first = gzip.compress(
        (
            "#shapemix_member\tGSM1_E11A.bed.gz\n"
            "chr1\t0\t50\tcellA\t2\n"
            "chr1\t400\t500\tunknown\t3\n"
        ).encode()
    )
    second = gzip.compress(
        (
            "#shapemix_member\tGSM2_E11B.bed.gz\n"
            "chr1\t100\t300\tcellB\t1\n"
        ).encode()
    )
    return first + second


def common_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    labels = tmp_path / "labels.tsv.gz"
    peaks = tmp_path / "peaks.tsv.gz"
    chrom_sizes = tmp_path / "chrom.sizes"
    write_gzip(
        labels,
        "cell_id\tcell_type_index\tround4_barcode\n"
        "cellA\t0\tE11A\n"
        "cellB\t1\tE11B\n",
    )
    write_gzip(
        peaks,
        "source_index\tpeak_id\tchrom\tstart\tend\n"
        "0\tchr1:0-100\tchr1\t0\t100\n"
        "1\tchr1:100-200\tchr1\t100\t200\n",
    )
    chrom_sizes.write_text("chr1\t1000\n")
    return labels, peaks, chrom_sizes


def parse_summary(path: Path) -> dict[str, str]:
    rows = [line.split("\t", 1) for line in path.read_text().splitlines()]
    assert rows[0] == ["key", "value"]
    return dict(rows[1:])


def test_statistics_mode_filters_and_counts_exactly(
    streamer: Path, tmp_path: Path
) -> None:
    labels, peaks, chrom_sizes = common_inputs(tmp_path)
    normalized = tmp_path / "retained.tsv.gz"
    statistics = tmp_path / "statistics.bin"
    cell_totals = tmp_path / "cell_totals.tsv"
    summary = tmp_path / "summary.tsv"
    subprocess.run(
        [
            str(streamer),
            "--mode",
            "statistics",
            "--labels",
            str(labels),
            "--peaks",
            str(peaks),
            "--chrom-sizes",
            str(chrom_sizes),
            "--output-fragments",
            str(normalized),
            "--output-statistics",
            str(statistics),
            "--output-cell-totals",
            str(cell_totals),
            "--output-summary",
            str(summary),
        ],
        input=source_payload(),
        check=True,
        cwd=ROOT,
    )

    with statistics.open("rb") as handle:
        assert handle.read(8) == b"SM216C01"
        types, features = struct.unpack("<QQ", handle.read(16))
        counts = np.frombuffer(handle.read(types * features * 8), dtype="<u8")
        coverage = np.frombuffer(handle.read(features), dtype="u1")
        assert handle.read() == b""
    np.testing.assert_array_equal(counts.reshape(types, features), [[2, 0], [0, 1]])
    np.testing.assert_array_equal(coverage, [1, 1])
    assert gzip.open(normalized, "rt").read() == (
        "chr1\t0\t50\tcellA\t2\n"
        "chr1\t100\t300\tcellB\t1\n"
    )
    assert cell_totals.read_text().splitlines() == [
        "cell_id\tbed_rows\tread_support_sum",
        "cellA\t1\t2",
        "cellB\t1\t1",
    ]
    values = parse_summary(summary)
    assert values["total_rows"] == "3"
    assert values["valid_rows"] == "3"
    assert values["unknown_barcodes"] == "1"
    assert values["retained_fragments"] == "2"
    assert values["assigned_cut_sites"] == "3"


def test_statistics_mode_rejects_retained_barcode_in_wrong_well(
    streamer: Path, tmp_path: Path
) -> None:
    labels, peaks, chrom_sizes = common_inputs(tmp_path)
    command = [
        str(streamer),
        "--mode",
        "statistics",
        "--labels",
        str(labels),
        "--peaks",
        str(peaks),
        "--chrom-sizes",
        str(chrom_sizes),
        "--output-fragments",
        str(tmp_path / "retained.tsv.gz"),
        "--output-statistics",
        str(tmp_path / "statistics.bin"),
        "--output-cell-totals",
        str(tmp_path / "cell_totals.tsv"),
        "--output-summary",
        str(tmp_path / "summary.tsv"),
    ]
    observed = subprocess.run(
        command,
        input=gzip.compress(
            (
                "#shapemix_member\tGSM2_E11B.bed.gz\n"
                "chr1\t0\t50\tcellA\t2\n"
            ).encode()
        ),
        cwd=ROOT,
        capture_output=True,
    )
    assert observed.returncode != 0
    assert b"does not match the active tar-member Round4 well" in observed.stderr


def test_shape_mode_emits_packed_layer_events(
    streamer: Path, tmp_path: Path
) -> None:
    labels, peaks, chrom_sizes = common_inputs(tmp_path)
    events = tmp_path / "events.bin"
    summary = tmp_path / "summary.tsv"
    subprocess.run(
        [
            str(streamer),
            "--mode",
            "shape",
            "--labels",
            str(labels),
            "--peaks",
            str(peaks),
            "--chrom-sizes",
            str(chrom_sizes),
            "--output-events",
            str(events),
            "--output-summary",
            str(summary),
        ],
        input=gzip.compress(
            (
                "chr1\t0\t50\tcellA\t2\n"
                "chr1\t100\t300\tcellB\t1\n"
            ).encode()
        ),
        check=True,
        cwd=ROOT,
    )

    observed = np.fromfile(events, dtype=EVENT_DTYPE)
    assert observed.tolist() == [(0, 0, 0), (0, 0, 0), (1, 1, 1)]
    values = parse_summary(summary)
    assert values["event_records"] == "3"
    assert values["cut_sites_per_bin.fragment_length_lt_100"] == "2"
    assert values["cut_sites_per_bin.fragment_length_100_249"] == "1"
