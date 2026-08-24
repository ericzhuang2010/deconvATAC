"""Build sparse parent-fragment-length layers from 10x ATAC fragments.

The counter in this module deliberately treats each five-column fragments row
as one deduplicated fragment.  ``readSupport`` is retained only as a QC total;
it never weights the two Tn5 cut-site observations emitted by a row.

``pysam`` is an optional ShapeMix dependency and is therefore imported only by
the tabix entry point, not when :mod:`deconvatac.pp` is imported.
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, MutableMapping, Optional, Sequence, Tuple, Union

import numpy as np
from scipy import sparse


FRAGMENT_SHAPE_SCHEMA_VERSION = 1
FRAGMENT_SHAPE_AXIS = "parent_fragment_length_bp"
FRAGMENT_SHAPE_COUNT_UNIT = "deduplicated_cut_sites"
FRAGMENT_SHAPE_READ_SUPPORT_POLICY = "ignore"
FRAGMENT_SHAPE_PEAK_ASSIGNMENT = "containing_nonoverlapping_peak"


@dataclass(frozen=True)
class FragmentLengthBin:
    """One half-open parent-fragment-length interval and its output layer."""

    name: str
    min_inclusive: int
    max_exclusive: Optional[int]
    layer: str

    def contains(self, length: int) -> bool:
        """Return whether ``length`` belongs to this bin."""
        return length >= self.min_inclusive and (
            self.max_exclusive is None or length < self.max_exclusive
        )


DEFAULT_FRAGMENT_LENGTH_BINS: Tuple[FragmentLengthBin, ...] = (
    FragmentLengthBin(
        name="short",
        min_inclusive=0,
        max_exclusive=100,
        layer="fragment_length_lt_100",
    ),
    FragmentLengthBin(
        name="mono",
        min_inclusive=100,
        max_exclusive=250,
        layer="fragment_length_100_249",
    ),
    FragmentLengthBin(
        name="long",
        min_inclusive=250,
        max_exclusive=None,
        layer="fragment_length_ge_250",
    ),
)
FRAGMENT_SHAPE_LAYER_NAMES: Tuple[str, ...] = tuple(bin_.layer for bin_ in DEFAULT_FRAGMENT_LENGTH_BINS)


def parse_fragment_length_bins(
    bins: Optional[Sequence[Union[FragmentLengthBin, Mapping[str, Any]]]] = None,
) -> Tuple[FragmentLengthBin, ...]:
    """Parse and validate a complete, ordered fragment-length partition.

    Bins must start at zero, be contiguous and non-overlapping, and end in one
    unbounded interval.  This makes every valid positive fragment length map to
    exactly one layer.
    """
    raw_bins: Sequence[Union[FragmentLengthBin, Mapping[str, Any]]]
    raw_bins = DEFAULT_FRAGMENT_LENGTH_BINS if bins is None else bins
    parsed = []
    for raw in raw_bins:
        if isinstance(raw, FragmentLengthBin):
            bin_ = raw
        elif isinstance(raw, Mapping):
            missing = {"name", "min_inclusive", "layer"}.difference(raw)
            if missing:
                raise ValueError(f"Fragment-length bin is missing fields: {sorted(missing)}")
            max_exclusive = raw.get("max_exclusive")
            bin_ = FragmentLengthBin(
                name=str(raw["name"]),
                min_inclusive=int(raw["min_inclusive"]),
                max_exclusive=None if max_exclusive is None else int(max_exclusive),
                layer=str(raw["layer"]),
            )
        else:
            raise TypeError("Each fragment-length bin must be a FragmentLengthBin or mapping.")
        parsed.append(bin_)

    if not parsed:
        raise ValueError("At least one fragment-length bin is required.")
    if parsed[0].min_inclusive != 0:
        raise ValueError("Fragment-length bins must begin at zero.")
    if len({bin_.name for bin_ in parsed}) != len(parsed):
        raise ValueError("Fragment-length bin names must be unique.")
    if len({bin_.layer for bin_ in parsed}) != len(parsed):
        raise ValueError("Fragment-length layer names must be unique.")

    for index, bin_ in enumerate(parsed):
        if not bin_.name or not bin_.layer:
            raise ValueError("Fragment-length bin names and layer names cannot be empty.")
        if bin_.min_inclusive < 0:
            raise ValueError("Fragment-length bin boundaries cannot be negative.")
        if bin_.max_exclusive is not None and bin_.max_exclusive <= bin_.min_inclusive:
            raise ValueError("Each bounded fragment-length bin must have positive width.")
        if index < len(parsed) - 1:
            if bin_.max_exclusive is None:
                raise ValueError("Only the final fragment-length bin may be unbounded.")
            if parsed[index + 1].min_inclusive != bin_.max_exclusive:
                raise ValueError("Fragment-length bins must be contiguous and ordered.")
        elif bin_.max_exclusive is not None:
            raise ValueError("The final fragment-length bin must be unbounded.")
    return tuple(parsed)


def fragment_length_bin(
    length: int,
    bins: Optional[Sequence[Union[FragmentLengthBin, Mapping[str, Any]]]] = None,
) -> int:
    """Return the index of the unique bin containing ``length``."""
    if isinstance(length, bool) or not isinstance(length, (int, np.integer)):
        raise TypeError("Fragment length must be an integer.")
    if length < 0:
        raise ValueError("Fragment length cannot be negative.")
    parsed = parse_fragment_length_bins(bins)
    return _fragment_length_bin_from_parsed(int(length), parsed)


def _fragment_length_bin_from_parsed(
    length: int,
    bins: Tuple[FragmentLengthBin, ...],
) -> int:
    """Return a bin index after the ordered partition has been validated once."""
    for index, bin_ in enumerate(bins):
        if bin_.contains(length):
            return index
    raise RuntimeError("Validated bins did not cover the fragment length.")


@dataclass(frozen=True)
class FragmentRecord:
    """One parsed five-column Cell Ranger ARC fragments record."""

    chrom: str
    start: int
    end: int
    barcode: str
    read_support: int

    @property
    def length(self) -> int:
        """Return the parent-fragment length in base pairs."""
        return self.end - self.start


class FragmentParseError(ValueError):
    """A fragments row that cannot satisfy the five-column data contract."""

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = category


def parse_fragment_line(line: Union[str, bytes]) -> Optional[FragmentRecord]:
    """Parse one five-column fragments line.

    Comment/header lines return ``None``.  Schema and coordinate failures raise
    :class:`FragmentParseError` with category ``schema`` or ``coordinates`` so
    streaming callers can count and skip them without losing the full run.
    """
    if isinstance(line, bytes):
        try:
            line = line.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise FragmentParseError("Fragments row is not valid UTF-8.", "schema") from exc
    if not isinstance(line, str):
        raise TypeError("Fragments lines must be str or bytes.")
    stripped = line.rstrip("\r\n")
    if stripped.startswith("#"):
        return None
    fields = stripped.split("\t")
    if len(fields) != 5 or any(field == "" for field in fields):
        raise FragmentParseError("Expected exactly five non-empty tab-separated fields.", "schema")

    chrom, start_text, end_text, barcode, support_text = fields
    try:
        start = int(start_text)
        end = int(end_text)
    except ValueError as exc:
        raise FragmentParseError("Fragment coordinates must be integers.", "coordinates") from exc
    try:
        read_support = int(support_text)
    except ValueError as exc:
        raise FragmentParseError("readSupport must be an integer.", "schema") from exc

    if start < 0 or end <= start:
        raise FragmentParseError("Fragment coordinates must satisfy 0 <= start < end.", "coordinates")
    if read_support < 1:
        raise FragmentParseError("readSupport must be a positive integer.", "schema")
    return FragmentRecord(chrom=chrom, start=start, end=end, barcode=barcode, read_support=read_support)


@dataclass(frozen=True)
class PeakInterval:
    """A half-open genomic peak interval in output feature order."""

    chrom: str
    start: int
    end: int
    name: str


def parse_peak_name(name: str) -> PeakInterval:
    """Parse a canonical ``chrom:start-end`` peak identifier."""
    try:
        chrom, coordinates = str(name).rsplit(":", 1)
        start_text, end_text = coordinates.split("-", 1)
        start = int(start_text)
        end = int(end_text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Peak name is not in 'chrom:start-end' form: {name!r}") from exc
    return PeakInterval(chrom=chrom, start=start, end=end, name=str(name))


def read_peak_bed(path: Union[str, Path]) -> Tuple[PeakInterval, ...]:
    """Read and validate a BED file while preserving its data-row order."""
    peaks = []
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED row {line_number} has fewer than three columns.")
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(f"BED row {line_number} has non-integer coordinates.") from exc
            name = f"{fields[0]}:{start}-{end}"
            peaks.append(PeakInterval(fields[0], start, end, name))
    # Coordinate sorting is internal to PeakIndex, so a caller may also use a
    # ranked BED.  The returned feature order remains exactly the file order.
    return PeakIndex(peaks).peaks


# Plural alias reads naturally at call sites and preserves the public name used
# by early Step 2 design notes.
read_peaks_bed = read_peak_bed


class PeakIndex:
    """Coordinate-sorted lookup that retains the caller's feature order."""

    def __init__(self, peaks: Sequence[PeakInterval]) -> None:
        if not peaks:
            raise ValueError("At least one peak is required.")
        normalized = tuple(peaks)
        names = [peak.name for peak in normalized]
        if len(set(names)) != len(names):
            raise ValueError("Peak names must be unique.")

        by_contig: MutableMapping[str, list[Tuple[int, int, int]]] = {}
        for feature_index, peak in enumerate(normalized):
            if not peak.chrom or not peak.name:
                raise ValueError("Peak chromosome and name cannot be empty.")
            if peak.start < 0 or peak.end <= peak.start:
                raise ValueError(f"Invalid half-open peak coordinates for {peak.name!r}.")
            by_contig.setdefault(peak.chrom, []).append((peak.start, peak.end, feature_index))

        self._lookup: dict[str, Tuple[Tuple[int, ...], Tuple[int, ...], Tuple[int, ...]]] = {}
        for chrom, intervals in by_contig.items():
            intervals.sort(key=lambda value: (value[0], value[1], value[2]))
            previous_end = -1
            for start, end, _ in intervals:
                if start < previous_end:
                    raise ValueError(f"Peaks on {chrom!r} overlap; unique peak assignment is impossible.")
                previous_end = end
            self._lookup[chrom] = (
                tuple(value[0] for value in intervals),
                tuple(value[1] for value in intervals),
                tuple(value[2] for value in intervals),
            )
        self.peaks = normalized
        self.contigs = frozenset(self._lookup)

    def assign(self, chrom: str, coordinate: int) -> Optional[int]:
        """Return the output feature index containing a cut-site coordinate."""
        lookup = self._lookup.get(chrom)
        if lookup is None:
            return None
        starts, ends, feature_indices = lookup
        interval_index = bisect_right(starts, coordinate) - 1
        if interval_index >= 0 and coordinate < ends[interval_index]:
            return feature_indices[interval_index]
        return None


def build_peak_index(
    peaks: Sequence[Union[PeakInterval, str, Sequence[Any], Mapping[str, Any]]],
) -> PeakIndex:
    """Normalize common peak representations into a validated lookup index."""
    normalized = []
    for raw in peaks:
        if isinstance(raw, PeakInterval):
            peak = raw
        elif isinstance(raw, str):
            peak = parse_peak_name(raw)
        elif isinstance(raw, Mapping):
            missing = {"chrom", "start", "end"}.difference(raw)
            if missing:
                raise ValueError(f"Peak mapping is missing fields: {sorted(missing)}")
            name = raw.get("name", raw.get("peak_id"))
            if name is None:
                name = f"{raw['chrom']}:{raw['start']}-{raw['end']}"
            peak = PeakInterval(str(raw["chrom"]), int(raw["start"]), int(raw["end"]), str(name))
        else:
            values = tuple(raw)
            if len(values) not in (3, 4):
                raise ValueError("Peak tuples must contain chrom, start, end, and optionally name.")
            chrom, start, end = values[:3]
            name = values[3] if len(values) == 4 else f"{chrom}:{start}-{end}"
            peak = PeakInterval(str(chrom), int(start), int(end), str(name))
        normalized.append(peak)
    return PeakIndex(normalized)


@dataclass
class FragmentShapeQC:
    """Streaming preprocessing counters for one fragments input scope."""

    total_rows: int = 0
    header_rows: int = 0
    invalid_schema_rows: int = 0
    invalid_coordinate_rows: int = 0
    unknown_barcodes: int = 0
    filtered_contigs: int = 0
    valid_rows: int = 0
    retained_fragments: int = 0
    fragments_with_assigned_cut_sites: int = 0
    cut_sites_outside_peaks: int = 0
    assigned_cut_sites: int = 0
    read_support_total: int = 0
    cut_sites_per_bin: dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return H5AD/YAML-safe Python scalars."""
        counters = {
            "total_rows": int(self.total_rows),
            "header_rows": int(self.header_rows),
            "invalid_schema_rows": int(self.invalid_schema_rows),
            "invalid_coordinate_rows": int(self.invalid_coordinate_rows),
            "unknown_barcodes": int(self.unknown_barcodes),
            "filtered_contigs": int(self.filtered_contigs),
            "valid_rows": int(self.valid_rows),
            "retained_fragments": int(self.retained_fragments),
            "fragments_with_assigned_cut_sites": int(self.fragments_with_assigned_cut_sites),
            "cut_sites_outside_peaks": int(self.cut_sites_outside_peaks),
            "assigned_cut_sites": int(self.assigned_cut_sites),
            "read_support_total": int(self.read_support_total),
        }
        # Keep the provenance mapping scalar-valued for the typed data contract
        # and AnnData 0.9.  The in-memory dataclass retains the convenient map.
        for layer, count in self.cut_sites_per_bin.items():
            counters[f"cut_sites_per_bin.{layer}"] = int(count)
        return counters


class _ChunkedSparseAccumulator:
    """Accumulate bounded COO buffers into log-structured CSR layers.

    A flush produces one level-zero CSR run for every represented layer.  Like
    binary addition, an occupied level is merged with the incoming run and
    carried upward until an empty level is found.  Each event therefore crosses
    at most logarithmically many CSR additions instead of being rescanned by
    every later flush.
    """

    def __init__(self, shape: Tuple[int, int], n_layers: int, chunk_size: int) -> None:
        if isinstance(chunk_size, bool) or not isinstance(chunk_size, (int, np.integer)) or chunk_size <= 0:
            raise ValueError("chunk_size must be a positive integer.")
        self.shape = shape
        self.n_layers = n_layers
        self.chunk_size = int(chunk_size)
        self._rows: list[int] = []
        self._columns: list[int] = []
        self._layers: list[int] = []
        self._levels: list[list[Optional[sparse.csr_matrix]]] = [[] for _ in range(n_layers)]
        self._merge_count = 0
        self._finished: Optional[Tuple[sparse.csr_matrix, ...]] = None

    @staticmethod
    def _merge(left: sparse.csr_matrix, right: sparse.csr_matrix) -> sparse.csr_matrix:
        """Add two CSR runs and canonicalize their sparse representation."""
        merged = (left + right).tocsr()
        merged.sum_duplicates()
        merged.eliminate_zeros()
        merged.sort_indices()
        return merged

    def _push_run(self, layer: int, run: sparse.csr_matrix) -> None:
        """Insert one CSR run using binary carries within a single layer."""
        levels = self._levels[layer]
        level = 0
        while True:
            if level == len(levels):
                levels.append(run)
                return
            existing = levels[level]
            if existing is None:
                levels[level] = run
                return
            levels[level] = None
            run = self._merge(existing, run)
            self._merge_count += 1
            level += 1

    @property
    def occupied_levels(self) -> Tuple[Tuple[int, ...], ...]:
        """Return occupied binary levels per layer for diagnostics/tests."""
        return tuple(
            tuple(index for index, run in enumerate(levels) if run is not None)
            for levels in self._levels
        )

    @property
    def merge_count(self) -> int:
        """Return the number of CSR additions performed so far."""
        return self._merge_count

    def add(self, row: int, column: int, layer: int) -> None:
        if self._finished is not None:
            raise RuntimeError("Cannot add events after the accumulator is finished.")
        self._rows.append(row)
        self._columns.append(column)
        self._layers.append(layer)
        if len(self._rows) >= self.chunk_size:
            self.flush()

    def flush(self) -> None:
        if self._finished is not None:
            raise RuntimeError("Cannot flush an accumulator after it is finished.")
        if not self._rows:
            return
        rows = np.asarray(self._rows, dtype=np.int64)
        columns = np.asarray(self._columns, dtype=np.int64)
        layer_indices = np.asarray(self._layers, dtype=np.int64)
        for layer in range(self.n_layers):
            selected = layer_indices == layer
            if not np.any(selected):
                continue
            chunk = sparse.coo_matrix(
                (np.ones(int(selected.sum()), dtype=np.int64), (rows[selected], columns[selected])),
                shape=self.shape,
                dtype=np.int64,
            ).tocsr()
            chunk.sum_duplicates()
            chunk.sort_indices()
            self._push_run(layer, chunk)
        self._rows.clear()
        self._columns.clear()
        self._layers.clear()

    def finish(self) -> Tuple[sparse.csr_matrix, ...]:
        if self._finished is not None:
            return self._finished
        self.flush()
        matrices = []
        for levels in self._levels:
            matrix: Optional[sparse.csr_matrix] = None
            # Low-to-high merging keeps the partial result no larger than the
            # next binary run, avoiding repeated scans of the largest run.
            for run in levels:
                if run is None:
                    continue
                if matrix is None:
                    matrix = run
                else:
                    matrix = self._merge(matrix, run)
                    self._merge_count += 1
            if matrix is None:
                matrix = sparse.csr_matrix(self.shape, dtype=np.int64)
            matrix.sum_duplicates()
            matrix.eliminate_zeros()
            matrix.sort_indices()
            matrices.append(matrix)
        self._finished = tuple(matrices)
        self._levels = [[] for _ in range(self.n_layers)]
        return self._finished


@dataclass
class FragmentShapeResult:
    """Sparse fragment-shape layers and the ordering/provenance needed to use them."""

    barcodes: Tuple[str, ...]
    peaks: Tuple[PeakInterval, ...]
    bins: Tuple[FragmentLengthBin, ...]
    layers: dict[str, sparse.csr_matrix]
    qc: FragmentShapeQC
    right_cut_offset: int

    @property
    def X(self) -> sparse.csr_matrix:
        """Return the exact CSR sum of all ordered shape layers."""
        total = sparse.csr_matrix((len(self.barcodes), len(self.peaks)), dtype=np.int64)
        for bin_ in self.bins:
            total = (total + self.layers[bin_.layer]).tocsr()
        total.sum_duplicates()
        total.eliminate_zeros()
        total.sort_indices()
        return total

    def fragment_shape_metadata(self, provenance: Optional[Mapping[str, Any]] = None) -> dict[str, Any]:
        """Return H5AD-safe fragment-shape metadata with a resolved cut convention.

        Bins use ordered numeric mapping keys instead of a list of dictionaries
        because AnnData 0.9 cannot round-trip the latter.  Every value carries
        its stable bin name, and the absent ``max_exclusive`` field on the final
        bin means unbounded.  Source-run preprocessing QC remains immutable,
        while ``matrix_counters`` records totals for the object as stored and
        can therefore be recomputed after safe matrix transformations.
        """
        # Import locally to keep preprocessing module initialization independent
        # of the higher-level data package and its AnnData loaders.
        from deconvatac.data import ordered_feature_sha256

        layer_totals = {
            bin_.layer: int(self.layers[bin_.layer].sum()) for bin_ in self.bins
        }
        matrix_counters = {
            "assigned_cut_sites": sum(layer_totals.values()),
            **{
                f"cut_sites_per_bin.{layer_name}": count
                for layer_name, count in layer_totals.items()
            },
        }
        metadata: dict[str, Any] = {
            "schema_version": FRAGMENT_SHAPE_SCHEMA_VERSION,
            "axis": FRAGMENT_SHAPE_AXIS,
            "count_unit": FRAGMENT_SHAPE_COUNT_UNIT,
            "read_support_policy": FRAGMENT_SHAPE_READ_SUPPORT_POLICY,
            "peak_assignment": FRAGMENT_SHAPE_PEAK_ASSIGNMENT,
            "left_cut_offset": 0,
            "right_cut_offset": int(self.right_cut_offset),
            "feature_sha256": ordered_feature_sha256(peak.name for peak in self.peaks),
            "bins": {},
            "preprocessing_counters": self.qc.to_dict(),
            "matrix_counters": matrix_counters,
        }
        for order, bin_ in enumerate(self.bins):
            bin_metadata: dict[str, Any] = {
                "name": bin_.name,
                "order": order,
                "min_inclusive": int(bin_.min_inclusive),
                "layer": bin_.layer,
            }
            if bin_.max_exclusive is not None:
                bin_metadata["max_exclusive"] = int(bin_.max_exclusive)
            metadata["bins"][str(order)] = bin_metadata

        if provenance:
            reserved = set(metadata).intersection(provenance)
            if reserved:
                raise ValueError(f"Provenance cannot override fragment-shape fields: {sorted(reserved)}")
            metadata.update(dict(provenance))
        return metadata


def _validate_right_cut_offset(right_cut_offset: int) -> int:
    if isinstance(right_cut_offset, bool) or not isinstance(right_cut_offset, (int, np.integer)):
        raise TypeError("right_cut_offset must be the integer 0 or -1.")
    offset = int(right_cut_offset)
    if offset not in (-1, 0):
        raise ValueError("right_cut_offset must be 0 or -1.")
    return offset


def _count_records(
    records: Iterable[Union[str, bytes, FragmentRecord]],
    barcodes: Sequence[str],
    peak_index: PeakIndex,
    bins: Tuple[FragmentLengthBin, ...],
    right_cut_offset: int,
    chunk_size: int,
    initial_header_rows: int = 0,
) -> FragmentShapeResult:
    ordered_barcodes = tuple(str(barcode) for barcode in barcodes)
    if not ordered_barcodes:
        raise ValueError("At least one barcode is required.")
    if len(set(ordered_barcodes)) != len(ordered_barcodes):
        raise ValueError("Barcodes must be unique.")
    barcode_indices = {barcode: index for index, barcode in enumerate(ordered_barcodes)}

    qc = FragmentShapeQC(
        header_rows=int(initial_header_rows),
        cut_sites_per_bin={bin_.layer: 0 for bin_ in bins},
    )
    accumulator = _ChunkedSparseAccumulator(
        shape=(len(ordered_barcodes), len(peak_index.peaks)),
        n_layers=len(bins),
        chunk_size=chunk_size,
    )

    for raw_record in records:
        if isinstance(raw_record, FragmentRecord):
            record = raw_record
            qc.total_rows += 1
            if any(
                isinstance(value, (bool, np.bool_))
                or not isinstance(value, (int, np.integer))
                for value in (record.start, record.end)
            ):
                qc.invalid_coordinate_rows += 1
                continue
            if record.start < 0 or record.end <= record.start:
                qc.invalid_coordinate_rows += 1
                continue
            if (
                not isinstance(record.chrom, str)
                or not record.chrom
                or not isinstance(record.barcode, str)
                or not record.barcode
                or isinstance(record.read_support, (bool, np.bool_))
                or not isinstance(record.read_support, (int, np.integer))
                or record.read_support < 1
            ):
                qc.invalid_schema_rows += 1
                continue
        else:
            if isinstance(raw_record, bytes):
                is_header = raw_record.startswith(b"#")
            elif isinstance(raw_record, str):
                is_header = raw_record.startswith("#")
            else:
                raise TypeError("Records must be FragmentRecord, str, or bytes.")
            if is_header:
                qc.header_rows += 1
                continue
            qc.total_rows += 1
            try:
                record = parse_fragment_line(raw_record)
            except FragmentParseError as exc:
                if exc.category == "coordinates":
                    qc.invalid_coordinate_rows += 1
                else:
                    qc.invalid_schema_rows += 1
                continue
            if record is None:  # Defensive for unusual strings with a decoded header.
                qc.header_rows += 1
                qc.total_rows -= 1
                continue

        qc.valid_rows += 1
        qc.read_support_total += record.read_support
        unknown_barcode = record.barcode not in barcode_indices
        filtered_contig = record.chrom not in peak_index.contigs
        if unknown_barcode:
            qc.unknown_barcodes += 1
        if filtered_contig:
            qc.filtered_contigs += 1
        if unknown_barcode or filtered_contig:
            continue

        qc.retained_fragments += 1
        bin_index = _fragment_length_bin_from_parsed(record.length, bins)
        row_index = barcode_indices[record.barcode]
        assigned_for_fragment = 0
        for coordinate in (record.start, record.end + right_cut_offset):
            feature_index = peak_index.assign(record.chrom, coordinate)
            if feature_index is None:
                qc.cut_sites_outside_peaks += 1
                continue
            accumulator.add(row_index, feature_index, bin_index)
            qc.assigned_cut_sites += 1
            qc.cut_sites_per_bin[bins[bin_index].layer] += 1
            assigned_for_fragment += 1
        if assigned_for_fragment:
            qc.fragments_with_assigned_cut_sites += 1

    matrices = accumulator.finish()
    layers = {bin_.layer: matrix for bin_, matrix in zip(bins, matrices)}
    return FragmentShapeResult(
        barcodes=ordered_barcodes,
        peaks=peak_index.peaks,
        bins=bins,
        layers=layers,
        qc=qc,
        right_cut_offset=right_cut_offset,
    )


def count_fragment_shapes_from_records(
    records: Iterable[Union[str, bytes, FragmentRecord]],
    barcodes: Sequence[str],
    peaks: Union[PeakIndex, Sequence[Union[PeakInterval, str, Sequence[Any], Mapping[str, Any]]]],
    *,
    right_cut_offset: int,
    bins: Optional[Sequence[Union[FragmentLengthBin, Mapping[str, Any]]]] = None,
    chunk_size: int = 1_000_000,
) -> FragmentShapeResult:
    """Count an iterable of fragments into bounded sparse length-bin layers."""
    parsed_bins = parse_fragment_length_bins(bins)
    peak_index = peaks if isinstance(peaks, PeakIndex) else build_peak_index(peaks)
    offset = _validate_right_cut_offset(right_cut_offset)
    return _count_records(
        records=records,
        barcodes=barcodes,
        peak_index=peak_index,
        bins=parsed_bins,
        right_cut_offset=offset,
        chunk_size=chunk_size,
    )


def _tabix_lines(tabix_file: Any, contigs: Optional[Sequence[str]]) -> Iterator[str]:
    if contigs is None:
        yield from tabix_file.fetch()
        return
    seen = set()
    for contig in contigs:
        if contig in seen:
            raise ValueError(f"Tabix contigs must be unique; repeated {contig!r}.")
        seen.add(contig)
        yield from tabix_file.fetch(str(contig))


def count_fragment_shapes(
    fragments_path: Union[str, Path],
    barcodes: Sequence[str],
    peaks: Union[PeakIndex, Sequence[Union[PeakInterval, str, Sequence[Any], Mapping[str, Any]]]],
    *,
    right_cut_offset: int,
    bins: Optional[Sequence[Union[FragmentLengthBin, Mapping[str, Any]]]] = None,
    chunk_size: int = 1_000_000,
    contigs: Optional[Sequence[str]] = None,
) -> FragmentShapeResult:
    """Stream a BGZF/tabix fragments file into sparse shape layers.

    By default every indexed row is fetched, which permits complete unknown-
    barcode and filtered-contig QC.  ``contigs`` can restrict the input scope
    for diagnostics or shards; counters then describe only the fetched scope.
    """
    try:
        import pysam  # type: ignore[import-not-found]
    except ImportError as exc:
        raise ImportError(
            "count_fragment_shapes requires optional dependency 'pysam'; "
            "install deconvATAC[shapemix]."
        ) from exc

    parsed_bins = parse_fragment_length_bins(bins)
    peak_index = peaks if isinstance(peaks, PeakIndex) else build_peak_index(peaks)
    offset = _validate_right_cut_offset(right_cut_offset)
    with pysam.TabixFile(str(fragments_path)) as tabix_file:
        header_rows = len(tuple(tabix_file.header))
        return _count_records(
            records=_tabix_lines(tabix_file, contigs),
            barcodes=barcodes,
            peak_index=peak_index,
            bins=parsed_bins,
            right_cut_offset=offset,
            chunk_size=chunk_size,
            initial_header_rows=header_rows,
        )


def build_fragment_shape_anndata(
    result: FragmentShapeResult,
    *,
    obs: Optional[Any] = None,
    var: Optional[Any] = None,
    provenance: Optional[Mapping[str, Any]] = None,
) -> Any:
    """Create an AnnData object whose CSR ``.X`` is the exact layer sum."""
    import anndata as ad
    import pandas as pd

    if obs is None:
        obs = pd.DataFrame(index=pd.Index(result.barcodes, name="barcode"))
    else:
        obs = obs.copy()
        if tuple(str(value) for value in obs.index) != result.barcodes:
            raise ValueError("obs index must exactly match result barcode order.")
    peak_names = tuple(peak.name for peak in result.peaks)
    if var is None:
        var = pd.DataFrame(
            {
                "chrom": [peak.chrom for peak in result.peaks],
                "start": [peak.start for peak in result.peaks],
                "end": [peak.end for peak in result.peaks],
            },
            index=pd.Index(peak_names, name="peak"),
        )
    else:
        var = var.copy()
        if tuple(str(value) for value in var.index) != peak_names:
            raise ValueError("var index must exactly match result peak order.")

    adata = ad.AnnData(X=result.X.copy(), obs=obs, var=var)
    for bin_ in result.bins:
        layer = result.layers[bin_.layer]
        if not sparse.isspmatrix_csr(layer):
            raise TypeError(f"Layer {bin_.layer!r} must be CSR.")
        adata.layers[bin_.layer] = layer.copy()
    adata.uns["fragment_shape"] = result.fragment_shape_metadata(provenance=provenance)
    return adata
