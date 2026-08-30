// Single-threaded, fail-closed fragment streamer for the GSE216371 reference.
//
// The Python materializer owns source locking, tar-member validation, ontology
// freezing, atomic publication, and AnnData construction. This helper handles
// the row-scale work without making the multi-billion-row source a Python loop.

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <unistd.h>
#include <zlib.h>

namespace {

using Clock = std::chrono::steady_clock;
constexpr std::uint32_t kMissingCell = std::numeric_limits<std::uint32_t>::max();
constexpr std::uint32_t kMissingFeature = std::numeric_limits<std::uint32_t>::max();
constexpr std::size_t kCoverageThreshold = 10;
constexpr std::uint64_t kProgressEvery = 10'000'000;
static_assert(
    std::endian::native == std::endian::little,
    "The frozen ShapeMix binary cache format requires a little-endian host");

struct TransparentHash {
  using is_transparent = void;
  std::size_t operator()(std::string_view value) const noexcept {
    return std::hash<std::string_view>{}(value);
  }
  std::size_t operator()(const std::string& value) const noexcept {
    return std::hash<std::string_view>{}(value);
  }
};

struct Cell {
  std::uint32_t index;
  std::uint32_t type;
  std::uint32_t well;
};

struct Interval {
  std::uint64_t start;
  std::uint64_t end;
  std::uint32_t feature;
};

struct ChromIndex {
  std::vector<std::uint64_t> starts;
  std::vector<std::uint64_t> ends;
  std::vector<std::uint32_t> features;
};

struct PeakIndex {
  std::unordered_map<std::string, ChromIndex, TransparentHash, std::equal_to<>> by_chrom;
  std::uint64_t features = 0;

  std::uint32_t assign(std::string_view chrom, std::uint64_t coordinate) const {
    const auto found = by_chrom.find(chrom);
    if (found == by_chrom.end()) {
      return kMissingFeature;
    }
    const auto& index = found->second;
    const auto upper = std::upper_bound(index.starts.begin(), index.starts.end(), coordinate);
    if (upper == index.starts.begin()) {
      return kMissingFeature;
    }
    const auto position = static_cast<std::size_t>(upper - index.starts.begin() - 1);
    if (coordinate < index.ends[position]) {
      return index.features[position];
    }
    return kMissingFeature;
  }
};

struct Counters {
  std::uint64_t total_rows = 0;
  std::uint64_t header_rows = 0;
  std::uint64_t valid_rows = 0;
  std::uint64_t unknown_barcodes = 0;
  std::uint64_t retained_fragments = 0;
  std::uint64_t fragments_with_assigned_cut_sites = 0;
  std::uint64_t cut_sites_outside_peaks = 0;
  std::uint64_t assigned_cut_sites = 0;
  std::uint64_t read_support_total = 0;
  std::array<std::uint64_t, 3> cut_sites_per_bin{0, 0, 0};
};

struct Arguments {
  std::string mode;
  std::string labels;
  std::string peaks;
  std::string chrom_sizes;
  std::string output_fragments;
  std::string output_statistics;
  std::string output_cell_totals;
  std::string output_events;
  std::string output_summary;
};

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

std::string required_value(int argc, char** argv, int& index) {
  if (index + 1 >= argc) {
    fail(std::string("Missing value for ") + argv[index]);
  }
  ++index;
  return argv[index];
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string key = argv[index];
    if (key == "--mode") {
      result.mode = required_value(argc, argv, index);
    } else if (key == "--labels") {
      result.labels = required_value(argc, argv, index);
    } else if (key == "--peaks") {
      result.peaks = required_value(argc, argv, index);
    } else if (key == "--chrom-sizes") {
      result.chrom_sizes = required_value(argc, argv, index);
    } else if (key == "--output-fragments") {
      result.output_fragments = required_value(argc, argv, index);
    } else if (key == "--output-statistics") {
      result.output_statistics = required_value(argc, argv, index);
    } else if (key == "--output-cell-totals") {
      result.output_cell_totals = required_value(argc, argv, index);
    } else if (key == "--output-events") {
      result.output_events = required_value(argc, argv, index);
    } else if (key == "--output-summary") {
      result.output_summary = required_value(argc, argv, index);
    } else {
      fail("Unknown argument: " + key);
    }
  }
  if (result.mode != "statistics" && result.mode != "shape") {
    fail("--mode must be statistics or shape");
  }
  if (result.labels.empty() || result.peaks.empty() || result.chrom_sizes.empty() ||
      result.output_summary.empty()) {
    fail("--labels, --peaks, --chrom-sizes, and --output-summary are required");
  }
  if (result.mode == "statistics" &&
      (result.output_fragments.empty() || result.output_statistics.empty() ||
       result.output_cell_totals.empty())) {
    fail("statistics mode requires fragment, statistics, and cell-total outputs");
  }
  if (result.mode == "shape" && result.output_events.empty()) {
    fail("shape mode requires --output-events");
  }
  return result;
}

bool gz_getline(gzFile handle, std::string& output) {
  output.clear();
  std::array<char, 1 << 16> buffer{};
  while (true) {
    char* value = gzgets(handle, buffer.data(), static_cast<int>(buffer.size()));
    if (value == nullptr) {
      if (gzeof(handle)) {
        return !output.empty();
      }
      int error_number = Z_OK;
      const char* message = gzerror(handle, &error_number);
      fail(std::string("gzip read failed: ") + (message == nullptr ? "unknown" : message));
    }
    output.append(value);
    if (!output.empty() && output.back() == '\n') {
      output.pop_back();
      if (!output.empty() && output.back() == '\r') {
        output.pop_back();
      }
      return true;
    }
    if (gzeof(handle)) {
      return true;
    }
  }
}

template <typename Function>
void for_each_gzip_line(const std::string& path, Function function) {
  gzFile handle = gzopen(path.c_str(), "rb");
  if (handle == nullptr) {
    fail("Cannot open gzip input: " + path);
  }
  gzbuffer(handle, 1 << 20);
  try {
    std::string line;
    while (gz_getline(handle, line)) {
      function(line);
    }
  } catch (...) {
    gzclose(handle);
    throw;
  }
  if (gzclose(handle) != Z_OK) {
    fail("Failed to close gzip input: " + path);
  }
}

std::vector<std::string_view> split_tabs(std::string_view line) {
  std::vector<std::string_view> fields;
  std::size_t start = 0;
  while (true) {
    const std::size_t tab = line.find('\t', start);
    if (tab == std::string_view::npos) {
      fields.push_back(line.substr(start));
      break;
    }
    fields.push_back(line.substr(start, tab - start));
    start = tab + 1;
  }
  return fields;
}

std::uint64_t parse_uint64(std::string_view value, const std::string& field) {
  if (value.empty()) {
    fail("Empty integer field: " + field);
  }
  std::uint64_t result = 0;
  const char* first = value.data();
  const char* last = value.data() + value.size();
  const auto parsed = std::from_chars(first, last, result);
  if (parsed.ec != std::errc{} || parsed.ptr != last) {
    fail("Invalid nonnegative integer in " + field + ": " + std::string(value));
  }
  return result;
}

std::unordered_map<std::string, std::uint64_t, TransparentHash, std::equal_to<>>
load_chrom_sizes(const std::string& path) {
  std::ifstream handle(path);
  if (!handle) {
    fail("Cannot open chromosome sizes: " + path);
  }
  std::unordered_map<std::string, std::uint64_t, TransparentHash, std::equal_to<>> result;
  std::string chrom;
  std::uint64_t size = 0;
  while (handle >> chrom >> size) {
    if (chrom.empty() || size == 0 || !result.emplace(chrom, size).second) {
      fail("Invalid or duplicate chromosome-size row for " + chrom);
    }
  }
  if (!handle.eof() || result.empty()) {
    fail("Malformed or empty chromosome-size file: " + path);
  }
  return result;
}

struct LabelData {
  std::unordered_map<std::string, Cell, TransparentHash, std::equal_to<>> cells;
  std::vector<std::string> ordered_ids;
  std::unordered_map<std::string, std::uint32_t, TransparentHash, std::equal_to<>> wells;
  std::uint32_t types = 0;
};

LabelData load_labels(const std::string& path) {
  LabelData result;
  bool header_seen = false;
  std::size_t cell_column = 0;
  std::size_t type_column = 0;
  std::size_t well_column = 0;
  for_each_gzip_line(path, [&](const std::string& line) {
    const auto fields = split_tabs(line);
    if (!header_seen) {
      for (std::size_t index = 0; index < fields.size(); ++index) {
        if (fields[index] == "cell_id") {
          cell_column = index;
        }
        if (fields[index] == "cell_type_index") {
          type_column = index;
        }
        if (fields[index] == "round4_barcode") {
          well_column = index;
        }
      }
      if (fields.empty() || fields[cell_column] != "cell_id" ||
          fields[type_column] != "cell_type_index" ||
          fields[well_column] != "round4_barcode") {
        fail("Resolved label table lacks cell_id, cell_type_index, or round4_barcode");
      }
      header_seen = true;
      return;
    }
    if (fields.size() <= std::max({cell_column, type_column, well_column})) {
      fail("Resolved label row has too few columns");
    }
    const std::string cell_id(fields[cell_column]);
    const std::string well(fields[well_column]);
    const auto type_value = parse_uint64(fields[type_column], "cell_type_index");
    if (cell_id.empty() || well.empty() ||
        type_value > std::numeric_limits<std::uint32_t>::max()) {
      fail("Invalid resolved label row");
    }
    if (result.ordered_ids.size() >= kMissingCell) {
      fail("Too many resolved cells for uint32 indices");
    }
    const auto well_inserted = result.wells.emplace(
        well, static_cast<std::uint32_t>(result.wells.size()));
    const Cell cell{static_cast<std::uint32_t>(result.ordered_ids.size()),
                    static_cast<std::uint32_t>(type_value),
                    well_inserted.first->second};
    if (!result.cells.emplace(cell_id, cell).second) {
      fail("Duplicate resolved cell ID: " + cell_id);
    }
    result.ordered_ids.push_back(cell_id);
    result.types = std::max(result.types, static_cast<std::uint32_t>(type_value + 1));
  });
  if (!header_seen || result.cells.empty() || result.types == 0) {
    fail("Resolved label table is empty");
  }
  std::vector<bool> observed(result.types, false);
  for (const auto& value : result.cells) {
    observed[value.second.type] = true;
  }
  if (std::find(observed.begin(), observed.end(), false) != observed.end()) {
    fail("Resolved cell-type indices are not contiguous from zero");
  }
  return result;
}

PeakIndex load_peaks(
    const std::string& path,
    const std::unordered_map<std::string, std::uint64_t, TransparentHash, std::equal_to<>>&
        chrom_sizes) {
  PeakIndex result;
  std::unordered_map<std::string, std::vector<Interval>, TransparentHash, std::equal_to<>>
      intervals;
  bool header_seen = false;
  std::size_t chrom_column = 0;
  std::size_t start_column = 0;
  std::size_t end_column = 0;
  for_each_gzip_line(path, [&](const std::string& line) {
    const auto fields = split_tabs(line);
    if (!header_seen) {
      for (std::size_t index = 0; index < fields.size(); ++index) {
        if (fields[index] == "chrom") chrom_column = index;
        if (fields[index] == "start") start_column = index;
        if (fields[index] == "end") end_column = index;
      }
      if (fields.empty() || fields[chrom_column] != "chrom" ||
          fields[start_column] != "start" || fields[end_column] != "end") {
        fail("Peak table lacks chrom, start, or end");
      }
      header_seen = true;
      return;
    }
    if (fields.size() <= std::max({chrom_column, start_column, end_column})) {
      fail("Peak row has too few columns");
    }
    if (result.features >= kMissingFeature) {
      fail("Too many peak features for uint32 indices");
    }
    const std::string chrom(fields[chrom_column]);
    const auto found_size = chrom_sizes.find(chrom);
    const auto start = parse_uint64(fields[start_column], "peak start");
    const auto end = parse_uint64(fields[end_column], "peak end");
    if (found_size == chrom_sizes.end() || start >= end || end > found_size->second) {
      fail("Peak is outside the declared genome: " + chrom);
    }
    intervals[chrom].push_back(
        Interval{start, end, static_cast<std::uint32_t>(result.features)});
    ++result.features;
  });
  if (!header_seen || result.features == 0) {
    fail("Peak table is empty");
  }
  for (auto& entry : intervals) {
    auto& values = entry.second;
    std::sort(values.begin(), values.end(), [](const Interval& left, const Interval& right) {
      if (left.start != right.start) return left.start < right.start;
      if (left.end != right.end) return left.end < right.end;
      return left.feature < right.feature;
    });
    ChromIndex index;
    index.starts.reserve(values.size());
    index.ends.reserve(values.size());
    index.features.reserve(values.size());
    std::uint64_t previous_end = 0;
    bool first = true;
    for (const auto& interval : values) {
      if (!first && interval.start < previous_end) {
        fail("Peak intervals overlap on " + entry.first);
      }
      first = false;
      previous_end = interval.end;
      index.starts.push_back(interval.start);
      index.ends.push_back(interval.end);
      index.features.push_back(interval.feature);
    }
    result.by_chrom.emplace(entry.first, std::move(index));
  }
  return result;
}

std::size_t fragment_bin(std::uint64_t length) {
  if (length < 100) return 0;
  if (length < 250) return 1;
  return 2;
}

void write_summary(const std::string& path, const std::string& mode,
                   const Counters& counters, std::uint64_t cells,
                   std::uint64_t types, std::uint64_t features,
                   std::uint64_t event_records) {
  std::ofstream handle(path);
  if (!handle) {
    fail("Cannot create summary: " + path);
  }
  handle << "key\tvalue\n"
         << "mode\t" << mode << "\n"
         << "cells\t" << cells << "\n"
         << "cell_types\t" << types << "\n"
         << "features\t" << features << "\n"
         << "total_rows\t" << counters.total_rows << "\n"
         << "header_rows\t" << counters.header_rows << "\n"
         << "valid_rows\t" << counters.valid_rows << "\n"
         << "unknown_barcodes\t" << counters.unknown_barcodes << "\n"
         << "retained_fragments\t" << counters.retained_fragments << "\n"
         << "fragments_with_assigned_cut_sites\t"
         << counters.fragments_with_assigned_cut_sites << "\n"
         << "cut_sites_outside_peaks\t" << counters.cut_sites_outside_peaks << "\n"
         << "assigned_cut_sites\t" << counters.assigned_cut_sites << "\n"
         << "read_support_total\t" << counters.read_support_total << "\n"
         << "cut_sites_per_bin.fragment_length_lt_100\t"
         << counters.cut_sites_per_bin[0] << "\n"
         << "cut_sites_per_bin.fragment_length_100_249\t"
         << counters.cut_sites_per_bin[1] << "\n"
         << "cut_sites_per_bin.fragment_length_ge_250\t"
         << counters.cut_sites_per_bin[2] << "\n"
         << "event_records\t" << event_records << "\n";
  if (!handle) {
    fail("Failed while writing summary: " + path);
  }
}

void write_statistics(const std::string& path, std::uint64_t types,
                      std::uint64_t features, const std::vector<std::uint64_t>& counts,
                      const std::vector<std::uint8_t>& coverage) {
  std::ofstream handle(path, std::ios::binary);
  if (!handle) {
    fail("Cannot create feature statistics: " + path);
  }
  const std::array<char, 8> magic{'S', 'M', '2', '1', '6', 'C', '0', '1'};
  handle.write(magic.data(), magic.size());
  handle.write(reinterpret_cast<const char*>(&types), sizeof(types));
  handle.write(reinterpret_cast<const char*>(&features), sizeof(features));
  handle.write(reinterpret_cast<const char*>(counts.data()),
               static_cast<std::streamsize>(counts.size() * sizeof(counts[0])));
  handle.write(reinterpret_cast<const char*>(coverage.data()),
               static_cast<std::streamsize>(coverage.size() * sizeof(coverage[0])));
  if (!handle) {
    fail("Failed while writing feature statistics: " + path);
  }
}

void write_cell_totals(const std::string& path, const LabelData& labels,
                       const std::vector<std::uint64_t>& rows,
                       const std::vector<std::uint64_t>& supports) {
  std::ofstream handle(path);
  if (!handle) {
    fail("Cannot create cell totals: " + path);
  }
  handle << "cell_id\tbed_rows\tread_support_sum\n";
  for (std::size_t index = 0; index < labels.ordered_ids.size(); ++index) {
    handle << labels.ordered_ids[index] << '\t' << rows[index] << '\t'
           << supports[index] << '\n';
  }
  if (!handle) {
    fail("Failed while writing cell totals: " + path);
  }
}

struct EventWriter {
  std::ofstream handle;
  std::uint64_t records = 0;

  explicit EventWriter(const std::string& path) : handle(path, std::ios::binary) {
    if (!handle) {
      fail("Cannot create shape-event output: " + path);
    }
  }

  void write(std::uint32_t cell, std::uint32_t feature, std::uint8_t layer) {
    handle.write(reinterpret_cast<const char*>(&cell), sizeof(cell));
    handle.write(reinterpret_cast<const char*>(&feature), sizeof(feature));
    handle.write(reinterpret_cast<const char*>(&layer), sizeof(layer));
    if (!handle) {
      fail("Failed while writing shape-event output");
    }
    ++records;
  }
};

void process(const Arguments& arguments) {
  const auto chrom_sizes = load_chrom_sizes(arguments.chrom_sizes);
  const auto labels = load_labels(arguments.labels);
  const auto peaks = load_peaks(arguments.peaks, chrom_sizes);

  const bool statistics_mode = arguments.mode == "statistics";
  std::vector<std::uint64_t> type_counts;
  std::vector<std::uint8_t> coverage;
  std::vector<std::uint32_t> coverage_cells;
  std::vector<std::uint64_t> cell_rows(labels.ordered_ids.size(), 0);
  std::vector<std::uint64_t> cell_supports(labels.ordered_ids.size(), 0);
  if (statistics_mode) {
    if (peaks.features > std::numeric_limits<std::size_t>::max() / labels.types) {
      fail("Feature-statistics matrix size overflows size_t");
    }
    type_counts.assign(static_cast<std::size_t>(peaks.features) * labels.types, 0);
    coverage.assign(static_cast<std::size_t>(peaks.features), 0);
    coverage_cells.assign(
        static_cast<std::size_t>(peaks.features) * kCoverageThreshold, kMissingCell);
  }

  gzFile normalized = nullptr;
  if (statistics_mode) {
    normalized = gzopen(arguments.output_fragments.c_str(), "wb6");
    if (normalized == nullptr) {
      fail("Cannot create normalized fragment output: " + arguments.output_fragments);
    }
    gzbuffer(normalized, 1 << 20);
  }
  EventWriter events(statistics_mode ? "/dev/null" : arguments.output_events);

  const int duplicated_stdin = dup(STDIN_FILENO);
  if (duplicated_stdin < 0) {
    fail("Cannot duplicate stdin");
  }
  gzFile input = gzdopen(duplicated_stdin, "rb");
  if (input == nullptr) {
    close(duplicated_stdin);
    fail("Cannot open concatenated gzip stream on stdin");
  }
  gzbuffer(input, 1 << 20);

  Counters counters;
  const auto started = Clock::now();
  std::uint32_t current_well = kMissingCell;
  bool member_marker_seen = false;
  try {
    std::string line;
    while (gz_getline(input, line)) {
      if (!line.empty() && line[0] == '#') {
        ++counters.header_rows;
        constexpr std::string_view prefix = "#shapemix_member\t";
        if (line.rfind(prefix, 0) == 0) {
          const std::string_view member =
              std::string_view(line).substr(prefix.size());
          constexpr std::string_view suffix = ".bed.gz";
          const auto underscore = member.rfind('_');
          if (underscore == std::string_view::npos ||
              member.size() <= underscore + 1 + suffix.size() ||
              member.substr(member.size() - suffix.size()) != suffix) {
            fail("Invalid ShapeMix member marker: " + std::string(member));
          }
          const std::string_view well = member.substr(
              underscore + 1,
              member.size() - underscore - 1 - suffix.size());
          const auto found_well = labels.wells.find(well);
          current_well =
              found_well == labels.wells.end() ? kMissingCell : found_well->second;
          member_marker_seen = true;
        }
        continue;
      }
      ++counters.total_rows;
      const auto fields = split_tabs(line);
      if (fields.size() != 5 ||
          std::any_of(fields.begin(), fields.end(), [](std::string_view value) {
            return value.empty();
          })) {
        fail("Fragment row does not have exactly five non-empty fields at row " +
             std::to_string(counters.total_rows));
      }
      const auto found_chrom = chrom_sizes.find(fields[0]);
      const auto start = parse_uint64(fields[1], "fragment start");
      const auto end = parse_uint64(fields[2], "fragment end");
      const auto support = parse_uint64(fields[4], "read support");
      if (found_chrom == chrom_sizes.end() || start >= end ||
          end > found_chrom->second || support == 0) {
        fail("Fragment coordinate/support gate failed at row " +
             std::to_string(counters.total_rows));
      }
      ++counters.valid_rows;

      const auto found_cell = labels.cells.find(fields[3]);
      if (found_cell == labels.cells.end()) {
        ++counters.unknown_barcodes;
      } else {
        const Cell cell = found_cell->second;
        if (statistics_mode &&
            (!member_marker_seen || current_well == kMissingCell ||
             cell.well != current_well)) {
          fail(
              "Retained barcode does not match the active tar-member Round4 well at row " +
              std::to_string(counters.total_rows));
        }
        ++counters.retained_fragments;
        ++cell_rows[cell.index];
        if (std::numeric_limits<std::uint64_t>::max() - cell_supports[cell.index] <
            support) {
          fail("Per-cell read-support total overflow");
        }
        cell_supports[cell.index] += support;
        if (std::numeric_limits<std::uint64_t>::max() - counters.read_support_total <
            support) {
          fail("Global read-support total overflow");
        }
        counters.read_support_total += support;
        if (statistics_mode) {
          if (gzwrite(normalized, line.data(), static_cast<unsigned int>(line.size())) !=
                  static_cast<int>(line.size()) ||
              gzputc(normalized, '\n') == -1) {
            fail("Failed while writing normalized fragments");
          }
        }

        const std::size_t layer = fragment_bin(end - start);
        std::uint64_t assigned_for_fragment = 0;
        for (const std::uint64_t coordinate : {start, end}) {
          const auto feature = peaks.assign(fields[0], coordinate);
          if (feature == kMissingFeature) {
            ++counters.cut_sites_outside_peaks;
            continue;
          }
          ++assigned_for_fragment;
          ++counters.assigned_cut_sites;
          ++counters.cut_sites_per_bin[layer];
          if (statistics_mode) {
            auto& count =
                type_counts[static_cast<std::size_t>(cell.type) * peaks.features + feature];
            if (count == std::numeric_limits<std::uint64_t>::max()) {
              fail("Cell-type feature count overflow");
            }
            ++count;
            auto& covered = coverage[feature];
            bool seen = false;
            for (std::size_t slot = 0; slot < covered; ++slot) {
              if (coverage_cells[static_cast<std::size_t>(feature) *
                                     kCoverageThreshold +
                                 slot] == cell.index) {
                seen = true;
                break;
              }
            }
            if (!seen && covered < kCoverageThreshold) {
              coverage_cells[static_cast<std::size_t>(feature) *
                                 kCoverageThreshold +
                             covered] = cell.index;
              ++covered;
            }
          } else {
            events.write(cell.index, feature, static_cast<std::uint8_t>(layer));
          }
        }
        if (assigned_for_fragment > 0) {
          ++counters.fragments_with_assigned_cut_sites;
        }
      }

      if (counters.total_rows % kProgressEvery == 0) {
        const double seconds =
            std::chrono::duration<double>(Clock::now() - started).count();
        std::cerr << "gse216371_stream mode=" << arguments.mode
                  << " rows=" << counters.total_rows
                  << " retained=" << counters.retained_fragments
                  << " million_rows_per_minute="
                  << (static_cast<double>(counters.total_rows) / 1.0e6) /
                         (seconds / 60.0)
                  << std::endl;
      }
    }
  } catch (...) {
    gzclose(input);
    if (normalized != nullptr) gzclose(normalized);
    throw;
  }
  if (gzclose(input) != Z_OK) {
    fail("Input gzip stream failed its complete close/CRC gate");
  }
  if (normalized != nullptr && gzclose(normalized) != Z_OK) {
    fail("Normalized fragment gzip output failed to close");
  }

  if (statistics_mode) {
    write_statistics(arguments.output_statistics, labels.types, peaks.features,
                     type_counts, coverage);
    write_cell_totals(arguments.output_cell_totals, labels, cell_rows, cell_supports);
  }
  write_summary(arguments.output_summary, arguments.mode, counters,
                labels.ordered_ids.size(), labels.types, peaks.features,
                events.records);

  const double seconds = std::chrono::duration<double>(Clock::now() - started).count();
  std::cerr << "gse216371_stream mode=" << arguments.mode
            << " rows=" << counters.total_rows
            << " retained=" << counters.retained_fragments
            << " assigned_cut_sites=" << counters.assigned_cut_sites
            << " elapsed_seconds=" << seconds << " status=complete" << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    process(parse_arguments(argc, argv));
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "gse216371_stream status=failed error=" << error.what() << std::endl;
    return 1;
  }
}
