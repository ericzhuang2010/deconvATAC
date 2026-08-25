# GSE205055 and GSE263333 acquisition and preprocessing record

Status: complete on 2026-08-24 through source-ready preprocessing; labeled-reference and ShapeMix feature-axis gates remain intentionally pending

## Scope and source locks

The acquisition uses each GEO family's complete author-processed parent archive, not a selected-file subset. Raw SRA/FASTQ reads were not downloaded because the deposits already contain the fragment-level intervals, spatial assets, and orthogonal matrices required for this stage. Reprocessing raw reads is a separate protocol-reconstruction experiment and is not needed to preserve fragment lengths here.

| Family | Complete scope | Parent archive bytes | Parent archive SHA-256 | Tracked lock |
|---|---:|---:|---|---|
| [GSE205055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205055) | 22 samples; 38 supplementary files; parent plus seven SubSeries metadata files | 8,177,080,320 | `30c4fc6354cefab86f83fb1dee7b545f77233321f798f462bf636fdee08d2f0d` | `configs/data_sources/shapemix_gse205055_lock.yaml` |
| [GSE263333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263333) | 12 samples; 32 supplementary files; complete Series metadata | 3,306,485,760 | `a303eb357e3b3cebce83fe65efe0e366981511d975cbb7c245eabb2a290cd7a4` | `configs/data_sources/shapemix_gse263333_lock.yaml` |

GSE205055 includes GSE205051, GSE205052, GSE205054, GSE205180, GSE205181, GSE217091, and GSE218593. The six ATAC-bearing groups cover P21/P22 mouse brain, E13 mouse embryo, and adult human hippocampus. GSE263333 contains two ATAC-bearing groups: E13 mouse embryo and five-month-old EAE mouse brain. All related RNA, protein, histone/CUT&Tag, and spatial files in the deposits were retained and preprocessed as appropriate.

GEO calls GSM8494157-GSM8494159 `tissue: Embryo`. The [primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC11906265/) and the deposited `5M` names identify the group as five-month-old EAE mouse brain. The tracked GSE263333 manifest uses the publication-supported identity and preserves the GEO discrepancy explicitly.

## Storage and commands

Immutable provider files are under:

```text
data/raw/sources/ncbi_geo/GSE205055/
  series_metadata/<GSE>_family.soft.gz
  processed_downloads/GSE205055_RAW.tar

data/raw/sources/ncbi_geo/GSE263333/
  series_metadata/GSE263333_family.soft.gz
  processed_downloads/GSE263333_RAW.tar
```

Derived products are under:

```text
data/processed/shapemix/gse205055_spatial/
data/processed/shapemix/gse263333_spatial_mux/
```

Each processed tree contains `source_audit/`, `extracted_payload/`, `normalized_atac_fragments/`, `validation_modalities/`, `spatial_coordinates/`, `cross_modality_alignment/`, and `manifests/`. Normalized ATAC fragments are never mixed with histone/CUT&Tag fragments. RNA, protein, and histone matrices are validation evidence, not composition truth.

The reproducible entry points are:

```bash
.venv/bin/python scripts/download_shapemix_spatial.py \
  --config configs/data_sources/shapemix_gse205055.yaml --workers 4
.venv/bin/python scripts/download_shapemix_spatial.py \
  --config configs/data_sources/shapemix_gse263333.yaml --workers 4

.venv/bin/python scripts/preprocess_shapemix_spatial.py all \
  --config configs/data_sources/shapemix_gse205055.yaml
.venv/bin/python scripts/preprocess_shapemix_spatial.py all \
  --config configs/data_sources/shapemix_gse263333.yaml
```

Retained source and derived storage after deleting transfer-only segments is:

| Family | Immutable raw bytes | Processed bytes |
|---|---:|---:|
| GSE205055 | 8,177,103,162 | 16,441,652,194 |
| GSE263333 | 3,306,490,696 | 6,854,129,602 |
| Total | 11,483,593,858 | 23,295,781,796 |

The retained total is 34,779,375,654 bytes, approximately 32.4 GiB. `data/work/downloads/gse205055/` and `data/work/downloads/gse263333/` contain no remaining transfer payloads.

## Validation performed

Acquisition and extraction fail closed on official content length, frozen expected archive size, full tar traversal, unsafe paths or member types, full metadata gzip streams, and SHA-256. Extraction then compares the unique archive filenames exactly with every `!Sample_supplementary_file` URL in the parent/SubSeries SOFT files.

All 70 deposited payloads matched their metadata exactly. Every gzip or nested tar stream was fully read. All 25 fragment files use five columns, have positive intervals/support, are coordinate sorted, contain no invalid rows, and contain no adjacent duplicate `(chromosome, start, end, barcode)` rows. Each derived fragment file was written as BGZF, tabix-indexed, hashed, and successfully reopened through `pysam`. All 17 H5AD outputs were also reopened and matched their recorded dimensions.

| Family | Fragment files | ATAC files | Total fragment rows | ATAC rows | Validation-epigenome rows | Spatial archives | Sparse matrices |
|---|---:|---:|---:|---:|---:|---:|---:|
| GSE205055 | 11 | 6 | 765,673,669 | 481,802,842 | 283,870,827 | 16 | 11: 6 RNA and 5 histone |
| GSE263333 | 14 | 2 | 336,439,641 | 133,668,961 | 202,770,680 | 12 | 6: 4 RNA and 2 protein |

The ATAC inputs are:

| Family/sample | Pixels with ATAC fragments | Fragment rows | `<100 bp` | `100-249 bp` | `>=250 bp` |
|---|---:|---:|---:|---:|---:|
| GSE205055 GSM6204623, mouse brain P21 | 2,500 | 33,554,132 | 11,059,919 | 10,410,519 | 12,083,694 |
| GSE205055 GSM6204624, mouse embryo E13 25 um | 9,966 | 156,016,351 | 49,992,672 | 51,667,914 | 54,355,765 |
| GSE205055 GSM6206884, human hippocampus | 2,500 | 44,940,844 | 12,464,058 | 15,496,742 | 16,980,044 |
| GSE205055 GSM6758284, mouse brain P21 replicate | 2,500 | 48,112,158 | 16,883,546 | 16,169,670 | 15,058,942 |
| GSE205055 GSM6758285, mouse brain P22 100-barcode replicate | 10,000 | 154,191,592 | 55,410,616 | 49,666,111 | 49,114,865 |
| GSE205055 GSM6801813, mouse embryo E13 50 um | 2,500 | 44,987,765 | 10,011,508 | 13,737,039 | 21,239,218 |
| GSE263333 GSM8189706, mouse embryo E13 | 2,500 | 108,080,932 | 77,796,173 | 30,275,591 | 9,168 |
| GSE263333 GSM8494157, five-month EAE mouse brain | 10,000 | 25,588,029 | 7,711,321 | 8,033,674 | 9,843,034 |

## Barcode alignment

GSE205055 fragment barcodes universally have a terminal `-1`, while coordinate and matrix barcodes omit it. The manifest now declares a collision-checked canonical mapping that removes this suffix for comparisons while leaving deposited and normalized fragment files unchanged. The rule held for all 64,966 fragment-barcode entries across 11 files and was injective within every file.

After canonicalization:

- human hippocampus ATAC, RNA, and both coordinate sets match exactly at 2,500 pixels;
- the primary P21 mouse brain group has complete fragment/coordinate grids, with 2,373 RNA and 2,387 histone matrix pixels;
- the E13 25-um group has 9,966 ATAC-positive and 8,886 RNA-positive pixels on a 10,000 grid;
- the P21 replicate group has 2,500 ATAC pixels, 2,498 RNA pixels, and 2,499 histone-matrix pixels;
- the P22 group has 10,000 ATAC pixels and 9,215 RNA-matrix pixels, with 9,370-9,752 pixels in the three histone matrices; and
- the E13 50-um group has 2,500 ATAC pixels and 2,187 RNA-matrix pixels.

GSE263333 identifiers require no suffix mapping. The E13 ATAC/RNA/spatial group matches exactly at 2,500 pixels, and the five-month ATAC/RNA/protein/spatial group matches exactly at 10,000 pixels. Validation-only omissions are retained: one E13 protein matrix has 2,497 pixels, and the five-month H3K27ac fragment file has 9,994 barcodes.

## Remaining gate

These families are source-ready, not yet runnable ShapeMix datasets. They intentionally were not registered under `data/processed/datasets/`. A compatible labeled, fragment-level single-cell ATAC reference must still be selected and audited for each species/tissue/stage, followed by reference-only peak selection, signature estimation, spatial shape-layer construction, and `.X == sum(layers)` validation.

Current reference candidates are recorded as audit-only in the source manifests. GSE246791 is preferred for adult mouse brain, GSE216371 is the preferred E10.5-E13.5 whole-embryo reference, and GSE244618 is a large human-brain candidate requiring a frozen hippocampal subset. Each source still requires an explicit fragment, barcode, label, and coordinate gate before use.

Real spatial samples have no exact cell-proportion truth. RNA, protein, histone marks, tissue images, marker accessibility, anatomy, and replicate consistency must therefore remain orthogonal qualitative or proxy validation—not be written as `truth/proportions.csv`.
