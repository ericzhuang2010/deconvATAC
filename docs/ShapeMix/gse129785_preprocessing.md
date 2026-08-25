# GSE129785 acquisition and preprocessing record

Status: acquisition and preprocessing completed on 2026-08-24

## Scope

The frozen ShapeMix acquisition uses 30 GSE129785 scATAC fragment files:

- 14 physical dilution samples: seven CD4-memory/CD8-naive ratios and seven monocyte/T-cell ratios;
- 9 author-sorted immune populations for the single-cell reference;
- PBMC Rep1–Rep4 for independent unsorted-PBMC evaluation; and
- fresh, frozen-sorted, and frozen-unsorted PBMC samples for preparation-mismatch analysis.

Tumor samples, bone marrow, CD34 progenitors, human/mouse cell-line mixtures, and redundant sorted-cell replicates are excluded. The exact GSM, SRX, title, role, component, and nominal-ratio records are frozen in `configs/data_sources/shapemix_gse129785.yaml`.

## Storage

Immutable provider bytes are stored under:

```text
data/raw/sources/ncbi_geo/GSE129785/
  series_metadata/
  samples/<GSM>/source_files/<GSM>_<title>_fragments.tsv.gz
```

Incomplete transfers remain under `data/work/downloads/gse129785/` with a `.part` suffix. A file is promoted only after its byte size matches the provider `Content-Length`, payload integrity passes, its SHA-256 is computed, and fragment files pass a five-column schema check.

Reusable products are stored under:

```text
data/processed/shapemix/gse129785_immune/
  source_audit/
  normalized_fragments/hg19/<GSM>/
  labels/
  feature_axes/
  fragment_shape_cache/samples/
```

Twenty-nine provider fragment files are BGZF and are hard-linked into `normalized_fragments/` before tabix indexing. GSM3722026 is ordinary gzip; it is atomically recompressed to a separate derived BGZF file while the immutable raw gzip remains unchanged. Generated indexes stay outside immutable source directories. Removing `data/raw/` later does not remove retained hard-linked or recompressed processed fragments.

The standardized reference and runnable evaluation datasets are written to:

```text
data/processed/references/gse129785_immune/atac/reference.h5ad
data/processed/datasets/gse129785_shapemix_physical_dilution_*/
data/processed/datasets/gse129785_shapemix_pbmc_replicates/
data/processed/datasets/gse129785_shapemix_preparation_comparison/
```

Every external descriptor declares the ordered nine-type reference universe and fixed ShapeMix outer/inner seeds. This makes unlabelled PBMC and preparation cohorts runnable without pretending they have quantitative truth.

## Reproducible commands

Acquisition:

```bash
.venv/bin/python scripts/download_gse129785.py --workers 4
```

Preprocessing can run as one restartable pipeline:

```bash
.venv/bin/python scripts/preprocess_gse129785.py all
```

or by stage:

```bash
.venv/bin/python scripts/preprocess_gse129785.py audit
.venv/bin/python scripts/preprocess_gse129785.py select-peaks
.venv/bin/python scripts/preprocess_gse129785.py validate-coordinates
.venv/bin/python scripts/preprocess_gse129785.py build-shapes
.venv/bin/python scripts/preprocess_gse129785.py materialize
```

Existing audit and shape-cache outputs are reused unless `--overwrite` is supplied.

## Preprocessing decisions

1. The first within-study benchmark remains in hg19. Fragment intervals are not lifted over.
2. The nine sorted-population labels come from sample identities, not unsupervised cluster interpretation.
3. Sorted-reference cells, PBMC Rep1–Rep4, fresh PBMC, and frozen-unsorted PBMC use author-published filtered barcodes.
4. GEO does not publish usable raw barcode calls for the 14 dilution files, and the frozen-sorted barcode table contains `NA` in its raw-barcode field. These samples use a declared fallback of at least 1,000 fragment rows per barcode. Complete barcode-count tables and call sources are retained.
5. Author peak labels are one-based closed genomic ranges and are normalized to zero-based half-open intervals as `[start - 1, end)` before tabix queries or feature construction.
6. Reference peak selection reads only the nine sorted reference populations. The large Matrix Market file is streamed into per-type counts and per-peak cell coverage without loading the complete matrix into memory.
7. Exactly 5,000 peaks are ranked with ShapeMix protocol-v1 scoring: variance across log-normalized type aggregates, at least 10 nonzero reference cells, then deterministic coverage/count/identifier tie-breaking.
8. Each deduplicated fragment contributes two cut sites grouped by parent-fragment length: `[0,100)`, `[100,250)`, and `[250,infinity)` base pairs. `readSupport` is ignored as a weight and retained only for audit.
9. Tabix queries begin at `peakStart - 1` so a fragment whose reported end equals a peak start is available for the explicit right-cut convention test.

The coordinate audit used 64 author-called GSM3722032 cells and 256 selected peaks (9,528 published counts). After normalizing GEO peak labels from one-based closed to zero-based half-open intervals, `right_cut_offset = -1` reproduced all 9,528 counts with zero mismatches and zero error. Offset 0 produced 14 mismatched entries and absolute error 14. Weighting fragment rows by `readSupport` overcounted by more than 17,000 counts, confirming the declared one-row-per-deduplicated-fragment policy.

## Interpretation of physical truth

The CD4-memory/CD8-naive series is quantitative external validation because both components map directly to sorted reference populations. Its ratios are still nominal sample-level inputs, not exact recovered-cell proportions; sorting error and differential cell recovery must be reported as uncertainty.

The monocyte/T-cell series has only broad T-cell truth while the reference contains several T-cell subtypes. It is materialized with broad nominal proportions under `validation/`, but marked exploratory rather than exact subtype-level truth.

## Validation gates

Completion requires:

- all 43 resources (30 fragments and 13 author metadata/matrix files) to have exact sizes and SHA-256 hashes;
- full gzip/BGZF integrity for compressed sources and header/size validation for GEO matrices served as plain text despite `.gz` suffixes;
- zero invalid fragment rows and unambiguous called-cell namespaces;
- exact author matrix dimensions and ordered cell/peak axes;
- exact reconstruction of an author count-matrix subset with `right_cut_offset = -1` while recording the `0` candidate;
- identical ordered 5,000-peak axes in every reference and mixture;
- nonnegative integer CSR layers whose elementwise sum equals `.X`;
- a valid standardized reference label universe; and
- runner-readable dataset descriptors with nominal-evidence limitations stated explicitly
  and no nominal ratios declared through the exact-truth contract.

## Final observed results

- All 43 frozen resources are present with no extras or missing files. The exact downloaded payload is 40,071,277,814 bytes (37.319 GiB); sizes and SHA-256 values are frozen in `configs/data_sources/shapemix_gse129785_lock.yaml`.
- The 30 fragment files contain 3,112,258,770 rows and retain 497,297 sample barcodes/cells. By role: physical dilutions retain 442,270, sorted references 14,688, PBMC replicates 21,570, and preparation samples 18,769. All full scans report zero invalid rows.
- Twenty-nine fragment sources are retained as hard links and GSM3722026 is retained as an atomically recompressed derived BGZF. All 30 normalized fragments have tabix indexes and all 30 samples have validated 5,000-peak ShapeMix caches.
- The standardized reference contains 14,688 author-called cells across nine ordered immune types and 5,000 reference-only peaks. The coordinate audit exactly reproduced 9,528 published counts with offset `-1`.
- Sixteen runner-readable datasets were materialized and validated through `load_deconvolution_input`: 14 physical dilutions, one four-replicate PBMC cohort, and one three-preparation cohort. All 14 physical dilutions carry nominal sample-level ratios under `validation/nominal_broad_proportions.csv`; none declares exact truth. The two cohort datasets also remain explicitly no-truth.
- Retained processed products occupy 33,270,269,157 bytes (30.985 GiB): 33,208,206,654 bytes in the reusable family tree plus 62,062,503 bytes for the standardized reference and runnable datasets. Because 29 raw fragments are hard-linked, processed products add only 0.968 GiB while raw paths remain. Removing raw paths later would free about 7.302 GiB and leave the 30.985 GiB retained processed set; no raw file was removed in this run.
