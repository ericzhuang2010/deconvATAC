# GSE194122 acquisition and preprocessing record

Status: acquisition and source-ready preprocessing completed on 2026-08-24

## Scope and purpose

[GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122) contains healthy human bone-marrow mononuclear cells profiled at four sites. The Multiome arm has 69,249 author-filtered and annotated cells in 13 site/donor samples. Ten donors occur in the Multiome arm; donor 1 occurs at all four sites and donor 3 occurs at two sites.

This family is retained for quantitative leave-one-donor-out ShapeMix evaluation. This preprocessing stage makes the author counts, annotations, fragments, and donor/site split indices source-ready. It does not select fold-specific peaks, estimate signatures, generate pseudo-spots, or run deconvolution. Those operations must remain downstream so the held-out donor cannot leak into feature selection or signature estimation.

## Source resolution

GEO provides the processed `GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz` object and links the 13 ATAC libraries to SRA. NCBI retains the original Cell Ranger ARC ATAC BAM for every run, but those BAMs total 531,234,571,942 bytes.

The public Open Problems post-competition bucket contains the corresponding 13 author-generated Cell Ranger ARC fragment files, tabix indexes, and per-barcode metrics files. The primary research workflow explicitly acquires this prefix as the BMMC Multiome fragment source: [epi-immune workflow](https://github.com/timoast/epi-immune/blob/master/Snakefile). The complete 39-object fragment suite is 31,687,758,209 bytes, so it is both more authoritative for fragment semantics and about 16.8 times smaller than retaining the BAM fallback.

The tracked manifest freezes all 13 GSM, BioSample, SRX, SRR, donor/site, original BAM URL/size/MD5, and fragment-suite object size/ETag identities:

```text
configs/data_sources/shapemix_gse194122.yaml
configs/data_sources/shapemix_gse194122_lock.yaml
```

## Storage

Immutable provider bytes are stored under:

```text
data/raw/sources/ncbi_geo/GSE194122/
  series_metadata/GSE194122_family.soft.gz
  processed_downloads/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz
  samples/<GSM>/source_files/
    atac_fragments.tsv.gz
    atac_fragments.tsv.gz.tbi
    per_barcode_metrics.csv
```

Incomplete, resumable transfers use `data/work/downloads/gse194122/` and are never promoted into `data/raw/` until the frozen byte size and remote identity match and payload validation passes.

Reusable products are stored under:

```text
data/processed/shapemix/gse194122_bmmc/
  source_audit/source_objects/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad
  source_audit/
  normalized_fragments/GRCh38/<sample_key>/
  labels/source_broad7_v1/cells.tsv.gz
  feature_axes/source_axis_v1/features.tsv.gz
  splits/broad7_lodo_v1/donor_<N>/cells.tsv.gz
```

The normalized fragment, index, and metrics paths are hard links to the verified immutable downloads. This preserves exact provider bytes and keeps the processed products usable if raw path names are later removed, without consuming a second 31.7 GB of disk blocks.

The completed acquisition contains 41 validated source files totaling 34,604,880,348 bytes: two GEO core files plus 13 fragment/index/metrics trios. On the current filesystem, the immutable source tree occupies about 33 GB and the additional processed products about 3.0 GB; the hard-linked fragment suite is not stored twice. No incomplete transfer files remain in the work directory.

## Reproducible commands

Acquire the GEO core and all author fragment products:

```bash
.venv/bin/python scripts/download_gse194122.py --include-fragments --workers 4
```

On a server that throttles single HTTP connections, pass an installed `aria2c` executable with `--aria2c /path/to/aria2c`; each object then uses four resumable ranges without changing any identity or promotion check.

For a gated pilot, `--fragment-sample s3d6` limits fragment acquisition to that sample; repeat the option to request several declared sample keys.

Run the restartable preprocessing stages:

```bash
.venv/bin/python scripts/preprocess_gse194122.py metadata
.venv/bin/python scripts/preprocess_gse194122.py fragments
.venv/bin/python scripts/preprocess_gse194122.py matrix-audit
```

During a gated pilot, `fragments --sample-key s3d6` prepares only the declared pilot. The later unfiltered `fragments` run is still required and fails unless every H5AD cell has a barcode bridge.

After acquisition, the three stages can also be run in order with:

```bash
.venv/bin/python scripts/preprocess_gse194122.py all
```

`bam-fallback` exists only to reconstruct a pilot from a frozen NCBI BAM if the author fragment suite ever fails validation. It is not part of the successful primary route.

## H5AD audit

The decompressed source object is 3,116,868,580 bytes and has shape 69,249 cells by 129,921 features. It contains:

- 116,490 ATAC peaks and 13,431 GEX features;
- the unmodified author `cell_type` labels (22 types);
- 13 `batch`/`Samplename` values across four sites and ten Multiome donors;
- the raw sparse count matrix in `layers["counts"]` with 325,339,757 stored entries; and
- author ATAC and GEX QC, embeddings, pseudotime fields, donor metadata, and gene activity.

The feature table retains source order, original feature identifiers, zero-based peak starts/ends parsed from the deposited `contig-start-end` names, and primary or alternate GRCh38 contigs. The label table retains source cell order, the original merged barcode, canonical 16-base barcode sequence, site/donor identifiers, author label, and an initially identity-preserving harmonized label. No biological label was inferred or collapsed during source preprocessing.

## Barcode reconciliation

The H5AD identifiers use a merged common-barcode namespace and contain mixed suffix forms such as `-1-s1d1`, `-2-s1d2`, or no numeric suffix before the sample tag. Cell Ranger ARC fragment column 4 uses the metrics `barcode` namespace with a per-sample `-1` suffix, so it still cannot be joined to the merged H5AD identifier by string equality. The separate metrics `atac_barcode` field records the raw ATAC-library barcode and is retained for provenance, but it is not the barcode written to the author fragment file.

For every sample, preprocessing uses the author `per_barcode_metrics.csv` bridge:

```text
H5AD cell ID
  -> 16-base common/GEX barcode sequence
  -> metrics.barcode where is_cell == 1
  -> fragment column 4
```

The mapping must be one-to-one within sample, every processed H5AD cell must be an author-called Cell Ranger cell, and no two processed cells may map to the same common fragment barcode. The fragment barcode, separate ATAC-library barcode, and author fragment/cut-site metrics are retained in `labels/source_broad7_v1/cells.tsv.gz`.

## Fragment and count semantics

Each source fragment file must be readable with its author tabix index and must pass a five-column prefix audit: contig, start, end, common fragment barcode, and positive `readSupport`. The source headers identify Cell Ranger ARC 2.0.0 and `refdata-cellranger-arc-GRCh38-2020-A-2.0.0`. The [10x fragment specification](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/fragments-file) defines these as zero-based BED-like coordinates and `readSupport` as duplicate-read support for one unique fragment.

The pilot uses GSM5828489 (`s3d6`), the smallest sample. A deterministic subset of 64 author-annotated cells and 256 nonzero chr1 peaks is compared between the retained full per-sample fragments and the H5AD raw count layer. The merged H5AD has Cell Ranger ARC `aggr` GEM-well suffixes and depth-normalized counts; [10x documents that `aggr` subsamples higher-depth libraries](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.1/analysis/running-pipelines/aggregating-multiple-gem-wells-aggr). Therefore exact equality to the pre-aggregation fragment file is neither expected nor required. The gate instead requires every depth-normalized H5AD entry to be less than or equal to its full-fragment reconstruction.

Each deduplicated fragment row contributes one left cut and one right cut; `readSupport` is not used as a weight. Both `chromEnd` and `chromEnd - 1` are tested. The selected convention must satisfy entrywise containment and uniquely minimize the discrepancy at endpoint-boundary cases. The observed pilot selects `chromEnd` with no offset. Across 16,384 audited entries, the H5AD contains 1,992 depth-normalized cuts and the full fragments contain 3,277 cuts (retained fraction 0.607873); no H5AD entry exceeds the full-fragment count. The `chromEnd - 1` alternative adds one unsupported boundary cut and therefore has the uniquely larger error.

## Leave-one-donor-out indices

Ten fold tables are written under `splits/broad7_lodo_v1/`. Each preserves the complete H5AD cell axis and marks each cell as `training` or `heldout` for that donor. Donor 1's four sites and donor 3's two sites always stay on the same side of the split. These are cell-membership indices only; no peaks or signatures are learned at this stage.

## Completion gates

Completion requires:

- exact frozen sizes and SHA-256 hashes for the GEO object and all 39 fragment-suite objects;
- full gzip integrity for the H5AD and every fragment file;
- readable author tabix indexes and valid five-column fragment prefixes;
- all 69,249 H5AD cells mapped one-to-one through called Cell Ranger metrics to common fragment barcodes;
- pilot proof that the depth-normalized H5AD counts are an entrywise subsample of full-fragment cut counts under the uniquely supported `chromEnd` convention;
- source-order-preserving cell and feature tables;
- ten donor-held-out index tables covering the exact cell axis; and
- passing focused and repository-wide tests.

All gates above passed on 2026-08-24. The final fragment audit reports all 13 samples complete and all 69,249 processed cells bridged exactly once. The frozen acquisition lock records `fragment_suite_included: true`, `fragment_samples: all`, 41 files, and 34,604,880,348 validated bytes.

## Interpretation boundary

The author labels and donor/site structure make this dataset suitable for the planned external generalization experiment, but source preprocessing alone is not a deconvolution result. Fold-specific peak selection, minimum-support rules for rare labels, pseudo-spot generation, ShapeMix/count-only fits, and donor-level effect summaries remain separate experiments.
