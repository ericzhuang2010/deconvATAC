# ShapeMix additional dataset acquisition and preprocessing plan

Status: GSE129785, GSE194122, GSE205055, and GSE263333 acquisition and source-ready preprocessing completed on 2026-08-24; the real-spatial reference gates and downstream external experiments remain planned

This document specifies the additional data needed after the completed one-donor PBMC benchmark. It defines the scientific role of each source, where its raw and processed files should live, and the preprocessing required before it can be used by ShapeMix. It does not change, replace, or reinterpret the completed protocol-v1 PBMC result.

## 1. Objectives

The additional data serve four distinct purposes:

1. Create controlled stress datasets from the existing PBMC cells.
2. Test rare-cell recovery and preparation mismatch in independent PBMC/immune data.
3. Measure generalization across biological donors and collection sites.
4. Apply a frozen model to real non-blood spatial ATAC data with RNA, protein, histology, or anatomical evidence for qualitative validation.

The planned sources are:

| Family | Organism and tissue | Main role | Exact composition truth? | Raw fragment availability | Preprocessing required? |
|---|---|---|---|---|---|
| Existing 10x PBMC sensitivity datasets | Human peripheral blood | Controlled depth, spot-size, rare-cell, subtype, feature, and bin stress tests | Yes, from recorded source cells | Already available locally | Yes; derived datasets only |
| [GSE129785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785) | Human PBMC, sorted immune cells, bone marrow, and tumor samples | Independent immune data, physical dilution series, and fresh/frozen mismatch | No exact truth for physical dilutions; nominal ratios are validation evidence. Exact truth is available only for newly recorded pseudo-spots | Per-sample fragment TSV files and SRA reads | Yes |
| [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122) | Human bone marrow mononuclear cells from 12 healthy donors and four sites | Donor- and site-level quantitative generalization | Yes for pseudo-spots built from annotated held-out donor cells | Raw reads in SRA; processed Multiome H5AD on GEO | Yes; substantial |
| [GSE205055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205055), especially ATAC SubSeries [GSE205052](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205052) | Mouse embryo, mouse brain, and human brain spatial epigenome/transcriptome samples | Real spatial ATAC validation | No exact proportions | Spatial fragment TSVs and SRA reads | Yes |
| [GSE263333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263333) | Mouse spatial multi-omics tissue sections | Real spatial ATAC validation with RNA, protein, and histone-mark evidence | No exact proportions | Spatial ATAC fragment TSVs and SRA reads | Yes |

The existing Heart and Russell objects are not additional ShapeMix inputs. They contain collapsed peak matrices and cannot recover parent-fragment lengths. They can become ShapeMix datasets only if their original raw fragment-level data, labels, and donor metadata are obtained and separately audited.

## 2. Directory and tracking policy

### 2.1 Raw external sources

Preserve every downloaded file byte-for-byte under a provider/accession directory:

```text
data/raw/sources/ncbi_geo/<GSE_ACCESSION>/
  series_metadata/
  processed_downloads/
  samples/<GSM_ACCESSION>/
    source_files/
    sra/
```

Examples:

```text
data/raw/sources/ncbi_geo/GSE129785/samples/GSM3722011/source_files/
data/raw/sources/ncbi_geo/GSE194122/processed_downloads/
data/raw/sources/ncbi_geo/GSE205055/processed_downloads/GSE205055_RAW.tar
data/raw/sources/ncbi_geo/GSE263333/processed_downloads/GSE263333_RAW.tar
```

Do not modify, recompress, rename, lift over, sort, or index files inside `source_files/`. Any normalized BGZF copy, tabix index, aligned BAM, reconstructed fragment table, or genome-converted artifact is derived data and belongs under `data/processed/shapemix/`.

### 2.2 Tracked source manifests

Before a large download begins, add one tracked manifest per source family:

```text
configs/data_sources/shapemix_gse129785.yaml
configs/data_sources/shapemix_gse194122.yaml
configs/data_sources/shapemix_gse205055.yaml
configs/data_sources/shapemix_gse263333.yaml
```

Each manifest must record:

- GEO, GSM, BioProject, SRA experiment, and SRA run accessions;
- source page and direct file URL where stable;
- organism, tissue, assay, genome build, donor, site, and sample role;
- expected byte size when published or observed;
- SHA-256 after download;
- original compression and fragment-file schema;
- download date and license or reuse statement;
- whether a file is raw, author-processed, or repository-derived;
- the planned reference/test role; and
- any exclusion and its reason.

The local ignored manifest beside each raw source may record resolved filesystem paths and download state. It must not be the only record of source identity.

### 2.3 Processed ShapeMix products

Use one family-level preprocessing tree:

```text
data/processed/shapemix/<family>/
  source_audit/
  normalized_fragments/
  fragment_shape_cache/
  labels/
  feature_axes/
  splits/
  manifests/
```

All standardized, repository-generated references use the processed reference tree:

```text
data/processed/references/<reference_id>/atac/reference.h5ad
data/processed/references/<reference_id>/reference.yaml
```

Final mixture or spatial datasets use:

```text
data/processed/datasets/<dataset_id>/
  atac/spatial.h5ad
  atac/features/*.txt
  truth/proportions.csv                 # only when exact truth exists
  simulation/source_cells_by_spot.jsonl # only for simulated mixtures
  validation/                            # orthogonal real-spatial evidence
  dataset.yaml
```

Register a dataset in `data/registry/datasets.yaml` only after all source, shape-layer, axis, truth, and split validations pass. Never reuse or overwrite the 21 existing ShapeMix dataset IDs.

## 3. Common preprocessing contract

Every external source requires a dataset-specific adapter. Do not assume that a file named `fragments.tsv.gz` uses the exact Cell Ranger ARC 2.0 coordinate and column contract already frozen for the current PBMC data.

### 3.1 Source acquisition and integrity

1. Export the GEO sample table and SRA RunInfo before downloading sequence data.
2. Resolve exact GSM-to-SRR, donor, site, modality, and paired-modality mappings.
3. Estimate required download and temporary disk space from the resolved files.
4. Download resumably into the raw source tree.
5. Record byte sizes and SHA-256 hashes.
6. Validate gzip/BGZF, tar, HDF5, BAM, or SRA container integrity as applicable.
7. Retain author-provided metadata and README files with the data.

### 3.2 Fragment schema and coordinate audit

For each technology and source version:

1. Determine the number and meaning of columns, coordinate system, barcode convention, duplicate policy, and support-count semantics.
2. Confirm that fragment lengths are positive and reflect paired Tn5 ends rather than peak widths or single-read alignments.
3. Establish the left- and right-cut coordinate convention against an author-provided peak matrix or a small independently reconstructed matrix.
4. Decide whether the primary count unit remains deduplicated cut sites grouped by parent-fragment length.
5. Preserve the source assembly. Do not lift fragment intervals between assemblies as a substitute for a coordinate-valid fragment reconstruction.
6. If random access is needed, write a sorted BGZF copy and tabix index under `data/processed/shapemix/<family>/normalized_fragments/`.

The current PBMC finding that `right_cut = chromEnd` is not automatically transferable to an older Cell Ranger ATAC file, a BED-like spatial file, sci-ATAC, or a custom spatial-Mux file.

### 3.3 Labels, feature axes, and shape layers

1. Harmonize author labels to an explicit, versioned cell-type ontology without silently merging types after results are inspected.
2. Require adequate reference and test support for every declared cell type.
3. Select peaks from training-reference cells only for quantitative simulations.
4. For real spatial data, select peaks from the external single-cell reference; do not use spatial outcomes to optimize the feature axis.
5. Count the frozen fragment-length layers and require `.X` to equal their exact sum.
6. Record ordered feature, cell-type, source-file, split, and preprocessing hashes.
7. Audit shape coverage, entropy, between-type divergence, split/donor reproducibility, and correlations with depth and QC.

### 3.4 Evaluation roles

- Simulated pseudo-spots have exact truth and may be evaluated with `rmse_v1`, `jsd_v2`, per-type error, and rare-cell detection metrics.
- Physical dilution samples have nominal preparation ratios at the sample level. These ratios must not be represented as exact per-pixel or per-cell truth.
- Real spatial samples generally lack exact composition truth. RNA-, protein-, histology-, marker-, or anatomy-derived scores are orthogonal validation evidence, not exact ground truth.

## 4. Existing PBMC stress and sensitivity datasets

### 4.1 Description and role

These datasets reuse the already checksum-pinned 10x PBMC fragments, the five reference/held-out splits, and the held-out source-cell shape layers. They require no new external raw download. Their purpose is to determine when fragment length becomes useful or harmful as depth, mixture complexity, rare-cell abundance, feature count, and binning change.

### 4.2 Planned locations

```text
data/processed/shapemix/pbmc_granulocyte_sorted_10k/sensitivity_v1/
  depth_thinning/
  cells_per_spot/
  rare_cell_abundance/
  subtype_challenge/
  feature_count/
  fragment_length_bins/
  background/
  exposure_normalization/

data/processed/datasets/
  pbmc_granulocyte_sorted_10k_shapemix_sensitivity_<factor>_<level>_split_<seed>_mix_<seed>/
```

### 4.3 Required preprocessing

Yes, but only derived preprocessing:

- generate every condition from held-out cells layer by layer;
- use new dataset IDs and leave the primary objects unchanged;
- record exact source cells and any fragment-thinning mask;
- recompute `.X` from the new layers;
- independently select 5,000, 10,000, and 20,000 reference-only peaks where requested;
- create separately versioned two-, three-, and five-bin objects rather than relabeling the primary layers;
- freeze factor levels and seeds before inspecting results; and
- validate exact truth and source-cell disjointness.

Run one factor at a time first. A small, predeclared depth-by-rare-cell interaction can follow; a complete factorial grid would be unnecessarily large and difficult to interpret.

## 5. GSE129785 immune-cell and dilution data

Status: the frozen 30-sample immune scope was acquired and preprocessed on 2026-08-24. See the [GSE129785 preprocessing record](gse129785_preprocessing.md) and tracked size/hash lock `configs/data_sources/shapemix_gse129785_lock.yaml`.

### 5.1 Description and scientific use

[GSE129785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785) is a human scATAC-seq study containing PBMC replicates, sorted immune populations, preparation comparisons, and controlled two-population dilution series. The study used Cell Ranger ATAC 1.0 and reports hg19 coordinates. Individual GEO samples provide fragment TSV files; for example, [GSM3722011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3722011) provides the `0.1% CD4 memory / 99.9% CD8 naive` mixture and its fragment table.

The ShapeMix acquisition scope is the immune portion of the study, not all tumor samples.

| Sample group | Planned use |
|---|---|
| Paired CD4-memory/CD8-naive dilution series at 0.1/99.9, 0.5/99.5, 1/99, 50/50, 99/1, 99.5/0.5, and 99.9/0.1 | Physical rare-cell and mixture-ratio validation |
| Corresponding monocyte/T-cell dilution series | A second physical mixture family after its exact T-cell definition is resolved from author metadata |
| Sorted dendritic cells, monocytes, B cells, Treg, naive/memory CD4 cells, NK cells, and naive/memory CD8 cells | Reference construction and held-out pseudo-spots |
| PBMC Rep1–Rep4 | Independent unsorted PBMC evaluation and reference-stability audit |
| Fresh, frozen-sorted, and frozen-unsorted PBMC samples | Protocol/preparation mismatch analysis |

The exact GSM list and every physical ratio must be frozen in the tracked source manifest before results are generated. The series-level GEO page currently reports a 66.7 GB complete supplementary archive, so downloading selected per-sample fragment files is preferable to acquiring unrelated tumor samples.

### 5.2 Planned locations

```text
data/raw/sources/ncbi_geo/GSE129785/
  series_metadata/
  samples/<GSM>/source_files/<author_fragment_file>.tsv.gz

configs/data_sources/shapemix_gse129785.yaml

data/processed/shapemix/gse129785_immune/
  source_audit/
  normalized_fragments/hg19/
  fragment_shape_cache/
  labels/
  feature_axes/
  splits/

data/processed/references/gse129785_immune/atac/reference.h5ad

data/processed/datasets/
  gse129785_shapemix_pseudospot_<condition>_split_<seed>_mix_<seed>/
  gse129785_shapemix_physical_dilution_<mixture_family>_<ratio>/
  gse129785_shapemix_preparation_<reference_prep>_to_<test_prep>/
```

### 5.3 Required preprocessing

Yes:

1. Inspect the Cell Ranger ATAC 1.0 fragment schema and duplicate/support semantics.
2. Validate hg19 fragment coordinates against an author peak matrix.
3. Keep the first within-study benchmark in hg19. If cross-study transfer from the current GRCh38 PBMC reference is required, reprocess raw reads against GRCh38 with a pinned pipeline; do not lift over individual fragment intervals and assume their lengths remain equivalent.
4. Normalize barcode identifiers across samples without creating collisions.
5. Build sorted-population label tables and retain preparation/sample identity.
6. Construct cell-disjoint or sample-disjoint reference and pseudo-spot pools where the design permits.
7. Treat physical dilution ratios as nominal sample-level truth and report binomial sampling and sorting uncertainty separately.
8. Run a fresh/frozen/frozen-sorted shape-signal audit before using cross-preparation signatures.

This source should be the first external acquisition pilot because author-provided fragment files avoid a large raw-read reconstruction step and the dilution samples directly target the rare-cell hypothesis.

## 6. GSE194122 multi-donor BMMC Multiome data

Status: complete through source-ready preprocessing. See the [GSE194122 preprocessing record](gse194122_preprocessing.md) and tracked lock `configs/data_sources/shapemix_gse194122_lock.yaml`.

### 6.1 Description and scientific use

[GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122) contains single-cell multi-omics data from 12 healthy human bone-marrow donors collected at four sites. Its Multiome arm contains paired RNA and ATAC measurements and author-provided cell annotations. Donor 1 appears at multiple sites, which permits a useful separation of donor effects from site effects.

The primary role is quantitative leave-one-donor-out evaluation:

- select peaks and estimate signatures using training donors only;
- create pseudo-spots exclusively from the held-out donor;
- retain inner mixture seeds within donor;
- summarize effects first by donor; and
- never treat spots or mixture seeds as biological replicates.

The acquisition includes the processed `GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz` file and all 13 Multiome ATAC site/donor samples listed on GEO. GEO exposes the author Cell Ranger ARC 2.0 BAMs through SRA rather than as supplementary fragment files. The public Open Problems post-competition bucket, referenced by a public primary-analysis workflow, provides the 13 corresponding author Cell Ranger ARC fragment files, tabix indexes, and per-barcode metrics files. Those 39 objects total 31,687,758,209 bytes, compared with 531,234,571,942 bytes for the frozen NCBI BAM fallback.

### 6.2 Planned locations

```text
data/raw/sources/ncbi_geo/GSE194122/
  series_metadata/
  processed_downloads/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz
  samples/<GSM>/source_files/
    atac_fragments.tsv.gz
    atac_fragments.tsv.gz.tbi
    per_barcode_metrics.csv

configs/data_sources/shapemix_gse194122.yaml

data/processed/shapemix/gse194122_bmmc/
  source_audit/
  source_audit/source_objects/
  normalized_fragments/GRCh38/<sample_key>/
  fragment_shape_cache/
  labels/source_broad7_v1/
  feature_axes/
  splits/broad7_lodo_v1/<heldout_donor>/

data/processed/references/gse194122_bmmc_lodo_<heldout_donor>/atac/reference.h5ad

data/processed/datasets/
  gse194122_bmmc_shapemix_lodo_<heldout_donor>_<condition>_mix_<seed>/
```

### 6.3 Required preprocessing

Yes. Source-ready preprocessing now does the following:

1. Freeze all 13 Multiome ATAC GSM/SRR/BioSample, donor/site, original BAM, and author fragment-product identities.
2. Audit the 69,249-cell by 129,921-feature H5AD, including 22 author cell types, ten donors, four sites, 116,490 ATAC peaks, 13,431 genes, and the raw count layer.
3. Retain the 13 author fragments/indexes/metrics trios as verified hard links under the processed tree.
4. Map every processed cell by its 16-base sequence to the called Cell Ranger common `barcode` written in fragment column 4; retain the distinct raw ATAC-library barcode from `per_barcode_metrics.csv` only as provenance.
5. Compare a deterministic GSM5828489 count-matrix subset to the full fragments, proving that the depth-normalized aggregate H5AD counts are contained entrywise and freezing `chromEnd` as the right-cut coordinate with one count per deduplicated fragment row.
6. Preserve author labels without biological reinterpretation and write ten exact donor-held-out membership tables.

The remaining experiment-specific work is deliberately downstream: apply a predeclared training-support rule to the 22 labels, select peaks independently within each training fold, build fragment-shape caches and signatures on training donors only, construct held-out-donor pseudo-spots, and report within-site, cross-site, and cross-donor effects separately.

The complete SRA/BAM download is no longer necessary. The BAM identities remain frozen as a reconstruction fallback, but the authoritative author fragment suite is the primary route.

## 7. GSE205055 spatial ATAC–RNA data

Status: complete through source-ready preprocessing. See the [GSE205055/GSE263333 preprocessing record](gse205055_gse263333_preprocessing.md) and tracked lock `configs/data_sources/shapemix_gse205055_lock.yaml`.

### 7.1 Description and scientific use

[GSE205055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205055) is a spatial epigenome/transcriptome SuperSeries. The complete acquired family includes mouse CUT&Tag GSE205051, mouse ATAC [GSE205052](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205052), mouse RNA GSE205054, human ATAC GSE205180, human RNA GSE205181, mouse spatial CUT&Tag replicates GSE217091, and mouse RNA replicates GSE218593. Together these SubSeries contain mouse embryo, mouse brain, and human hippocampus sections at multiple spatial resolutions, with spatial coordinates and image assets.

The immutable acquisition is the complete parent processed-data archive rather than a hand-picked sample subset: 38 unique supplementary files from 22 samples. Six matched ATAC/RNA groups are available for later ShapeMix analysis; the remaining deposited modalities are retained as orthogonal validation evidence. [GSM6801813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6801813), for example, is an E13 mouse embryo 50-um ATAC sample with a deposited BED-like fragment file and spatial archive.

### 7.2 Planned locations

```text
data/raw/sources/ncbi_geo/GSE205055/
  series_metadata/<GSE>_family.soft.gz
  processed_downloads/GSE205055_RAW.tar

configs/data_sources/shapemix_gse205055.yaml

data/processed/shapemix/gse205055_spatial/
  source_audit/
  extracted_payload/<GSM>/
  normalized_atac_fragments/<GSM>/
  validation_modalities/{epigenome,rna}/<GSM>/
  spatial_coordinates/<GSM>/
  cross_modality_alignment/<group>.yaml
  feature_axes/<reference_id>/
  fragment_shape_spatial/<sample>/

data/processed/datasets/
  gse205055_<species>_<tissue>_<resolution>_<sample>_shapemix_real_spatial/
```

Each final real-spatial dataset should contain:

```text
validation/
  rna_cell_type_scores.csv
  marker_accessibility_scores.csv
  anatomical_region_annotations.csv
  cross_modality_alignment.yaml
```

These files are validation evidence, not `truth/proportions.csv`.

### 7.3 Required preprocessing

Yes:

1. Match ATAC, RNA, image, spatial-coordinate, resolution, species, tissue, and replicate metadata.
2. Parse and validate the custom BED-like fragment schema and spatial barcode format.
3. Determine whether rows are already deduplicated and whether the deposited interval endpoints are adjusted Tn5 coordinates.
4. Normalize fragment files into an indexed internal representation without altering the raw downloads.
5. Count spatial spots onto a feature axis selected from a compatible single-cell reference.
6. Align ATAC and RNA pixels exactly; document any missing or non-overlapping pixels.
7. Derive RNA-based cell-type scores with a method frozen independently of the ShapeMix comparison.
8. Evaluate spatial marker agreement, anatomy, replicate consistency, boundary behavior, and shape-versus-count map stability.

### 7.4 Required single-cell references

The spatial series supplies mixtures but not labeled single cells for fixed ShapeMix signatures. A compatible fragment-level scATAC or Multiome reference is mandatory before a section can be deconvolved.

Reference selection is a preprocessing gate, not an afterthought:

- match species, tissue/region, developmental stage when possible, and genome assembly;
- require fragment-level data and cell-type annotations;
- record protocol mismatch explicitly;
- freeze label harmonization and the shared cell-type universe before inspecting spatial maps; and
- select peaks from the single-cell reference only.

The preferred adult mouse-brain candidate is [GSE246791](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246791), whose processed per-sample H5AD objects retain raw mm10 fragments; its labels and parent-fragment semantics still require an audit. [GSE111586](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111586) remains a fallback, but its processed matrices are mm9 and lack fragment lengths. The human brain candidate is [GSE244618](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244618), frozen to hippocampal regions and donors before download. The embryo candidate is [GSE216371](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216371), an E10.5-E13.5 whole-embryo mm10 SPATAC atlas with cell-level annotations and deposited fragment BED files. These candidates still require fragment, barcode, label, and coordinate gates.

## 8. GSE263333 spatial-Mux-seq data

Status: complete through source-ready preprocessing. See the [GSE205055/GSE263333 preprocessing record](gse205055_gse263333_preprocessing.md) and tracked lock `configs/data_sources/shapemix_gse263333_lock.yaml`.

### 8.1 Description and scientific use

[GSE263333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263333) contains mouse spatial multi-omics measurements. It combines spatial ATAC with RNA and, for selected sections, proteins and histone modifications. This makes it a valuable orthogonal validation family after the simpler GSE205055 pilot.

The immutable acquisition contains the complete 12-sample, 32-file processed deposit. The eventual ShapeMix scope contains the two sample groups with explicit ATAC fragments:

| ATAC sample | Paired evidence | Notes |
|---|---|---|
| GSM8189706, `ME13_50um_3_ATAC_H3K4me3_H3K27me3` | GSM8189707 RNA plus H3K4me3/H3K27me3 and spatial archive | The ATAC fragment file is approximately 941 MB |
| GSM8494157, `5M_20um_ATAC_H3K27ac_H3K27me3` | GSM8494158 RNA, GSM8494159 protein, histone marks, and spatial archive | The ATAC fragment file is approximately 266 MB |

Only files explicitly identified as ATAC fragments may populate ShapeMix fragment-length layers. CUT&Tag fragment lengths measure a different assay and must remain separate validation modalities.

GEO labels GSM8494157-GSM8494159 as `tissue: Embryo`, but the primary publication and the deposited `5M` names identify this group as five-month-old EAE mouse brain. The preprocessing manifest uses the publication-supported identity and permanently records the GEO discrepancy; it must not be silently discarded.

### 8.2 Planned locations

```text
data/raw/sources/ncbi_geo/GSE263333/
  series_metadata/GSE263333_family.soft.gz
  processed_downloads/GSE263333_RAW.tar

configs/data_sources/shapemix_gse263333.yaml

data/processed/shapemix/gse263333_spatial_mux/
  source_audit/
  extracted_payload/<GSM>/
  normalized_atac_fragments/<sample>/
  spatial_coordinates/<sample>/
  cross_modality_alignment/<sample>/
  validation_modalities/{epigenome,rna,protein}/<sample>/
  fragment_shape_spatial/<sample>/

data/processed/datasets/
  gse263333_me13_50um_3_shapemix_real_spatial/
  gse263333_5m_20um_shapemix_real_spatial/
```

### 8.3 Required preprocessing

Yes:

1. Validate the custom spatial fragment and barcode schema for each section.
2. Separate ATAC, H3K27ac, H3K27me3, and H3K4me3 files by assay in both metadata and paths.
3. Recover exact spatial coordinates from the spatial archives.
4. Confirm that ATAC/RNA/protein pixels refer to the same grid before cross-modality comparison.
5. Build the ShapeMix layers from ATAC fragments only.
6. Select and preprocess a compatible labeled single-cell fragment reference using the same gate described for GSE205055.
7. Store RNA-, protein-, histone-, marker-, and anatomy-based validation results separately from truth-based metrics.

## 9. Acquisition order and gates

Acquire and preprocess in this order:

1. Generate the existing-PBMC stress datasets; no external download is required.
2. Download and validate one GSE129785 dilution fragment file and its matching pure/reference samples.
3. Freeze the full ShapeMix-relevant GSE129785 sample list, then acquire that immune subset.
4. Acquire the GSE194122 processed H5AD and complete 13-sample author fragment/index/metrics suite. Completed 2026-08-24; the 531 GB NCBI BAM inventory is retained only as a frozen fallback.
5. Audit the GSE194122 barcode bridge and raw count semantics, then freeze donor-held-out membership indices. Completed 2026-08-24; fold-specific peaks, signatures, pseudo-spots, and fits remain downstream experiments.
6. Acquire and source-preprocess the complete GSE205055 parent archive and all related SubSeries metadata. Completed 2026-08-24; reference selection remains gated.
7. Acquire and source-preprocess the complete GSE263333 archive, including both ATAC groups and all orthogonal modalities. Completed 2026-08-24; reference selection remains gated.
8. Select and audit compatible fragment-level references before constructing any runnable real-spatial ShapeMix dataset.
9. Expand real-spatial validation to additional sections or species only after the first reference-aligned section passes.

The following gate must be satisfied before a source advances:

| Gate | Required evidence |
|---|---|
| Download | Source accessions, URLs, sizes, and hashes are recorded |
| Fragment semantics | Column, duplicate, coordinate, length, and cut-site conventions are validated |
| Barcode mapping | Every retained fragment barcode maps unambiguously to a declared cell or spatial pixel |
| Feature axis | Reference and mixture/spatial objects use the same ordered peaks |
| Shape contract | Nonnegative integer CSR layers sum exactly to `.X` |
| Quantitative truth | Source-cell records reproduce every pseudo-spot truth row exactly |
| Donor separation | Held-out donor cells never enter peaks, signatures, smoothing, or tuning |
| Spatial alignment | ATAC and orthogonal modality pixels/coordinates are reconciled and missingness is reported |
| Provenance | Raw, code, configuration, label, split, and output hashes are complete |

## 10. Expected research products

After preprocessing, the additional data should support three separate result classes:

1. **Controlled sensitivity results:** exact paired ShapeMix-versus-count-only effects across depth, spot size, rare abundance, subtype similarity, feature count, bins, background, and exposure normalization.
2. **External quantitative results:** physical dilution performance in GSE129785 and donor-level effects in GSE194122.
3. **Real-spatial qualitative results:** spatial maps and orthogonal concordance in GSE205055 and GSE263333, explicitly without claiming exact composition accuracy.

The completed one-donor PBMC result remains protocol version 1. Any model tuning motivated by these new sources requires a new protocol version, development-only tuning data, and frozen external evaluation sets.
