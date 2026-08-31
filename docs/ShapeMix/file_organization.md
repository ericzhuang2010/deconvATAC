# ShapeMix data and results organization

## 1. Purpose and decisions

This document defines where downloaded data, derived data, runnable datasets, experiment outputs, and temporary files belong. It is the canonical filesystem and tracking policy for new ShapeMix work.

Two repository-wide decisions are frozen:

1. Standardized reference objects are derived data and belong under `data/processed/references/`, not `data/raw/references/`.
2. The entire `results/` directory is exposed to Git. No primary, sensitivity, external-validation, or real-spatial result subtree is intentionally ignored.

The large `data/` tree remains ignored. Reproducibility is provided by tracked source manifests, preparation scripts, experiment configurations, validation records, and the committed `results/` tree.

## 2. Canonical directory tree

```text
configs/
  data_sources/<source_id>.yaml
  datasets/<dataset_id>.yaml
  experiments/<campaign_id>.yaml
  methods/<method_or_variant>.yaml

data/
  work/
    downloads/<source_id>/
    preprocessing/<family>/<job_id>/

  raw/
    sources/<provider>/<accession_or_release>/

  processed/
    shapemix/<family>/
      source_audit/
      normalized_fragments/
      fragment_shape_cache/
      labels/
      feature_axes/
      splits/
      manifests/

    shapemix/real_spatial_validation/
      marker_features_v1/<reference_id>/
        marker_features.tsv
        selected_features.txt
        manifest.yaml

    references/<reference_id>/
      reference.yaml
      atac/reference.h5ad
      rna/reference.h5ad                 # when available

    datasets/<dataset_id>/
      dataset.yaml
      atac/spatial.h5ad
      atac/features/*.txt
      truth/proportions.csv              # only when exact truth exists
      simulation/source_cells_by_spot.jsonl
      simulation/manifest.yaml
      validation/                        # orthogonal real-spatial evidence
        alignment.yaml
        rna/rna.h5ad                    # when available
        protein/protein.h5ad            # when available
        epigenome/*.h5ad                # when available

  registry/
    datasets.yaml

results/
  development/<campaign_id>/
  primary/<campaign_id>/
  sensitivity/<campaign_id>/
  external_validation/<campaign_id>/
  real_spatial/<campaign_id>/

  external_validation/shapemix_full_evaluation_v1/
    evidence_summary.yaml
    evidence_table.tsv
    resource_table.tsv

/tmp/deconvatac-<task>-<unique_id>/
```

Existing Step 6 campaign directories remain valid historical layouts. New campaigns should place shards, baselines, summaries, and logs under one stable `<campaign_id>` directory when the runner configuration permits it.

## 3. File lifecycle

### 3.1 Download staging

Incomplete downloads, extraction directories, SRA conversion intermediates, and retry state go under:

```text
data/work/downloads/<source_id>/
```

A download may enter `data/raw/sources/` only after its expected size, checksum, compression integrity, and basic schema have been validated. A completed file should be moved atomically when possible.

### 3.2 Immutable external sources

Files obtained from a provider are stored byte-for-byte under:

```text
data/raw/sources/<provider>/<accession_or_release>/
```

This includes author-published processed matrices or fragment files when they are preserved without modification. Do not sort, recompress, rename, lift over, index, or edit a finalized source file in place. Our normalized or indexed copy belongs under `data/processed/shapemix/<family>/`.

Every source family must have a tracked manifest under `configs/data_sources/` containing accessions, URLs, source roles, observed sizes, SHA-256 hashes, genome build, assay, organism, tissue, donor metadata, download date, and validation status.

### 3.3 Restartable preprocessing work

Large intermediate files that may be needed after a job restart go under:

```text
data/work/preprocessing/<family>/<job_id>/
```

They are not benchmark inputs and may be removed after validated final products have been written. Truly disposable per-process scratch belongs under a unique `/tmp/deconvatac-*` directory.

### 3.4 Reusable ShapeMix products

Normalized fragments, tabix indices made by this project, fragment-shape caches, label harmonization tables, feature axes, donor/split assignments, and preprocessing manifests go under:

```text
data/processed/shapemix/<family>/
```

These products may be expensive to rebuild but remain derived data. Each stage must record its input hashes, parameters, code revision, software versions, random seeds, dimensions, and validation results.

Reference-only marker panels used to score orthogonal real-spatial evidence are reusable derived products, not runnable datasets or experiment results. They therefore belong under `data/processed/shapemix/real_spatial_validation/marker_features_v1/<reference_id>/`. The panel manifest must record the reference hash and selection parameters, and feature selection must never inspect spatial outcomes.

The large GSE216371 embryo-reference stream uses these exact locations:

```text
data/processed/shapemix/gse216371_embryo_reference/
  labels/author_e13_5_v1/                 # author-label audit
  labels/major_types_v1/                  # frozen broad ontology
  feature_axes/author_ccres_v1/            # audited 830,873-cCRE source axis
  feature_axes/major_types_v1/             # selected 5,000-cCRE axis
  normalized_fragments/major_types_v1/     # retained E13.5 fragments and statistics
  manifests/embryo_fragment_coordinate_audit.yaml

data/work/preprocessing/gse216371_reference/
  bin/                                     # compiled streaming helper
  shape_events/major_types_v1/             # disposable packed sparse events

data/processed/references/gse216371_mouse_embryo_e13_major_types_v1/
  reference.yaml
  atac/reference.h5ad
```

The normalized-fragment cache retains its member hashes, per-cell
fragment-total audit, and reference-only aggregate feature statistics because
they are expensive reusable derivations. The compiled executable and packed
event stream are restartable work and are not authoritative inputs; immutable
manifests record their source/binary hashes and final counters without
depending on a surviving `data/work/` path.

### 3.5 Standardized references

Reference H5ADs produced or standardized by repository scripts go under:

```text
data/processed/references/<reference_id>/
```

The `reference.yaml` file declares modality paths, label keys, source identity, construction parameters, and source hashes. Dataset configurations must point to these processed paths. A reference ID is immutable after use in a frozen experiment; a semantic change requires a new reference ID or explicit version.

### 3.6 Runnable datasets

Only validated, runner-compatible data belong under:

```text
data/processed/datasets/<dataset_id>/
```

Family caches and incomplete preprocessing outputs must not be placed here. Register a dataset in `data/registry/datasets.yaml` only after reference/spatial feature alignment, shape-layer semantics, declared cell types, truth, split disjointness, and required provenance have passed validation. Never overwrite an ID used by a completed campaign.

For real-spatial datasets, the ShapeMix model input remains under `atac/`. Copied RNA, protein, histone, anatomical-alignment, or other cross-modality objects belong under `validation/` and must be declared as orthogonal evidence rather than exact composition truth. The dataset manifest records their hashes and alignment checks.

### 3.7 Experiment results

Development, primary, sensitivity, external-validation, and real-spatial outputs are separated by scope under `results/`. Each completed run should retain:

```text
<run_id>/
  run.yaml
  inputs.yaml
  environment.txt
  output_sha256.yaml
  results/
    proportions.csv
    abundance.csv                    # when produced
    truth.csv                        # when exact truth exists
    diagnostics.json
    raw_method_output/
```

A campaign should also retain its resolved protocol, batch manifest, run table, comparison table, failures, summary tables, logs needed to interpret failures, and provenance hashes. Completed result directories are append-free and must not be silently overwritten; reruns use a new campaign ID or a clearly versioned campaign directory.

The cross-family synthesis for the frozen full evaluation belongs at `results/external_validation/shapemix_full_evaluation_v1/`. It links to, but does not copy or pool, the primary, sensitivity, external-validation, and real-spatial campaign evidence classes.

## 4. Naming rules

- `<source_id>` identifies one provider release or accession family.
- `<family>` identifies one biological/source family across reusable preprocessing stages.
- `<reference_id>` identifies one immutable standardized reference and its label universe.
- `<dataset_id>` identifies one exact runner input, including split, condition, and seed when those change its contents.
- `<campaign_id>` identifies a frozen experimental protocol and includes a version such as `_v1` when results will be cited.
- Seeds belong in dataset or manifest metadata and in IDs when they distinguish materialized datasets. Timestamps are used for temporary job IDs, not as substitutes for protocol versions.

All authoritative paths stored in configs and manifests should be project-relative. An absolute execution path may be retained as informational provenance, but it must not be the only locator.

## 5. Git tracking policy

Tracked:

- code, tests, notebooks, and documentation;
- source, dataset-template, method, and experiment configurations;
- tracked source manifests and checksums;
- the complete `results/` directory, including all campaign scopes and per-run provenance.

Ignored:

- `data/raw/`, `data/work/`, `data/processed/`, and `data/registry/`;
- environment directories, package caches, test caches, and operating-system metadata;
- disposable `/tmp` work.

Because data are ignored, a clean clone must be able to reconstruct required data using tracked manifests and scripts. Because results are exposed, adding new result artifacts can materially increase repository size; campaign outputs should therefore contain evidence needed for review and reproduction, while download caches and preprocessing scratch must stay outside `results/`.

## 6. Reference migration completed in this workspace

The standardized reference tree was moved without changing the H5AD bytes:

```text
data/raw/references/human_cardiac_niches/
  -> data/processed/references/human_cardiac_niches/

data/raw/references/pbmc_granulocyte_sorted_10k_multiome/
  -> data/processed/references/pbmc_granulocyte_sorted_10k_multiome/

data/raw/references/russell_250/
  -> data/processed/references/russell_250/
```

Live configs, preparation scripts, current recreation documentation, standardized dataset descriptors, ShapeMix split manifests, and reference manifests use the processed paths.

Historical records are not rewritten solely to conceal an old path. In particular, completed-run `inputs.yaml` files record the configuration used when those runs executed, and `docs/migration/data_files_*.tsv` records an earlier data snapshot. Those files may therefore retain `data/raw/references/...` as historical provenance; they are not active path consumers.

## 7. Step 1-6 mapping

The completed work maps onto this policy as follows:

| Step | Persistent artifact class | Location |
|---|---|---|
| 1 | Immutable 10x source files | `data/raw/sources/10x_genomics/.../cellranger_arc_2.0.0/` |
| 1 | Tracked source identity and hashes | `configs/data_sources/pbmc_granulocyte_sorted_10k_cellranger_arc_2.0.0.yaml` |
| 2 | Shape counter and data contract | `src/`, `tests/`, `configs/`, and `docs/` |
| 3 | Reusable ShapeMix caches and splits | `data/processed/shapemix/pbmc_granulocyte_sorted_10k/` |
| 3 | Twenty-one runnable ShapeMix datasets | `data/processed/datasets/*shapemix*` |
| 4 | Fixed-signature MAP implementation | `src/deconvatac/shapemix/` and tests |
| 5 | Runner smoke evidence | `results/development/shapemix_smoke_development/` |
| 6 | Primary shards, baseline runs, summaries, and controls | `results/primary/` and `results/development/shapemix_negative_controls_v1/` |

The one Step 6 resource pilot under `/private/tmp/deconvatac-step6-profile-20260822/` remains outside the primary results because it was a temporary sizing run, not benchmark evidence.
