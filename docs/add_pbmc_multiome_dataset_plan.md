# Plan to Add 10x PBMC Multiome

## Scope

This document records the PBMC Multiome addition plan and the current execution state.

Add the 10x Genomics PBMC Multiome dataset:

```text
PBMC from a Healthy Donor - Granulocytes Removed Through Cell Sorting (10k)
```

10x describes this as an Epi Multiome dataset analyzed with Cell Ranger ARC 2.0.0. The dataset page reports:

```text
estimated cells: 11,898
ATAC peaks: 143,887
ATAC median high-quality fragments per cell: 17,306
GEX median genes per cell: 1,826
GEX median UMI counts per cell: 3,778
license: CC BY 4.0
```

This is a **Case 2** dataset for our benchmark design:

```text
source cells have ATAC/RNA profiles and labels after annotation
source cells do not have real per-cell tissue coordinates
```

So it should be used to generate simulated spatial ATAC spots with synthetic regions and exact ground truth.

## Candidate Dataset ID

Use:

```text
pbmc_granulocyte_sorted_10k_multiome
```

For simulated benchmark datasets derived from it, use explicit IDs such as:

```text
pbmc_granulocyte_sorted_10k_sim_equal_celltype
pbmc_granulocyte_sorted_10k_sim_observed_abundance
```

## Execution Steps

Execute this plan one step at a time. Do not register the PBMC dataset or add it to benchmark configs until the validation steps pass.

Current execution status as of 2026-07-04:

```text
Step 1: complete
Step 2: complete
Step 3: complete
Step 4: complete with local SnapATAC2/GET RNA label file stored under data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad
Step 5: complete with SnapATAC2/GET rna.h5ad labels
Step 6: complete with regenerated feature sets
Step 7: complete for equal_celltype and observed_abundance with real labels
Step 8: complete for loader/truth/shape validation with real labels
Step 9: complete
Step 10: complete
Step 11: complete and verified for all_methods_all_atac_datasets.yaml
Step 12: optional follow-up, deferred
```

Important Step 4 update:

```text
The preferred SnapATAC2/GET label source could not be downloaded earlier because
the SnapATAC2 helper points to a Mendeley draft file URL that returned HTTP 403.
The user has now provided the missing RNA label object, which is stored at
data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad.
Use that file as the canonical PBMC label source going forward.
Do not use CellTypist or any other fallback label source without explicit user
approval.
```

Canonical label source:

```text
data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad
  obs index: cells
  obs label column: cell_type
  observed rows in file: 9,631
  labeled rows in file: 9,631
  missing cell_type rows in file: 0
  cell_type categories: 19
  manifest: data/raw/sources/snapatac2/pbmc10k_multiome/manifest.yaml
```

When converting the original 10x matrix, keep only barcodes that are present in
the raw-source `rna.h5ad` and have a non-missing `obs["cell_type"]`. Ignore/drop original 10x
cells that do not have a cell-type label. With the currently inspected files,
this means retaining 9,627 labeled 10x cells and dropping 2,271 of the original
11,898 10x cells. The `rna.h5ad` label file itself has 9,631 labeled rows; four
of those barcodes are not present in the 10x filtered matrix.

Previous fallback label outputs, now superseded for the canonical PBMC build:

```text
data/raw/sources/celltypist/pbmc_granulocyte_sorted_10k/
  manifest.yaml
  cell_type_mapping.csv
  cell_type_summary.csv
  data/models/models.json
  data/models/Immune_All_High.pkl
  data/models/Immune_All_Low.pkl
```

Reference outputs:

```text
data/raw/references/pbmc_granulocyte_sorted_10k_multiome/
  reference.yaml
  atac/reference.h5ad
  rna/reference.h5ad
```

Validation summary:

```text
input cells: 11,898
labeled rows in rna.h5ad: 9,631
labeled cells retained after 10x intersection: 9,627
labels not found in 10x matrix: 4
dropped original 10x cells missing label: 2,271
ATAC features: 143,887
RNA features: 36,601
SnapATAC2/GET cell_type labels: 19
ATAC and RNA reference cells match: yes
missing cell_type labels: 0
missing ATAC peak coordinates: 0
```

Processed simulation outputs:

```text
data/processed/datasets/pbmc_granulocyte_sorted_10k_sim_equal_celltype/
data/processed/datasets/pbmc_granulocyte_sorted_10k_sim_observed_abundance/
```

Simulation validation summary:

```text
datasets generated: 2
spots per dataset: 1,024
grid shape: 32 x 32
mean cells per spot parameter: 10
selected SnapATAC2/GET cell types: 16 labels with at least 100 retained source cells
selected cell type threshold: at least 100 source cells
truth columns: all 19 SnapATAC2/GET labels, with zero proportions for unsampled rare labels
feature lists: 20,000 ATAC highly_variable, 20,000 ATAC highly_accessible, 20,000 RNA highly_variable
loader validation: passed for ATAC and RNA highly_variable on both datasets
source-cell provenance reconstructs truth: yes, max absolute difference about 8.3e-17
```

### Step 1: Verify Download URLs

Use header-only checks such as `wget --spider` or `curl -I` to verify the 10x CDN URLs before downloading large files.

Expected result:

```text
confirmed URL list for required and recommended 10x files
decision recorded for whether fragments and cCREs are deferred
```

Proceed only if the required matrix URL exists:

```text
pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
```

### Step 2: Download Required And Recommended 10x Files

Download the first-pass files into the raw source directory.

Expected result:

```text
data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/
  pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
  pbmc_granulocyte_sorted_10k_atac_peaks.bed
  pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv
  pbmc_granulocyte_sorted_10k_web_summary.html
```

Do not download fragments or cCREs in the first pass unless we explicitly decide to build a harmonized cCRE feature universe immediately.

### Step 3: Write Raw Data Manifest And Checksums

Create a manifest for the downloaded 10x files.

Expected result:

```text
data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/manifest.yaml
```

The manifest should record:

```text
source dataset
source URL for each file
download date
file size
checksum
license
Cell Ranger ARC version
```

### Step 4: Prepare Cell-Type Mapping

Use the local SnapATAC2/GET PBMC10k Multiome RNA object now stored under
`data/raw/`:

```text
data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad
```

This is the real label source intended by the GET Foundation `prepare_pbmc.ipynb`
workflow. It contains a barcode index in `.obs["cells"]` and cell-type labels in
`.obs["cell_type"]`.

The canonical raw-source location should be:

Expected result:

```text
data/raw/sources/snapatac2/pbmc10k_multiome/
  manifest.yaml
  rna.h5ad
  cell_type_mapping.csv
  cell_type_summary.csv
```

Expected columns in `cell_type_mapping.csv`:

```text
barcode
cell_type
cell_type_source
```

Use:

```text
cell_type_source = snapatac2_pbmc10k_multiome_get_prepare_pbmc
```

Drop any barcode whose `cell_type` is missing. Also drop original 10x barcodes
that are not present in `rna.h5ad`, because they do not have labels from this
source. Record retained and dropped counts in the manifest.

Current local inspection:

```text
original 10x cells: 11,898
rna.h5ad rows: 9,631
rna.h5ad rows with non-missing cell_type: 9,631
rna.h5ad labels not found in original 10x matrix: 4
original 10x cells retained after label intersection: 9,627
original 10x cells dropped because no label is available: 2,271
cell_type categories: 19
source manifest: data/raw/sources/snapatac2/pbmc10k_multiome/manifest.yaml
```

### Step 5: Convert 10x Matrix To Labeled References

Read the 10x filtered feature-barcode matrix, keep labeled barcodes, attach cell-type labels, and split RNA and ATAC features.

Expected result:

```text
data/raw/references/pbmc_granulocyte_sorted_10k_multiome/
  reference.yaml
  atac/reference.h5ad
  rna/reference.h5ad
```

The ATAC and RNA references must share identical `.obs_names` and cell-type labels.

### Step 6: Compute Feature Sets

Compute first-pass feature sets from the Cell Ranger ARC peak matrix.

Expected result after simulation datasets are created:

```text
atac/features/highly_variable.txt
atac/features/highly_accessible.txt
rna/features/highly_variable.txt
```

Do not create a cCRE feature set in this step.

### Step 7: Simulate PBMC Spatial Datasets

Use the labeled PBMC references to create synthetic spatial datasets with exact truth.

First target datasets:

```text
pbmc_granulocyte_sorted_10k_sim_equal_celltype
pbmc_granulocyte_sorted_10k_sim_observed_abundance
```

Expected result:

```text
data/processed/datasets/<pbmc_sim_dataset_id>/
  dataset.yaml
  atac/spatial.h5ad
  rna/spatial.h5ad
  truth/proportions.csv
  simulation/source_cells_by_spot.jsonl
```

The regenerated real-label implementation uses selected labels from
`rna.h5ad`. With the current `min_source_cells >= 100` rule, the selected labels
are:

```text
CD14 Mono
CD4 Naive
CD8 Naive
CD4 TCM
CD16 Mono
NK
CD8 TEM_1
CD8 TEM_2
Intermediate B
Memory B
CD4 TEM
cDC
Treg
gdT
MAIT
Naive B
```

Labels with fewer than 100 source cells are kept in the reference and truth
columns but are not sampled into first-pass synthetic spots unless the threshold
is lowered intentionally.

Sampling designs:

```text
equal_celltype: selected real labels are sampled with equal probability
observed_abundance: selected real labels are sampled according to observed abundance after renormalizing within the selected set
```

### Step 8: Validate Dataset Contract

Validate the processed dataset before registry changes.

Required checks:

```text
ATAC and RNA references have identical cells
cell_type exists and has no missing values
retained and dropped cell counts are recorded
ATAC features have genomic coordinates
truth/proportions.csv rows sum to 1
simulation/source_cells_by_spot.jsonl reconstructs truth exactly
load_deconvolution_input() works for ATAC highly_variable
```

Do not continue if any of these checks fail.

### Step 9: Register Dataset

Add the validated simulated dataset to:

```text
data/registry/datasets.yaml
```

Expected entry:

```yaml
pbmc_granulocyte_sorted_10k_sim_equal_celltype:
  config: data/processed/datasets/pbmc_granulocyte_sorted_10k_sim_equal_celltype/dataset.yaml
pbmc_granulocyte_sorted_10k_sim_observed_abundance:
  config: data/processed/datasets/pbmc_granulocyte_sorted_10k_sim_observed_abundance/dataset.yaml
```

### Step 10: Run Smoke Benchmark

Run a small smoke benchmark before adding PBMC to full experiment configs.

Minimum methods:

```text
nnls
tangram smoke
```

Expected result:

```text
results/<pbmc_smoke_run>/
  runs.csv
  comparison.csv
```

Execution result as of 2026-07-04:

```text
temporary config: /private/tmp/deconvatac_pbmc_step10_smoke.yaml
temporary output: /private/tmp/deconvatac_pbmc_step10_smoke_results/pbmc_step10_smoke/
cleanup: temporary config and output directory removed after inspection
jobs: 4
failures: 0
```

Smoke metrics:

```text
pbmc_granulocyte_sorted_10k_sim_equal_celltype / nnls:
  rmse = 0.03868826852068771
  jsd  = 0.2952651247505604

pbmc_granulocyte_sorted_10k_sim_equal_celltype / tangram smoke:
  rmse = 0.09622688718319712
  jsd  = 0.6646955279972313

pbmc_granulocyte_sorted_10k_sim_observed_abundance / nnls:
  rmse = 0.038031266336587
  jsd  = 0.27666656300328507

pbmc_granulocyte_sorted_10k_sim_observed_abundance / tangram smoke:
  rmse = 0.06866858438156663
  jsd  = 0.5052711575786042
```

### Step 11: Add PBMC To Experiment Configs

Add validated PBMC simulation datasets to benchmark configs such as:

```text
configs/experiments/all_methods_all_atac_datasets.yaml
```

Execution result as of 2026-07-04:

```text
config: configs/experiments/all_methods_all_atac_datasets.yaml
PBMC datasets included:
  pbmc_granulocyte_sorted_10k_sim_equal_celltype
  pbmc_granulocyte_sorted_10k_sim_observed_abundance
validated job expansion:
  total jobs: 42
  datasets: 7
  methods: 6
  PBMC jobs: 12
```

### Step 12: Optional cCRE/Fragments Follow-Up

Only after the Cell Ranger peak-based PBMC benchmark works end to end, decide whether to download:

```text
pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
cCRE_hg38.tsv.gz
```

Use these only if we need a harmonized cCRE feature universe or fragment-level re-counting.

## Download URLs to Verify

The 10x dataset page confirms the dataset, assay, Cell Ranger ARC version, web summary link, and key metrics. The following CDN URLs follow the 10x filename pattern for that dataset and should be verified with `wget --spider`, `curl -I`, or an equivalent header-only check before downloading:

```text
https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peaks.bed
```

Also verify whether these additional 10x output files exist:

```text
https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_web_summary.html
https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv
```

The proposed cCRE URL is:

```text
http://catlas.org/catlas_downloads/humantissues/cCRE_hg38.tsv.gz
```

I could not verify that cCRE URL from the browser. Treat it as optional until checked with a download-header or `wget --spider` step.

Verification status:

```text
10x dataset page: confirmed
10x web summary link on dataset page: confirmed
filtered_feature_bc_matrix.h5 URL: confirmed HTTP 200, downloaded
atac_peaks.bed URL: confirmed HTTP 200, downloaded
per_barcode_metrics.csv URL: confirmed HTTP 200, downloaded
web_summary.html URL: confirmed HTTP 200, downloaded
atac_fragments.tsv.gz URL: confirmed HTTP 200, deferred
atac_fragments.tsv.gz.tbi URL: confirmed HTTP 200, deferred
cCRE_hg38.tsv.gz URL: HTTP 404, deferred
```

## What Is Exactly Needed

### Required For First Pass

```text
pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
```

This is the main required file. It should contain both:

```text
Gene Expression
Peaks
```

Cell Ranger ARC documentation says the feature-barcode matrix contains both gene-expression features and peak features. For peaks, matrix values are cut-site counts for a feature/barcode pair.

This file is enough to build:

```text
RNA reference matrix
ATAC reference matrix
shared cell/barcode IDs
```

### Required Annotation Step

The 10x output matrix does **not** by itself solve the cell-type label requirement.

We still need:

```text
cell barcode -> cell type label
```

Use the GET Foundation `prepare_pbmc.ipynb` label-mapping workflow with the local
RNA label object now available in this repo.

That notebook loads the 10x filtered feature-barcode matrix, then loads a
preprocessed SnapATAC2 PBMC10k Multiome RNA object:

```python
rna = snap.read(snap.datasets.pbmc10k_multiome(modality="RNA"), backed=None)
```

The raw-source file `data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad` is
that RNA label object. It already contains `.obs["cell_type"]`. The conversion
now should:

1. load the 10x filtered feature-barcode matrix
2. load or inspect `data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad`
3. keep only 10x barcodes present in that label file
4. drop any barcode with missing `obs["cell_type"]`
5. create a barcode-to-cell-type dictionary from the retained RNA object
6. map those labels back onto the 10x matrix by barcode
7. split the labeled 10x matrix into RNA and ATAC views

Conceptually:

```python
ad = sc.read_10x_h5(
    "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
    gex_only=False,
)
rna = sc.read_h5ad("data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad")

rna = rna[~rna.obs["cell_type"].isna()].copy()
ad = ad[ad.obs_names.isin(rna.obs_names)].copy()
barcode_to_celltype = rna.obs["cell_type"].astype(str).to_dict()
ad.obs["cell_type"] = ad.obs_names.map(barcode_to_celltype)
ad = ad[~ad.obs["cell_type"].isna()].copy()
```

The local file has 9,631 cells and 19 `cell_type` categories. It currently has
no missing labels within the file, but the original 10x matrix has 11,898 cells.
Four labels in `rna.h5ad` are not present in the 10x filtered matrix. Therefore,
the conversion retains 9,627 labeled 10x cells and drops 2,271 original 10x
cells that do not have usable labels from this source.

Implementation note: the current environment's installed `anndata` may fail to
read this file because of a newer/null field under `uns/log1p/base`. If that
happens, either read the label columns directly with `h5py` from
`/obs/cells` and `/obs/cell_type`, or clean that `uns` field with a compatible
AnnData version before using `scanpy.read_h5ad`.

Fallback approach only if the canonical raw-source `rna.h5ad` cannot be used,
and only after explicit user approval:

1. Ask the user for approval to use a fallback label source.
2. Use the paired RNA modality from the same filtered feature-barcode matrix.
3. Annotate cells using a PBMC reference method, for example CellTypist, Azimuth/Seurat, or a curated PBMC marker workflow.
4. Store both broad and fine labels in `.obs`.

Suggested `.obs` columns:

```text
cell_type
cell_type_broad
cell_type_fine
cell_type_source
cell_type_confidence
```

For the SnapATAC2-derived labels, set:

```text
cell_type_source = snapatac2_pbmc10k_multiome_get_prepare_pbmc
cell_type_confidence = not_reported
```

### Recommended Provenance/QC Files

```text
pbmc_granulocyte_sorted_10k_atac_peaks.bed
pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv
pbmc_granulocyte_sorted_10k_web_summary.html
```

Why:

- `atac_peaks.bed`: explicit peak intervals and peak provenance.
- `per_barcode_metrics.csv`: barcode-level ATAC/GEX QC metrics and cell-calling status.
- `web_summary.html`: human-readable Cell Ranger ARC QC report.

Cell Ranger ARC documentation says `per_barcode_metrics.csv` contains paired ATAC/GEX barcode sequences, QC metrics, and cell-associated partition status.

### Optional Files

```text
pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
cCRE_hg38.tsv.gz
```

The fragments file and tabix index are useful if we want to:

- re-call peaks
- re-count fragments into another feature universe
- intersect with cCREs
- compute additional ATAC QC
- build a harmonized cCRE-based benchmark

They are **not required** for the first pass if we accept the Cell Ranger ARC peak matrix.

The cCRE file is also **not required** for the first pass. It is only needed if we decide to move from Cell Ranger peaks to a common cCRE feature universe.

## Proposed Raw Layout

Store immutable downloaded files under:

```text
data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/
  manifest.yaml
  pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
  pbmc_granulocyte_sorted_10k_atac_peaks.bed
  pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv
  pbmc_granulocyte_sorted_10k_web_summary.html
  pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz        # optional
  pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi    # optional
```

Canonical SnapATAC2/GET label source:

```text
data/raw/sources/snapatac2/pbmc10k_multiome/
  manifest.yaml
  rna.h5ad
  cell_type_mapping.csv
  cell_type_summary.csv
```

The user-provided source file has been moved into the canonical raw-source
directory. The mapping should be built from `.obs["cells"]` and
`.obs["cell_type"]`.

Previous CellTypist fallback output, retained only for provenance until the
PBMC files are regenerated:

```text
data/raw/sources/celltypist/pbmc_granulocyte_sorted_10k/
  manifest.yaml
  cell_type_mapping.csv
  cell_type_summary.csv
  data/models/models.json
  data/models/Immune_All_High.pkl
  data/models/Immune_All_Low.pkl
```

Expected columns:

```text
barcode
cell_type
cell_type_source
```

If cCREs are used later:

```text
data/raw/feature_universes/hg38/catlas_cCRE_hg38/
  cCRE_hg38.tsv.gz
  source.yaml
```

## Proposed Reference Layout

After converting and annotating the 10x matrix, write shared labeled references under:

```text
data/raw/references/pbmc_granulocyte_sorted_10k_multiome/
  reference.yaml
  atac/reference.h5ad
  rna/reference.h5ad
```

The ATAC and RNA references should share the same cells and `.obs` labels.

Minimum `.obs` columns:

```text
cell_type
cell_type_source
donor_id
organism
tissue
assay
source_dataset_id
```

Optional `.obs` columns if a separate broad/fine mapping is added later:

```text
cell_type_broad
cell_type_fine
```

Minimum `.var` columns for ATAC:

```text
feature_type = Peaks
chrom
start
end
```

Minimum `.var` columns for RNA:

```text
feature_type = Gene Expression
gene_id
gene_name
```

## Conversion Plan

1. Download and checksum the required 10x files.
2. Export the barcode-to-cell-type mapping from the local `rna.h5ad` file.
3. Read the filtered feature-barcode matrix with `scanpy.read_10x_h5(..., gex_only=False)`.
4. Keep only barcodes present in `rna.h5ad` with non-missing `cell_type`; ignore/drop cells with no label.
5. Attach label columns to `.obs`.
6. Split features:

```text
var["feature_types"] == "Peaks"            -> ATAC AnnData
var["feature_types"] == "Gene Expression"  -> RNA AnnData
```

7. Normalize or preserve matrices carefully:
   - store raw counts in `.X`
   - optionally add normalized layers later
8. Attach QC metrics from `per_barcode_metrics.csv` if available.
9. Copy final label columns to both ATAC and RNA references.
10. Write `reference.h5ad` files under `data/raw/references/`.

The next canonical PBMC rebuild should record:

```text
label_source = data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad
cell_type_source = snapatac2_pbmc10k_multiome_get_prepare_pbmc
input_10x_cells = 11898
retained_labeled_cells = 9631
dropped_unlabeled_or_unmatched_cells = 2267
```

## Simulation Plan

Because there are no real per-cell tissue coordinates, simulate spots using the Heart-style random sampling algorithm:

1. Choose target cell-type groups.
2. Define artificial spatial regions.
3. For each synthetic spot, sample source cells from the allowed cell types.
4. Sum ATAC profiles to create `atac/spatial.h5ad`.
5. Optionally sum RNA profiles to create `rna/spatial.h5ad`.
6. Compute exact truth from sampled source-cell labels.
7. Store source-cell provenance.

Expected processed dataset layout:

```text
data/processed/datasets/<pbmc_sim_dataset_id>/
  dataset.yaml
  atac/
    spatial.h5ad
    features/
      highly_variable.txt
      highly_accessible.txt
  rna/
    spatial.h5ad
    features/
      highly_variable.txt
  truth/
    proportions.csv
  simulation/
    source_cells_by_spot.jsonl
```

For PBMC, useful first simulations:

```text
equal_celltype: selected real PBMC cell types sampled with equal probability
observed_abundance: selected real PBMC cell types sampled by observed reference abundance
immune_regions: future synthetic regions enriched for T/NK, B/plasma, myeloid, granulocyte-like cells
rare_cell_stress: future simulation including rare labels to test whether methods recover low-abundance types
```

The retained reference label categories and observed counts after intersecting
`rna.h5ad` labels with the 10x filtered matrix are:

```text
CD14 Mono         2551
CD4 Naive         1382
CD8 Naive         1353
CD4 TCM           1113
CD16 Mono          442
NK                 403
CD8 TEM_1          322
CD8 TEM_2          315
Intermediate B     300
Memory B           298
CD4 TEM            286
cDC                180
Treg               157
gdT                143
MAIT               130
Naive B            125
pDC                 98
HSPC                17
Plasma              12
```

If the existing minimum source-cell threshold of 100 is kept, `pDC`, `HSPC`,
and `Plasma` will not be sampled in the first two simulations. They can remain
as zero-valued truth columns, or a future rare-cell stress simulation can lower
the threshold intentionally.

## Feature-Set Plan

For the first pass, keep Cell Ranger ARC peaks as the ATAC feature universe.

This matches the practical route in `prepare_pbmc.ipynb`: the notebook notes that a union peak/cCRE set can reduce domain shift for zero-shot work, but keeps the Cell Ranger peak set for convenience. For this repository, the Cell Ranger peak set is the right first pass because it avoids requiring fragment re-counting before the benchmark is usable.

Compute:

```text
highly_variable.txt
highly_accessible.txt
```

Do not switch to cCREs until the Cell Ranger peak workflow is working end to end.

Later, if cCRE harmonization is needed:

1. Verify and download `cCRE_hg38.tsv.gz`.
2. Intersect fragments or peaks with cCRE intervals.
3. Create a second processed dataset or feature set named explicitly, for example:

```text
cCRE_hg38
```

## Registry Plan

Register simulated datasets only after references and processed simulation files exist.

Add entries to:

```text
data/registry/datasets.yaml
```

Example:

```yaml
pbmc_granulocyte_sorted_10k_sim_equal_celltype:
  config: data/processed/datasets/pbmc_granulocyte_sorted_10k_sim_equal_celltype/dataset.yaml
pbmc_granulocyte_sorted_10k_sim_observed_abundance:
  config: data/processed/datasets/pbmc_granulocyte_sorted_10k_sim_observed_abundance/dataset.yaml
```

Do not register the raw 10x download directly unless it has been converted to the unified format.

## Validation Plan

Before registering:

1. Confirm ATAC and RNA references have identical `.obs_names`.
2. Confirm `obs["cell_type"]` exists and has no missing values for retained cells.
3. Confirm the number of retained labeled cells is recorded. With the current local `rna.h5ad` and 10x matrix, this is 9,627.
4. Confirm dropped cells are explained. With the current local `rna.h5ad` and 10x matrix, this is 2,271 original 10x cells without usable labels; four `rna.h5ad` labels are not present in the 10x matrix.
5. Confirm ATAC features have genomic coordinates.
6. Confirm simulated `truth/proportions.csv` rows sum to 1.
7. Confirm `simulation/source_cells_by_spot.jsonl` reconstructs truth exactly.
8. Confirm `load_deconvolution_input()` works for:

```text
dataset_id = pbmc_granulocyte_sorted_10k_sim_equal_celltype
modality = atac
feature_set = highly_variable
```

Current validated dataset IDs:

```text
pbmc_granulocyte_sorted_10k_sim_equal_celltype
pbmc_granulocyte_sorted_10k_sim_observed_abundance
```

9. Run smoke methods:

```text
nnls
tangram smoke
```

10. Add the dataset to an experiment config after loader validation and, ideally, smoke validation.

## Download Script Plan

Create a script later:

```text
scripts/download_pbmc_multiome_10x.sh
```

The script should:

1. create the raw source directory
2. download only the selected required/recommended files
3. write checksums
4. write `manifest.yaml`
5. not overwrite existing files unless `--overwrite` is passed

Use `wget --spider` or equivalent header checks first to verify exact URLs before real download.

The label-preparation script has been implemented:

```text
scripts/prepare_pbmc_multiome_labels.py
```

The script now implements the SnapATAC2/GET `rna.h5ad` mapping path by default
and accepts:

```text
--rna-h5ad data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad
--output-dir data/raw/sources/snapatac2/pbmc10k_multiome
```

Implemented behavior:

1. read barcode IDs from `.obs["cells"]` or `.obs_names`
2. read labels from `.obs["cell_type"]`
3. drop rows with missing `cell_type`
4. write `data/raw/sources/snapatac2/pbmc10k_multiome/cell_type_mapping.csv`
5. write `cell_type_summary.csv`
6. write `manifest.yaml` with retained and dropped cell counts

Keep the CellTypist path as a fallback mode only, and require explicit user
approval before using it. Keep this separate from
`scripts/download_pbmc_multiome_10x.sh` because label preparation depends on the
SnapATAC2/GET RNA object rather than only the 10x CDN files.

The current fallback script enforces this policy by requiring:

```text
--allow-celltypist-fallback
```

That flag should be used only after explicit user approval.

Implemented helper scripts:

```text
scripts/prepare_pbmc_multiome_labels.py
scripts/prepare_pbmc_multiome_references.py
scripts/regenerate_pbmc_simulations.py
```

## Open Questions

- Should a future rare-cell PBMC simulation lower the `min_source_cells >= 100` threshold to include `pDC` with 98 cells?
- Do we need fragments immediately, or can we defer them until cCRE harmonization?
- The CATlas cCRE URL checked here returned HTTP 404. Is there a current CATlas cCRE URL, or should a different cCRE source be used?

## Recommendation

The first implementation pass downloaded:

```text
filtered_feature_bc_matrix.h5
atac_peaks.bed
per_barcode_metrics.csv
web_summary.html
```

Then initially converted, annotated, simulated, and validated the two first PBMC
datasets using temporary CellTypist fallback labels.

The canonical rebuild has now been completed from the intended SnapATAC2/GET RNA
label object at `data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad`. The
reference rebuild keeps only barcodes with non-missing `obs["cell_type"]` that
also appear in the 10x filtered matrix, retaining 9,627 labeled cells and
dropping 2,271 original 10x cells without usable labels. Do not use fallback
labels without explicit user approval.

Defer:

```text
fragments.tsv.gz
fragments.tsv.gz.tbi
cCRE_hg38.tsv.gz
```

until the Cell Ranger peak-based benchmark is working end to end.

## Sources Checked

- 10x Genomics dataset page for PBMC granulocyte-sorted 10k Multiome:
  `https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0`
- 10x Genomics Cell Ranger ARC Feature-Barcode Matrices documentation:
  `https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/feature-barcode-matrices`
- 10x Genomics Cell Ranger ARC ATAC Fragments documentation:
  `https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/fragments-file`
- 10x Genomics Cell Ranger ARC ATAC Peaks documentation:
  `https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/atac-peak-file`
- 10x Genomics Cell Ranger ARC Per Barcode QC Metrics documentation:
  `https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/per-barcode-qc-metrics`
- 10x Genomics Cell Ranger ARC Cell Type Annotation Outputs documentation:
  `https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/cr-arc-cell-annotation-outputs`
- GET Foundation `prepare_pbmc.ipynb`, especially the SnapATAC2-derived barcode-to-cell-type mapping:
  `https://github.com/GET-Foundation/get_model/blob/master/tutorials/prepare_pbmc.ipynb`
