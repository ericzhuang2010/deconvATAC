# Dataset Search Requirements

## Goal

Find additional datasets that can support spatial ATAC deconvolution method development and benchmarking.

The strongest benchmark datasets are not necessarily real spatial ATAC datasets. The strongest benchmark inputs are annotated single-cell or single-nucleus ATAC or multiome datasets that let us simulate spatial spots with known source cells.

Ground truth requires knowing which source cells created each simulated spot.

## Required Information to Record

For every candidate dataset, record:

```text
dataset_name
source_url_or_accession
publication_or_preprint
organism
tissue
assay_type
file_format
genome_build
number_of_cells
number_of_features
cell_type_label_column
whether_cell_coordinates_exist
whether_raw_counts_exist
whether_peak_coordinates_exist
download_size
license_or_access_restrictions
notes
```

## Case 1: Source Cells Have Coordinates

This is the Russell-like case.

Minimum requirements:

- Single-cell or single-nucleus ATAC, or multiome with ATAC.
- Cell-by-peak count matrix.
- Cell-type labels for each source cell.
- Per-cell spatial coordinates.
- Stable cell IDs connecting the matrix, labels, and coordinates.
- Raw or count-like accessibility values, not only normalized embeddings.
- Peak identifiers or genomic intervals.

Useful optional fields:

- Paired RNA modality.
- Fine-grained and broad cell-type labels.
- Tissue section or sample ID.
- Donor or condition metadata.
- UMAP or clustering metadata.

How it will be used:

1. Bin cells by their real coordinates.
2. Sum cells in each bin to create synthetic spatial ATAC spots.
3. Compute ground-truth cell-type proportions from the cells in each bin.
4. Store source-cell provenance in:

This case preserves real source-cell geometry and gives exact benchmark truth.

## Case 2: Source Cells Do Not Have Coordinates

This is the Heart-like case.

Minimum requirements:

- Single-cell or single-nucleus ATAC, or multiome with ATAC.
- Cell-by-peak count matrix.
- Cell-type labels for each source cell.
- Stable cell IDs connecting the matrix and labels.
- Raw or count-like accessibility values.
- Peak identifiers or genomic intervals.
- Enough cells per cell type to support repeated sampling.

Useful optional fields:

- Paired RNA modality.
- Broad and fine cell-type annotations.
- Donor, tissue region, disease, or condition metadata.
- Published marker annotations or expected tissue composition.

How it will be used:

1. Create an artificial spatial grid.
2. Define artificial spatial regions.
3. Assign allowed cell types or mixture preferences to each region.
4. Randomly sample source cells from the allowed cell types.
5. Sum sampled cells to create synthetic spatial ATAC spots.
6. Compute ground-truth cell-type proportions from the sampled cell labels.
7. Store source-cell provenance in:


This case gives exact benchmark truth, but the spatial geometry is synthetic.

## Real Spatial ATAC Datasets

Real spatial ATAC datasets are useful, but they usually do not provide exact cell-type proportion truth.

They are still useful for:

- running methods on real spatial ATAC inputs
- visualizing predicted cell-type maps
- checking predictions against known tissue regions
- comparing predictions to marker peaks or external annotations

They are not sufficient for quantitative RMSE or JSD benchmarking unless they include an independent ground-truth source.

Minimum requirements for real spatial ATAC use:

- Spatial spot-by-peak count matrix.
- Spot coordinates.
- Peak identifiers or genomic intervals.
- Matching or compatible single-cell ATAC reference with cell-type labels.
- Genome build and peak naming compatible with the reference, or enough information to harmonize peaks.

## Strong Candidate Criteria

Prioritize datasets that have:

- public download access
- clear license or reuse terms
- raw count matrices
- reliable cell-type labels
- enough cells per major cell type
- documented genome build
- peak coordinates as genomic intervals
- h5ad, h5mu, h5, mtx, fragment, or other standard files
- publication or dataset documentation

## Weak Candidate Criteria

Avoid or deprioritize datasets that only provide:

- processed embeddings without count matrices
- clusters without biological cell-type labels
- inaccessible or restricted files
- unclear genome build
- no peak identifiers
- very small numbers of cells per cell type
- only bulk ATAC
- only spatial data without a compatible labeled reference

