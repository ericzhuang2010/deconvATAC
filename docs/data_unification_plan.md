# Data Unification Plan

## Short Answer

For Heart: RNA exists and is registered, but it is not needed for the main spatial ATAC method-development workflow. Keep RNA as an optional modality for paper reproduction and cross-modality comparison, but make ATAC the default path.

The regenerated Heart simulation dataset IDs should use explicit, self-describing names:

```text
human_cardiac_niches_sim_1zone_3ct_low_density
human_cardiac_niches_sim_1zone_10ct
human_cardiac_niches_sim_4zone_stripes
human_cardiac_niches_sim_4zone_circles
```

The shared single-cell reference atlas should remain separate:

```text
human_cardiac_niches
```

## Current End State

This plan has been executed in the current workspace.

The registered method-ready datasets are:

```text
data/processed/datasets/russell_250/
data/processed/datasets/human_cardiac_niches_sim_1zone_3ct_low_density/
data/processed/datasets/human_cardiac_niches_sim_1zone_10ct/
data/processed/datasets/human_cardiac_niches_sim_4zone_stripes/
data/processed/datasets/human_cardiac_niches_sim_4zone_circles/
```

The Heart simulation datasets now use the canonical per-modality layout:

```text
data/processed/datasets/<heart_dataset_id>/
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

The Russell dataset also uses a file-based canonical truth table:

```text
data/raw/references/russell_250/
  reference.yaml
  atac/reference.h5ad

data/processed/datasets/russell_250/
  dataset.yaml
  atac/
    spatial.h5ad
  truth/
    proportions.csv
```

The shared Human cardiac niches references are stored once under:

```text
data/raw/references/human_cardiac_niches/
  reference.yaml
  atac/reference.h5ad
  rna/reference.h5ad
```

The old `heart_1` through `heart_4` processed dataset directories were removed after validation. The old shared `data/processed/feature_sets/` directory was also removed because feature lists now live inside each dataset's modality directory.

The central registry is `data/registry/datasets.yaml`. It is the source of truth for runnable datasets. The old short Heart names are no longer registered, but `scripts/regenerate_heart_simulations.py` still accepts aliases such as `heart_1`, `heart_2`, `heart_3`, and `heart_4` and maps them to the long dataset IDs.

Validation completed after migration:

- all registered dataset paths exist
- all four Heart ATAC datasets load through the unified loader
- one Heart RNA dataset loads through the unified loader
- Russell truth loads from `data/processed/datasets/russell_250/truth/proportions.csv`
- unit tests pass
- a Heart ATAC `nnls` smoke run completed and was evaluated

## Recommended Unified Format

Use one canonical processed dataset layout:

```text
data/processed/datasets/<dataset_id>/
  dataset.yaml
  <modality>/
    spatial.h5ad
    features/
      highly_variable.txt
      highly_accessible.txt   # ATAC only, when available
  truth/
    proportions.csv
  simulation/
    source_cells_by_spot.jsonl       # simulations only
```

References should stay in `data/raw/` if they are shared and large. Use a separate shared-reference namespace:

```text
data/raw/references/<reference_id>/
  reference.yaml
  <modality>/
    reference.h5ad
```

For Heart, the target layout should be:

```text
data/raw/references/human_cardiac_niches/
  reference.yaml
  atac/
    reference.h5ad
  rna/
    reference.h5ad
```

The pre-migration files mapped to that target as:

```text
data/raw/human_cardiac_niches/Adult_Peaks.h5ad -> data/raw/references/human_cardiac_niches/atac/reference.h5ad
data/raw/human_cardiac_niches/Global_raw.h5ad -> data/raw/references/human_cardiac_niches/rna/reference.h5ad
```

The `dataset.yaml` should point to these raw references instead of copying them into every Heart dataset.

## Data Concepts

### Reference Data

Reference data is the labeled single-cell atlas used to define what each cell type looks like. It is stored as cells by features and must have cell-type labels, usually in `obs["cell_type"]`.

Examples:

```text
data/raw/references/human_cardiac_niches/atac/reference.h5ad
data/raw/references/human_cardiac_niches/rna/reference.h5ad
data/raw/references/russell_250/atac/reference.h5ad
```

Keep labeled references under `data/raw/references/<reference_id>/` and point to them from each dataset config.

### Spatial Data

Spatial data is the target dataset to deconvolve. Each row is a spatial spot, not a known single cell. The deconvolution method estimates the cell-type proportions for each spatial spot using the reference data.

Examples:

```text
data/processed/datasets/russell_250/atac/spatial.h5ad
data/processed/datasets/human_cardiac_niches_sim_4zone_circles/atac/spatial.h5ad
```

### Modality

`<modality>` means the measurement type. The current maintained modalities are:

```text
atac
rna
```

For this project, `atac` is the primary modality. `rna` is optional and mainly retained for paper reproduction and cross-modality comparison.

### Truth

`truth/proportions.csv` stores the ground-truth cell-type proportions for each simulated spatial spot.

For Heart simulations, ATAC and RNA share the same truth because both modalities are generated from the same source cells per spot. The ATAC and RNA matrices are different measurements, but the underlying cell-type mixture is identical.

For Russell, truth is now exported to `data/processed/datasets/russell_250/truth/proportions.csv` and `dataset.yaml` points to that file. The original `spatial.h5ad` may still contain `obsm["proportions"]` as a convenience copy, but the CSV is the canonical truth source.

### Feature Sets

Feature-set files list which features should be used for a run.

For ATAC:

- `highly_variable.txt`: peaks that vary most across cell types/clusters; usually the preferred deconvolution feature set.
- `highly_accessible.txt`: peaks with high overall accessibility/count signal; useful as a benchmark condition but not necessarily cell-type-specific.

For RNA:

- `highly_variable.txt`: genes that vary most across cell types/clusters.

## Naming Convention

Use lowercase `snake_case` for IDs and directory names. Do not include method names, feature-set names, or run IDs in dataset names.

### Dataset IDs

Dataset IDs identify the spatial dataset to deconvolve:

```text
russell_250
human_cardiac_niches_sim_1zone_3ct_low_density
human_cardiac_niches_sim_1zone_10ct
human_cardiac_niches_sim_4zone_stripes
human_cardiac_niches_sim_4zone_circles
```

Rules:

- Use `dataset_id` for spatial datasets under `data/processed/datasets/<dataset_id>/`.
- Keep dataset IDs stable after registration.
- Do not encode modality in the dataset ID. Use modality subdirectories instead.
- Do not encode feature set in the dataset ID. Use `features/<feature_set>.txt` instead.

Heart rename mapping:

```text
heart_1 -> human_cardiac_niches_sim_1zone_3ct_low_density
heart_2 -> human_cardiac_niches_sim_1zone_10ct
heart_3 -> human_cardiac_niches_sim_4zone_stripes
heart_4 -> human_cardiac_niches_sim_4zone_circles
```

### Reference IDs

Reference IDs identify shared labeled single-cell references:

```text
human_cardiac_niches
russell_250
```

Rules:

- Use `reference_id` for shared references under `data/raw/references/<reference_id>/`.
- Use the same `reference_id` for ATAC and RNA when they come from the same biological atlas.
- Put modality-specific reference files under `<modality>/reference.h5ad`.
- Preserve source file names and download metadata in `reference.yaml`, not in the canonical filename.

### Modality Names

Use only these modality directory names unless the loader is extended:

```text
atac
rna
```

### Feature Set Names

Use concise biological/statistical names:

```text
all
highly_variable
highly_accessible
```

Feature-list files should be named:

```text
features/highly_variable.txt
features/highly_accessible.txt
```

### Standard File Names

Use the same file names across datasets:

```text
dataset.yaml
<modality>/spatial.h5ad
<modality>/features/<feature_set>.txt
truth/proportions.csv
simulation/source_cells_by_spot.jsonl
```

For references:

```text
data/raw/references/<reference_id>/<modality>/reference.h5ad
data/raw/references/<reference_id>/reference.yaml
```

### Run IDs

For standardized outputs, use a separator that makes fields unambiguous:

```text
<dataset_id>__<modality>__<feature_set>__<method>__<tag>
```

Examples:

```text
russell_250__atac__highly_variable__nnls__smoke
human_cardiac_niches_sim_4zone_circles__atac__highly_variable__rctd
human_cardiac_niches_sim_4zone_circles__rna__highly_variable__cell2location
```

## Plan

1. Document the canonical processed-data contract in `docs/` and `data/processed/README.md`.

2. Convert Heart spatial files from `.h5mu` to per-modality `.h5ad`:

```text
human_cardiac_niches_sim_4zone_circles/
  atac/spatial.h5ad
  rna/spatial.h5ad
```

3. Move feature lists into each dataset directory so `data/processed/feature_sets/` can go away:

```text
human_cardiac_niches_sim_4zone_circles/atac/features/highly_variable.txt
human_cardiac_niches_sim_4zone_circles/atac/features/highly_accessible.txt
human_cardiac_niches_sim_4zone_circles/rna/features/highly_variable.txt
```

4. Rename regenerated Heart dataset IDs and directories:

```text
heart_1 -> human_cardiac_niches_sim_1zone_3ct_low_density
heart_2 -> human_cardiac_niches_sim_1zone_10ct
heart_3 -> human_cardiac_niches_sim_4zone_stripes
heart_4 -> human_cardiac_niches_sim_4zone_circles
```

5. Update all Heart `dataset.yaml` files to use the same style as Russell: explicit `reference`, `spatial`, `truth`, and `feature_sets` paths.

6. Update `scripts/regenerate_heart_simulations.py` so regenerated Heart datasets write this layout directly.

7. Validate:

- load all renamed Heart simulation datasets for ATAC
- load RNA only for reproduction workflows
- run tests
- run one small `nnls` smoke run

8. After validation, remove old Heart `spatial.h5mu` files and remove `data/processed/feature_sets/`.

## Migration Status

Completed in this workspace:

- Shared Human cardiac niches references were moved to `data/raw/references/human_cardiac_niches/`.
- Regenerated Heart simulations now use explicit long dataset IDs.
- Heart spatial data were split from `.h5mu` into per-modality `.h5ad` files.
- Heart feature-list files were moved into each dataset's `{modality}/features/` directory.
- The central registry now points to the canonical renamed datasets.
- The old `data/processed/feature_sets/` directory and old `heart_1` through `heart_4` processed dataset directories were removed after validation.
- Registered dataset paths, unit tests, and a Heart ATAC NNLS smoke run were validated.

## RNA Decision

Keep RNA, but mark it as optional/non-primary.

- Use ATAC for new spatial ATAC method development.
- Use ATAC for default `cell2location` / `rctd` comparison.
- Keep RNA because `configs/experiments/reproduce_paper.yaml` includes RNA and the original benchmark tables include RNA results.
