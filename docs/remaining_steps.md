# Remaining Steps

The original restructuring plan now includes some completed historical checklist items. This file lists the current practical remaining work.

## Highest Priority

1. Register missing benchmark datasets.
   - `brain_1`, `brain_2`, `brain_3`, `brain_4`
   - This is blocked until their input `.h5mu` files are present locally or downloaded/restored.

## Medium Priority

2. Full paper reproduction.
   - Run all registered datasets and methods through the unified workflow.
   - Generate standardized comparison tables and figures from `results/`.

## Completed

- Core folder structure.
- Data/result contracts.
- Registry/loader.
- Cell2location and RCTD adapters.
- Tangram, DestVI, and SpatialDWLS adapters.
- NNLS baseline method.
- Metrics.
- Unified runner/evaluator.
- Standardized notebooks under `notebooks/reproduction/`, `notebooks/method_development/`, and `notebooks/exploratory/`.
- Archived legacy result copies under `data/archive/legacy_results/`.
- Russell feature annotation preprocessing.
- Russell `highly_variable` smoke comparison through Cell2location and RCTD.
  - Output: `results/comparisons/russell_250_atac_highly_variable_smoke.csv`
  - Note: Cell2location used `configs/methods/cell2location_smoke.yaml` with 5 epochs per training stage, so this verifies the pipeline but is not a final benchmark-quality Cell2location result.
- Project-local R dependencies for RCTD under `.r-lib`, including `spacexr`, `S4Vectors`, and `SingleCellExperiment`.
- Heart shared feature-list preprocessing.
  - Feature lists are now stored under each canonical Heart dataset's `{modality}/features/` directory.
  - Heart configs reuse shared raw references without duplicating the large reference `.h5ad` files.
- Heart simulation regeneration.
  - Script: `scripts/regenerate_heart_simulations.py`
  - Registry entries: `human_cardiac_niches_sim_1zone_3ct_low_density`, `human_cardiac_niches_sim_1zone_10ct`, `human_cardiac_niches_sim_4zone_stripes`, `human_cardiac_niches_sim_4zone_circles`
  - Outputs: per-modality `spatial.h5ad` files under `data/processed/datasets/<dataset_id>/{atac,rna}/`
  - Sampled-cell tables were written as JSONL because no parquet engine is installed locally.
- Data unification migration.
  - Shared Human cardiac niches references now live under `data/raw/references/human_cardiac_niches/`.
  - Regenerated Heart datasets now use explicit long dataset IDs.
  - Heart spatial inputs are split into per-modality `.h5ad` files.
- Heart NNLS smoke run through the unified runner.
  - Run: `results/human_cardiac_niches_sim_4zone_circles__atac__highly_variable__nnls__smoke`
  - Metrics: `results/comparisons/human_cardiac_niches_sim_4zone_circles__atac__highly_variable__nnls__smoke.csv`
- Tests: `15 passed, 1 skipped`.
