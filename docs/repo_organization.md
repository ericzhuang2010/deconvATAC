# Repository Organization

## Mental Model

This repository is organized around a unified deconvolution workflow:

```text
registered dataset -> standardized loader -> method adapter -> standardized run output -> evaluation
```

## Source Code

- `src/deconvatac/data/`: shared data/result contracts, dataset registry loading, input validation.
- `src/deconvatac/methods/`: unified adapters for `cell2location`, `rctd`, `tangram`, `destvi`, `spatialdwls`, and `nnls`.
- `src/deconvatac/metrics/`: standardized metrics such as RMSE and JSD.
- `src/deconvatac/pp/`: preprocessing helpers, including feature selection.
- `src/deconvatac/tl/`: older/legacy method wrappers that the new adapters call.

## Configs

- `configs/datasets/`: dataset configs for older/local example datasets.
- `configs/methods/`: method parameter configs, including smoke-test configs.
- `configs/experiments/`: executable experiment-level configs for batch runs.
- `data/registry/datasets.yaml`: central lookup table mapping dataset IDs like `russell_250` and `human_cardiac_niches_sim_4zone_circles` to their dataset configs.

## Data

- `data/raw/`: intended place for immutable downloaded raw data.
- `data/processed/references/human_cardiac_niches/`: large shared Heart reference files, with canonical ATAC/RNA `reference.h5ad` files plus `reference.yaml`.
- `data/processed/references/russell_250/`: Russell ATAC reference file plus `reference.yaml`.
- `data/processed/datasets/`: standardized method-ready inputs.
  - `russell_250/`
  - `human_cardiac_niches_sim_1zone_3ct_low_density/`
  - `human_cardiac_niches_sim_1zone_10ct/`
  - `human_cardiac_niches_sim_4zone_stripes/`
  - `human_cardiac_niches_sim_4zone_circles/`
- `data/processed/datasets/{dataset_id}/{modality}/features/`: selected feature-list files used by the runner.
- `data/archive/legacy_results/`: canonical archive of legacy paper-style outputs with a manifest.
- `results/`: new standardized outputs from `scripts/run_deconvolution.py`, including batch `comparison.csv` files.

## Scripts

- `scripts/run_deconvolution.py`: main unified runner for one registered dataset/method or an experiment config.
- `scripts/evaluate_runs.py`: evaluates standardized runs.
- `scripts/prepare_feature_annotations.py`: writes feature annotations into processed `.h5ad` files.
- `scripts/prepare_shared_feature_sets.py`: computes shared feature-list text files for large references.
- `scripts/regenerate_heart_simulations.py`: regenerates the Human cardiac niches simulation datasets.
- `scripts/migrate_legacy_results.py`: archives old result folders.
- `scripts/legacy/`: old root scripts retained for reference only.
- `scripts/{method}/`: older SEML/Sacred-style method-specific runners.

## Notebooks

- `notebooks/reproduction/`: paper reproduction notebooks.
- `notebooks/reproduction/legacy/`: old root notebooks retained for reference only.
- `notebooks/method_development/`: method inspection and development notebooks.
- `notebooks/exploratory/`: ad hoc exploration.
- `data/example_notebooks/`: large input data kept for compatibility with existing notebooks.

## Documents And Tests

- `docs/restructuring_plan_spatial_atac_deconvolution.md`: overall restructuring plan and design rationale.
- `docs/remaining_steps.md`: current remaining work.
- `docs/repo_organization.md`: this repository organization guide.
- `tests/`: unit tests for loaders, contracts, method registry, feature selection, and result standardization.

## Current State

Brain datasets are still not registered because the required input files are not present. Russell and the renamed Human cardiac niches simulations are usable through the unified registry/runner. Legacy outputs have been consolidated under `data/archive/legacy_results/`; new outputs should be written under top-level `results/`.
