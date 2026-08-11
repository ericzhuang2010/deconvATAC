# deconvATAC spatial ATAC deconvolution workspace

This repository is organized as a local research workspace for reproducing and extending the benchmark from:

> Spatial transcriptomics deconvolution methods generalize well to spatial chromatin accessibility data

The current goal of the workspace is to make spatial ATAC deconvolution experiments reproducible across multiple methods. The main path is:

```text
registered dataset -> standardized loader -> method adapter -> standardized run output -> evaluation
```

The maintained entry points are the dataset registry in `data/registry/datasets.yaml` and the unified runner in `scripts/run_deconvolution.py`.

## Repository Structure

```text
configs/
  datasets/              Dataset configs for older/local examples.
  experiments/           Experiment specs and benchmark plans.
  methods/               Method parameter configs for the unified runner.

data/
  registry/              Central dataset registry.
  raw/references/        Immutable shared single-cell references.
  processed/datasets/    Standardized method-ready datasets.
  archive/legacy_results/ Archived outputs from old scripts.
  example_notebooks/     Large compatibility inputs for reproduction notebooks.

results/                 Standardized outputs from new runs and comparisons.

src/deconvatac/
  data/                  Dataset loading, schemas, validation, and result writing.
  methods/               Unified method adapters.
  metrics/               Evaluation metrics.
  pp/                    Preprocessing and feature-selection helpers.
  tl/                    Older method/simulation code used by adapters and scripts.

scripts/
  run_deconvolution.py           Main unified method runner.
  evaluate_runs.py               Evaluates standardized run directories.
  regenerate_heart_simulations.py Regenerates Heart spatial ATAC/RNA simulations.
  prepare_feature_annotations.py Computes feature annotations into `.var`.
  prepare_shared_feature_sets.py Helper for recomputing reference-derived feature lists.
  legacy/                       Old root scripts retained for reference.
  cell2location/, rctd/, ...     Older method-specific experiment runners.

notebooks/
  reproduction/          Paper-reproduction notebooks and tables.
  method_development/    Method inspection and development notebooks.
  exploratory/           Ad hoc exploration.

docs/
  related papers/        Reference papers.
  notes/                 Project notes.
  ideas/                 Method-development notes.
  *.md                   Cleanup, restructuring, and workflow documentation.

tests/                   Unit tests for the maintained Python workflow.
```

Local environment folders may also exist:

- `.venv/`: local Python virtual environment.
- `.r-lib/`: local R package library used by `configs/methods/rctd_local.yaml`.

## Registered Inputs

Datasets are registered in `data/registry/datasets.yaml`. The currently registered IDs are:

- `russell_250`: processed Russell ATAC example.
- `human_cardiac_niches_sim_1zone_3ct_low_density`: regenerated Heart simulation, one zone, three selected cell types, lower cell count.
- `human_cardiac_niches_sim_1zone_10ct`: regenerated Heart simulation, one zone, ten selected cell types.
- `human_cardiac_niches_sim_4zone_stripes`: regenerated Heart simulation, four striped zones.
- `human_cardiac_niches_sim_4zone_circles`: regenerated Heart simulation, four circular zones.

The runner loads each dataset through its dataset config. A config defines:

- spatial data path
- reference data path
- modality, usually `atac` or `rna`
- cell-type label column
- spatial coordinate key
- truth proportions, when available
- feature sets such as `all`, `highly_variable`, and `highly_accessible`

Feature sets can come from `.var` columns in `.h5ad` files or from text files under each dataset's `<modality>/features/` directory.

## Setup

Use the local virtual environment if it already exists:

```bash
source .venv/bin/activate
export PYTHONPATH=src
```

For a fresh Python environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test,simulation]"
export PYTHONPATH=src
```

Install method-specific optional dependencies only when needed:

```bash
python -m pip install -e ".[cell2location]"
python -m pip install -e ".[destvi]"
python -m pip install tangram-sc==1.0.4
```

RCTD and SpatialDWLS require R packages. This workspace can use the project-local R library at `.r-lib/` for RCTD through `configs/methods/rctd_local.yaml`.

Prefer `pyproject.toml` extras for new installs when they work. Tangram can be installed directly as `tangram-sc==1.0.4`. `requirements.txt` is a pinned local environment snapshot and may contain machine-specific paths.

## Run Deconvolution Workflows

All maintained methods should be called through `scripts/run_deconvolution.py`.

### Run an Experiment Config

Use this mode for benchmark runs that span multiple datasets, methods, modalities, or feature sets.

```bash
PYTHONPATH=src .venv/bin/python scripts/run_deconvolution.py \
  --experiment-config configs/experiments/all_methods_all_atac_datasets.yaml \
  --overwrite
```

The experiment config defines:

- datasets
- modalities
- feature sets
- methods
- method config files
- metrics
- output location

The all-methods ATAC config runs:

```text
7 registered ATAC datasets x 6 methods x 1 feature set = 42 requested runs
```

Batch output is written under a timestamped folder:

```text
results/<timestamp>_all_methods_all_atac_datasets/
  runs.csv
  failures.csv                 only if any run fails
  comparison.csv               one table for cross-method comparison
  <dataset>__<modality>__<feature_set>__<method>/
    run.yaml
    inputs.yaml
    environment.txt
    results/
      proportions.csv
      diagnostics.json
      abundance.csv            optional
      truth.csv
      raw_method_output/       optional
```

`comparison.csv` contains RMSE and Jensen-Shannon divergence rows for successful runs. If `continue_on_error: true` is set in the experiment config, failed runs are also recorded with `status=failed` and an error message.

Available experiment configs:

- `configs/experiments/all_methods_all_atac_datasets.yaml`: all current ATAC datasets, all maintained methods, `highly_variable`.
- `configs/experiments/develop_spatial_atac_method.yaml`: smaller Russell ATAC development run.
- `configs/experiments/compare_cell2location_rctd.yaml`: Russell Cell2location versus RCTD comparison.
- `configs/experiments/reproduce_paper.yaml`: broader paper-reproduction style spec.

### Run One Method Manually

Basic pattern:

```bash
PYTHONPATH=src .venv/bin/python scripts/run_deconvolution.py \
  --dataset <dataset_id> \
  --modality <atac_or_rna> \
  --feature-set <feature_set> \
  --method <method> \
  --config configs/methods/<method>.yaml \
  --run-id <run_id> \
  --overwrite
```

Supported method names are:

- `cell2location`
- `rctd`
- `tangram`
- `destvi`
- `spatialdwls`
- `nnls`

Each run writes to:

```text
results/<run_id>/
  run.yaml
  inputs.yaml
  environment.txt
  results/proportions.csv
  results/diagnostics.json
  results/abundance.csv          optional
  results/raw_method_output/     optional
```

## Example Workflows

### NNLS Baseline on Russell ATAC

```bash
PYTHONPATH=src .venv/bin/python scripts/run_deconvolution.py \
  --dataset russell_250 \
  --modality atac \
  --feature-set highly_variable \
  --method nnls \
  --config configs/methods/nnls.yaml \
  --run-id russell_250_atac_highly_variable_nnls \
  --overwrite
```

### RCTD on Russell ATAC

```bash
PYTHONPATH=src .venv/bin/python scripts/run_deconvolution.py \
  --dataset russell_250 \
  --modality atac \
  --feature-set highly_variable \
  --method rctd \
  --config configs/methods/rctd_local.yaml \
  --run-id russell_250_atac_highly_variable_rctd \
  --overwrite
```

Use `configs/methods/rctd.yaml` instead of `rctd_local.yaml` if you want R to use its normal library paths instead of `.r-lib/`.

### Cell2location on Russell ATAC

```bash
PYTHONPATH=src .venv/bin/python scripts/run_deconvolution.py \
  --dataset russell_250 \
  --modality atac \
  --feature-set highly_variable \
  --method cell2location \
  --config configs/methods/cell2location.yaml \
  --run-id russell_250_atac_highly_variable_cell2location \
  --overwrite
```

Cell2location usually needs the `cell2location` optional dependency and may require GPU-related configuration depending on the machine.

### RCTD on a Heart Simulation

```bash
PYTHONPATH=src .venv/bin/python scripts/run_deconvolution.py \
  --dataset human_cardiac_niches_sim_4zone_circles \
  --modality atac \
  --feature-set highly_variable \
  --method rctd \
  --config configs/methods/rctd_local.yaml \
  --run-id human_cardiac_niches_sim_4zone_circles__atac__highly_variable__rctd \
  --overwrite
```

### Tangram Smoke-Size Run

```bash
PYTHONPATH=src .venv/bin/python scripts/run_deconvolution.py \
  --dataset russell_250 \
  --modality atac \
  --feature-set highly_variable \
  --method tangram \
  --config configs/methods/tangram_smoke.yaml \
  --run-id russell_250_atac_highly_variable_tangram_smoke \
  --overwrite
```

## Evaluate Runs

Experiment-config mode evaluates successful runs automatically and writes `comparison.csv` in the batch output directory.

Evaluate one or more standardized run directories:

```bash
PYTHONPATH=src .venv/bin/python scripts/evaluate_runs.py \
  --runs results/russell_250_atac_highly_variable_nnls \
         results/russell_250_atac_highly_variable_rctd \
  --output results/comparisons/russell_250_nnls_vs_rctd.csv
```

By default this computes RMSE and Jensen-Shannon divergence against the truth file copied into each run directory. You can override the truth file:

```bash
PYTHONPATH=src .venv/bin/python scripts/evaluate_runs.py \
  --runs results/<run_id> \
  --truth path/to/truth.csv \
  --output results/comparisons/<comparison_name>.csv
```

## Rebuild Inputs

Most runs should use existing processed inputs. Rebuild them only when raw inputs, feature-selection settings, or simulation parameters change.

### Regenerate Heart Simulations

Regenerate all Heart datasets:

```bash
PYTHONPATH=src .venv/bin/python scripts/regenerate_heart_simulations.py \
  --overwrite
```

Regenerate a subset:

```bash
PYTHONPATH=src .venv/bin/python scripts/regenerate_heart_simulations.py \
  --datasets human_cardiac_niches_sim_4zone_circles \
  --overwrite
```

The script writes per-modality spatial files, truth proportions, simulation provenance, dataset configs, and registry entries under `data/processed/datasets/<dataset_id>/`.

### Recompute `.var` Feature Annotations

```bash
PYTHONPATH=src .venv/bin/python scripts/prepare_feature_annotations.py \
  --dataset russell_250 \
  --modalities atac \
  --n-top-features 20000 \
  --overwrite
```

This computes `highly_variable` and, for ATAC, `highly_accessible` feature annotations and writes processed `.h5ad` inputs.

## Add New Inputs or Methods

To add a new dataset:

1. Place raw immutable files under `data/raw/` or another stable data directory.
2. Create a dataset config describing spatial/reference inputs, labels, truth, and feature sets.
3. Register the config in `data/registry/datasets.yaml`.
4. Confirm it loads with `scripts/run_deconvolution.py` using `--method nnls` first.

To add a new method:

1. Add a method adapter under `src/deconvatac/methods/`.
2. Make the adapter return a standardized `DeconvolutionResult`.
3. Register the method name in `src/deconvatac/methods/registry.py`.
4. Add a method config under `configs/methods/`.
5. Add it to an experiment config and run it through `scripts/run_deconvolution.py --experiment-config`.

## Tests

Run the maintained test suite with:

```bash
PYTHONPATH=src .venv/bin/python -m pytest
```

## Legacy Files

Old root scripts and notebooks were moved to:

- `scripts/legacy/`
- `notebooks/reproduction/legacy/`

They are kept for reference. New work should use the unified runner, dataset registry, and standardized output directories.

## Citation

Sarah Ouologuem, Laura D Martens, Anna C Schaar, Maiia Shulman, Julien Gagneur, Fabian J Theis.
Spatial transcriptomics deconvolution methods generalize well to spatial chromatin accessibility data.
Bioinformatics, Volume 41, Issue Supplement_1, July 2025, Pages i314-i322.
https://doi.org/10.1093/bioinformatics/btaf268
