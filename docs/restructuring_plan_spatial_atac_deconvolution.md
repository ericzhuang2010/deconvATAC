# Restructuring Plan for Spatial ATAC Deconvolution Development

> Historical note: this plan predates the data unification migration. The current canonical data naming and layout are documented in `docs/data_unification_plan.md`.

## Purpose

This repository currently mixes three roles:

1. Reproducing results from `Spatial transcriptomics deconvolution methods generalize well to spatial chromatin accessibility data`.
2. Running existing deconvolution methods, especially Cell2location and RCTD.
3. Developing a new spatial ATAC deconvolution method.

The restructuring should keep the reproduction workflow intact while making the repository easier to use as a method-development and benchmarking framework. The central design goal is that every deconvolution method receives the same standardized input object and returns the same standardized result object, even when the underlying method has different dependencies, output files, or assumptions.

## Implementation Status

Implemented in this repository:

- Standardized input and result contracts under `src/deconvatac/data/`.
- Dataset registry and local dataset configs under `data/registry/` and `configs/datasets/`.
- Unified method interface and method registry under `src/deconvatac/methods/`.
- Adapters for Cell2location and RCTD around the existing `src/deconvatac/tl/` wrappers.
- Adapters for Tangram, DestVI, and SpatialDWLS around the existing `src/deconvatac/tl/` wrappers.
- The plain NNLS baseline is registered as `nnls`.
- Standardized proportion metrics under `src/deconvatac/metrics/`.
- Unified run, evaluation, and legacy-result dry-run migration scripts under `scripts/`.
- Method and experiment configs under `configs/methods/` and `configs/experiments/`.
- Focused tests for data contracts, method registry behavior, result standardization, and the NNLS baseline.
- Feature-selection preprocessing via `scripts/prepare_feature_annotations.py`.
- Processed Russell ATAC inputs with real `highly_variable` and `highly_accessible` annotations under `data/processed/datasets/russell_250/`.
- Heart feature lists now live inside each canonical Heart dataset under `{modality}/features/`.
- Regenerated Human cardiac niches simulation inputs now use explicit dataset IDs under `data/processed/datasets/human_cardiac_niches_sim_*/`.
- Standardized notebooks under `notebooks/reproduction/`, `notebooks/method_development/`, and `notebooks/exploratory/`.
- Copied legacy result folders into `data/archive/legacy_results/` with a manifest.
- Russell `highly_variable` smoke comparison across Cell2location and RCTD under `results/`.

Not yet migrated:

- Full Brain benchmark input `h5mu` files are not present locally, so their dataset configs are not fully registered.
- Legacy outputs have been consolidated under `data/archive/legacy_results/`. New generated outputs should use standardized run folders under top-level `results/`.

Feature-selection status:

- `russell_250` now uses processed copies with `reference.var["highly_variable"]`, `reference.var["highly_accessible"]`, and matching spatial `.var` columns. Each selected feature set currently contains 20,000 ATAC peaks.
- `heart_homogeneous_1zone`, `heart_heterogeneous_4zones`, and regenerated `heart_1` through `heart_4` now use shared feature-list files instead of `mode: all` placeholders for `highly_variable` and `highly_accessible` feature sets.
- Completed: for Heart inputs, feature masks are computed once per shared Heart reference and reused across `heart_1`, `heart_2`, `heart_3`, and `heart_4`. The local Heart ATAC reference is `139,835 x 429,828` (`6.8G`) and the Heart RNA reference is `704,296 x 32,732` (`8.5G`), so the configs now reference text feature lists rather than duplicated processed reference files.
- Completed: Heart simulations were regenerated with:

```bash
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/private/tmp PYTHONPATH=src .venv/bin/python \
  scripts/regenerate_heart_simulations.py \
  --overwrite
```

  The regenerated datasets are registered as `heart_1`, `heart_2`, `heart_3`, and `heart_4`. Each has ATAC shape `961 x 429,828`, RNA shape `961 x 32,732`, truth proportions in `.obsm["proportions"]`, spatial coordinates in `.obsm["spatial"]`, and simulation provenance records under `simulation/source_cells_by_spot.jsonl`. The original bundled examples remain available as `heart_homogeneous_1zone` and `heart_heterogeneous_4zones`.
- The preprocessing command used for Russell was:

```bash
PYTHONPATH=src MPLCONFIGDIR=/private/tmp .venv/bin/python scripts/prepare_feature_annotations.py \
  --dataset russell_250 \
  --modalities atac \
  --n-top-features 20000 \
  --overwrite
```

## Current Repository Observations

- Existing method wrappers live under `src/deconvatac/tl/`:
  - `cell2location.py`
  - `rctd.py`
  - `tangram.py`
  - `destvi.py`
  - `spatialdwls.py`
  - `metrics.py`
  - `simulate.py`
- Existing benchmark runners live under `scripts/{method}/experiment_runner.py` with method-specific Sacred/SEML configs.
- Existing results are stored mainly under `data/deconvolution_results/{method}/{modality}/{dataset_feature_set}/`.
- There are also top-level result folders, including `cell2location_results/` and `rctd_results/`.
- Input data currently appears in several places:
  - `data/example_notebooks/`
  - `data/raw/references/human_cardiac_niches/`
  - notebooks and notebook-derived files
- The current configs already encode important benchmark axes:
  - dataset
  - modality: `atac` or `rna`
  - feature set: `highly_variable` or `highly_accessible`
  - reference dataset
  - spatial dataset
  - label key, usually `cell_type`

## Target Principles

1. Separate immutable input data from generated outputs.
2. Keep reproduction assets distinct from active method-development assets.
3. Store all reusable input datasets behind a small registry instead of hard-coded paths.
4. Give all method adapters the same public interface.
5. Standardize result files so metrics and plotting do not need method-specific logic.
6. Keep method-specific dependencies isolated because RCTD uses R/spacexr and Cell2location uses Python/scvi dependencies.
7. Preserve current scripts during migration, then convert them into thin wrappers around the unified interface.

## Proposed Top-Level Layout

```text
deconvATAC/
  README.md
  pyproject.toml

  documents/
    restructuring_plan_spatial_atac_deconvolution.md
    related papers/
    ...

  data/
    raw/
      external/
      paper_reproduction/
    processed/
      datasets/
        {dataset_id}/
          dataset.yaml
          reference.h5mu
          spatial.h5mu
          truth/
            proportions.csv
          features/
            highly_variable.txt
            highly_accessible.txt
    registry/
      datasets.yaml
      feature_sets.yaml
    archive/
      legacy_results/
      legacy_notebook_outputs/

  results/
    {run_id}/
      run.yaml
      inputs.yaml
      environment.txt
      logs/
      results/
        proportions.csv
        abundance.csv
        diagnostics.json
        model_artifacts/
    comparisons/

  configs/
    datasets/
      russell_250.yaml
      heart_1.yaml
      heart_2.yaml
      heart_3.yaml
      heart_4.yaml
      brain_1.yaml
      brain_2.yaml
      brain_3.yaml
      brain_4.yaml
    methods/
      cell2location.yaml
      rctd.yaml
      nnls.yaml
    experiments/
      reproduce_paper.yaml
      compare_cell2location_rctd.yaml
      develop_spatial_atac_method.yaml

  src/deconvatac/
    data/
      registry.py
      loaders.py
      validators.py
      schemas.py
    methods/
      base.py
      registry.py
      cell2location.py
      rctd.py
      nnls.py
    workflows/
      run_deconvolution.py
      evaluate.py
      compare.py
    metrics/
      proportions.py
    plotting/
      spatial.py
      benchmark.py
    tl/
      legacy wrappers retained during migration

  scripts/
    run_deconvolution.py
    evaluate_runs.py
    migrate_legacy_results.py

  notebooks/
    exploratory/
    reproduction/
    method_development/

  tests/
    test_data_contract.py
    test_method_interface.py
    test_result_standardization.py
```

## Data Storage Plan

### `data/raw/`

Use for immutable downloaded or externally generated data. Files here should not be modified in place.

Recommended subfolders:

```text
data/raw/
  external/
    zenodo/
    cellxgene/
  paper_reproduction/
    original_files_used_by_publication/
```

Expected rule: if a file came from Zenodo, CELLxGENE, a collaborator, or the paper's original reproduction workflow, store it here with its original filename plus a small metadata note.

### `data/processed/datasets/{dataset_id}/`

Use for standardized inputs that are ready to feed into any deconvolution method.

Each dataset folder should contain:

```text
data/processed/datasets/human_cardiac_niches_sim_1zone_3ct_low_density/
  dataset.yaml
  reference.h5mu
  spatial.h5mu
  truth/
    proportions.csv
  features/
    highly_variable.txt
    highly_accessible.txt
```

`dataset.yaml` should include:

```yaml
dataset_id: human_cardiac_niches_sim_1zone_3ct_low_density
source: paper_reproduction
reference_path: reference.h5mu
spatial_path: spatial.h5mu
modalities:
  - atac
  - rna
labels_key: cell_type
spatial_key: spatial
available_feature_sets:
  - highly_variable
  - highly_accessible
truth:
  proportions: truth/proportions.csv
notes: "Human cardiac niches simulation, zone 1."
```

### `data/registry/`

Use for global lookup tables. This avoids hard-coded paths in scripts.

Example `data/registry/datasets.yaml`:

```yaml
human_cardiac_niches_sim_1zone_3ct_low_density:
  path: data/processed/datasets/human_cardiac_niches_sim_1zone_3ct_low_density/dataset.yaml
human_cardiac_niches_sim_1zone_10ct:
  path: data/processed/datasets/human_cardiac_niches_sim_1zone_10ct/dataset.yaml
russell_250:
  path: data/processed/datasets/russell_250/dataset.yaml
```

### `results/{run_id}/`

Use for all new generated outputs. Do not write new results directly into method-specific top-level folders.

Recommended run ID format:

```text
{date}_{dataset_id}_{modality}_{feature_set}_{method}
```

Example:

```text
results/2026-07-03_heart_1_atac_highly_accessible_rctd/
```

Each run folder should contain:

```text
run.yaml
inputs.yaml
environment.txt
logs/
results/
  proportions.csv
  abundance.csv
  diagnostics.json
  raw_method_output/
```

`proportions.csv` should be the main comparison file for metrics. Rows should be spatial spots/locations and columns should be cell types. Values should sum to 1 across cell types when the method returns proportions.

`abundance.csv` should be optional and used when a method naturally returns abundance rather than normalized proportions, as Cell2location often does.

`diagnostics.json` should store method-specific quality information, training loss summaries, runtime, convergence flags, selected hyperparameters, and warnings.

## Standard Input Contract

Create a small data object used by all methods:

```python
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import anndata as ad
import pandas as pd


@dataclass
class DeconvolutionInput:
    dataset_id: str
    modality: str
    feature_set: str
    spatial: ad.AnnData
    reference: ad.AnnData
    labels_key: str
    spatial_key: str = "spatial"
    truth: Optional[pd.DataFrame] = None
    output_dir: Optional[Path] = None
```

All loaders should guarantee:

- `spatial.var_names` and `reference.var_names` are aligned.
- `reference.obs[labels_key]` exists.
- spatial coordinates are available in `spatial.obsm[spatial_key]`.
- raw count and normalized layers are documented.
- feature selection has already been applied.
- optional ground-truth proportions are aligned to `spatial.obs_names`.

## Standard Result Contract

Create a single result object:

```python
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import pandas as pd


@dataclass
class DeconvolutionResult:
    method: str
    dataset_id: str
    modality: str
    feature_set: str
    proportions: pd.DataFrame
    abundance: Optional[pd.DataFrame] = None
    diagnostics: dict[str, Any] = field(default_factory=dict)
    output_dir: Optional[Path] = None
```

All methods should return `DeconvolutionResult`. Method-specific files can still be saved, but downstream evaluation should read from this object or from standardized files written by this object.

## Unified Method Interface

Add a base class under `src/deconvatac/methods/base.py`:

```python
from abc import ABC, abstractmethod


class BaseDeconvolver(ABC):
    method_name: str

    def __init__(self, **kwargs):
        self.config = kwargs

    @abstractmethod
    def run(self, data):
        """Run deconvolution and return DeconvolutionResult."""
```

Add a method registry under `src/deconvatac/methods/registry.py`:

```python
METHODS = {}


def register_method(name, cls):
    METHODS[name] = cls


def get_method(name):
    return METHODS[name]
```

Expected usage:

```python
from deconvatac.data.loaders import load_deconvolution_input
from deconvatac.methods.registry import get_method

data = load_deconvolution_input(
    dataset_id="heart_1",
    modality="atac",
    feature_set="highly_accessible",
)

method = get_method("rctd")(doublet_mode="full")
result = method.run(data)
result.write("results/2026-07-03_heart_1_atac_highly_accessible_rctd")
```

## Method Adapter Plan

### Cell2location

Create `src/deconvatac/methods/cell2location.py`.

The adapter should:

- Accept `DeconvolutionInput`.
- Call the existing implementation in `src/deconvatac/tl/cell2location.py` during migration.
- Save method-native files under `results/raw_method_output/`.
- Load `means_cell_abundance_w_sf.csv` and `q05_cell_abundance_w_sf.csv`.
- Store the raw mean abundance in `DeconvolutionResult.abundance`.
- Create `DeconvolutionResult.proportions` by row-normalizing abundance, with clear handling for zero rows.
- Record `N_cells_per_location`, `detection_alpha`, GPU flag, epoch counts, and training diagnostics in `diagnostics.json`.

Cell2location-specific config:

```yaml
method: cell2location
params:
  N_cells_per_location: 8
  detection_alpha: 20
  use_gpu: true
  max_epochs_spatial: 30000
  max_epochs_ref: null
```

### RCTD

Create `src/deconvatac/methods/rctd.py`.

The adapter should:

- Accept `DeconvolutionInput`.
- Call the existing implementation in `src/deconvatac/tl/rctd.py` during migration.
- Save method-native files under `results/raw_method_output/`.
- Load `estimated_proportions.csv`.
- Standardize index and cell-type columns.
- Store output in `DeconvolutionResult.proportions`.
- Record R library path, `doublet_mode`, spacexr settings, and R session details in `diagnostics.json`.

RCTD-specific config:

```yaml
method: rctd
params:
  doublet_mode: full
  r_lib_path: null
  create_rctd_kwargs:
    CELL_MIN_INSTANCE: 0
    gene_cutoff: 0
    fc_cutoff: 0
    gene_cutoff_reg: 0
    fc_cutoff_reg: 0
    UMI_min: 0
```

Note: the method is conventionally spelled `RCTD`, not `RTCD`.

### Future Spatial ATAC Method

The placeholder prototype has been removed. Future spatial ATAC method development should add a clearly named adapter under `src/deconvatac/methods/`, register it in `src/deconvatac/methods/registry.py`, and provide a method config under `configs/methods/`.

The plain non-negative least-squares baseline remains available as `nnls` for smoke tests and input validation.

## CLI and Workflow Plan

Add one command-line entry point for all methods:

```bash
python scripts/run_deconvolution.py \
  --dataset heart_1 \
  --modality atac \
  --feature-set highly_accessible \
  --method rctd \
  --config configs/methods/rctd.yaml \
  --output-root results
```

The script should:

1. Resolve dataset paths through `data/registry/datasets.yaml`.
2. Load and validate `DeconvolutionInput`.
3. Instantiate the selected method adapter.
4. Run the method.
5. Write standardized results.
6. Write the exact resolved config and environment information.

Add a comparison command:

```bash
python scripts/evaluate_runs.py \
  --runs results/2026-07-03_heart_1_atac_highly_accessible_cell2location \
         results/2026-07-03_heart_1_atac_highly_accessible_rctd \
         results/2026-07-03_heart_1_atac_highly_accessible_nnls \
  --metrics rmse jsd \
  --output results/comparisons/heart_1_atac_highly_accessible.csv
```

## Reproduction Workflow Plan

Keep paper reproduction separate from active method development.

Recommended folders:

```text
notebooks/reproduction/
  01_prepare_inputs.ipynb
  02_run_paper_methods.ipynb
  03_recreate_tables.ipynb
  04_recreate_figures.ipynb

notebooks/method_development/
  01_baseline_model.ipynb
  02_method_prototype.ipynb
  03_compare_against_cell2location_rctd.ipynb
```

Existing notebooks should be moved gradually. Do not move everything at once unless paths are updated and notebooks are smoke-tested.

## Migration Plan

### Phase 1: Add structure without breaking existing code

- Add `data/registry/`.
- Add `configs/datasets/`, `configs/methods/`, and `configs/experiments/`.
- Add `src/deconvatac/data/` schemas, loaders, and validators.
- Add `src/deconvatac/methods/` base classes and registry.
- Keep existing `src/deconvatac/tl/` wrappers unchanged.
- Keep existing `scripts/{method}/experiment_runner.py` unchanged.

### Phase 2: Wrap Cell2location and RCTD

- Implement `Cell2LocationDeconvolver`.
- Implement `RCTDDeconvolver`.
- Add tests using small AnnData objects.
- Confirm both adapters write standardized `proportions.csv`.
- Confirm row/column alignment with ground truth.

### Phase 3: Centralize experiment execution

- Add `scripts/run_deconvolution.py`.
- Add method configs for Cell2location and RCTD.
- Add dataset configs for the current benchmark datasets.
- Re-run one small Cell2location example and one RCTD example through the unified runner.

### Phase 4: Standardize evaluation

- Move or wrap existing RMSE/JSD metrics into `src/deconvatac/metrics/`.
- Add `scripts/evaluate_runs.py`.
- Require evaluation to read standardized `proportions.csv`.
- Generate a comparison table with columns:
  - `run_id`
  - `dataset_id`
  - `modality`
  - `feature_set`
  - `method`
  - `metric`
  - `value`

### Phase 5: Develop the new method

- Start with a deterministic baseline.
- Add the new model behind the same method interface once it has a specific name.
- Compare first against Cell2location and RCTD on one small dataset.
- Scale to the full paper reproduction datasets only after the small comparison is stable.

### Phase 6: Archive legacy outputs

- Copy top-level `cell2location_results/`, top-level `rctd_results/`, and `data/deconvolution_results/` into `data/archive/legacy_results/`.
- Keep a manifest that records old path to archived path.
- Avoid deleting old outputs until all figures/tables have been reproduced from standardized run folders.

## Testing Plan

Add focused tests before moving large files:

- `test_data_contract.py`
  - verifies feature alignment
  - verifies labels key exists
  - verifies spatial coordinates exist
  - verifies ground truth aligns to spatial observations

- `test_method_interface.py`
  - verifies every registered method has `.run(data)`
  - verifies every method returns `DeconvolutionResult`
  - uses lightweight mock adapters where real dependencies are unavailable

- `test_result_standardization.py`
  - verifies `proportions.csv` row sums
  - verifies cell-type columns are stable
  - verifies Cell2location abundance can be normalized
  - verifies RCTD output parsing

## Immediate Next Steps

Historical checklist from the original restructuring plan. For current open work, use `documents/remaining_steps.md`.

1. Create empty target folders:
   - `configs/datasets/`
   - `configs/methods/`
   - `configs/experiments/`
   - `data/registry/`
   - `src/deconvatac/data/`
   - `src/deconvatac/methods/`

2. Add the first dataset registry entries for:
   - `russell_250`
   - `heart_1`
   - `heart_2`
   - `heart_3`
   - `heart_4`
   - `brain_1`
   - `brain_2`
   - `brain_3`
   - `brain_4`

3. Add method configs for:
   - `cell2location`
   - `rctd`
   - `nnls`

4. Implement the data/result dataclasses and validators.

5. Implement Cell2location and RCTD adapters around the existing wrappers.

6. Run one small end-to-end comparison:
   - dataset: `russell_250`
   - modality: `atac`
   - feature set: `highly_variable`
   - methods: `cell2location`, `rctd`

7. Only after that works, start moving legacy notebooks and results into the new folder structure.

## Decisions to Make Before Implementation

These decisions are accepted for the current restructure:

- Standardized processed inputs should prefer one `h5mu` per multi-modal dataset. Single-modality examples may remain as `h5ad` files. The loader supports both formats.
- Cell2location comparisons should use row-normalized proportions by default for fair comparison against RCTD and other proportion-returning methods. Raw Cell2location abundance should still be saved as `abundance.csv`.
- Run IDs should remain human-readable for now. A config hash can be added later when large parameter sweeps make collisions or provenance ambiguity more likely.
- SEML/Sacred should become optional behind the unified runner. Existing SEML/Sacred scripts should remain available for cluster execution, but new local workflows should use `scripts/run_deconvolution.py`.
- Any future new method should start as ATAC-only or single-modality deconvolution. Joint RNA/ATAC modeling should be deferred until the unified single-modality interface is stable and benchmarked.
