# Runtime Estimate

## Target Config

The estimate here is for:

```text
configs/experiments/all_methods_all_atac_datasets.yaml
```

This config expands to:

```text
42 jobs = 7 datasets x 6 methods x 1 feature set
```

Datasets:

```text
russell_250
human_cardiac_niches_sim_1zone_3ct_low_density
human_cardiac_niches_sim_1zone_10ct
human_cardiac_niches_sim_4zone_stripes
human_cardiac_niches_sim_4zone_circles
pbmc_granulocyte_sorted_10k_sim_equal_celltype
pbmc_granulocyte_sorted_10k_sim_observed_abundance
```

Methods:

```text
nnls
rctd
cell2location
tangram
destvi
spatialdwls
```

Feature set:

```text
highly_variable
```

## Data Size Context

Approximate local data sizes:

```text
data/raw/references/human_cardiac_niches      15G
data/raw/references/russell_250              94M
data/raw/references/pbmc_granulocyte_sorted_10k_multiome  848M
data/processed/datasets/russell_250          74M
Heart processed datasets                     748M to 1.5G each
PBMC processed datasets                      about 136M each
```

The Heart reference is the largest runtime driver.

## Timing Probes

Small timing probes were run and then cleaned up from:

```text
results/_runtime_probe
```

Measured results:

```text
Russell NNLS:                 about 8 seconds
Heart low-density NNLS:       about 27 seconds
Tangram smoke on Russell:     about 8 seconds after installing tangram-sc==1.0.4
RCTD on Russell:              still running after about 100 seconds, then stopped
```

The RCTD probe used:

```text
configs/methods/rctd_local.yaml
```

The initial Tangram smoke probe used:

```text
configs/methods/tangram_smoke.yaml
```

and failed with:

```text
module 'tangram' has no attribute 'pp_adatas'
```

The cause was that `tangram-sc` was not installed, and the legacy local directory `scripts/tangram/` was being imported instead. This has been fixed by:

- installing `tangram-sc==1.0.4`
- updating `src/deconvatac/tl/tangram.py` so it avoids legacy `scripts/tangram/` shadowing

After the fix, the Russell Tangram smoke run completed in about 8 seconds.

## Rough Estimate

For the current serial full config, budget:

```text
minimum: several hours
realistic local run: 8-24 hours
possible CPU-only worst case: 1-2 days
```

The estimate is rough because the heavy methods dominate runtime and depend strongly on hardware, installed backends, and package behavior.

## Main Bottlenecks

NNLS is not the bottleneck.

Likely bottlenecks:

- `cell2location`
  - Full config uses `max_epochs_spatial: 30000`.
  - Runtime can be long without a working GPU path.
- `destvi`
  - Full config uses `max_epochs_ref: 300` and `max_epochs_spatial: 2000`.
- `rctd`
  - RCTD was not finished after about 100 seconds even on the small Russell dataset.
  - Heart runs should be longer.
- `spatialdwls`
  - R-based and likely nontrivial on the larger Heart inputs.
- `tangram`
  - Smoke mode now runs after installing `tangram-sc==1.0.4`.
  - Full `num_epochs: 1000` runtime still needs to be measured.

## Practical Recommendation

Do not launch the full config casually.

First run a smoke config with reduced epochs and fewer methods, then run the full benchmark when it can run unattended overnight.

Suggested immediate checks:

1. Run `nnls` across all datasets to verify data loading and output writing.
2. Run `rctd` on Russell only and let it complete to get a real RCTD timing.
3. Measure a full Tangram run on Russell before launching Tangram across all Heart datasets.
4. Use smoke configs for `cell2location` and `destvi` before running full epoch counts.

## Command For Full Run

```bash
PYTHONPATH=src .venv/bin/python scripts/run_deconvolution.py \
  --experiment-config configs/experiments/all_methods_all_atac_datasets.yaml \
  --overwrite
```

Expected batch output:

```text
results/<timestamp>_all_methods_all_atac_datasets/
  runs.csv
  failures.csv        # only if methods fail
  comparison.csv
  <dataset>__<modality>__<feature_set>__<method>/
```
