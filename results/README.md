# Results

New deconvolution outputs should be written here by `scripts/run_deconvolution.py`.

Single-run output directories contain:

- `run.yaml`
- `inputs.yaml`
- `environment.txt`
- `results/proportions.csv`
- `results/diagnostics.json`
- optional `results/abundance.csv`
- optional `results/raw_method_output/`

Experiment-config runs create a batch directory:

```text
results/<run_group>/
  runs.csv
  comparison.csv
  failures.csv              optional
  <dataset>__<modality>__<feature_set>__<method>/
```

Older standalone comparison CSVs may be written under `results/comparisons/`, but new config-driven runs should use the batch-local `comparison.csv`.

## Tracking policy

The entire `results/` directory is intentionally exposed to Git. This includes development checks, primary campaigns, sensitivity analyses, external validation, real-spatial analyses, summaries, logs, failure records, and per-run provenance. Do not add nested ignore rules for result subdirectories.

Write partial downloads, preprocessing shards, package caches, and other disposable files under `data/work/` or `/tmp`, not under `results/`. A result directory should contain reviewable experiment evidence and should not be overwritten after its campaign is declared complete.
