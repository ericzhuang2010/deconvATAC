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
