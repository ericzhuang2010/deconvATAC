# What `run_deconvolution.py` Uses

`scripts/run_deconvolution.py` does not use the top-level `notebooks/` directory.

When you run the unified deconvolution workflow, the runner uses:

- `data/registry/datasets.yaml`
- registered `data/processed/datasets/*/dataset.yaml` files
- data files under `data/processed/datasets/`
- shared references under `data/raw/references/`
- method configs under `configs/methods/`
- experiment configs under `configs/experiments/` when `--experiment-config` is used
- Python code under `src/deconvatac/`
- output directory `results/`

It does not import, execute, or read any `.ipynb` files.

In single-run mode, the runner uses command-line flags such as `--dataset`, `--method`, and `--feature-set`.

In experiment-config mode, the runner uses a YAML file such as:

```text
configs/experiments/all_methods_all_atac_datasets.yaml
```

That config chooses the datasets, modalities, feature sets, methods, method configs, and metrics. The runner then writes a batch directory with `runs.csv`, optional `failures.csv`, and `comparison.csv`.

## Caveat

There are two legacy dataset configs under `configs/datasets/` that reference `data/example_notebooks/`, but those configs are not registered in `data/registry/datasets.yaml`, so they are not used by the normal runner.

Also, `data/example_notebooks/` is separate from top-level `notebooks/`.

## Practical Conclusion

Top-level `notebooks/` is safe to ignore for `scripts/run_deconvolution.py`.

Removing `notebooks/` would not break the unified runner, but it would remove reproduction and exploratory notebooks.
