# PBMC stress evaluation: support-safe version-2 amendment

## Status and reason for the amendment

The original `shapemix_pbmc_stress_v1` result tree is retained as failed-run
evidence and is not overwritten. The first evaluation seed exposed a data
contract problem after 51 durable jobs: ShapeMix requires every selected peak
to have positive pooled training-reference counts, but protocol v1 selected
the 5,000-peak axis before applying the subtype and reference-cell subsets.
The restricted references therefore contained unsupported peaks. NNLS could
run because it does not enforce the fragment-shape signature contract, but
those rows cannot make the incomplete campaign scientifically valid.

The audit found the following numbers of unsupported peaks on the frozen v1
5,000-peak axis:

| Reference variant | Unsupported v1 peaks |
|---|---:|
| all 16 types, all reference cells | 0 |
| broad 3 types, all reference cells | 0 |
| broad 3 types, 250 cells/type | 21 |
| broad 3 types, 100 cells/type | 80 |
| broad 3 types, 50 cells/type | 158 |
| CD4-related 3 types, all reference cells | 27 |

The union contains 174 v1 peaks. The failures were detected from input-contract
messages. No prediction values, truth values, accuracy metrics, or paired
effects were inspected when defining this repair.

## Frozen version-2 rule

Version 2 uses one common axis for every factor so the one-factor-at-a-time
comparisons remain interpretable:

1. Recompute the existing protocol-v1 ranking from the split-1103 training
   reference only and require its first 5,000 entries to equal the immutable
   v1 list.
2. For every candidate peak, determine whether its pooled count is positive in
   all six predeclared reference variants shown above. Reference-support caps
   use the same smallest-SHA-256-barcode rule as v1.
3. Retain candidates passing the original ten-nonzero-training-cell rule and
   this universal positive-support rule.
4. Take the first 5,000 retained candidates in the unchanged protocol-v1 rank
   order. This replaces the 174 unsupported v1 peaks.
5. Stream the raw fragments once on that axis using five bins. Derive the
   canonical three-bin and two-bin matrices by exact layer summation, validate
   the collapsed recount against the official accessibility matrix, and then
   materialize all references and held-out pseudo-spots.

This selection still uses training-reference counts only. Held-out counts,
pseudo-spot truth, predictions, and evaluation metrics do not enter the rule.

## File organization

Immutable and reusable v2 artifacts are separate from v1:

```text
data/processed/shapemix/pbmc_granulocyte_sorted_10k/sensitivity/stress_v2/
  feature_axis/
    selected_peaks.txt
    feature_support_audit.csv
    reference_cells.h5ad
    heldout_cells.h5ad
  fragment_shape_cache/
    bins_three_full_split_1103.h5ad
    bins_two_full_split_1103.h5ad
    bins_five_full_split_1103.h5ad
    manifest.yaml
  manifests/
    source_axis.yaml
    design.csv
    materialization.yaml
  splits/
    bins_two/heldout_cells.h5ad
    bins_five/heldout_cells.h5ad

data/processed/references/shapemix_pbmc_stress_v2/<variant>/
  reference.yaml
  atac/reference.h5ad

data/processed/datasets/pbmc_shapemix_stress_v2_<condition>_mix_<seed>/
  dataset.yaml
  atac/spatial.h5ad
  atac/features/highly_variable.txt
  truth/proportions.csv
  simulation/manifest.yaml
  simulation/source_cells_by_spot.jsonl

data/work/preprocessing/shapemix_pbmc_stress_v2/
  source_axis.<temporary-id>/

results/sensitivity/shapemix_pbmc_stress_v2/
  shapemix_pbmc_stress_protocol_v2_cuda/
```

Temporary recount/materialization state is confined to `data/work/`. The
completed `results/` tree remains exposed to Git, consistent with
`file_organization.md`.

## Execution and acceptance

The complete 40-dataset, 120-job PBMC design is rerun. Successful v1 jobs are
not mixed with v2 jobs because their feature axes differ. The v2 runner is
fail-closed (`continue_on_error: false`), and the summary is accepted only when
all 120 jobs and all 240 metric rows are present, unique, finite, and
successful.

Tracked entry points:

- `configs/datasets/shapemix_pbmc_stress_v2.yaml`
- `configs/experiments/shapemix_pbmc_stress_v2.yaml`
- `scripts/prepare_shapemix_pbmc_sensitivity_v2.py`
- `scripts/summarize_shapemix_pbmc_sensitivity_v2.py`
- `scripts/run_shapemix_full_evaluation_v2.sh`
- `configs/experiments/shapemix_full_evaluation_v2.yaml`
