# ShapeMix Step 6 execution report

Status: completed 2026-08-23

This is a non-normative execution report for the frozen [benchmark protocol](benchmark_protocol.md). It records what ran, what did not run, the result, and the limits of interpretation. It does not amend the preregistered endpoints after seeing results.

## Bottom line

The three-bin fragment-length term did **not** improve deconvolution in the frozen one-donor PBMC benchmark. The preregistered directional-support rule was false for both mixture conditions:

- In observed-abundance mixtures, RMSE was essentially mixed, while Jensen–Shannon divergence was worse with fragment length in all five outer splits.
- In equal-cell-type mixtures, both RMSE and Jensen–Shannon divergence were worse with fragment length in all five outer splits.
- The conditional primary negative-control campaign was therefore not triggered.

This is a valid negative result for the implemented fixed-signature, three-bin MAP model. It does not show that every possible fragment-shape model is unhelpful, and it does not support donor- or population-level generalization.

## Frozen design and integrity

The primary grid used five outer reference/held-out cell splits, two inner mixture seeds nested within each split, and two mixture conditions. Each of the resulting 20 registered datasets was fit with the shape-aware and count-only arms, giving 40 primary runs. Every dataset has 1,024 held-out pseudo-spatial spots, 16 declared cell types, and the same reference-only 5,000-peak axis for its outer split.

The effect is `shapemix_length - shapemix_count_only`. Negative values favor fragment length for the two lower-is-better primary endpoints. Inner-seed effects were averaged within each outer split before the five outer effects were summarized. The 95% intervals are percentile bootstraps of those five outer effects with 10,000 frozen replicates. Exact paired sign-flip tests are exploratory.

Before the summary was opened:

- all five deterministic shards completed with eight successful runs each;
- all 40 fits reported convergence;
- all 40 output manifests and required files rehashed successfully;
- the 20 paired dataset units were complete and disjoint across shards;
- batch/run cross-references, resolved assignments, and current config, protocol, registry, and recorded code hashes matched; and
- there were no failed runs or unavailable inner/outer metrics.

Frozen provenance includes:

- benchmark protocol SHA-256: `1c540b1f195af1d398b2ac55a80c06236f4d01aebd2e4daaeb3a97028f8667b6`;
- primary experiment source SHA-256: `0604074d476b125df006565e59a909fe6d5f0983ae6773455889acf8ab16b09b`;
- executed-code manifest SHA-256: `68896dbdb4f575e2b315f896c95bde1fb8a2231c500fe900d57e8d3e9bdc7c40`;
- recorded Git commit: `31128334f6f8d64aa3b7db767a1324c298dab457` with the Step 0–6 worktree recorded as dirty; and
- strict `summary.yaml` SHA-256: `457496f7a5c8f84a0b72f78c6e652cf73babdcfc74c2dcc9ec27366564b5b22c`.

## Primary endpoints

| Condition | Endpoint | Count-only | Length | Length minus count | Bootstrap 95% interval | Outer splits favoring length | Sign-flip p | Rule met |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Observed abundance | `rmse_v1` | 0.046732309 | 0.046757315 | +0.000025006 | [-0.000107341, +0.000126247] | 2/5 | 0.75 | No |
| Observed abundance | `jsd_v2` | 0.122298176 | 0.128065489 | +0.005767313 | [+0.005492311, +0.006042315] | 0/5 | 0.0625 | No |
| Equal cell type | `rmse_v1` | 0.052396156 | 0.057626206 | +0.005230051 | [+0.005017263, +0.005442838] | 0/5 | 0.0625 | No |
| Equal cell type | `jsd_v2` | 0.156500738 | 0.183554373 | +0.027053635 | [+0.026397299, +0.027640329] | 0/5 | 0.0625 | No |

The rule required both co-primary mean effects to be below zero and at least four of five outer effects to favor length for each endpoint within a condition. Neither condition was close to satisfying that joint rule. Relative to count-only, length increased observed-abundance JSD by 4.72%; in equal-cell-type mixtures it increased RMSE by 9.98% and JSD by 17.29%.

## Descriptive cell-type and rare-cell results

These analyses explain the primary result but do not replace it.

- Observed-abundance per-type RMSE improved on average for 6 of 16 types. Only CD14 Mono, cDC, CD4 TCM, and CD4 Naive improved in all five outer splits. Correlation decreased for 15 of 16 types; CD14 Mono was the exception.
- Equal-cell-type RMSE and correlation worsened on average for all 16 types.
- Among the five frozen rare types (`cDC`, `Treg`, `gdT`, `MAIT`, and `Naive B`), only cDC RMSE improved under observed abundance; it worsened under equal-cell-type mixtures. Treg, gdT, MAIT, and Naive B RMSE worsened in both conditions.
- Per-type rare-cell AUPRC decreased for every rare type in both conditions.

For the pooled rare-type detection summaries, positive changes would favor length because these scores are higher-is-better:

| Condition | Precision change | Recall change | F1 change | AUPRC change |
|---|---:|---:|---:|---:|
| Observed abundance | -0.008881 | -0.000340 | -0.011043 | -0.051427 |
| Equal cell type | -0.015535 | -0.002557 | -0.012879 | -0.029214 |

All rare-cell metrics were defined; no cell type was silently omitted.

## Reconstruction, convergence, and resources

The shape-aware model made very small improvements to shape-bin reconstruction RMSE in all three bins and all five outer splits: mean changes were between -0.00028 and -0.00039 under observed abundance and between -0.00059 and -0.00102 under equal-cell-type mixtures. Those gains did not improve composition accuracy. Total-count reconstruction RMSE instead increased by 0.00071 and 0.00139, respectively.

All 40 primary fits succeeded and converged.

| Condition | Count-only seconds/run | Length seconds/run | Paired runtime ratio | Count-only mean peak RSS | Length mean peak RSS |
|---|---:|---:|---:|---:|---:|
| Observed abundance | 2,038.8 | 3,352.2 | 1.631x | 1,620.0 MB | 1,466.3 MB |
| Equal cell type | 1,931.9 | 3,584.9 | 1.855x | 1,583.5 MB | 1,611.3 MB |

Across arms, the mean paired runtime ratio was 1.74x. Accumulated per-run time was 30.30 hours; sampled peak RSS ranged from 1,286.1 to 1,681.9 MB. The memory measurements are process-tree high-water diagnostics, not a claim that the length arm intrinsically uses less memory.

## Negative controls

The generated development controls produced hash-manifested local evidence:

- homogenized cell-type `omega` exactly matched count-only proportions (`max abs difference = 0`);
- one-bin shape exactly matched one-bin count-only, with shape log likelihood equal to zero;
- the factorized Poisson likelihood matched the independent-bin form to `1.78e-15`; and
- a deterministic, non-identity cell-type permutation of `omega` completed and was frozen as a diagnostic.

The homogenized control initially exposed that a data-only shape-likelihood constant could alter the optimizer's stopping scale. Before the primary launch, the likelihood gained an exact homogeneous-`omega` abundance-invariant path with zero `z` gradient, and the MAP stopping rule was centered on the initial shape likelihood. The development controls were then rerun successfully, and the frozen primary code manifest includes both corrections.

The primary protocol required dataset-level controls only after a positive directional result. Because neither condition triggered, no primary control outputs were generated or truth-scored. The permutation diagnostic is consequently not evidence for or against the negative primary result.

Before using this control after any future positive result, freeze a numeric definition of "the gain disappears or reverses," bind the attestation to the triggering strict-summary hash and supported conditions, and add campaign-level completeness/failure tracking. This prospective clarification is not needed to interpret the present negative result.

## NNLS and gated baselines

NNLS was the only activated external baseline. Its 20/20 runs succeeded, its output manifests revalidated, and it used byte-identical truth files and the same 1,024-spot by 16-type evaluation contract. Values below are means plus or minus SD across the five outer splits after averaging the two inner mixtures:

| Condition | Endpoint | NNLS | Count-only | Length |
|---|---:|---:|---:|---:|
| Observed abundance | RMSE | 0.055269 ± 0.000496 | 0.046732 ± 0.000659 | 0.046757 ± 0.000618 |
| Observed abundance | JSD | 0.146447 ± 0.002226 | 0.122298 ± 0.001467 | 0.128065 ± 0.001398 |
| Equal cell type | RMSE | 0.060405 ± 0.000884 | 0.052396 ± 0.000842 | 0.057626 ± 0.000692 |
| Equal cell type | JSD | 0.183971 ± 0.002168 | 0.156501 ± 0.002633 | 0.183554 ± 0.001948 |

NNLS was worse than count-only on all four endpoints in every outer split. It was also worse than length for three endpoints in every split; for equal-cell-type JSD, length was better in four of five splits. This descriptive comparison supports the functionality of the count-only implementation but does not alter the paired primary conclusion. NNLS averaged 3.16 seconds/run and 2,187 MB sampled peak RSS under the runner's measurement scope.

The following are explicit gates, not failed or missing result rows:

- Cell2location needs a frozen full-size device, time, and memory validation before its 30,000-epoch configuration can run.
- RCTD is gated because the `spacexr` R package is absent.
- SpatialDWLS is gated because Giotto dependencies are absent.

## Stress and sensitivity gates

No separately versioned stress datasets exist, so `configs/experiments/shapemix_stress_tests.yaml` deliberately contains no executable dataset jobs. Depth thinning, cells per spot, rare-cell enrichment/depletion, subtype challenge, 10,000/20,000 peak axes, alternate fragment-length cutoffs, background, smoothing, and signature-yield choices remain future, separately versioned analyses. No Step 6 claim is based on them.

## Blinded resource-pilot disclosure

Before the primary launch, the shape-aware `pbmc_granulocyte_sorted_10k_shapemix_equal_celltype_split_1103_mix_101` fit was used only as a blinded resource pilot. The inspected fields were dimensions, runtime and memory demand, convergence status, and output shape. Prediction values, truth values, accuracy metrics, and paired effects were not inspected. The retained diagnostic records 1,268.645 seconds of fit time; the older artifact does not record its thread environment. It remains under `/private/tmp/deconvatac-step6-profile-20260822/primary_profile_shapemix_length`, outside `results/primary`, and was excluded from the strict summary. It is a documented protocol departure used only to choose a safe five-shard, two-thread primary execution plan.

## Interpretation limits

All 20 datasets are nested resamples from one PBMC donor. The five outer splits quantify conditional variability under this preparation and simulation protocol; they are not biological replicates. This result therefore does not estimate donor-level uncertainty, establish population generalization, or validate ShapeMix on independent real spatial ATAC tissue. The MAP point estimates also have no calibrated posterior uncertainty.

The result directly supports only this statement: with the frozen three-bin, per-peak, fixed-signature model and these same-protocol held-out pseudo-spots, adding the conditional fragment-length likelihood did not improve the preregistered proportion endpoints and usually degraded them.

## Artifact map

Primary run groups:

```text
results/primary/shapemix_primary_ablation_protocol_v1__shard_00_of_05/
results/primary/shapemix_primary_ablation_protocol_v1__shard_01_of_05/
results/primary/shapemix_primary_ablation_protocol_v1__shard_02_of_05/
results/primary/shapemix_primary_ablation_protocol_v1__shard_03_of_05/
results/primary/shapemix_primary_ablation_protocol_v1__shard_04_of_05/
```

Strict summary:

```text
results/primary/shapemix_primary_ablation_protocol_v1_summary/
  summary.yaml
  run_metrics.csv
  paired_effects.csv
  outer_effects.csv
  primary_summary.csv
  cell_type_metrics.csv
  cell_type_paired_effects.csv
  rare_cell_metrics.csv
  rare_cell_paired_effects.csv
  performance.csv
  reconstruction.csv
  failures.csv
  provenance.csv
```

Baseline and development-control evidence:

```text
results/primary/shapemix_primary_baselines_protocol_v1/
  batch_manifest.yaml
  runs.csv
  comparison.csv

results/development/shapemix_negative_controls_v1/
  pbmc_granulocyte_sorted_10k_shapemix_equal_celltype_split_000_mix_000_smoke/
    control_evidence.yaml
    control_proportions.csv
    output_sha256.yaml
```

The complete `results/` tree, including these `results/primary` artifacts, is intentionally exposed to Git. Reproduce the ignored data first with the [data-directory recreation guide](<../recreate_data_directory (important).md>), then run the frozen Step 6 configurations and strict summarizer recorded in the [implementation plan](implementation_plan.md).

## Verification

The final repository suite passed with **230 tests and one existing skip** using the repository import path `PYTHONPATH=src:.`.
