# ShapeMix results summary

Status: Step 6 completed 2026-08-23

## Main result

The frozen three-bin fragment-length term did **not** improve spatial ATAC deconvolution in the one-donor PBMC benchmark. Effects are shape-aware minus count-only, so negative values favor ShapeMix.

| Condition | Endpoint | Mean effect | Bootstrap 95% interval | Outer splits favoring shape |
|---|---:|---:|---:|---:|
| Observed abundance | `rmse_v1` | +0.000025 | [-0.000107, +0.000126] | 2/5 |
| Observed abundance | `jsd_v2` | +0.005767 | [+0.005492, +0.006042] | 0/5 |
| Equal cell type | `rmse_v1` | +0.005230 | [+0.005017, +0.005443] | 0/5 |
| Equal cell type | `jsd_v2` | +0.027054 | [+0.026397, +0.027640] | 0/5 |

Neither condition met the preregistered directional-support rule. Fragment length worsened JSD in every outer split under both conditions and worsened equal-cell-type RMSE in every split. Observed-abundance RMSE was mixed and close to zero.

## Execution and validation

- The benchmark contained 20 datasets: five outer reference/held-out splits, two inner mixture seeds, and two mixture conditions.
- Both arms ran on every dataset, giving 40/40 successful and converged primary fits.
- All required files, output manifests, paired assignments, and frozen code/config/protocol hashes validated before results were opened.
- There were no failed runs or unavailable paired or outer-split metrics.
- The final repository suite passed with 230 tests and one existing skip.

## Secondary findings

- Equal-cell-type RMSE worsened for all 16 cell types. Under observed abundance, only 6 of 16 types improved on average.
- Pooled rare-cell precision, recall, F1, and AUPRC decreased in both conditions; per-type AUPRC decreased for every frozen rare type.
- Shape-bin reconstruction improved only slightly, while total-count reconstruction worsened slightly. These reconstruction changes did not improve proportion accuracy.
- The shape-aware arm required approximately 1.74 times the count-only runtime. All primary runs stayed within approximately 1,286–1,682 MB sampled peak RSS.

## Baselines and controls

- NNLS completed 20/20 same-data runs but was worse than count-only on both primary endpoints in every outer split.
- Development degeneracy controls passed their exact algebraic checks. The deterministic permuted-signature run remained a diagnostic rather than a truth-scored result.
- Primary negative controls were not triggered because neither condition showed preregistered positive directional support.
- Cell2location, RCTD, SpatialDWLS, and stress/sensitivity experiments remain explicit unexecuted gates.

## Interpretation limit

All datasets are nested resamples from one PBMC donor. The result supports conditional conclusions for this fixed-signature, three-bin model and simulation protocol only. It does not establish donor-level uncertainty, population generalization, or performance on independent real spatial ATAC tissue.

See the [detailed Step 6 execution report](step6_results.md) for exact arm means, intervals, rare-cell results, reconstruction diagnostics, runtime, provenance, the blinded resource-pilot disclosure, and artifact locations.
