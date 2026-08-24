# ShapeMix-ATAC benchmark protocol

Status: normative benchmark protocol, version 1, frozen 2026-08-22

This protocol predeclares the first ShapeMix experiment before primary-seed results are inspected. It tests the model in [model_specification.md](model_specification.md) using one 10x Genomics PBMC Multiome donor. Results quantify conditional resampling behavior in this dataset; they do not establish donor-level or population-level generalization.

## 1. Primary question and comparison

The primary question is whether parent-fragment length composition improves spot-level cell-type proportion recovery beyond total peak counts.

The primary paired comparison is:

```text
shapemix_length     use_shape: true
shapemix_count_only use_shape: false
```

The two variants must use the same raw fragments, cell labels, split, held-out source cells, spots, peaks, total counts, fixed `A`, fixed `phi_ref`, abundance prior, initialization policy, optimizer settings, seeds, stopping rule, and compute budget. The sole model difference is inclusion of the conditional multinomial shape term. A comparison is invalid if any other setting differs.

## 2. Dataset and admissible cells

### 2.1 Source

Use the 10x Genomics `pbmc_granulocyte_sorted_10k` Multiome output from Cell Ranger ARC 2.0.0, the version-matched peak BED and peak-barcode matrix, and the repository's SnapATAC2/GET `prepare_pbmc` barcode labels. Raw fragment processing follows the five-column and cut-site contract in the model specification.

This source contains one donor. RNA may support the fixed cell labels, but neither RNA counts nor held-out ATAC counts may enter ShapeMix signatures or peak selection.

The already resolved local inputs are pinned as follows:

| Input | SHA-256 |
|---|---|
| `data/raw/sources/snapatac2/pbmc10k_multiome/cell_type_mapping.csv` | `3aa94d5f636c01c0159324984930cb14e775e49611102cb4dd02a491b63bf298` |
| `data/raw/sources/snapatac2/pbmc10k_multiome/cell_type_summary.csv` | `809a3b285d996024daed95e0c2a8a17299ba93d1e5f0e884988dac9362150f4d` |
| `data/processed/references/pbmc_granulocyte_sorted_10k_multiome/atac/reference.h5ad` | `cdaefffbfd5dd3cb36318f68158d4bec0df1d2b2bf63562b14065f1055cf6ee6` |
| `data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/pbmc_granulocyte_sorted_10k_atac_peaks.bed` | `3975a4057f9caa3fb69ddaecc6ae9e530e77551717a1464c2d93ac9d73cb60ab` |
| `data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5` | `f6824171378787baab244f559b8b438f79db2eb39f78d17b2196f7ecd2c03549` |
| `data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz` | `5075e32a0e9c6dded35b060bf90d6144375b150e131ffb0be121a93e3b5e1e38` |
| `data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi` | `3a516291d0e6e5ddf9f651f470b6312b83eeb26f153ea12cfa9d082760a5e7f5` |

The complete versioned source record is [tracked here](../../configs/data_sources/pbmc_granulocyte_sorted_10k_cellranger_arc_2.0.0.yaml). A hash change is a new input version, not an in-place substitution.

### 2.2 Frozen cell-type universe

Intersect the prepared ATAC reference observations with the pinned labels, then include the 16 labels with at least 100 cells in that executable pre-split pool, in this exact order:

| Order | Cell type | Pre-split cells |
|---:|---|---:|
| 1 | CD14 Mono | 2,551 |
| 2 | CD4 Naive | 1,382 |
| 3 | CD8 Naive | 1,353 |
| 4 | CD4 TCM | 1,113 |
| 5 | CD16 Mono | 442 |
| 6 | NK | 403 |
| 7 | CD8 TEM_1 | 322 |
| 8 | CD8 TEM_2 | 315 |
| 9 | Intermediate B | 300 |
| 10 | Memory B | 298 |
| 11 | CD4 TEM | 286 |
| 12 | cDC | 180 |
| 13 | Treg | 157 |
| 14 | gdT | 143 |
| 15 | MAIT | 130 |
| 16 | Naive B | 125 |

The expected prepared reference has 9,627 labeled ATAC observations, of which these 16 types contribute exactly 9,500 cells. The raw label mapping contains 9,631 labeled barcodes before matrix intersection (with aggregate counts recorded in `cell_type_summary.csv`); three `CD14 Mono` barcodes and one `CD8 Naive` barcode are absent from the prepared ATAC reference and therefore cannot enter an executable split. Fail if the pinned hashes, reference dimensions, or post-intersection counts differ; do not silently recompute a new universe. `pDC`, `HSPC`, and `Plasma` are excluded from protocol version 1 because their post-intersection counts are below 100.

Every generated dataset YAML must declare this same ordered list at `truth.cell_types`. It must exactly equal the reference label universe and the ordered truth columns. A pseudo-spot may contain zero cells of a declared type; the type remains in the schema and metric denominator.

## 3. Split design and random-number contract

### 3.1 Random-number generator

Use NumPy `Generator(PCG64(SeedSequence(...)))`; record NumPy version, bit generator, full seed tuple, and output hashes. Never use process-global random state or Python's randomized hash order.

Use namespace integer `20260822` in every seed tuple. Sort input barcodes lexicographically before any seeded operation.

### 3.2 Development seeds

The following seeds may be used for code development, visualization, tuning, and acceptance tests:

```yaml
development:
  outer_split_seed: 0
  inner_mixture_seed: 0
```

Development results are never included in the primary summaries.

### 3.3 Frozen primary seeds

The five outer split seeds are:

```yaml
outer_split_seeds: [1103, 2203, 3301, 4409, 5501]
```

Within every outer split, use these two inner mixture seeds:

```yaml
inner_mixture_seeds: [101, 211]
```

Create the split RNG from `SeedSequence([20260822, outer_split_seed])`. Create each mixture RNG independently from `SeedSequence([20260822, outer_split_seed, inner_mixture_seed])`. A mixture stream is therefore nested in, and unique to, its split even though the two inner seed labels repeat.

Create the two-fold dispersion assignment from `SeedSequence([20260822, outer_split_seed, 17])`. Method restart streams use `SeedSequence([20260822, outer_split_seed, inner_mixture_seed, method_seed, restart_index])`; `method_seed` is 0 for both arms.

No primary configuration, threshold, peak rule, or exclusion rule may be changed after outputs from these primary pairs are inspected. A necessary correction requires a documented protocol version increment, invalidation of all prior primary runs, and complete reruns of both arms.

### 3.4 Stratified 70/30 split

For each outer seed:

1. Start from the exact 9,500-cell universe above.
2. Within each cell type, sort barcodes lexicographically.
3. Derive an independent type stream from `SeedSequence([20260822, outer_split_seed, cell_type_order])` and permute that type's sorted barcodes.
4. Put the first `floor(0.70 * n_c)` barcodes in the training-reference pool and the remainder in the held-out source pool.
5. Restore canonical cell-type order and lexicographic barcode order when writing files.

Fail if reference and held-out barcodes overlap, if a source barcode is duplicated, or if either pool lacks a declared type. Record split membership and its SHA-256 hash. The reference pool alone is visible to peak selection, signature estimation, smoothing, dispersion estimation, and model development for that outer split.

## 4. Deterministic reference-only peak selection

Select 5,000 peaks independently within each outer split. Use only the training-reference collapsed cut-site matrix and labels; do not read shape layers, held-out counts, held-out labels, pseudo-spots, or truth.

For every cell type `c` and candidate Cell Ranger peak `p`, let

```text
K[c,p] = sum of training-reference collapsed cut-site counts for type c at peak p.
```

For every type, fail if `sum_p K[c,p] == 0`. Define

```text
L[c,p] = log2(1 + 10^4 * K[c,p] / sum_p K[c,p]).
```

A peak is eligible when it is nonzero in at least 10 training-reference cells. Score each eligible peak with population variance across the 16 types:

```text
score[p] = variance_c(L[:,p], ddof=0).
```

Sort eligible peaks by this exact key:

1. descending `score[p]`;
2. descending number of nonzero training cells;
3. descending total training count;
4. ascending canonical peak ID using bytewise/ASCII order.

Take the first 5,000 and preserve this ranked order everywhere. Fail rather than changing `P` if fewer than 5,000 peaks are eligible. Record the candidate-universe hash, ordered selected-peak hash, scores, coverage, totals, and selection parameters.

The official full Cell Ranger peak set is used to validate cut-site coordinate semantics before applying this selector. Both ShapeMix arms and every baseline in a split receive the exact same selected peaks and order.

## 5. Pseudo-spatial mixtures

### 5.1 Primary conditions

Generate both conditions for every `(outer_split_seed, inner_mixture_seed)` pair:

- `equal_celltype`: all 16 declared types have probability `1/16`;
- `observed_abundance`: type probabilities equal the frozen post-intersection counts in Section 2.2 divided by 9,500.

For each condition:

- generate 1,024 spots on a 32 by 32 row-major grid;
- name them `spot_0000` through `spot_1023`;
- draw cells per spot as `max(1, Poisson(10))`;
- draw cell types from the condition's categorical probabilities;
- for each drawn type, sample uniformly from that type's held-out source barcodes;
- do not reuse a barcode within one spot when enough cells of that type are available;
- allow reuse across spots, because the held-out pool is smaller than all simulated nuclei combined;
- sum each source cell's three sparse layers to construct the spot;
- set spot `.X` to the exact layer sum;
- do not depth-thin the primary datasets;
- write exact source-cell provenance and truth proportions from sampled cell counts.

The two conditions must use separate deterministic substreams derived by appending condition index 0 (`equal_celltype`) or 1 (`observed_abundance`) to the mixture seed tuple. Source-cell sampling uses a further distinct substream from cell-type/count sampling.

The observed-abundance condition is the primary realism condition. Equal-celltype is a co-reported controlled condition that tests all types with similar expected representation. Neither may be dropped after results are known.

### 5.2 Smoke dataset

Use only development seeds for the smoke path: 32 spots, 100–500 top-ranked peaks, two or three declared types copied into a smoke-specific universe, and the same layer-conservation rules. Smoke output is for integration testing and is excluded from scientific comparisons.

### 5.3 Stress tests

Depth thinning, altered cells per spot, rare-cell enrichment/depletion, subtype challenges, alternate peak counts, alternate bins, backgrounds, and signature smoothing are secondary sensitivity analyses. They require separate dataset IDs/configurations. They cannot replace or be pooled with the two primary conditions.

## 6. Fixed signatures and fitting

Estimate `A`, `omega`, and cross-fitted `phi_ref` from the training-reference pool exactly as specified in the model specification. The two model arms consume the same serialized signature artifact or an artifact with the same content hash.

Use these model-version-1 values in both arms:

```yaml
signature_rate_pseudocount: 0.5
signature_shape_concentration: 1.0
total_likelihood: negative_binomial
conditional_shape_likelihood: multinomial
dispersion_mode: reference_crossfit_global_scaled_by_abundance
dispersion_crossfit_folds: 2
dispersion_alpha_floor: 1.0e-8
exposure_mode: absorbed_in_abundance
background_mode: none
abundance_prior: gamma
abundance_prior_shape: 2.0
abundance_prior_rate: 1.0
optimizer: adam
learning_rate: 0.03
max_steps: 2000
patience: 100
tolerance: 1.0e-5
restarts: 3
spot_batch_size: 64
peak_chunk_size: 512
seed: 0
device: cpu
dtype: float32
```

`shapemix_count_only` changes only `use_shape` from `true` to `false`. Runs with non-finite values, no finite converged restart, missing output types/spots, or invalid proportion rows fail and remain represented as failures; they are not silently excluded.

The single `signature_shape_concentration` value is applied at both hierarchy levels defined in the model specification: pooled peak distributions shrink toward the dataset-global distribution, and cell-type/peak distributions shrink toward those pooled peak targets.

## 7. Metric input contract

Metrics receive the fixed ordered universe from `truth.cell_types`, never the intersection of truth and prediction columns.

Before scoring:

1. Truth and prediction indices and columns must each be unique.
2. Prediction must contain exactly the truth spot set. Reordering is allowed only after set equality is established.
3. Missing declared prediction columns are inserted as zero so omission is penalized.
4. Prediction columns outside the declared universe are an error.
5. Truth columns and the reference label universe must exactly equal the declared universe and order.
6. All entries must be numeric, finite, and nonnegative, with no all-zero rows.
7. Every truth and prediction row must sum to one within absolute tolerance `1e-6`. Do not silently normalize or drop invalid rows.

All metric outputs must be finite. A metric failure fails evaluation and remains in the run manifest.

## 8. Primary endpoints

### 8.1 RMSE

Over all `S × C` entries in the fixed universe,

```text
rmse_v1 = sqrt(mean_{s,c}((truth[s,c] - prediction[s,c])^2)).
```

Lower is better.

### 8.2 Base-2 Jensen–Shannon divergence

For each spot, let `M_s = (truth_s + prediction_s) / 2`. With `0 * log(0/q) = 0`, define

```text
JSD_2(truth_s, prediction_s)
  = 0.5 * KL_2(truth_s || M_s)
    + 0.5 * KL_2(prediction_s || M_s).

jsd_v2 = mean_s JSD_2(truth_s, prediction_s).
```

Equivalently, this is

```text
mean(jensenshannon(truth, prediction, base=2)^2).
```

`jsd_v2` has range `[0,1]`; lower is better. The repository's historical `jsd` function computes an unsquared Jensen–Shannon distance and must be reported, if needed, under the distinct name `js_distance_v1`. It is not interchangeable with `jsd_v2`.

RMSE and `jsd_v2` are co-primary and are reported separately for observed-abundance and equal-celltype datasets. No metric may be selected after viewing results.

## 9. Secondary endpoints and frozen rare-cell rules

Report per-cell-type MAE and RMSE, Pearson and Spearman correlation across spots, runtime, peak resident memory, convergence, and reconstruction diagnostics.

For protocol version 1, a **rare reference type** is a declared type with pre-split prevalence below 2% of the 9,500-cell executable universe. The frozen rare-type set is:

```text
cDC, Treg, gdT, MAIT, Naive B
```

For every `(spot, rare_type)` pair:

- truth presence is `truth proportion > 0` because truth comes from exact source-cell provenance;
- predicted presence is `predicted proportion >= 0.01`;
- a value exactly `0.01` is positive;
- the threshold is never tuned on primary seeds.

Report pooled micro precision, recall, and F1 across all rare-type pairs; per-type precision/recall/F1; macro averages over types with both evaluable classes; and pooled and per-type AUPRC using continuous predicted proportions. If a metric is undefined because only one class occurs, report it as undefined with the reason rather than deleting the type.

These rare-cell metrics are secondary. In the equal-celltype condition, the labels remain “rare reference types” even though their simulated spot abundance is not rare.

## 10. Paired analysis across nested seeds

For metric `M`, condition `q`, outer split `o`, and inner mixture `j`, define the paired effect

```text
delta[o,j,q] = M[shapemix_length,o,j,q]
               - M[shapemix_count_only,o,j,q].
```

For error metrics, negative values favor ShapeMix. First average the two inner effects within each outer split:

```text
delta_outer[o,q] = mean_j delta[o,j,q].
```

The five `delta_outer` values, not spots or the ten inner datasets, are the resampling units. For each condition and endpoint, report:

- both inner-pair effects;
- each outer-split mean effect;
- mean and median outer effect;
- standard deviation, minimum, and maximum across outer splits;
- a percentile 95% bootstrap interval formed by resampling the five outer effects with replacement using seed tuple `[20260822, 9001]` and 10,000 replicates;
- the number of outer splits with improvement (`delta_outer < 0`);
- an exact two-sided paired sign-flip test across the five outer effects, labeled exploratory because its resolution is low.

Do not bootstrap spots or treat inner mixtures as independent biological replicates. Intervals and tests describe conditional resampling variability within one donor.

The preregistered directional support rule is met for a condition only when both co-primary mean effects are below zero and at least four of five outer-split effects favor ShapeMix for each endpoint. Report effect sizes regardless of whether this descriptive rule is met. Do not translate the rule or exploratory p-values into a donor-generalization claim.

## 11. Required negative controls

Before interpreting a positive shape result, run these controls with development seeds and then the same frozen primary pairs if the development behavior is correct:

1. Replace every cell type's `omega[c,p,:]` with the pooled `u[p,:]`; ShapeMix and peak-only estimates should agree within numerical tolerance.
2. Permute complete cell-type shape signatures with one deterministic permutation per outer split while leaving `A` fixed; a legitimate gain should disappear or reverse.
3. Collapse to one bin; the conditional term must be zero and estimates must agree.
4. Verify Poisson total factorization against independent Poisson bin likelihood on toy data.

Label all other changes—alternative bins, peak counts, Poisson totals, fixed background, smoothing, exposure adjustment, or grouped shape signatures—as sensitivities rather than primary analyses.

## 12. Exclusions, failures, and reporting

- Do not exclude a spot, type, split, seed pair, method run, or metric after seeing its effect.
- Input-contract violations fail before fitting.
- Optimizer failures and resource-limit failures remain in `failures.csv` and benchmark summaries.
- If one arm fails for a pair, the paired effect is unavailable and the failure is reported; do not substitute an unpaired run.
- Never omit a declared prediction type to improve a score.
- Record dataset, split, feature, signature, config, code revision, environment, and output hashes.
- Report count and shape likelihood components, count reconstruction, bin-composition reconstruction, convergence, runtime, and memory alongside accuracy.
- State prominently that all primary resampling is from one donor.

## 13. Protocol completion gate

The benchmark can begin only after:

- raw fragments and tabix index pass pinned hashes;
- the right-cut coordinate convention reconstructs the official peak matrix to the documented acceptance tolerance;
- all layer/data-contract invariants pass;
- split disjointness and peak-selection determinism pass;
- model golden tests and negative controls pass;
- both variants complete the development smoke dataset;
- this protocol and all primary configuration files are committed or otherwise content-hashed before primary outputs are opened.

Any departure must be recorded with its rationale and classified as a protocol amendment, new protocol version, sensitivity analysis, or exploratory analysis.
