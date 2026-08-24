# ShapeMix-ATAC MVP model specification

Status: normative MVP specification, version 1, frozen 2026-08-22

This document defines the data semantics and statistical model that an implementation named `shapemix` must follow. Changes that alter the count unit, bin boundaries, likelihood, depth interpretation, dispersion estimator, or fixed-signature construction require a schema or model version change and a corresponding update to the benchmark protocol.

## 1. Scientific estimand

The primary estimand is the cell-type proportion vector in each mixed spatial ATAC spot. ShapeMix tests whether the **parent-fragment length composition of observed Tn5 cut sites** adds deconvolution information beyond ordinary peak-level cut-site totals.

The primary comparison is nested:

- **Peak-only:** fit total peak counts.
- **ShapeMix:** fit those identical total peak counts and, conditionally on each total, their allocation across parent-fragment length bins.

The two arms differ only by the conditional shape likelihood. They do not use different count units or independent implementations.

## 2. Raw input semantics

### 2.1 Cell Ranger ARC 2.0 fragments schema

The canonical PBMC raw input is a BGZF-compressed, tabix-indexed fragments file with five tab-separated fields:

1. `chrom` — chromosome or contig name;
2. `chromStart` — zero-based adjusted fragment start;
3. `chromEnd` — adjusted half-open fragment end;
4. `barcode` — cell barcode;
5. `readSupport` — number of read pairs supporting the deduplicated fragment.

There is no strand column. Header/comment lines begin with `#` and are not records.

Every valid row is one deduplicated fragment. The primary analysis sets `read_support_policy: ignore`: it counts the row once and does not weight it by `readSupport`, because weighting would reintroduce PCR support into the biological count.

Reject and count records with a non-integer coordinate, `chromEnd <= chromStart`, or a missing field. Ignore and count valid records whose barcode is outside the requested cell set or whose contig is outside the canonical peak set.

### 2.2 Length and cut-site observations

For a valid fragment row, define parent-fragment length

```text
L = chromEnd - chromStart
```

Assign `L` to exactly one bin:

| Bin index | Biological label | Interval (bp) | AnnData layer |
|---:|---|---|---|
| 0 | short | `[0,100)` | `fragment_length_lt_100` |
| 1 | mono | `[100,250)` | `fragment_length_100_249` |
| 2 | long | `[250,∞)` | `fragment_length_ge_250` |

Thus, lengths 99, 100, 249, and 250 map to bins 0, 1, 1, and 2 respectively.

Each fragment emits a left and right Tn5 cut-site observation. Both observations retain the same parent-length bin. Assign each cut site independently to the unique non-overlapping Cell Ranger peak containing its coordinate. A fragment can therefore contribute:

- two counts to one peak;
- one count to each of two different peaks;
- one count if only one end is in a selected peak; or
- zero counts if neither end is in a selected peak.

The right-end coordinate convention is a data-validation decision, not a free model option. Step 2 compared the two plausible mappings (`chromEnd` and `chromEnd - 1`) by reconstructing the official Cell Ranger peak matrix. On a deterministic sample of 2,048 peaks and 512 barcodes (1,048,576 entries across 24 contigs), `chromEnd` reproduced all 158,472 official counts with zero mismatches. `chromEnd - 1` produced 58 mismatched entries and total absolute error 59. An independent all-peak chr21 audit also matched exactly only with `chromEnd`. Therefore the canonical mapping is `left_cut = chromStart`, `right_cut = chromEnd`, and `right_cut_offset = 0`, with half-open peak membership `peakStart <= cut < peakEnd`. The full audit statistics are pinned in the [source manifest](../../configs/data_sources/pbmc_granulocyte_sorted_10k_cellranger_arc_2.0.0.yaml). This offset 0 mapping is frozen for the primary Cell Ranger ARC 2.0 benchmark. Each external source must independently reconstruct its official matrix and may select offset 0 or -1; reference and mixture objects within one benchmark must use the same validated mapping.

When querying tabix peak by peak, begin the query at `max(0, peakStart - 1)`. Otherwise, a fragment ending exactly at `peakStart` is not returned even though its `chromEnd` cut site belongs to that peak. Contig-streaming implementations do not require this query extension.

This gives the primary count unit:

```text
deduplicated Tn5 cut sites, grouped by parent-fragment length
```

A unique-fragment tensor is not interchangeable with this tensor. It would require a separately versioned dataset and a matched fragment-count baseline.

### 2.3 Peak assignment and conservation

The canonical Cell Ranger peaks must be sorted and non-overlapping. A cut site is counted only if it belongs to exactly one selected peak. Ambiguous overlaps are a validation error rather than a double count.

Let `R[i,p,b]` be reference-cell counts and `Y[s,p,b]` be mixed-spot counts. Define collapsed counts:

```text
T[i,p] = sum_b R[i,p,b]
N[s,p] = sum_b Y[s,p,b]
```

The collapsed matrices must use the same cut-site definition as the bin layers. It is invalid to compare shape-layer cut sites with a fragment-count baseline.

## 3. Shape-aware AnnData contract

Reference and spatial objects use cells/spots as observations and selected peaks as variables. Each object must contain the three ordered CSR layers above. Required invariants are:

- layer entries are finite, nonnegative integers;
- every layer has exactly the `.X` observation and variable axes;
- `.X` is a finite, nonnegative integer matrix;
- `.X` equals the exact elementwise sum of the three layers;
- reference and spatial objects have identical ordered peak IDs;
- reference observations have non-missing values at the declared `labels_key`;
- spatial observations have the declared spatial coordinates and complete truth rows;
- the declared cell-type order equals the reference label universe and truth columns;
- slicing features preserves layer order and the conservation identity.

Both objects store equivalent metadata under `.uns["fragment_shape"]`:

```yaml
schema_version: 1
axis: parent_fragment_length_bp
count_unit: deduplicated_cut_sites
read_support_policy: ignore
peak_assignment: containing_nonoverlapping_peak
left_cut_offset: 0
right_cut_offset: 0
bins:
  - name: short
    min_inclusive: 0
    max_exclusive: 100
    layer: fragment_length_lt_100
  - name: mono
    min_inclusive: 100
    max_exclusive: 250
    layer: fragment_length_100_249
  - name: long
    min_inclusive: 250
    max_exclusive: null
    layer: fragment_length_ge_250
source_sha256: null
feature_sha256: null
split_sha256: null
coordinate_validation: null
software_versions: null
preprocessing_counters: null
matrix_counters: null
```

The `null` provenance values above are schema illustrations only. Canonical data must replace them with resolved values. `source_sha256` is a mapping that includes `fragments` and `tabix_index`; `split_sha256` identifies the cell split. `feature_sha256` is bound to the current ordered `.var_names` axis by hashing each UTF-8 feature name after prefixing it with its unsigned eight-byte big-endian byte length. This avoids ambiguous concatenations and makes reordering detectable.

`coordinate_validation` must record a `selected_right_cut_offset` of 0 or -1 that equals the stored `right_cut_offset`, plus `matrix_match: exact`, `mismatched_entries: 0`, and `absolute_error: 0`. `software_versions` must identify the software used to build the object. `preprocessing_counters` are immutable source-run QC: row validity, filtering, fragment assignment, read support, and per-bin totals from the raw-fragment counting run. They intentionally remain attached when an object is safely feature-sliced. `matrix_counters` instead describe the object as currently stored: they contain `assigned_cut_sites` and one `cut_sites_per_bin.<layer>` value per declared layer. Validators compare those values exactly with `.X.sum()` and each layer sum; the loader recomputes them after feature slicing.

The source-run counters obey these conservation identities:

```text
valid_rows + invalid_schema_rows + invalid_coordinate_rows = total_rows
assigned_cut_sites + cut_sites_outside_peaks = 2 * retained_fragments
sum(cut_sites_per_bin) = assigned_cut_sites
```

## 4. Indices and observed data

Let:

- `s = 1,...,S` index spots;
- `p = 1,...,P` index selected peaks;
- `b = 1,...,B` index the three length bins;
- `c = 1,...,C` index declared cell types;
- `Y[s,p,b]` be observed mixed-spot cut-site counts;
- `N[s,p] = sum_b Y[s,p,b]` be ordinary peak totals;
- `R[i,p,b]` be training-reference counts for cell `i`;
- `T[i,p] = sum_b R[i,p,b]` be training-reference peak totals.

All signatures and nuisance estimates are functions of the training-reference split only. Held-out source cells and pseudo-spots are not used to select peaks, estimate signatures, set smoothing, estimate dispersion, initialize global parameters, or select hyperparameters.

## 5. Fixed reference signatures

### 5.1 Total accessibility signature

For cell type `c`, let `I_c` be its training-reference cells and `n_c = |I_c|`. The primary same-protocol exposure is `h[i] = 1`. Define the pooled per-cell peak target

```text
g[p] = sum_i T[i,p] / sum_i h[i].
```

With `alpha_A = 0.5`, estimate

```text
A[c,p] = (sum_{i in I_c} T[i,p] + alpha_A * g[p])
         / (sum_{i in I_c} h[i] + alpha_A).
```

`A[c,p]` is the smoothed mean cut-site rate at peak `p` per training-reference cell of type `c`. It absorbs the tutorial's accessibility, cell-type yield, and peak-detectability factors; those factors are not separately inferred in the MVP.

### 5.2 Conditional shape signature

Let the dataset-global distribution be

```text
u_global[b] = sum_{i,p} R[i,p,b] / sum_{i,p,b'} R[i,p,b'].
```

Preprocessing fails if the global denominator is zero or any global bin has zero probability: that would mean the frozen three-bin representation is unsupported by the training reference. Use the single frozen `signature_shape_concentration` value `alpha_omega = 1.0` at both levels of a two-stage hierarchy. First smooth every pooled peak toward the global distribution:

```text
u[p,b]
  = (sum_i R[i,p,b] + alpha_omega * u_global[b])
    / (sum_{i,b'} R[i,p,b'] + alpha_omega).
```

Then smooth every cell-type/peak distribution toward that strictly positive pooled-peak target:

```text
omega[c,p,b]
  = (sum_{i in I_c} R[i,p,b] + alpha_omega * u[p,b])
    / (sum_{i in I_c,b'} R[i,p,b'] + alpha_omega).
```

For every `(c,p)`, `omega[c,p,:]` must be finite, strictly positive after smoothing, and sum to one within numerical tolerance. The layer and cell-type orders are never inferred from set iteration; they come from the dataset contract.

## 6. Cross-fitted reference dispersion

The MVP uses one global per-reference-cell inverse-dispersion. Estimate it without fitting each cell to a signature containing itself:

1. Within every cell type, deterministically divide training-reference cells into two folds using the benchmark's dispersion-fold seed.
2. For each reference cell `i`, estimate its held-out expected total `m[i,p]` from the opposite fold's cell-type mean with the same `alpha_A = 0.5` smoothing rule.
3. Compute

```text
numerator = sum_{i,p} ((T[i,p] - m[i,p])^2 - m[i,p])
denominator = sum_{i,p} m[i,p]^2

alpha_ref_raw = numerator / denominator
alpha_ref = max(alpha_ref_raw, 1e-8)
phi_ref = 1 / alpha_ref.
```

The calculation fails if a type cannot populate both folds, if `denominator <= 0`, or if any result is non-finite. Record the fold membership hash, seed, numerator, denominator, `alpha_ref_raw`, `alpha_ref`, and `phi_ref`. Re-estimate production `A` from all training-reference cells after cross-fitting.

In the mean/inverse-dispersion parameterization used here,

```text
E[X] = mu
Var[X] = mu + mu^2 / phi.
```

The spot-level inverse-dispersion is scaled by effective abundance as specified below. Peak-specific or learned dispersion is not part of model version 1.

## 7. Generative model

### 7.1 Parameters and deterministic quantities

For each spot and cell type, infer positive effective abundance `z[s,c]`. Define

```text
e[s] = sum_c z[s,c].
```

The primary benchmark has no background, so `beta[p] = 0`. The general fixed-background notation is retained only to make the factorization explicit:

```text
v[s,p,b] = beta[p] * u_background[p,b]
           + sum_c z[s,c] * A[c,p] * omega[c,p,b]

mu[s,p] = sum_b v[s,p,b]
        = beta[p] + sum_c z[s,c] * A[c,p]

rho[s,p,b] = v[s,p,b] / mu[s,p].
```

For model version 1, set `beta[p] = 0` and do not infer `u_background`. Guard `mu` with a documented numerical epsilon only inside logarithms/division; do not change the reported mean. Peaks with zero expected rate for all types are invalid selected features.

### 7.2 Factorized likelihood

The total-count likelihood is

```text
N[s,p] ~ NegativeBinomial(
    mean = mu[s,p],
    inverse_dispersion = max(e[s], epsilon) * phi_ref
).
```

Therefore,

```text
E[N[s,p]] = mu[s,p]
Var[N[s,p]] = mu[s,p]
              + mu[s,p]^2 / (max(e[s], epsilon) * phi_ref).
```

Conditionally on the observed total,

```text
Y[s,p,1:B] | N[s,p]
  ~ Multinomial(total_count=N[s,p], probabilities=rho[s,p,1:B]).
```

Implement its log probability directly:

```text
log p(Y | N, rho)
  = lgamma(N + 1)
    - sum_b lgamma(Y[b] + 1)
    + sum_{b:Y[b] > 0} Y[b] * log(rho[b]).
```

Compute `log(rho)` through a stable log-sum-exp normalization. A zero-total `(s,p)` row contributes exactly zero to the conditional term. A zero observed bin contributes exactly zero, even if its stabilized log probability is large in magnitude. Non-finite objectives or gradients are hard failures.

Independent negative-binomial likelihoods for the three bins are not the MVP. They generally do not collapse to the same negative-binomial total and would therefore confound the count model with the shape representation.

### 7.3 Prior and MAP objective

Use the shape/rate Gamma convention:

```text
z[s,c] ~ Gamma(shape=2.0, rate=1.0).
```

Parameterize optimization variables as

```text
z[s,c] = softplus(raw_z[s,c]) + epsilon.
```

With all signatures and `phi_ref` fixed, the two objectives are:

```text
peak-only objective
  = sum_{s,p} log p(N[s,p] | z,A,phi_ref)
    + sum_{s,c} log p(z[s,c])

ShapeMix objective
  = peak-only objective
    + sum_{s,p} log p(Y[s,p,:] | N[s,p],z,A,omega).
```

Both arms use the exact same `N`, `A`, prior, parameter tensor, initialization, restarts, optimizer, stopping rule, numeric type, and compute budget. `use_shape: false` removes only the last term. It does not collapse or regenerate data through a separate path.

The initial `shape=2, rate=1` choice is frozen for protocol version 1. A later change based on toy/development runs requires a protocol version bump before final-seed results are inspected; it can never be tuned independently for the two arms.

## 8. Depth and output interpretation

Do not include an observed spatial library-size offset in model version 1. Total spot depth relative to the average training-reference cell is absorbed into `z`.

Consequently:

- `abundance.csv` stores effective reference-cell-equivalent, depth-scaled abundance;
- it is not a calibrated nucleus count and need not sum to the simulated number of cells;
- the standardized scientific output is

```text
pi[s,c] = z[s,c] / sum_c z[s,c];
```

- each row of `proportions.csv` must be finite, nonnegative, include every declared type, and sum to one within `1e-6` absolute tolerance.

No credible intervals are reported by the MAP MVP.

## 9. Frozen optimization behavior

Protocol-version-1 runs use:

```yaml
optimizer: adam
learning_rate: 0.03
max_steps: 2000
patience: 100
tolerance: 1.0e-5
restarts: 3
spot_batch_size: 64
peak_chunk_size: 512
device: cpu
dtype: float32
```

Initialize each restart from collapsed-count nonnegative least squares with a uniform positive fallback and deterministic restart perturbations. Select the finite converged restart with the largest complete objective for that arm. Record restart seeds, loss components, step count, stopping reason, convergence state, and non-finite events. Do not materialize the complete dense `S × P × B` expected tensor.

These values are the first frozen benchmark settings, not universal defaults. Any revision must be made jointly for both arms from synthetic unit cases or development seeds and must increment the protocol version before primary results are opened.

## 10. Required model invariants

An implementation is not conforming until these tests pass:

1. **Conservation:** layer sums equal collapsed counts exactly.
2. **Poisson equivalence:** with a Poisson total, the factorized model equals independent Poisson bin counts up to numerical tolerance.
3. **Identical shapes:** if all cell types have identical `omega`, the conditional likelihood contains no information about `z`, and both arms return the same estimate within tolerance.
4. **One bin:** the conditional term is constant zero and both arms agree.
5. **Shape-identifiable toy:** equal total signatures with opposing shape signatures are recoverable only by ShapeMix.
6. **Count-identifiable toy:** distinct total signatures are recoverable by both arms.
7. **Permutation negative control:** permuting shape signatures across cell types removes or reverses any legitimate shape benefit.
8. **Zero handling:** zero-total rows and zero observed bins never produce NaNs.
9. **Chunk parity:** chunked and unchunked objectives and gradients agree on a small problem.
10. **Order invariance:** reordering input cells before deterministic sorting does not change signatures, dispersion, or fitted outputs.

## 11. Versioned departures

The following changes are new models or sensitivity analyses and must not silently replace version 1:

- weighting fragments by `readSupport`;
- counting unique fragments instead of cut sites;
- changing length bins;
- independent-bin negative binomials;
- Dirichlet-multinomial shape overdispersion;
- observed depth offsets or library-adjusted signatures;
- fixed or learned background;
- peak-specific or learned dispersion;
- learned accessibility/shape signatures;
- spatial smoothing or variational inference.
