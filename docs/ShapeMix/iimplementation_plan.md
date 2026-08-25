# ShapeMix-ATAC implementation plan

Status: active implementation roadmap; Steps 0 through 6 completed 2026-08-23; all four planned external GEO source families acquired and source-preprocessed by 2026-08-24; full cross-dataset evaluation is planned in Section 13

This document turns the ShapeMix concept into a repository-specific engineering and research plan and records the execution status of each step. The implementation preserves all current datasets and methods, creates new ShapeMix-specific data products, and makes the shape-aware versus peak-only comparison a controlled, reproducible ablation.

## 1. Executive recommendation

The minimum credible ShapeMix implementation is:

1. Download the raw 10x PBMC ATAC fragments file and its tabix index.
2. Split labeled PBMC cells into disjoint reference and held-out test pools.
3. Select peaks using only the reference pool.
4. Count deduplicated Tn5 cut sites by peak and by the length of their parent fragment, using three length bins.
5. Store those counts as three sparse `AnnData` layers, with `.X` equal to their sum.
6. Aggregate only held-out test cells into pseudo-spatial spots.
7. Learn fixed cell-type accessibility and shape signatures from reference cells only.
8. Fit positive spot-by-cell-type abundances with a non-spatial PyTorch MAP model.
9. Compare the shape-aware model with the exact same model after omitting only the conditional shape likelihood.
10. Run paired simulations across nested split and mixture seeds and compare RMSE, Jensen–Shannon divergence, rare-cell recovery, runtime, and diagnostic checks.

Do not begin with variational inference, a spatial prior, motif footprints, learned reference signatures, or a large set of nuisance parameters. Those are gated extensions after the three-bin MAP model demonstrates that fragment length contributes useful information.

No existing tracked file or current dataset should be deleted. In particular, preserve the current PBMC simulations even though they are unsuitable as the primary ShapeMix benchmark because their source-cell pool is also used as the deconvolution reference.

The immediate priority after acquisition is now the full evaluation campaign in Section 13. It keeps the completed protocol-v1 PBMC result immutable, runs the already materialized GSE129785 evaluation first, builds leakage-free donor folds for GSE194122, and advances GSE205055/GSE263333 only after compatible labeled fragment-level references pass their own gates. GPU acceleration is an engineering prerequisite for that campaign, not a change to the statistical model.

## 2. Source documents and authority

The plan reconciles these local documents:

- [ShapeMix research proposal](../research_class/ShapeMix_ATAC_proposal_draft.md): latest experimental scope and primary hypothesis.
- [Current Bayesian tutorial source](tutorials/ShapeMix_ATAC_Bayesian_Model_Tutorial.tex): current formal model and priors.
- [Beginner tutorial](tutorials/ShapeMix_ATAC_Tutorial.pdf): biological intuition and broader roadmap.
- [High-level idea](ShapeMix_high_level.docx): original concept.
- `tutorials/ShapeMix_ATAC_Bayesian_Model_Tutorial_deprecated.docx`: historical rationale only; it should not override the current proposal or TeX source.

The latest proposal and TeX source agree on a three-bin, reference-guided, non-spatial MAP first version. Where the documents disagree, this plan uses the following precedence:

1. The latest research proposal defines the primary experiment.
2. The current Bayesian TeX defines model terminology.
3. This implementation plan resolves computational and statistical ambiguities required to make the experiment valid.
4. Older tutorials remain useful background but are not the executable specification.

Before Step 0, the source documents described `Y` as a fragment count. Step 0 reconciled them on 2026-08-22 around the primary count unit: **deduplicated Tn5 cut sites, each annotated with its parent fragment's length bin**. Summing the bins therefore reconstructs the same kind of peak count used by Cell Ranger and the count-only baselines. The scientific claim is: *does the parent-fragment length composition of observed insertions improve deconvolution beyond their peak totals?* A unique-fragment-count formulation would be a separately versioned sensitivity analysis with a matched fragment-count baseline, not an interchangeable input.

## 3. Repository constraints and resolved gaps

ShapeMix could not be added only as another method adapter. The work also required upstream data, evaluation, and execution changes:

- `DeconvolutionInput` now supports an optional validated fragment-shape contract while remaining backward-compatible with aligned two-dimensional inputs. The legacy PBMC objects still contain collapsed peak counts only; Step 3 added separately versioned split-specific objects with canonical fragment-length layers.
- Before Step 1, the required 10x fragments file, approximately 2.4 GB compressed, and its `.tbi` index were deferred. They are now downloaded, checksum-pinned, and tabix-validated. Step 2 added the reusable fragment counter and optional shape-layer contract, and Step 3 used them to produce the canonical split-specific PBMC layers.
- The legacy PBMC simulator samples cells from the same full object later used as the reference. Those datasets remain preserved for provenance, while Step 3 resolved the problem for ShapeMix by generating cell-disjoint reference and held-out pools.
- Heart and Russell data contain only collapsed peak matrices. Fragment length cannot be recovered from them, so they cannot test ShapeMix's core claim.
- Step 5 resolved the runner limitation by adding named method variants, so the shape-aware and count-only arms now run together with distinct IDs and fully resolved configuration provenance.
- Step 6 replaced intersection-based scoring for the primary benchmark with a strict declared-cell-type and exact-spot contract. Missing declared prediction types are penalized as zero, while unknown types, missing/extra spots, invalid values, and invalid row sums fail evaluation.
- Step 6 added canonical `rmse_v1` and squared base-2 Jensen–Shannon divergence `jsd_v2`, while retaining the historical unsquared value as `js_distance_v1`. It also added rare-cell metrics, nested-seed summaries, failure propagation, runtime/memory reporting, and hash-verified provenance. Calibrated posterior uncertainty remains outside the MAP MVP.

Steps 2–6 resolved the implementable gaps without replacing legacy data IDs or changing existing-method behavior. The one-donor source and lack of recoverable fragment length in collapsed-matrix datasets remain scientific constraints rather than software defects.

## 4. Decisions to freeze before coding

| Topic | MVP decision | Reason |
|---|---|---|
| Shape axis | Parent-fragment length | It is directly available and is the least sparse shape signal proposed. |
| Length definition | `chromEnd - chromStart` from the adjusted, half-open ARC 2.0 fragment interval | This is deterministic and keeps bin assignment independent of which end overlaps a peak. |
| Length bins | `[0,100)`, `[100,250)`, `[250,∞)` bp | This is the latest proposal/current Bayesian tutorial definition. Four or more bins are sensitivity analyses. |
| Count unit | Deduplicated Tn5 cut sites grouped by parent-fragment length | Cell Ranger ARC peak matrices count cut sites, so the bin layers can collapse to the same count unit used by existing baselines. |
| 10x `readSupport` | Ignore in the primary analysis | Each fragments-file row is one deduplicated fragment; using read support would reintroduce PCR duplicates. Retain it only as a sensitivity option. |
| Peak assignment | Assign each cut site to the unique non-overlapping Cell Ranger peak containing it | Cell Ranger peaks are sorted and non-overlapping. Validate the right-end coordinate convention against the vendor matrix before freezing canonical data. |
| Peak selection | Deterministic reference-only accessibility-variance ranking; 5,000 peaks for the first benchmark | The scoring formula, coverage filter, and tie-breaking are frozen below. It never uses shape bins or held-out cells. Add 10,000/20,000-peak sensitivity runs later. |
| Reference/test split | Stratified 70/30 cell split with outer split seeds | The available PBMC data contain one donor, so donor-level separation is impossible. Both sets must retain enough cells per included type. Mixture seeds are nested within split seeds. |
| Storage | Sparse `AnnData` layers, one cell/spot-by-peak CSR matrix per bin; `.X = sum(layers)` | Existing feature slicing and count-only methods remain compatible. |
| Primary likelihood | Total-count negative binomial plus conditional multinomial shape likelihood | This creates an exactly nested shape-versus-count ablation. |
| Inference | Fixed signatures; infer only positive abundance with PyTorch MAP | It is debuggable and avoids the major identifiability problems in the full model. |
| Depth and abundance | No observed spot-depth offset; total depth is absorbed into `z` | `abundance.csv` is effective reference-cell-equivalent, depth-scaled abundance, not guaranteed nucleus counts. Row-normalized `z` is the proportion output. |
| Dispersion | One reference-cross-fitted global inverse-dispersion, scaled by each spot's effective abundance | This is an exact, reproducible MVP estimator and approximates aggregation of independent reference-depth cells. Peak-specific dispersion is a later sensitivity. |
| Background | Zero for the primary same-dataset simulation; fixed small background as a robustness check | A learned background can absorb cell-type signal and obscure the shape test. |
| Spatial smoothing | Excluded | It would confound the first test of fragment shape. |
| Uncertainty | No credible intervals in the MAP MVP | Report point estimates and convergence diagnostics honestly; add VI only after the MAP result is credible. |
| Metric universe | Ordered cell types declared by the dataset; mean base-2 Jensen–Shannon **divergence** | Omitted declared types are filled with zero, unknown extras are rejected, and the fixed denominator is identical across methods. |

### Length-bin interpretation and boundary tests

The three bins describe the length of the parent ATAC fragment—the distance between its two Tn5 cut sites—not the length of a peak or a separate measurement at each cut site:

- **Short, `[0,100)` bp:** fragments below 100 bp often arise from relatively exposed or nucleosome-free DNA, where Tn5 can cut nearby positions.
- **Nucleosome-scale, `[100,250)` bp:** fragments in this broad interval are compatible with DNA spanning approximately one nucleosome. About 147 bp of DNA wraps around the histone core; linker DNA, cleavage positions, and assay variation broaden the observed fragment length beyond exactly 147 bp.
- **Long, `[250,∞)` bp:** these fragments are more likely to span multiple nucleosomes or other larger protected chromatin structures.

These names are biological interpretations, not definitive classifications of individual fragments. A short fragment can reflect transcription-factor protection or technical effects, and a long fragment does not prove a particular nucleosome count. The exact 100 bp and 250 bp thresholds come from the ShapeMix proposal and current Bayesian tutorial; they are frozen rather than estimated from this PBMC dataset so held-out evaluation outcomes cannot influence their definition. Alternative thresholds or additional bins belong in explicitly versioned sensitivity analyses.

The acceptance-test lengths `99`, `100`, `249`, and `250` bp are chosen specifically to test the half-open boundaries. They must map to `short`, `nucleosome-scale`, `nucleosome-scale`, and `long`, respectively. This pairwise testing catches an inclusive-versus-exclusive or other off-by-one implementation error at both thresholds.

The motivation for retaining these bins is that equal total accessibility can conceal different fragment-length compositions. For example, two cell types might each contribute 100 cut sites at a peak, while one contributes `80/15/5` and the other `20/50/30` across the short/nucleosome-scale/long bins. A peak-only model sees the same total of 100; ShapeMix can use the conditional bin composition as additional evidence.

The version-matched Cell Ranger ARC 2.0 documentation says that each fragments-file row is a unique fragment, the fifth column is read-pair support after duplicate collapsing, and the five-column file is BGZF/tabix indexed: [Cell Ranger ARC 2.0 ATAC fragments file](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/analysis/fragments-file). It also states that peak-barcode matrix entries count cut sites in peaks: [Cell Ranger ARC 2.0 feature-barcode matrices](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/analysis/feature-barcode-matrices).

### Coordinate validation result

A fragment has an adjusted left and right cut site. Step 2 reconstructed the official peak matrix under both plausible right-cut mappings. On the deterministic primary sample of 2,048 peaks and 512 barcodes, `start + end` reproduced all 158,472 official counts with zero mismatches; `start + end - 1` produced 58 mismatched entries and absolute error 59. An independent all-peak chr21 audit again matched exactly only with `end`. The frozen convention is therefore `left = chromStart`, `right = chromEnd`, `right_cut_offset = 0`, with half-open peak membership. The complete comparison statistics and raw-file hashes are recorded in the tracked source manifest.

## 5. Canonical MVP model

### 5.1 Data objects

Let:

- `s = 1,...,S` index spots;
- `p = 1,...,P` index selected peaks;
- `b = 1,...,B` index fragment-length bins;
- `c = 1,...,C` index cell types;
- `Y[s,p,b]` be observed cut-site counts whose parent fragments fall in bin `b`;
- `N[s,p] = sum_b Y[s,p,b]` be the ordinary peak count;
- `A[c,p]` be the fixed mean total cut-site rate per reference cell;
- `omega[c,p,b]` be the fixed conditional length-bin distribution, summing to one over `b`;
- `z[s,c] > 0` be the effective abundance to infer.
- `e[s] = sum_c z[s,c]` be the spot's effective number of reference-depth cells;
- `phi_ref` be the global per-reference-cell negative-binomial inverse-dispersion estimated only from training-reference cells.

For an optional fixed background, let `beta[p]` be its total rate and `u[p,b]` its shape distribution. The primary simulated benchmark uses `beta = 0`.

### 5.2 Factorized likelihood

Define:

```text
v[s,p,b] = beta[p] * u[p,b]
             + sum_c z[s,c] * A[c,p] * omega[c,p,b]

mu[s,p]  = sum_b v[s,p,b]
           = beta[p] + sum_c z[s,c] * A[c,p]

rho[s,p,b] = v[s,p,b] / mu[s,p]
```

Fit:

```text
N[s,p] ~ NegativeBinomial(
             mean=mu[s,p],
             dispersion=max(e[s], epsilon) * phi_ref
         )

Y[s,p,1:B] | N[s,p]
    ~ Multinomial(total_count=N[s,p], probabilities=rho[s,p,1:B])
```

The negative-binomial parameterization is:

```text
E[N]   = mu
Var[N] = mu + mu^2 / (e * phi_ref)
```

The primary objectives are therefore:

```text
peak-only:
    log p(N | z, A, phi_ref) + log p(z)

ShapeMix:
    log p(N | z, A, phi_ref)
  + log p(Y | N, z, A, omega)
  + log p(z)
```

Everything except the conditional shape term must be identical between the two arms: data split, peaks, total counts, signatures, parameters, priors, initialization, optimizer, seed, stopping rule, and compute budget.

Implement the conditional shape term directly in log space rather than constructing a `torch.distributions.Multinomial` object for every spot/peak with a different total:

```text
log p(Y | N, rho)
  = lgamma(N + 1) - sum_b lgamma(Y[b] + 1)
    + sum_{b:Y[b]>0} Y[b] * log(rho[b])
```

Compute `log(rho)` with a stable log-sum-exp normalization. Entries with `Y=0` contribute exactly zero, and rows with `N=0` are explicitly masked to a zero conditional contribution. Tests must cover zero totals and extremely small probabilities without NaNs.

This factorization is preferable to the tutorial's independent negative-binomial likelihood for every bin. Independent negative binomials generally do not collapse into the same negative-binomial total-count model, so changing from one bin to three would change both the information and the stochastic weighting. Retain independent-bin NB only as a labeled sensitivity model.

Two valuable invariants follow:

- In the zero-background MVP, if all cell types have the same `omega[c,p,:]`, the shape likelihood contains no information about `z`; ShapeMix and peak-only estimates should agree. With fixed background, the same invariant requires its `u[p,:]` to equal that shared shape.
- With a Poisson total model, the factorization is algebraically equivalent to independent Poisson counts per bin. Use that equivalence as a golden implementation test before enabling NB totals.

PyTorch supplies the differentiable tensor operations and negative-binomial distribution needed for MAP; implement the conditional term manually as specified above: [PyTorch distributions](https://docs.pytorch.org/docs/stable/distributions.html).

### 5.3 Fixed reference signatures

Let `R[i,p,b]` be the training-reference shape count for cell `i`, and `T[i,p] = sum_b R[i,p,b]`. Estimate signatures only from reference-split cells.

For cell type `c`, use smoothed empirical estimates:

```text
A[c,p] = (sum_{i in c} T[i,p] + alpha_A * g[p])
         / (sum_{i in c} h[i] + alpha_A)

omega[c,p,b]
    = (sum_{i in c} R[i,p,b] + alpha_omega * u[p,b])
      / (sum_{i in c} T[i,p] + alpha_omega)
```

Here `g[p]` is a pooled peak-rate target, `u[p,b]` is the pooled training-reference shape distribution for the peak after smoothing it toward the dataset-global bin distribution with the same `alpha_omega`, and `h[i]` is an exposure factor. Require every global bin to have positive training-reference support. Use `h[i] = 1` for the primary same-protocol pseudo-spot benchmark so `A` is an average per reference cell; evaluate library-size-adjusted exposures as a sensitivity analysis.

When a peak is too sparse to estimate a cell-type-specific `omega`, shrink toward the pooled peak distribution, then toward the dataset-global distribution. Do not use unscaled raw reference counts as a later Dirichlet concentration: abundant cell types would become artificially overconfident. Record effective smoothing concentrations explicitly.

The MVP absorbs the tutorial's `q[c]`, `d[p]`, and `theta[c,p]` into fixed `A[c,p]`. Estimating all of `z`, `q`, `d`, and `theta` jointly is scale-nonidentifiable and is not needed to test the shape hypothesis.

### 5.4 Dispersion, depth, and abundance semantics

Use a deterministic, cell-type-stratified two-fold cross-fit inside the training-reference pool to estimate the one MVP overdispersion value. For each reference cell `i`, estimate its held-out mean `m[i,p]` from the opposite fold's cell-type signature, then compute:

```text
alpha_ref_raw
  = sum_{i,p} ((T[i,p] - m[i,p])^2 - m[i,p])
    / sum_{i,p} m[i,p]^2

alpha_ref = max(alpha_ref_raw, 1e-8)
phi_ref   = 1 / alpha_ref
```

Fail if the denominator or result is non-finite, and record the fold seed, numerator, denominator, `alpha_ref`, and `phi_ref`. The production signature `A` is then re-estimated from all training-reference cells. At fit time use spot dispersion `e[s] * phi_ref`, where `e[s] = sum_c z[s,c]`. This is a global MVP approximation; per-peak or learned dispersion is a later, explicitly labeled model.

Do not supply observed spatial library size as an offset in the primary model. `z` absorbs total peak depth relative to the average training-reference cell as well as mixture abundance. Consequently, `abundance.csv` contains **effective reference-cell-equivalent, depth-scaled abundance** and must not be interpreted as a calibrated nucleus count. The standardized scientific output is `pi[s,c] = z[s,c] / sum_c z[s,c]`.

### 5.5 MAP inference

- Parameterize `z = softplus(raw_z) + epsilon`.
- Initialize with collapsed-count NNLS and a uniform fallback.
- Use multiple deterministic restarts.
- Optimize with Adam; optionally finish the best restart with LBFGS after tests establish stability.
- Fix `phi_ref` with the cross-fitted estimator above; do not tune it on held-out mixtures. Peak-specific or learned dispersion is later work.
- Process spots and peaks in chunks; never construct the full dense `S × P × B` expected tensor.
- Fail on non-finite loss, gradient, abundance, or proportions rather than allowing normalization to turn failures into zeros.
- Return `pi[s,c] = z[s,c] / sum_c z[s,c]` as the standardized proportion estimate.

Step 1 added a tested ShapeMix optional extra in `pyproject.toml` and validated PyTorch 2.11 and pysam 0.24 in the project environment. Pyro may remain recorded in the complete `requirements.txt` environment snapshot, but it is not an MVP ShapeMix dependency; add it to the ShapeMix extra only when variational inference is implemented. Pyro provides automatic guides for that later stage: [Pyro automatic guides](https://docs.pyro.ai/en/stable/infer.autoguide.html).

## 6. Data and file flow

```text
10x fragments + index + labels + Cell Ranger peaks
                         |
                         v
        stratified split and reference-only peak selection
                         |
              +----------+----------+
              |                     |
              v                     v
  reference cell shape H5AD   held-out cell shape H5AD
              |                     |
              |                     v
              |           layer-wise pseudo-spot sums
              |                     |
              +----------+----------+
                         v
         registered ShapeMix dataset (.X plus 3 layers)
                         |
              +----------+----------+
              |                     |
              v                     v
        peak-only model       ShapeMix model
              |                     |
              +----------+----------+
                         v
          paired metrics, diagnostics, and ablations
```

### Shape-aware `AnnData` contract

Suggested layer names:

```text
fragment_length_lt_100
fragment_length_100_249
fragment_length_ge_250
```

Required invariants:

- all layers are nonnegative integer CSR matrices;
- all layers have the same observation and peak axes as `.X`;
- `.X` equals their elementwise sum exactly;
- reference and spatial objects have identical ordered peaks;
- reference observations have the declared cell-type label;
- spatial observations have coordinates and all truth rows;
- bin edges and layer order agree between YAML metadata and both H5AD objects;
- `.uns["fragment_shape"]` records the count unit, bin edges, cut-site convention, read-support policy, source hashes, feature hash, and split hash.

Proposed modality metadata:

```yaml
fragment_shape:
  schema_version: 1
  axis: parent_fragment_length_bp
  count_unit: deduplicated_cut_sites
  read_support_policy: ignore
  peak_assignment: containing_nonoverlapping_peak
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
```

Sparse H5AD layers are appropriate because AnnData stores CSR/CSC matrices on disk and retains aligned layers when observations or variables are sliced: [AnnData sparse on-disk format](https://anndata.readthedocs.io/en/stable/fileformat-prose.html).

## 7. Phased implementation

### Step 0 — Freeze the scientific contract

Goal: create one executable specification before model code or canonical data are produced.

Execution status: **completed 2026-08-22**. The canonical README, model specification, and benchmark protocol were created; the proposal and Bayesian TeX were reconciled; and the 18-page PDF was regenerated successfully from the updated TeX.

| Action | File | Planned change |
|---|---|---|
| Add | `docs/ShapeMix/README.md` | Point readers to the canonical specification and implementation plan; explicitly label the high-level DOCX and beginner/deprecated tutorials as historical, non-normative context. |
| Add | `docs/ShapeMix/model_specification.md` | Record the cut-site-by-parent-length tensor semantics, likelihood, cross-fitted dispersion, priors, signature construction, bin boundaries, and effective-abundance interpretation. |
| Add | `docs/ShapeMix/benchmark_protocol.md` | Predeclare datasets, nested split/mixture seeds, endpoints, metric versions, rare-cell definitions, and statistical comparisons. |
| Modify | `docs/ShapeMix/tutorials/ShapeMix_ATAC_Bayesian_Model_Tutorial.tex` | Change `Y` from ambiguous fragment counts to cut-site counts annotated by parent-fragment length, state the resulting scientific claim, and replace or qualify independent-bin NB with the factorized likelihood. |
| Modify | `docs/ShapeMix/tutorials/ShapeMix_ATAC_Bayesian_Model_Tutorial.pdf` | Regenerate from the updated TeX; do not edit the PDF manually. |
| Modify | `docs/research_class/ShapeMix_ATAC_proposal_draft.md` | Reconcile fragment-versus-cut-site terminology, state the parent-length-composition claim, clarify the nested likelihood, and record the one-donor split limitation. |
| Delete | None | Retain all earlier documents as provenance; label the deprecated document rather than removing it. |

Acceptance criteria:

- Three-bin semantics, cut-site count unit, parent-length annotation, duplicate policy, reference/test split, factorized likelihood, and primary endpoints are unambiguous.
- The peak-only arm is defined as the same model with only the conditional shape term disabled.
- The benchmark protocol is frozen before inspecting results from the final outer split and inner mixture seed lists.
- The README identifies the proposal, current TeX/PDF, model specification, and benchmark protocol as canonical and all older binary tutorials as historical/non-normative.

### Step 1 — Acquire and pin raw fragment inputs

Goal: make fragment-level data reproducible on any machine.

Execution status: **completed 2026-08-22**. The 2,403,785,496-byte fragments file (`SHA-256 5075e32a...e3b5e1e38`) and 1,089,534-byte tabix index (`SHA-256 3a516291...60a5e7f5`) were downloaded from the versioned 10x URLs. Full gzip integrity passed; pysam 0.24.0 opened all 39 indexed contigs and validated records from the first, primary sex-chromosome, and terminal contigs. Exact hashes are stored in the tracked and local manifests.

| Action | File | Planned change |
|---|---|---|
| Produce, git-ignored | `data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz` | Download the official approximately 2.4 GB BGZF file. |
| Produce, git-ignored | Same directory, `pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi` | Download the official tabix index. |
| Modify, git-ignored | Same directory, `manifest.yaml` | Move both files from deferred to required ShapeMix inputs and record byte sizes, SHA-256 hashes, URLs, and download date. |
| Add | `configs/data_sources/pbmc_granulocyte_sorted_10k_cellranger_arc_2.0.0.yaml` | Track the versioned official URLs, filenames, expected byte sizes, SHA-256 hashes, and five-column schema in Git; the ignored local manifest records resolved paths and download state. |
| Modify | `docs/recreate_data_directory (important).md` | Add resumable download, checksum, and validation instructions for the two ShapeMix inputs. |
| Modify | `pyproject.toml` | Add a `shapemix` optional dependency group with tested PyTorch and `pysam` bounds. Do not add Pyro until the VI phase. |
| Modify later | `requirements.txt` | Regenerate the complete tested environment snapshot and include `pysam`. Pyro may remain in that snapshot but is not an MVP ShapeMix dependency. |
| Delete | None | The tracked source manifest and recreation guide preserve source history; the existing ignored local manifest is updated in place. |

Acceptance criteria:

- Both files pass recorded checksums.
- `pysam.TabixFile` opens the file with the adjacent index and reports expected contigs.
- The Cell Ranger ARC 2.0 five-column schema—chromosome, start, end, barcode, and `readSupport`—and its coordinate semantics are recorded exactly; no strand column is expected for this input version.
- Disk-space guidance in the recreation document includes raw and derived ShapeMix data.

### Step 2 — Add the fragment-shape data contract and counter

Goal: turn raw fragments into validated sparse layers without making old datasets invalid.

Execution status: **completed 2026-08-22**. The streaming counter, typed bins, ranked-output/non-overlapping peak index, bounded COO accumulation with balanced CSR merges, H5AD builder, optional loader/schema integration, provenance validation, and legacy-compatible opt-in checks are implemented. Ordered feature hashes are bound to the current feature axis. Immutable source-run QC is kept separate from matrix totals, which are checked against the stored layers and recomputed after loader feature slicing. The focused and full suites pass, and a real-data 200-barcode by 200-peak integration check reconstructed all 14,415 official counts with zero mismatches. Canonical split-specific H5AD production remains Step 3.

| Action | File | Planned change |
|---|---|---|
| Add | `src/deconvatac/pp/fragment_shapes.py` | Typed bin parsing, tabix streaming, cut-site-to-peak assignment, chunked sparse accumulation, provenance counters, and H5AD helpers. |
| Modify | `src/deconvatac/pp/__init__.py` | Export reusable fragment-shape preprocessing functions. Do not alter `reads_to_fragments.py`; it cannot recover length and remains a legacy helper. |
| Modify | `src/deconvatac/data/schemas.py` | Add optional backward-compatible `FragmentShapeSpec` and ordered `cell_types` metadata to `DeconvolutionInput`, both defaulting to `None` for old datasets. |
| Modify | `src/deconvatac/data/loaders.py` | Parse optional `fragment_shape` metadata and an ordered `truth.cell_types` universe; normal datasets remain unchanged. |
| Modify | `src/deconvatac/data/validators.py` | Add shape-layer, integer, alignment, bin, ordered-feature-hash, provenance-counter, `.X == sum(layers)`, matrix-total, and declared-cell-type validation. Invoke shape checks only when shape metadata are declared. |
| Modify | `src/deconvatac/data/__init__.py` | Export the new optional contract and validator. |
| Add | `tests/test_fragment_shapes.py` | Test parsing, boundaries, cut-site coordinates, barcode filtering, duplicate policy, sparse accumulation, and determinism. |
| Add | `tests/test_shape_data_contract.py` | Test valid layers, missing/mismatched layers, negative/noninteger values, metadata mismatches, ordered cell types, feature slicing, and old-input compatibility. |
| Delete | None | No current loader or data schema is removed. |

The counter must emit explicit counts for total rows, header rows, invalid coordinates, unknown barcodes, filtered contigs, cut sites outside selected peaks, cut sites per bin, and read-support totals. These are immutable source-run QC. A separate matrix-counter block must remain synchronized with the currently stored object and equal its layer and collapsed-matrix totals. The implementation accumulates bounded COO chunks and merges CSR runs in balanced binary levels rather than retaining the complete fragments table in memory or repeatedly adding every chunk to one growing matrix.

Import `pysam` inside the fragment-reading functions rather than at module import time. The top-level package imports `deconvatac.pp`, so an eager optional import would otherwise break every method in environments without the ShapeMix extra.

Acceptance criteria:

- Lengths 99, 100, 249, and 250 bp enter the expected bins.
- Unknown barcodes and invalid fragments are counted and excluded deterministically.
- Summed bin layers equal `.X` exactly.
- Reconstructed `.X` agrees with the official selected peak matrix after the cut-site coordinate convention is fixed. Any remaining discrepancy is quantified and explained.
- Peak/cell ordering and output hashes are stable across repeated runs.

### Step 3 — Build the first PBMC ShapeMix benchmark with held-out evaluation cells

Goal: create reference signatures and pseudo-spots from cell-disjoint reference and held-out evaluation pools.

This step is deliberately PBMC-specific because it is the first source for which the repository has a pinned raw fragments file, a matching tabix index, an official peak matrix, and an executable cell-label mapping. Existing Heart and Russell peak matrices cannot recover parent-fragment lengths. They need their own raw fragment-level inputs, coordinate audit, labels, donor-aware split, and versioned benchmark protocol before they can become ShapeMix datasets; they must not be manufactured by copying the PBMC assumptions.

| Action | File | Planned change |
|---|---|---|
| Add | `scripts/prepare_shapemix_pbmc.py` | Create the stratified split, select reference-only peaks, build reference/held-out shape-layer H5ADs, and write manifests atomically. |
| Add | `scripts/regenerate_shapemix_pbmc_simulations.py` | Sum held-out cells layer by layer into the frozen equal-celltype and observed-abundance pseudo-spots; keep rare-cell, subtype, and depth variants as separately identified later sensitivities. |
| Add | `scripts/audit_shapemix_signal.py` | Report coverage, shape entropy, split-half fingerprint reproducibility, between-cell-type divergence, and technical confounding. |
| Add | `scripts/shapemix_provenance.py` | Share code/protocol hashes, Git state, package versions, shape semantics, dimensions, and nonzero-count summaries across generated manifests. |
| Add | `scripts/validate_shapemix_step3.py` | Recheck registered artifacts, hashes, axes, source-cell membership, exact truth, and layer-wise aggregation after production. |
| Modify | `src/deconvatac/pp/feature_selection.py` | Add a deterministic reference-only peak ranker with explicit coverage, score, tie-breaking, and ordered output; preserve current selectors' behavior. |
| Modify | `tests/test_feature_selection.py` | Test the exact score, minimum coverage, stable ties, input-order invariance, training-only use, and ordered output. |
| Add | `tests/test_shapemix_pbmc_preparation.py` | Test stratification, no overlap, reference-only feature selection, manifests, hashes, and deterministic output ordering. |
| Add | `tests/test_shapemix_simulation.py` | Test layer-wise aggregation, exact truth, source-cell provenance, depth thinning, and seed independence. |
| Produce, git-ignored | `data/processed/shapemix/pbmc_granulocyte_sorted_10k/split_<split_seed>/` | Split-specific reference/test H5ADs, split table, selected peaks, QC summaries, signal audit, and manifest without cross-seed collisions. |
| Produce, git-ignored | `data/processed/datasets/pbmc_granulocyte_sorted_10k_shapemix_<condition>_split_<split_seed>_mix_<mixture_seed>/` | Registered spatial H5AD, truth, source-cell JSONL, dataset YAML, and simulation manifest per nested seed pair. |
| Modify, git-ignored | `data/registry/datasets.yaml` | Register generated ShapeMix dataset IDs. |
| Delete | None | Do not overwrite the existing PBMC reference or either existing PBMC simulation. |

Freeze this peak selector in `benchmark_protocol.md` and implement it as a new function without changing the legacy selectors:

1. Use only the training-reference collapsed matrix `.X` and its cell-type labels—never shape layers or held-out cells.
2. Require a peak to be nonzero in at least 10 training cells (`min_reference_cells: 10`).
3. For every cell type `c`, compute `L[c,p] = log2(1 + 10^4 * K[c,p] / sum_p K[c,p])`, where `K[c,p]` is the summed training count. Fail if a cell type has zero total count.
4. Score each eligible peak as `variance_c(L[:,p])` using population variance (`ddof=0`).
5. Sort by descending score, then descending number of nonzero training cells, then descending total training count, then ascending peak ID. Take the first `P` and preserve that ranked order in the feature file.

The dataset YAML must declare one ordered `truth.cell_types` list. It must exactly match the reference label universe and truth columns; pseudo-spots may have zero abundance for a declared type, but the type cannot disappear from the schema.

Suggested intermediate tree:

```text
data/processed/shapemix/pbmc_granulocyte_sorted_10k/split_000/
  reference_cells.h5ad
  heldout_test_cells.h5ad
  split.csv
  selected_peaks.txt
  peak_selection.csv
  signal_audit.csv
  signal_audit_summary.yaml
  manifest.yaml
```

Suggested dataset tree:

```text
data/processed/datasets/pbmc_granulocyte_sorted_10k_shapemix_equal_split_000_mix_000/
  atac/spatial.h5ad
  atac/features/highly_variable.txt
  truth/proportions.csv
  simulation/source_cells_by_spot.jsonl
  simulation/manifest.yaml
  dataset.yaml
```

The implemented smoke dataset uses 32 spots, 200 peaks, three cell types, and one split/mixture seed pair. Each primary dataset uses 1,024 spots, approximately 10 cells per spot, 5,000 peaks, all 16 retained cell types, and either equal-cell-type or observed-abundance sampling. The frozen grid has five outer split seeds with two independently derived mixture seeds inside each split. Broader depth, peak-count, and cells-per-spot grids remain separately versioned, post-primary gates.

Acceptance criteria:

- Reference and held-out barcode sets are disjoint.
- Peak rankings are identical under input-row reordering and use only training-reference totals/labels.
- Every pseudo-spot source barcode belongs to the held-out pool.
- Feature selection and signature audit never read held-out labels/counts except for final evaluation.
- `.X` and every layer of a pseudo-spot equal the corresponding sums of recorded source cells, adjusted only by a declared depth-thinning operation.
- The one-donor limitation is present in the dataset manifest and benchmark report.
- Split and mixture seeds are stored separately; their results are described as conditional resampling variability, not donor or biological generalization uncertainty.
- The signal audit is reviewed before expensive model development. If per-peak shape is too sparse, add hierarchical shrinkage or peak-group shape signatures; do not silently change the primary representation.

Execution status: **completed 2026-08-22**. The 9,500-cell raw-fragment pass produced one validated 9,500 × 5,754 official-order union cache. Five primary splits (`1103`, `2203`, `3301`, `4409`, and `5501`) each contain 6,644 reference cells, 2,856 held-out cells, and a separately ranked 5,000-peak axis. The development smoke split uses seed 0, the first three canonical types, 3,699 reference cells, and 200 independently ranked peaks. One 32-spot smoke dataset and all 20 frozen primary datasets (five outer seeds × two inner seeds × two conditions) were generated and registered; no existing PBMC dataset was overwritten or deleted.

All 21 dataset configurations round-trip through the maintained loader with their declared axes and truth universes. The read-only campaign validator hashed 3.89 GiB of unique artifacts and exactly reconstructed every stored layer, collapsed matrix, and truth count from 204,670 held-out source-cell assignments. The generated ShapeMix cache/split tree occupies about 348 MiB and the 21 dataset trees occupy about 541 MiB, approximately 890 MiB together. The reference-only signal review supports retaining the three-bin per-peak representation with the already planned hierarchical smoothing: across the five primary splits, selected peaks have a median of 2,973 cut sites, split-half shape Spearman correlation has median 0.747 and range 0.641–0.799, and median per-peak between-type generalized Jensen–Shannon divergence is 0.077–0.079 bits. Depth-versus-bin-fraction correlations are moderate (about 0.22–0.35 in absolute value), so technical confounding remains an explicit diagnostic and external-data priority rather than a claim of causal nucleosome occupancy. Rare-cell, subtype, depth, alternate-peak-count, Heart, Russell, and external-donor datasets were not built in this step; they remain separately versioned sensitivity or generalization work.

### Step 4 — Implement and test the fixed-signature MAP model

Goal: implement one core code path that supports the nested peak-only and shape-aware objectives.

| Action | File | Planned change |
|---|---|---|
| Add | `src/deconvatac/shapemix/__init__.py` | Public ShapeMix core API. |
| Add | `src/deconvatac/shapemix/config.py` | Typed configuration and rejection of unknown/misspelled parameters. |
| Add | `src/deconvatac/shapemix/signatures.py` | Estimate smoothed `A` and `omega`, cross-fit the global reference dispersion, preserve stable axis order, and hash the result. |
| Add | `src/deconvatac/shapemix/likelihood.py` | Poisson golden likelihood, abundance-scaled NB total likelihood, stable manual conditional-multinomial log likelihood, priors, and component diagnostics. |
| Add | `src/deconvatac/shapemix/map.py` | Positive parameterization, initialization, deterministic restarts, chunked optimization, convergence handling, and abundance/proportion extraction. |
| Add | `src/deconvatac/shapemix/diagnostics.py` | Reconstruction, residual, coverage, convergence, and multiple-start summaries without saving the full expected tensor. |
| Add | `tests/test_shapemix_signatures.py` | Aggregation, smoothing, cross-fitted dispersion formula, stable cell-type order, finite zero-count behavior, and probability normalization. |
| Add | `tests/test_shapemix_likelihood.py` | NB mean/variance and abundance scaling, factorization, zero totals, tiny probabilities, identical-shape invariant, permutation invariance, and Poisson equivalence. |
| Add | `tests/test_shapemix_map.py` | Known-mixture recovery, nonnegativity, finite gradients, reproducibility, chunk parity, and failure handling. |
| Delete | None | Do not create a second independent peak-only implementation; both variants must share the same model code. |

Essential synthetic tests:

1. Two cell types have equal peak totals but opposite length distributions. ShapeMix should recover the mixture; peak-only should be unidentifiable.
2. Cell types have distinct peak totals. Both variants should recover the mixture.
3. All cell types have identical `omega`. ShapeMix and peak-only should return the same estimates within tolerance.
4. Shape signatures are permuted across cell types. Any gain should disappear or reverse.
5. One bin produces the peak-only likelihood plus a constant conditional term.
6. Chunked and unchunked objectives/gradients agree on a small dense example.
7. All-zero peaks and spots either follow a documented policy or fail clearly; they never create NaNs.

The implementation freezes the remaining optimizer details before primary runs. Every restart starts from collapsed-count NNLS with each abundance clipped to `0.05`; an all-zero or invalid NNLS result uses a uniform positive fallback. Restart `r` receives the multiplicative perturbation `exp(N(0, 0.20))` from `PCG64(SeedSequence([20260822, outer_split_seed, inner_mixture_seed, 0, r]))`. Convergence is reached when the complete objective fails to improve its anchor by `tolerance * max(1, abs(anchor))` for `patience` consecutive steps, or when the full-gradient L2 norm is at most `tolerance`. Reaching `max_steps` alone is not convergence. Each restart retains its best finite state, and the fit selects the largest complete objective only among finite, converged restarts.

Acceptance criteria:

- All invariants above pass on CPU in a small test suite.
- Fixed seeds produce estimates within declared numerical tolerances.
- Loss components are reported separately as count likelihood, shape likelihood, abundance prior, and total objective.
- The optimizer records convergence state, stopping reason, restart selection, and all non-finite events.
- No full dense `S × P × B` tensor is materialized.

Execution status: **completed 2026-08-22**. The strict model configuration, immutable fixed-signature estimator, opposite-fold-only dispersion cross-fit, stable factorized likelihood, compact diagnostics, and streamed deterministic MAP fitter are implemented in one shared core. Dense toy arrays and ordered sparse bin layers use the same fitting path; only the current spot/peak chunk is densified, and expected bin rates are never materialized for the complete dataset. Both arms use identical reference hashes, NNLS initialization, restart streams, optimizer behavior, and total-count likelihood; `use_shape: false` removes only the conditional multinomial term.

The 41 focused Step 4 tests cover formula-level signature aggregation, hierarchical smoothing, fold and input-order determinism, shared ablation hashes, NB parameterization, Poisson factorization, exact zero/tiny-probability handling, one-bin and identical-shape controls, chunked objective/gradient parity, shape- and count-identifiable recovery, permuted-signature failure, sparse inputs, deterministic restarts, zero-spot fallback, and hard optimizer failures. The complete repository suite at the end of Step 4 passed with 138 tests and one pre-existing skip. A read-only real-data check on the 3,699-cell by 200-peak development reference produced `A` with shape `3 × 200`, `omega` with shape `3 × 200 × 3`, `phi_ref = 0.8008838701125459`, fold hash `5ab4bd516d9940cd0b6bd826610e05d15a5c7ca8a79c062b61377d4698895c0c`, and signature hash `3f642a07f02ddc60d6e25daf7e2745c63b51853b08b224a9ffd55d7b4d9d80a7`. Both frozen arms then converged through the sparse-layer core on the 32-spot smoke dataset. This was an implementation check, not a primary benchmark result; maintained-runner registration and standardized outputs were completed in Step 5.

### Step 5 — Integrate ShapeMix with the maintained runner

Goal: make ShapeMix a normal registered method with standardized outputs and named ablation variants.

| Action | File | Planned change |
|---|---|---|
| Add | `src/deconvatac/methods/shapemix.py` | ATAC-only adapter; validate shape input, call the MAP core, write compact native diagnostics, and return `DeconvolutionResult`. |
| Modify | `src/deconvatac/methods/registry.py` | Register only `shapemix`; keep optional imports lazy. |
| Add | `configs/methods/shapemix.yaml` | Shape-aware MAP configuration with `use_shape: true`. |
| Add | `configs/methods/shapemix_count_only.yaml` | Identical configuration with `use_shape: false`. |
| Add | `configs/experiments/shapemix_smoke.yaml` | Small shape-ready dataset and both named method runs. |
| Modify | `scripts/run_deconvolution.py` | Add the fully specified backward-compatible `method_runs` schema below, include the variant in jobs/run IDs/manifests, reject collisions before execution, and record seed, Torch, pysam, device, and determinism metadata. |
| Modify | `tests/test_method_interface.py` | Assert ShapeMix registration and standardized return types. |
| Modify | `tests/test_run_deconvolution_experiment.py` | Test two named configs for one method, unique run directories, metadata, and compatibility with the old `methods:` list syntax. |
| Add | `tests/test_shapemix_adapter.py` | Test valid outputs, missing layers, RNA rejection, non-finite rejection, and absence of accidental large artifacts. |
| Delete | None | Do not add `shapemix_count_only` as a second registry method solely to work around the runner schema. |

Proposed experiment syntax:

```yaml
method_runs:
  - id: shapemix_length
    method: shapemix
    config: configs/methods/shapemix.yaml
  - id: shapemix_count_only
    method: shapemix
    config: configs/methods/shapemix_count_only.yaml
```

Runner semantics must be explicit and tested:

- An experiment defines exactly one of `methods` or `method_runs`. Defining both, neither, duplicate variant IDs, or duplicate resolved run IDs is an error before any job starts.
- Every `method_runs` entry has a unique `id`, a registry `method`, and either a `config` path/inline mapping or the method's default config. Validate that the config's `method` matches the registry method.
- Legacy `methods: [nnls, ...]` is internally expanded to method runs whose `id == method` and keeps the existing `method_configs`/default-config resolution.
- Every job dictionary contains both `method` and `method_run_id`. The default run template becomes `{dataset}__{modality}__{feature_set}__{method_run_id}`; because legacy IDs equal method names, legacy default run IDs remain unchanged. Custom templates may use either placeholder.
- `run.yaml`, `runs.csv`, `failures.csv`, and `comparison.csv` record `method` as the registry implementation and `method_run_id` as the variant. The fully resolved config and its source path/hash are also recorded.

Proposed primary method configuration:

```yaml
method: shapemix
params:
  use_shape: true
  total_likelihood: negative_binomial
  conditional_shape_likelihood: multinomial
  signature_rate_pseudocount: 0.5
  signature_shape_concentration: 1.0
  exposure_mode: absorbed_in_abundance
  dispersion_mode: reference_crossfit_global_scaled_by_abundance
  dispersion_crossfit_folds: 2
  dispersion_alpha_floor: 1.0e-8
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

`shapemix_count_only.yaml` must change only `use_shape` to `false`. Fragment bins and layer names belong to the dataset contract, not the method configuration.

The numeric values above are smoke-test starting points, not prevalidated scientific defaults. Tune or freeze them using toy cases and development seeds only; final test seeds must not influence the primary configuration.

Do not add ShapeMix to `configs/experiments/all_methods_all_atac_datasets.yaml`: most current datasets have no shape layers. Use ShapeMix-specific experiment configs until the runner has generic method-capability filtering.

Standard outputs remain:

```text
results/proportions.csv
results/abundance.csv
results/diagnostics.json
run.yaml
inputs.yaml
environment.txt
```

ShapeMix-specific native outputs should be compact:

```text
results/raw_method_output/
  training_history.csv
  restart_summary.csv
  reconstruction_summary.csv
  residuals_by_peak_and_bin.csv.gz
  signature_summary.yaml
```

Do not save the full fitted `S × P × B` expectation tensor by default.

Acceptance criteria:

- Both variants complete through `scripts/run_deconvolution.py` on the smoke dataset.
- They have identical count-likelihood configuration and differ only in `use_shape`.
- Legacy `methods` experiments retain their existing default run IDs; mixed schemas, duplicate IDs, and resolved run-ID collisions fail before execution.
- Both the `shapemix` registry method and its variant ID are present in run metadata and comparison rows.
- Every expected spot and declared reference cell type is present in output.
- Proportion rows sum to one within tolerance.
- Running another method does not import PyTorch/pysam through ShapeMix.
- Repeating a CPU run with the same seed yields the same result within tolerance.

Execution status: **completed 2026-08-22**. One lazily registered `shapemix` adapter now serves both nested arms, validates the complete ATAC fragment-shape contract, estimates fixed signatures from the reference only, returns standardized abundance and proportions, and streams reconstruction diagnostics without retaining the full expected tensor. The runner now accepts named `method_runs` while preserving the legacy `methods` syntax and default run IDs. It preflights schema, config-method agreement, duplicate variant IDs, resolved run-ID collisions, existing paths, and path escapes before creating the batch directory; every run and manifest surface records the registry method, variant ID, fully resolved config, canonical and source-file hashes, seed, device, dtype, determinism policy, and conditional Torch/pysam versions.

The 31 focused Step 5 adapter, registry, configuration, and runner tests pass, and the complete repository suite at that checkpoint passed with 166 tests and one pre-existing skip. Both frozen configurations completed through `scripts/run_deconvolution.py` on the 32-spot by 200-peak development dataset, produced finite `32 × 3` outputs whose rows sum to one, wrote all standard files plus the five compact native files, and used the identical reference-signature hash `3f642a07f02ddc60d6e25daf7e2745c63b51853b08b224a9ffd55d7b4d9d80a7`. A second CPU run in an independent output root reproduced proportions, abundances, selected restarts, objectives, histories, reconstruction summaries, compressed peak/bin residuals, signature summaries, and stable metadata exactly. Those historical smoke metrics remain development-only execution diagnostics; the strict Step 6 evaluator produced the first interpretable paired benchmark result.

### Step 6 — Fix evaluation and execute the primary ablation

Goal: produce a statistically valid, paired assessment of whether shape contributes information.

| Action | File | Implemented change or execution status |
|---|---|---|
| Modify | `src/deconvatac/metrics/proportions.py` | Score against the dataset-declared ordered cell-type universe, fill omitted declared predictions with zero, reject unknown extras, strictly validate values/rows, and implement versioned base-2 Jensen–Shannon divergence. |
| Modify | `src/deconvatac/metrics/__init__.py` | Export added per-type, correlation, and presence/detection metrics. |
| Add | `tests/test_proportion_metrics.py` | Test fixed-denominator penalties, unknown/duplicate types, exact spot sets, NaN/infinity/negative/all-zero failures, row sums, and the divergence formula. |
| Add | `scripts/summarize_shapemix_benchmark.py` | Combine run groups, compute paired nested-seed effects, resampling intervals, per-type results, rare-cell metrics, runtime, and memory. |
| Add | `tests/test_shapemix_benchmark_summary.py` | Test pairing, outer/inner seed grouping, metric versions, thresholds, intervals, and failure propagation. |
| Add | `configs/experiments/shapemix_primary_ablation.yaml` | Shape-aware and peak-only runs on equal/imbalanced datasets across frozen outer split and inner mixture seeds. |
| Add | `configs/experiments/shapemix_baselines.yaml` | Activate NNLS on all 20 primary datasets and declare Cell2location, RCTD, and SpatialDWLS behind explicit resource/dependency gates. NNLS completed 20/20; the gated methods were not executed. |
| Add | `configs/experiments/shapemix_stress_tests.yaml` | Declare depth, cells-per-spot, rare-cell, subtype, feature-count, and cutoff sensitivities behind a non-executable no-versioned-datasets gate. No stress evidence was produced. |
| Add | `configs/experiments/shapemix_negative_controls.yaml` and `scripts/run_shapemix_negative_controls.py` | Execute checksum-pinned development degeneracy controls and require an explicit positive-primary attestation before primary controls can run. |
| Add | `tests/test_shapemix_negative_controls.py`, `tests/test_shapemix_step6_configs.py`, and `tests/test_evaluate_runs.py` | Test control identities, execution gates, Step 6 configurations, strict evaluation, and failure behavior. |
| Modify | `src/deconvatac/shapemix/likelihood.py` and `src/deconvatac/shapemix/map.py` | Make cell-type-homogeneous shape signatures exactly abundance-invariant with zero `z` gradient, and center the stopping criterion on the initial shape likelihood so a data-only constant cannot change convergence timing. |
| Modify | `scripts/evaluate_runs.py` and `scripts/run_deconvolution.py` | Use one shared versioned metric registry, pass the declared universe, record full provenance and resource use, write per-run hash manifests, support deterministic paired sharding, and make resume fail closed. |
| Produce, tracked | `results/primary/shapemix_primary_ablation_protocol_v1__shard_*` and `results/primary/shapemix_primary_ablation_protocol_v1_summary/` | Write 40 primary runs, five shard manifests, comparisons, paired/outer effects, cell-type and rare-cell metrics, reconstruction, performance, failure, provenance, and strict summary artifacts. No figure was generated. |
| Add | `docs/ShapeMix/step6_results.md` | Record the non-normative execution report, results, gates, limitations, artifact map, and resource-pilot disclosure without modifying the frozen benchmark protocol. |
| Delete | None | Failed runs remain represented in manifests rather than being removed from summaries. |

Primary comparison:

- ShapeMix versus ShapeMix peak-only.
- Pair methods within the exact `(outer_split_seed, inner_mixture_seed, condition)` dataset. Average inner-mixture effects within each outer split before reporting the distribution across outer splits; do not treat spots or inner mixtures as independent biological replicates.
- These intervals quantify conditional resampling variability within one donor. They are not donor-level uncertainty or evidence of biological generalization.
- Primary metrics: RMSE and mean base-2 Jensen–Shannon divergence.
- Secondary metrics: per-cell-type MAE/RMSE, Pearson/Spearman correlation, rare-cell precision/recall/F1/AUPRC, runtime, and peak memory.
- Freeze the definition of a rare cell type, true presence, and prediction-positive threshold in `benchmark_protocol.md` using development seeds only.

Metric input and alignment contract:

1. `truth.cell_types` in the dataset YAML is the one ordered universe and must exactly equal the truth columns and reference-label universe.
2. Truth and prediction indices and columns must be unique. Prediction must have exactly the same spot set as truth; reorder only after checking set equality.
3. A missing declared prediction column is inserted as zero so omission is penalized. A prediction column outside the declared universe is an error rather than changing the denominator.
4. Truth and prediction values must be numeric, finite, nonnegative, and have no all-zero rows. Since these are standardized proportions, every row must sum to one within a frozen tolerance (for example `atol=1e-6`); invalid rows fail the run and are never filtered or silently normalized.
5. RMSE uses all `number_of_spots × number_of_declared_types` entries.
6. Define `jsd_v2` as the mean per-spot base-2 Jensen–Shannon **divergence**, `mean(jensenshannon(truth, prediction, base=2)^2)`, with range `[0,1]`. The current unsquared distance is a different metric; preserve its old results as `js_distance_v1` if historical comparison is needed.
7. Record the metric name/version and declared universe in comparison outputs. Non-finite metric output fails evaluation; never drop a row to make a mean finite.

Development controls completed before primary interpretation:

- Cell-type-homogenized `omega` exactly reproduced count-only proportions (`max abs difference = 0`).
- A one-bin shape model exactly reproduced its one-bin count-only partner and had a zero shape log likelihood.
- The Poisson factorization identity agreed to `1.78e-15` absolute error.
- A deterministic non-identity cell-type permutation of `omega` completed as a frozen development diagnostic. It was not truth-scored, and no same-pair primary permutation runs were launched because the preregistered positive-direction trigger was false.

Primary negative controls were conditional. For a condition, both co-primary mean effects had to be below zero and at least four of five outer effects had to favor ShapeMix for each endpoint. Neither condition met that rule, so no primary-control campaign was run. Before any future positive-result campaign, the permutation control's phrase "disappear or reverse" should be given a frozen numeric scoring rule and bound to the triggering summary hash.

The following sensitivities remain separately gated and unexecuted: two/three/five-bin cutoffs, zero versus fixed background, reference smoothing, 5,000/10,000/20,000 peaks, raw versus library-adjusted yield, grouped/global shape signatures, and separately versioned stress datasets for depth, cells per spot, rare-cell composition, and subtype difficulty. They must not be described as Step 6 evidence.

The first result need not be positive. If ShapeMix does not improve the predeclared endpoints, report that result and use the diagnostics to distinguish inadequate shape signal, excessive sparsity, reference/test instability, optimization failure, or a genuinely unhelpful feature.

Acceptance criteria:

- The evaluator cannot improve a score by omitting a truth cell type or spot.
- Every primary comparison has matched data, nested seeds, initialization policy, and optimization budget.
- Results include paired effects across nested resamples and label their uncertainty as one-donor conditional resampling variability.
- Count and shape reconstruction diagnostics accompany accuracy metrics.
- Baseline failures and resource limits are reported, not silently dropped.

Execution status: **completed 2026-08-23**. The frozen primary grid contained 20 registered pseudo-spatial datasets: five outer reference/held-out splits, two inner mixture seeds within each split, and two mixture conditions. Both arms ran on every dataset for 40/40 successful, converged fits. The five deterministic shards contained eight runs each; all output manifests, required files, run/batch cross-references, resolved campaign assignments, and current config/protocol/registry/code hashes revalidated before summarization. There were no failed runs, unavailable paired metrics, or unavailable outer summaries. The complete suite passes with **230 tests and one existing skip**.

The preregistered effect is shape-aware minus count-only, so negative values favor ShapeMix. Results across the five outer splits were:

| Condition | Endpoint | Count-only mean | Shape-aware mean | Mean paired effect | Bootstrap 95% interval | Outer splits favoring shape | Directional support |
|---|---:|---:|---:|---:|---:|---:|---:|
| Observed abundance | `rmse_v1` | 0.046732 | 0.046757 | +0.000025 | [-0.000107, +0.000126] | 2/5 | No |
| Observed abundance | `jsd_v2` | 0.122298 | 0.128065 | +0.005767 | [+0.005492, +0.006042] | 0/5 | No |
| Equal cell type | `rmse_v1` | 0.052396 | 0.057626 | +0.005230 | [+0.005017, +0.005443] | 0/5 | No |
| Equal cell type | `jsd_v2` | 0.156501 | 0.183554 | +0.027054 | [+0.026397, +0.027640] | 0/5 | No |

Fragment length therefore did **not** improve either condition under the frozen rule. It degraded JSD in all five outer splits for both conditions and degraded equal-cell-type RMSE in all five; observed-abundance RMSE was mixed and near zero. Descriptive diagnostics agreed with the negative primary result: all 16 cell types had worse mean RMSE in equal-cell-type mixtures, rare-cell pooled precision/recall/F1/AUPRC declined in both conditions, and tiny shape-bin reconstruction gains did not translate into more accurate proportions. The shape-aware arm took 1.74 times the count-only runtime on average; all 40 fits converged and sampled peak RSS stayed between 1,286 and 1,682 MB.

NNLS was the only activated external baseline and completed 20/20 datasets. It was worse than count-only on both endpoints in every outer split, so the negative ablation result is not explained by a generally nonfunctional count-only implementation. Cell2location still requires a frozen full-size device/resource configuration; RCTD lacks `spacexr`; SpatialDWLS lacks Giotto dependencies. These declarations are gates, not missing result rows.

One protocol departure is disclosed rather than hidden: before the primary launch, the shape-aware `pbmc_granulocyte_sorted_10k_shapemix_equal_celltype_split_1103_mix_101` fit was used only as a blinded resource pilot. Only dimensions, runtime/resource demand, convergence status, and output shape were inspected; prediction values, truth values, accuracy metrics, and paired effects were not. Its retained diagnostic records 1,268.645 seconds of fit time but does not record the thread environment. The pilot stayed under `/private/tmp/deconvatac-step6-profile-20260822/`, outside `results/primary`, and was excluded from every strict summary. Its sole purpose was selecting a safe five-shard, two-thread primary execution plan.

The result supports only conditional resampling conclusions within one PBMC donor. The 20 datasets are nested resamples, not 20 donors or biological replicates; this benchmark does not establish donor-level uncertainty, population generalization, or performance on independent real spatial ATAC tissue. See [the Step 6 execution report](step6_results.md) for detailed secondary results, provenance, gates, and artifact locations.

### Step 7 — Gated extensions, not MVP requirements

Steps 0–6 now provide a reliable independent-spot result, and that result is negative for the three-bin length term. Any extension below must be a separately versioned follow-up motivated by the recorded diagnostics; it cannot redefine or replace the completed primary analysis.

#### 7A. Shape overdispersion

| Action | File | Planned change |
|---|---|---|
| Add | `src/deconvatac/shapemix/dirichlet_multinomial.py` | Replace the multinomial conditional with a Dirichlet-multinomial using concentration `kappa[p] * rho[s,p,:]`. |
| Add | `configs/methods/shapemix_dm.yaml` | Fixed or empirical concentration configuration. |
| Add | `tests/test_shapemix_dirichlet_multinomial.py` | Limit to multinomial, mean/variance checks, finite gradients, and posterior-predictive dispersion. |
| Delete | None | Keep the multinomial as the simpler primary reference model. |

#### 7B. Variational uncertainty

| Action | File | Planned change |
|---|---|---|
| Add | `src/deconvatac/shapemix/vi.py` | Pyro SVI for uncertainty in `z` first; signature uncertainty comes later. |
| Add | `configs/methods/shapemix_vi.yaml` | Guide, ELBO, learning-rate, sample, and seed settings. |
| Add | `tests/test_shapemix_uncertainty.py` | Quantile ordering, coverage on toy data, ELBO stability, and deterministic posterior summaries. |
| Modify | `src/deconvatac/methods/shapemix.py` | Select MAP or VI backend lazily. |
| Modify | `src/deconvatac/data/schemas.py` | Add generic result-table support only if uncertainty should be standardized across methods. |
| Delete | None | MAP remains available as the reproducible baseline. |

Expected later outputs:

```text
results/proportion_q05.csv
results/proportion_q50.csv
results/proportion_q95.csv
results/posterior_predictive_summary.csv
results/raw_method_output/training_history.csv
```

Validate interval coverage, calibration, width versus depth, simulation-based calibration, multiple guide seeds, and posterior predictive zero rates/length ratios before calling the intervals trustworthy.

#### 7C. Background, mismatch, and latent signatures

Only add learned background, peak-specific dispersion, cross-protocol scaling, or latent `theta/omega/q/d` when residuals show they are necessary. Add one component at a time with a negative control and identifiability test. Do not infer spot depth, abundance scale, cell-type yield, peak detectability, accessibility, and shape simultaneously.

| Action | File | Planned change |
|---|---|---|
| Add | `src/deconvatac/shapemix/background.py` | Add fixed-background support first and learned background only if residual diagnostics motivate it. |
| Add | `configs/methods/shapemix_background.yaml` | Isolate the background choice from all other model settings. |
| Add | `tests/test_shapemix_background.py` | Test zero-background equivalence, recovery with known background, rare-type preservation, and scale identifiability. |
| Modify | `src/deconvatac/shapemix/likelihood.py`, `src/deconvatac/shapemix/signatures.py`, and `src/deconvatac/methods/shapemix.py` | Introduce only the single diagnostic-motivated component and expose its contribution separately. |
| Delete | None | Keep the fixed-signature, zero-background MVP as the primary benchmark. |

#### 7D. Richer shape features

Add cut-position, TSS, or motif-footprint layers only after fragment length succeeds or clearly fails for an understood reason. New axes require versioned data schemas and new dataset IDs; do not reinterpret existing length layers in place.

| Action | File | Planned change |
|---|---|---|
| Add | `scripts/prepare_shapemix_richer_shapes.py` | Generate a separately versioned position-, TSS-, or motif-aware dataset from raw fragments. |
| Add | `tests/test_shapemix_richer_shapes.py` | Test axis-specific boundaries, conservation rules, provenance, and isolation from fragment-length datasets. |
| Modify | `src/deconvatac/pp/fragment_shapes.py` and the data contract files | Generalize only the reusable sparse-layer machinery; require an explicit axis/schema version. |
| Produce, git-ignored | `data/processed/shapemix/<new_shape_axis>/` and new dataset IDs | Preserve the original fragment-length inputs and make cross-axis comparisons reproducible. |
| Delete | None | Never reinterpret or overwrite the primary length-bin layers. |

#### 7E. Spatial prior

| Action | File | Planned change |
|---|---|---|
| Add | `src/deconvatac/shapemix/spatial.py` | Sparse graph-Laplacian/GMRF penalty; do not materialize the tutorial's dense inverse covariance. |
| Add | `configs/methods/shapemix_spatial.yaml` | Neighbor graph and smoothing parameters. |
| Add | `tests/test_shapemix_spatial.py` | No-smoothing equivalence, disconnected graphs, boundary preservation, and oversmoothing detection. |
| Delete | None | Independent-spot ShapeMix remains the primary shape claim. |

#### 7F. External validation and hardening

The external source acquisition originally anticipated here is complete for GSE129785, GSE194122, GSE205055, and GSE263333. Source acquisition is not identical to experiment readiness: GSE194122 still needs fold-specific shape layers and pseudo-spots, while the real-spatial families still require compatible labeled fragment-level references. Section 13 supersedes this subsection with the complete execution order, comparison matrix, GPU work, resource policy, and acceptance gates. Real spatial ATAC validation remains qualitative unless independent cell-composition ground truth exists.

| Action | File | Planned change |
|---|---|---|
| Add | `scripts/profile_shapemix.py` | Profile time, CPU/GPU memory, convergence, and chunk parity across representative dimensions. |
| Add | `docs/ShapeMix/user_guide.md` and `docs/ShapeMix/results_interpretation.md` | Document stable preparation, execution, diagnostics, limitations, and result interpretation. |
| Add | A versioned external dataset config and preparation script, paths chosen when the dataset is selected | Apply donor/dataset-level validation without weakening provenance or the held-out design. |
| Modify | `README.md` and the ShapeMix tutorial sources | Link the stable workflow and clearly separate simulated, qualitative real-data, and externally validated evidence. |
| Modify later | `requirements.txt` and CI configuration | Regenerate and test the supported environment after the API and dependency set stabilize. |
| Delete | None | Retain the PBMC benchmark and all negative results as provenance. |

## 8. Complete tracked-file impact summary

### Implemented additions through Step 6

```text
docs/ShapeMix/iimplementation_plan.md         # active roadmap and execution record
docs/ShapeMix/README.md
docs/ShapeMix/model_specification.md
docs/ShapeMix/benchmark_protocol.md
docs/ShapeMix/step6_results.md

configs/data_sources/pbmc_granulocyte_sorted_10k_cellranger_arc_2.0.0.yaml

src/deconvatac/pp/fragment_shapes.py
src/deconvatac/shapemix/__init__.py
src/deconvatac/shapemix/config.py
src/deconvatac/shapemix/signatures.py
src/deconvatac/shapemix/likelihood.py
src/deconvatac/shapemix/map.py
src/deconvatac/shapemix/diagnostics.py
src/deconvatac/methods/shapemix.py

scripts/prepare_shapemix_pbmc.py
scripts/regenerate_shapemix_pbmc_simulations.py
scripts/audit_shapemix_signal.py
scripts/shapemix_provenance.py
scripts/summarize_shapemix_benchmark.py
scripts/run_shapemix_negative_controls.py
scripts/validate_shapemix_step3.py

configs/methods/shapemix.yaml
configs/methods/shapemix_count_only.yaml
configs/experiments/shapemix_smoke.yaml
configs/experiments/shapemix_primary_ablation.yaml
configs/experiments/shapemix_baselines.yaml
configs/experiments/shapemix_stress_tests.yaml
configs/experiments/shapemix_negative_controls.yaml

tests/test_fragment_shapes.py
tests/test_shape_data_contract.py
tests/test_shapemix_pbmc_preparation.py
tests/test_shapemix_simulation.py
tests/test_shapemix_signatures.py
tests/test_shapemix_likelihood.py
tests/test_shapemix_map.py
tests/test_shapemix_adapter.py
tests/test_proportion_metrics.py
tests/test_shapemix_benchmark_summary.py
tests/test_evaluate_runs.py
tests/test_shapemix_negative_controls.py
tests/test_shapemix_signal_audit.py
tests/test_shapemix_step5_configs.py
tests/test_shapemix_step6_configs.py
```

Gated extension files are excluded from the MVP list and should be added only when their step begins.

### Implemented modifications through Step 6

```text
.gitignore
configs/data_sources/pbmc_granulocyte_sorted_10k_cellranger_arc_2.0.0.yaml
docs/ShapeMix/README.md
docs/ShapeMix/iimplementation_plan.md
docs/ShapeMix/model_specification.md
docs/recreate_data_directory (important).md

src/deconvatac/pp/__init__.py
src/deconvatac/pp/feature_selection.py
src/deconvatac/data/schemas.py
src/deconvatac/data/loaders.py
src/deconvatac/data/validators.py
src/deconvatac/data/__init__.py
src/deconvatac/methods/registry.py
src/deconvatac/metrics/proportions.py
src/deconvatac/metrics/__init__.py

scripts/run_deconvolution.py
scripts/evaluate_runs.py

tests/test_method_interface.py
tests/test_run_deconvolution_experiment.py
tests/test_feature_selection.py
```

### Deletions through Step 6

None.

Do not delete or silently replace:

- existing PBMC reference or simulation H5ADs;
- current Heart/Russell data;
- current method adapters;
- `src/deconvatac/pp/reads_to_fragments.py`;
- historical ShapeMix tutorials;
- legacy scripts or archived results.

Temporary preprocessing shards may be removed only after final outputs are atomically written, checksummed, and validated. They are runtime scratch data, not tracked project files.

## 9. Generated artifacts and tracking policy

The repository must be reproducible even though `data/` is ignored. The entire `results/` tree, including development, primary, sensitivity, external-validation, and real-spatial campaigns, is intentionally exposed to Git so completed outputs and their provenance can be reviewed and committed. The scripts, tracked configs, hashes, and recreation guide remain part of the implementation, not optional documentation.

```text
data/raw/sources/10x_genomics/.../
  pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
  pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
  manifest.yaml

data/processed/shapemix/pbmc_granulocyte_sorted_10k/split_<split_seed>/
  reference_cells.h5ad
  heldout_test_cells.h5ad
  split.csv
  selected_peaks.txt
  signal_audit.csv
  manifest.yaml

data/processed/datasets/pbmc_granulocyte_sorted_10k_shapemix_<condition>_split_<split_seed>_mix_<mixture_seed>/
  atac/spatial.h5ad
  atac/features/*.txt
  truth/proportions.csv
  simulation/source_cells_by_spot.jsonl
  simulation/manifest.yaml
  dataset.yaml

results/<run_group>/<run_id>/
  results/proportions.csv
  results/abundance.csv
  results/diagnostics.json
  results/raw_method_output/*
  run.yaml
  inputs.yaml
  environment.txt

results/primary/shapemix_primary_ablation_protocol_v1__shard_<index>_of_05/
  batch_manifest.yaml
  runs.csv
  comparison.csv
  <eight method-run directories forming four pairs>/output_sha256.yaml

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

results/primary/shapemix_primary_baselines_protocol_v1/
  batch_manifest.yaml
  runs.csv
  comparison.csv
  <twenty NNLS run directories>/output_sha256.yaml

results/development/shapemix_negative_controls_v1/<dataset_id>/
  control_evidence.yaml
  control_proportions.csv
  output_sha256.yaml
```

Every generated manifest should include source hashes, code commit, package versions, full parameters, separate outer split and inner mixture seeds, bin edges, count semantics, split/feature hashes, declared cell types, dimensions, and nonzero counts.

## 10. Test and review gates

| Gate | Must be true before continuing |
|---|---|
| Data semantics | The cut-site-by-parent-fragment-length convention reproduces or explains the vendor peak matrix; bin layers sum exactly to `.X`; all canonical docs use the same terminology while older binary tutorials are labeled non-normative. |
| Reference/held-out separation | Reference cells, held-out evaluation cells, deterministic reference-only peak ranking, and pseudo-spot provenance pass disjointness checks. |
| Shape signal | Reference-only audit finds reproducible length distributions at sufficient coverage, or documents the need for grouped shrinkage. |
| Likelihood | Poisson equivalence, abundance-scaled NB parameterization, cross-fitted dispersion, stable zero-total shape term, identical-shape negative control, and chunk parity tests pass. |
| Optimization | Known mixtures are recovered, failures are explicit, and deterministic CPU runs agree. |
| Integration | Both variants run through the maintained runner, registry method and variant IDs are unambiguous, and old method-list experiments retain their run IDs. |
| Evaluation | The declared cell-type/spot universe is fixed, invalid output fails, `jsd_v2` is versioned, and paired nested-seed summaries/failure records are complete. |
| Scientific claim | Shape-aware improvement, no difference, or degradation is reported according to the frozen protocol without post-hoc endpoint changes. |
| Extensions | VI, background, spatial smoothing, and richer features begin only when a specific diagnostic motivates them. |

## 11. Principal risks and mitigations

| Risk | Consequence | Mitigation |
|---|---|---|
| Per-peak shape counts are sparse | `omega` is unstable and shape adds noise | Start with three bins and 5,000 covered peaks; shrink toward pooled shapes; test grouped signatures. |
| Technical depth drives length patterns | Apparent cell-type signal is confounded | Audit length versus QC/depth, use stratified splits, test exposure normalization, and add external data later. |
| One donor is reused across pools | Result does not establish donor generalization | State the limitation; require a second donor/dataset before generalization claims. |
| Full model is nonidentifiable | Optimizer finds arbitrary decompositions | Fix signatures/dispersion/background in MVP and infer only abundance. |
| Shape arm gets a different likelihood advantage | Ablation does not isolate shape | Use total NB in both arms and add only the conditional shape term. |
| Background absorbs rare types | False negatives or misleading fit | Use zero background in primary synthetic data; add fixed/learned background only as a labeled robustness analysis. |
| Dense tensor causes memory failure | 1,024 × 20,000 × 3 becomes expensive | Sparse H5AD layers, reference aggregation, spot batches, and peak chunks; never save full fitted means. |
| Evaluator changes its denominator or drops invalid rows | Metrics are falsely favorable | Use the dataset-declared cell-type universe and exact spot set, fill only omitted declared types with zero, reject extras/invalid rows, and version JSD. |
| One split is treated as replication | Intervals overstate generalization | Nest mixture seeds within multiple outer cell-split seeds and label all intervals as one-donor conditional resampling variability. |
| Optional dependencies break other methods | Existing workflows regress | Lazy imports, optional extras, CPU smoke tests, and no ShapeMix entry in incompatible all-dataset configs. |
| Post-hoc choices inflate a positive result | Scientific conclusion is unreliable | Freeze primary bins, endpoints, thresholds, seeds, and tests in `benchmark_protocol.md`. |

## 12. Definition of the minimum publishable implementation

ShapeMix is ready for a first scientific conclusion when all of the following are true:

- Raw fragment inputs and checksums are reproducible on another machine.
- Deterministic cell-disjoint reference/held-out splits, reference-only peak ranks, and the fragment counter are validated.
- Shape-aware reference and pseudo-spot H5ADs satisfy the sparse-layer contract.
- The fixed-signature total-NB plus conditional-multinomial MAP model passes toy, negative-control, chunking, and reproducibility tests.
- Shape-aware and peak-only variants run through the unified interface with one shared implementation.
- The primary comparison uses the same data, total-count model, parameters, initialization, and compute budget.
- Paired nested resamples report RMSE, versioned base-2 Jensen–Shannon divergence, rare-cell results, runtime, and model diagnostics without implying biological replication.
- Existing methods and current datasets still pass their tests unchanged.
- Limitations include the one-donor reference, same-protocol pseudo-spots, and lack of calibrated posterior uncertainty.
- Results are reported even if fragment length provides no improvement.

Execution status: **satisfied for the first one-donor conditional conclusion on 2026-08-23**. The conclusion is negative: the frozen fragment-length term did not improve the co-primary endpoints and generally degraded them. This milestone is not external validation and does not satisfy the separate requirement for a second donor or real-tissue generalization.

Variational inference, spatial smoothing, motif footprints, and real-tissue maps are valuable follow-up work, but they are not required to answer the first ShapeMix research question.

## 13. Full cross-dataset evaluation after source acquisition

This section is the executable roadmap after completion of the four external acquisitions. It extends the completed one-donor protocol-v1 result; it does not replace, rerun, or reinterpret that frozen result. New configurations, datasets, and output groups require new IDs and hashes.

The canonical path and lifecycle authority for every new artifact is [the ShapeMix file-organization specification](file_organization.md). If this section and that specification ever conflict, `file_organization.md` wins and the campaign must stop until the plan or implementation is corrected.

### 13.1 Readiness is not the same as acquisition

All planned GEO source families are present locally with pinned identities and source-ready preprocessing. They are at different analysis-readiness levels:

| Family | Current state on 2026-08-24 | Experiment class | Remaining gate before fitting |
|---|---|---|---|
| Existing 10x PBMC | 20 primary datasets evaluated; fragment cache and five reference splits retained | Controlled mechanism and sensitivity | Materialize separately versioned stress datasets; do not alter protocol-v1 objects |
| GSE129785 | Standardized nine-type reference plus 16 runner-readable datasets are registered; all physical ratios are validation-only and no descriptor declares exact truth | Seven quantitative nominal CD4-memory/CD8-naive dilutions; seven broad monocyte/T-cell dilutions; PBMC replicate and preparation robustness | Freeze the external experiment configuration and complete the all-descriptor runner preflight |
| GSE194122 | All 13 fragment/index/metrics trios, 69,249 annotated cells, ten donor membership tables, and the processed Multiome object are validated | Quantitative leave-one-donor-out pseudo-spot evaluation | Freeze label support/harmonization, select peaks within each training fold, build shape layers and held-out-donor pseudo-spots, then register datasets |
| GSE205055 | Complete source archive and related SubSeries are extracted; six matched ATAC/RNA spatial groups are aligned | Real-spatial validation without exact composition truth | Acquire and validate species/tissue/stage-compatible labeled scATAC or Multiome references and build shared feature axes |
| GSE263333 | Complete source archive is extracted; two ATAC sections are aligned to RNA and available protein/histone evidence | Real-spatial multi-omic validation without exact composition truth | Acquire and validate compatible mouse embryo and adult/EAE-brain labeled fragment references and build shared feature axes |

The last two rows are deliberately blocked from deconvolution until their reference gates pass. A spatial fragment file supplies the mixture but cannot by itself supply the fixed labeled cell-type signatures required by ShapeMix.

### 13.2 Frozen comparison and evidence rules

Every runnable dataset receives the same core three-method comparison:

1. `shapemix_length`: total-count negative binomial plus the conditional fragment-length term;
2. `shapemix_count_only`: the exact nested model with only the conditional length term disabled; and
3. `nnls`: a fast collapsed-count baseline on the same reference and ordered feature axis.

The two ShapeMix arms must use identical device, dtype, optimizer, initialization, restarts, convergence budget, signatures, spot ordering, and feature ordering. Only `use_shape` may differ. Cell2location, RCTD, and SpatialDWLS remain secondary baseline tracks: activate them only after their device/dependency pilots and frozen configurations pass. A missing optional baseline is reported as a gated limitation, never as an omitted or silently failed result.

Evidence classes must remain separate:

- Exact pseudo-spot truth supports `rmse_v1`, `jsd_v2`, per-type errors, correlations, and rare-cell detection metrics.
- GSE129785 physical dilution ratios are nominal sample-level inputs. Report mixture-ratio error, monotonicity, calibration, and rare-component detection, but retain sorting and differential-recovery uncertainty.
- Unlabelled PBMC/preparation cohorts support replicate stability and preparation sensitivity, not truth-based accuracy.
- Real spatial sections support RNA/protein/histone/marker/anatomy concordance and map stability, not RMSE or JSD against a manufactured truth table.
- The biological replication unit is the donor for GSE194122, the independently prepared sample/section where available for physical or spatial data, and never the spot, pixel, or repeated optimizer seed.

### 13.3 Minimum core campaign matrix

| Campaign | Materialized evaluation units | Frozen core design | Core jobs before optional baselines |
|---|---:|---|---:|
| Completed PBMC protocol v1 | 20 pseudo-spot datasets | Five outer cell splits × two inner mixtures × two abundance conditions | Already complete: 40 ShapeMix fits plus 20 NNLS fits |
| GSE129785 external immune | 16 registered descriptors | Run both ShapeMix arms and NNLS on every descriptor; compare seven CD4-memory/CD8-naive samples with nominal validation evidence outside the runner truth contract | 48 |
| GSE194122 BMMC | 40 planned pseudo-spot datasets | Ten held-out donors × two inner mixture seeds (`101`, `211`) × two conditions (`observed_abundance`, `equal_celltype`) | 120 |
| GSE205055 real spatial | Six matched ATAC/RNA groups | One frozen compatible reference per species/tissue/stage track; run both ShapeMix arms and NNLS per section | 18 |
| GSE263333 real spatial | Two ATAC-bearing multi-omic sections | Embryo and adult/EAE-brain reference tracks; run both ShapeMix arms and NNLS per section | 6 |

The minimum new core is therefore 192 jobs, of which NNLS jobs should be short. This count excludes separately versioned PBMC stress datasets and optional Cell2location/RCTD/SpatialDWLS runs. Do not launch the 192-job campaign as one undifferentiated batch: each stage below has a fail-closed gate and its own summary.

### 13.4 Ordered execution stages

#### Evaluation stage E0 — Freeze the external protocol and preflight inputs

- Add `docs/ShapeMix/full_evaluation_protocol.md` with immutable dataset roles, label universes, thresholds, endpoints, seeds, and analysis units.
- Add tracked dataset reconstruction templates under `configs/datasets/` and separate experiment configurations under `configs/experiments/` for GSE129785, GSE194122, GSE205055, and GSE263333; do not place truth-bearing and qualitative datasets in one summary group.
- Add a layout validator that enforces `file_organization.md` before promotion, registry mutation, campaign launch, and final summarization.
- **Completed 2026-08-24:** migrated the seven GSE129785 CD4-memory/CD8-naive nominal ratio files out of the reserved `truth/` contract, removed their runner `truth` declarations, and updated the materializer so all 14 physical dilutions use `validation/nominal_broad_proportions.csv`.
- Validate every registered dataset through `load_deconvolution_input`, exact shape-layer conservation, ordered-axis hashes, reference/test separation, and output-path collision checks.
- Require project-relative authoritative paths; absolute paths may appear only as supplemental execution provenance.
- Freeze harmonized label tables before fitting. Any later merge or exclusion creates a new protocol version and new dataset IDs.
- Record a campaign inventory with planned, runnable, gated, successful, and failed units. A source-ready but non-runnable family must be shown as gated rather than absent.
- Verify that the entire `results/` tree remains exposed to Git and that no data cache, download, or scratch artifact enters it.

#### Evaluation stage E1 — Implement and qualify the CUDA backend

Implement the GPU work in Sections 13.6–13.9 and pass CPU/CUDA parity, determinism, memory, and performance gates before external results are inspected. The development smoke dataset and one representative full-size protocol-v1 dataset are resource pilots only. They must be written outside the new external result roots and excluded from scientific summaries.

#### Evaluation stage E2 — Run GSE129785

This is first because all 16 descriptors and the standardized reference already exist.

- Primary nominal-validation family: seven CD4-memory/CD8-naive ratios from 0.1/99.9 through 99.9/0.1. Read ratios from `validation/nominal_broad_proportions.csv`, never from the exact-truth loader contract. Report nominal-fraction RMSE/JSD descriptively, absolute error of the rare component, detection at 0.1%, 0.5%, and 1%, calibration slope/intercept, rank monotonicity across the series, and predicted off-target mass.
- Exploratory physical family: seven monocyte/T-cell ratios. Score the broad monocyte-versus-total-T contrast only; do not present the broad ratios as subtype-level truth.
- Unlabelled family: four PBMC replicates and three preparation conditions. Report replicate dispersion, preparation shifts, off-target mass, shape/count differences, reconstruction, and convergence without inventing truth.
- Analyze one prepared sample as one sample. Barcodes aggregated within that preparation are not biological replicates.
- Summarize the paired `length - count_only` effect separately for the quantitative-nominal, broad-nominal, and no-truth families.

#### Evaluation stage E3 — Build and run GSE194122 donor-held-out evaluation

For each of the ten Multiome donors:

1. Exclude every occurrence of that donor across all sites before peak selection, signature estimation, smoothing, or tuning.
2. Apply a frozen support rule to the 22 author labels. Produce one shared broad ontology for cross-donor summaries and retain original fine labels for secondary reporting. Fail a fold if a declared type lacks the predeclared training support; do not remove it after viewing predictions.
3. Select 5,000 peaks from training donors only with the frozen count-only selector.
4. Count three fragment-length layers for training and held-out cells on that exact axis and require their sum to equal collapsed counts.
5. Create 1,024 held-out-donor pseudo-spots for each of two mixture seeds and each of the observed-abundance/equal-cell-type conditions. Record every source cell.
6. Keep repeated sites for donor 1 and donor 3 on the held-out side together. Preserve site identity so within-donor site sensitivity can be reported without treating sites as extra donors.
7. Run the three core methods, then summarize inner mixtures within donor before computing donor-level effects and intervals.

Primary endpoints remain `rmse_v1` and `jsd_v2`; secondary endpoints include per-type error, rare-cell recovery, reconstruction, convergence, runtime, and site-stratified stability. This is the first campaign that can support donor-level biological generalization claims.

#### Evaluation stage E4 — Run diagnostic PBMC stress tests

Because protocol v1 was negative, stress tests are explanatory rather than a search for a favorable subset. Build one factor at a time from held-out cells, with the three-bin/5,000-peak/ten-cells-per-spot setting as the anchor:

**Frozen 2026-08-24, pending guarded execution:** the exact levels, two
evaluation seeds, alternate-bin edges, and controls are declared in
`full_evaluation_protocol.md`. The design contains 40 unique datasets and 120
core jobs. Its tracked reconstruction template, experiment config, builder,
summarizer, and focused tests are present; no sensitivity prediction was read
before the design was frozen.

- fragment/cut-site depth retention;
- cells per spot;
- controlled rare-cell fraction;
- closely related subtype mixtures;
- reference training-cell support;
- feasible peak-count levels after a source-universe audit; and
- two-, three-, and five-bin representations.

Run development seeds to validate each generator, then freeze the final factor levels and evaluation seeds before viewing paired accuracy. Do not run a complete factorial grid. Only a predeclared depth-by-rare-fraction interaction may follow the one-factor campaign.

#### Evaluation stage E5 — Obtain references and run real-spatial sections

Create three separately gated reference tracks rather than one universal atlas:

- mouse embryo for the E13 sections;
- adult mouse brain, with disease/protocol mismatch explicitly recorded for the five-month EAE section; and
- adult human hippocampus for the human section.

The frozen candidates are GSE246791 for adult mouse brain, a nine-sample GSE244618 hippocampal donor/region subset for human brain, and GSE216371 for E13-compatible whole mouse embryo. Candidate status is not approval: each source still requires acquisition, genome-build verification, labeled-cell support, parent-fragment semantics, and coordinate validation. Peak selection uses the reference only.

For each section, freeze the reference, label universe, feature axis, and cross-modality alignment before fitting. Report:

- ShapeMix-versus-count-only abundance-map correlation and absolute difference;
- concordance with independently constructed RNA cell-type scores;
- protein concordance where available;
- histone/marker accessibility and anatomical-region enrichment;
- replicate consistency for matched sections;
- local spatial continuity and boundary preservation, reported descriptively because no spatial prior is used; and
- off-reference mass or residual warnings rather than forcing every spatial signal into a known type.

Do not tune labels, peaks, smoothing, or model parameters to maximize RNA/protein concordance on the same section later used for reporting.

#### Evaluation stage E6 — Cross-family synthesis

- Produce one summary per evidence class and a final cross-family table; never pool exact pseudo-truth, nominal physical truth, and qualitative spatial concordance into one effect estimate.
- Treat GSE194122 donors as the population-level quantitative units. Present GSE129785 ratios and spatial sections individually because their small sample counts do not justify asymptotic uncertainty claims.
- Relate gains or failures to depth, reference support, shape entropy/divergence, protocol mismatch, and residual diagnostics using predeclared analyses.
- Report all failures, convergence exclusions, gated references, and unavailable baselines.
- Preserve the negative protocol-v1 result and label all follow-up model changes as new versions.

### 13.5 Canonical file organization for the new campaigns

All generated files must follow [the canonical ShapeMix file-organization specification](file_organization.md). The paths below are normative for this campaign; scripts must not introduce parallel roots or source-specific synonyms. Authoritative paths written into YAML/JSON manifests are project-relative.

#### Configuration and source identity

```text
configs/data_sources/<source_id>.yaml
configs/datasets/<dataset_id>.yaml
configs/experiments/<campaign_id>.yaml
configs/methods/<method_or_variant>.yaml
```

Source identities for any newly acquired spatial reference remain under `configs/data_sources/`. Every generated runnable dataset also gets a tracked reconstruction template under `configs/datasets/`; the ignored realized descriptor remains `data/processed/datasets/<dataset_id>/dataset.yaml`. Experiment and method configs never belong inside `data/` or `results/` as their only copy.

#### Staging, immutable source, and disposable work

```text
data/work/downloads/<source_id>/
data/work/preprocessing/<family>/<job_id>/
data/raw/sources/<provider>/<accession_or_release>/
/tmp/deconvatac-<task>-<unique_id>/
```

Downloads enter `data/raw/sources/` only after identity, size, checksum, compression, and schema validation. Restartable intermediates stay in `data/work/preprocessing/`; per-process scratch stays under a unique `/tmp/deconvatac-*` directory. Neither class may be used as a registered benchmark input or copied into `results/`.

#### Reusable preprocessing, references, and runnable datasets

```text
data/processed/shapemix/<family>/
  source_audit/
  normalized_fragments/
  fragment_shape_cache/
  labels/
  feature_axes/
  splits/
  manifests/

data/processed/references/<reference_id>/
  reference.yaml
  atac/reference.h5ad
  rna/reference.h5ad                 # only when available

data/processed/datasets/<dataset_id>/
  dataset.yaml
  atac/spatial.h5ad
  atac/features/*.txt
  truth/proportions.csv              # exact pseudo-spot truth only
  simulation/source_cells_by_spot.jsonl
  simulation/manifest.yaml
  validation/                        # nominal or orthogonal evidence

data/registry/datasets.yaml
```

GSE194122 fold products must use `data/processed/shapemix/gse194122_bmmc/{fragment_shape_cache,labels,feature_axes,splits,manifests}/`. Each held-out-donor reference gets an immutable versioned ID under `data/processed/references/`, and each split/condition/seed combination gets a distinct dataset ID under `data/processed/datasets/`. GSE205055/GSE263333 reference-source preprocessing follows the same reusable tree before standardized reference H5ADs are promoted to `data/processed/references/`.

Only exact source-cell pseudo-spot composition belongs in `truth/proportions.csv`. GSE129785 physical ratios live at `validation/nominal_broad_proportions.csv` and are referenced as nominal validation evidence, including the seven CD4-memory/CD8-naive samples. The legacy `truth` declarations and directories were removed on 2026-08-24, and `scripts/preprocess_gse129785.py` was updated so regeneration cannot recreate them. RNA, protein, histone, marker, image, and anatomical evidence for real spatial data also stays under `validation/`.

A dataset may enter `data/registry/datasets.yaml` only after its exact reference/spatial axes, fragment-shape contract, declared labels, evidence class, split disjointness, and provenance pass. IDs used by completed campaigns are append-free and never overwritten; a semantic change receives a new versioned reference, dataset, or campaign ID.

#### Tracked result scopes

```text
results/development/shapemix_gpu_qualification_v1/
results/sensitivity/shapemix_pbmc_stress_v1/
results/external_validation/shapemix_gse129785_v1/
results/external_validation/shapemix_gse194122_lodo_v1/
results/real_spatial/shapemix_gse205055_v1/
results/real_spatial/shapemix_gse263333_v1/
results/external_validation/shapemix_full_evaluation_v1/
```

There is no `results/external/` or `results/summary/` root. Cross-family synthesis belongs in the versioned full-evaluation campaign above and records the hashes of the contributing campaign summaries. Every run directory retains `run.yaml`, `inputs.yaml`, `environment.txt`, `output_sha256.yaml`, standardized results, diagnostics, and raw method output. Each campaign root retains its resolved protocol, batch manifest, run/comparison/failure/resource tables, logs required to interpret failures, provenance hashes, and summary. Completed result directories are append-free and fully exposed to Git.

Add `scripts/validate_shapemix_file_layout.py` and focused tests before materialization begins. The validator must reject noncanonical result scopes, absolute authoritative paths, raw-source mutation, references outside `data/processed/references/`, family caches inside runnable dataset roots, nominal/qualitative evidence under `truth/`, premature registry entries, missing per-run provenance, and writes into completed campaign IDs. Run it before promotion, before campaign launch, and during final result revalidation. Resume remains fail-closed: an existing run may be reused only when every declared input, configuration, code, and output hash matches.

### 13.6 Why GPU acceleration is appropriate

The observed host has an Intel i7-10700K (8 physical cores/16 threads) and one NVIDIA GeForce RTX 3080 with 10,353,442,816 bytes of VRAM. The installed environment is already CUDA-capable: PyTorch `2.11.0+cu130`, CUDA runtime 13.0, compute capability 8.6, and `torch.cuda.is_available() == True`.

The current fitter does not use that device. `ShapeMixConfig` accepts only `device: cpu`; `map.py` rejects non-CPU dense tensors; fixed signatures, abundance parameters, and every count chunk are explicitly created on CPU. Each optimization step loops over spot/peak chunks twice—gradient plus complete-objective evaluation—and repeatedly converts SciPy slices to dense NumPy/Torch arrays. The historical full-size resource pilot took 1,268.645 seconds for one shape-aware fit, while the primary campaign required 40 ShapeMix fits. The streamed likelihood and gradient are therefore the first GPU target.

GPU acceleration does not help every stage:

| Workload | Device policy | Reason |
|---|---|---|
| Tabix fragment parsing, coordinate validation, BGZF/index work, H5AD I/O | CPU, bounded workers | Branch-heavy and I/O-bound; no useful tensor kernel |
| Sparse reference aggregation, label/axis audits, SciPy NNLS initialization | CPU, one or two math threads | One-time work; preserve the validated algorithms first |
| ShapeMix likelihood, gradient, Adam state, fixed signatures, abundance parameters | GPU | Repeated dense matrix products and elementwise `lgamma`/log operations dominate fit time |
| Reconstruction summaries and CSV/H5AD output | CPU initially | One pass only; move later only if profiling shows it matters |
| Cell2location, if activated | Same single GPU, separate scheduled job | It already has a GPU path but must not contend with ShapeMix |
| RCTD/SpatialDWLS, if activated | CPU with explicit thread limits | Their maintained adapters are not CUDA workloads |

### 13.7 CUDA implementation design

Implement one backend, selected by an explicit device, rather than a separate scientific model:

| File | Required change |
|---|---|
| `src/deconvatac/shapemix/config.py` | Accept `cpu` and explicit `cuda:<index>` devices; reject unavailable/malformed devices. Separate statistical-model hashing from execution hashing so device/chunk/cache choices do not masquerade as different likelihoods, while retaining the full resolved hash. Keep final campaigns explicit; `auto` is development-only if added at all. |
| `src/deconvatac/shapemix/map.py` | Move `A`, `omega`, `raw_z`, and Adam state to the selected device. Add a count cache that materializes sparse bin layers blockwise into one dense `float32 [S,P,B]` CUDA tensor when it fits, with a bounded host-to-device streaming fallback. Remove hard-coded CPU devices. |
| `src/deconvatac/shapemix/map.py` | Accumulate objective scalars on device and synchronize once per complete objective, not once per peak chunk. Defer finite-status synchronization to a safe step boundary. Preserve the same mathematical sum, restart seeds, and convergence rule. |
| `src/deconvatac/methods/shapemix.py` | Record device identity, GPU name, compute capability, CUDA/runtime versions, peak allocated/reserved VRAM, cache mode, and CPU thread policy in diagnostics. |
| `scripts/run_deconvolution.py` | Fail closed when a requested CUDA device is absent or already over the frozen memory threshold. Preserve one-process-per-GPU execution and campaign resume semantics. |
| `scripts/profile_shapemix.py` | Benchmark CPU and CUDA on synthetic shapes, the smoke dataset, one 1,024 × 5,000 dataset, and one 10,000 × 5,000 spatial-shaped input; synchronize around timing and record peak VRAM/RSS. |
| `configs/methods/shapemix_cuda.yaml` and `configs/methods/shapemix_count_only_cuda.yaml` | Freeze identical CUDA execution settings and differ only in `use_shape`. |
| `tests/test_shapemix_gpu.py` | Cover availability failure, objective/gradient/output parity, chunk/cache parity, deterministic repeats on one device, OOM-safe fallback, diagnostics, and CPU backward compatibility. Skip hardware tests only when CUDA is genuinely unavailable. |

For 5,000 selected peaks, even the largest current 10,000-pixel section needs about 600,000,000 bytes for the dense three-bin `float32` count cache. This is small enough for the 10 GB RTX 3080 while leaving room for likelihood intermediates, signatures, parameters, and the desktop. The implementation must calculate the requirement from `S × P × B × dtype_bytes` and actual free VRAM; it must not assume that every future dataset fits.

Use this cache policy:

1. Query currently free VRAM immediately before a fit.
2. Reserve at least 25% of free VRAM for autograd workspaces, allocator fragmentation, and the display process.
3. Use full GPU count caching only when the estimated cache plus a measured workspace margin fits below the remaining 75% budget.
4. Otherwise cache one or more spot blocks and transfer each block once per restart, not once per optimization step.
5. On CUDA OOM, cleanly discard the partial restart and retry the predeclared smaller cache/chunk level. Record the fallback; never silently change dtype or model parameters.

Start with `float32` and `torch.set_float32_matmul_precision("highest")`. Do not enable autocast, `float16`, `bfloat16`, or TF32 in the scientific campaign until a separately versioned numerical-sensitivity study passes. The likelihood contains `lgamma`, low-count probabilities, and convergence comparisons for which reduced precision can change optimization behavior.

PyTorch documents that CUDA execution is asynchronous, so timing must use `torch.cuda.synchronize()` or CUDA events; it also exposes `max_memory_allocated()` for peak device-memory reporting. Host-to-device copies may use `non_blocking=True` after measurement shows a benefit, but manual pinning is not assumed to be faster. See the official [CUDA semantics](https://docs.pytorch.org/docs/stable/notes/cuda.html), [CUDA peak-memory API](https://docs.pytorch.org/docs/stable/generated/torch.cuda.max_memory_allocated.html), and [transfer benchmark guidance](https://docs.pytorch.org/tutorials/intermediate/pinmem_nonblock.html).

### 13.8 Reproducibility and CPU/CUDA qualification gates

PyTorch explicitly warns that bitwise equality is not guaranteed between CPU and GPU or across releases/platforms. The campaign must therefore freeze the exact device/software stack and demonstrate numerical equivalence rather than claim cross-device byte identity. Follow the official [PyTorch reproducibility guidance](https://docs.pytorch.org/docs/stable/notes/randomness.html).

Before CUDA becomes the production backend:

1. Enable `torch.use_deterministic_algorithms(True)` with `warn_only=False`, retain the existing NumPy seed tuples, set the Torch seed, and record the deterministic policy. Set any cuBLAS workspace requirement before process start.
2. Require CPU/CUDA agreement on toy objective components and gradients with predeclared absolute/relative tolerances.
3. Require CPU/CUDA smoke proportions to agree within `1e-5` maximum absolute error and full-size development proportions within `1e-4`; require `rmse_v1` and `jsd_v2` to agree within `1e-5`. Any threshold change must occur before external predictions are inspected.
4. Repeat the CUDA smoke and full-size pilot twice on the same hardware/software stack and require identical convergence status plus numerical agreement within `1e-7`.
5. Verify cached and streamed CUDA modes agree within the same tolerances and produce the same declared spot/type axes.
6. Benchmark synchronized end-to-end fit time, peak VRAM, peak RSS, average CPU utilization, and host-to-device time. Adopt CUDA for the full campaign only if it passes parity and is at least twice as fast on the representative full-size fit. Otherwise retain CPU and use the resource limits below while profiling the bottleneck.
7. Use the same qualified backend for both arms of every pair. Never compare a CPU count-only fit with a CUDA shape-aware fit in the primary paired analysis.

### 13.9 CPU-load containment and scheduling

Co-tenant constraint, added 2026-08-24: another workload on this host may consume four to five physical cores. Treat that workload as higher priority and preserve at least one physical core of headroom. This constraint remains active until the user explicitly removes it.

Only one deconvATAC process may run at a time, and only one fit process may own the RTX 3080. Do not overlap a GPU fit, NNLS baseline, preprocessing job, checksum pass, download, test suite, or summary job from this repository. Do not launch if `nvidia-smi` shows an unrelated compute process using the GPU; GPU sharing requires explicit approval.

Before every job, record the one-minute load average, available RAM, and GPU process/memory state. On this eight-physical-core/16-logical-CPU host, do not start a new deconvATAC task while the one-minute load average is `>= 6.0`; leave it queued and recheck later. This is a conservative launch gate, not permission to interrupt either workload. A running job may finish, but no additional repository process starts while the gate is closed.

All new fits and other material work must enter through the guarded launcher:

```bash
scripts/run_shapemix_low_impact.sh \
  .venv/bin/python scripts/run_deconvolution.py --experiment-config <config>
```

The launcher holds a host-wide deconvATAC lock, enforces the load/GPU/worker gates, lowers CPU and I/O priority, and caps the common host math libraries before Python imports them. Bypassing it is not permitted while the co-tenant constraint is active.

The runner must also call `torch.set_num_threads(1)` and `torch.set_num_interop_threads(1)` before fitting. GPU and CPU-only fits, including NNLS, are restricted to one host thread. The earlier two-thread GPU-input exception is disabled while the co-tenant constraint is active.

Fragment counting, preprocessing, and downloads may use at most two workers (`--workers 2`) and must also run at lowered priority. Checksum scans use one process. If a stage cannot obey these limits, keep it gated rather than increasing parallelism. Never automatically fall back from CUDA to an unconstrained CPU fit.

Use one sequential GPU shard per campaign and resume by completed, hash-verified run. Record the launch-gate measurements, priority, worker/thread caps, wall time, CPU time, peak RSS, peak allocated/reserved VRAM, device-utilization samples, and any cache/chunk fallback. A load, temperature, memory-pressure, or GPU-occupancy monitor pauses launching the next run but does not terminate or rewrite a completed result.

### 13.10 Full-evaluation acceptance criteria

The full campaign is complete only when:

- every source family is classified as completed or explicitly gated, with no silent omission;
- all core runnable units have paired ShapeMix arms on the same qualified device and an NNLS baseline;
- GSE194122 has no donor leakage and is summarized at the donor level;
- all GSE129785 conclusions respect nominal versus absent truth, with no nominal or qualitative evidence stored under `truth/`;
- GSE205055/GSE263333 use audited compatible references and make no exact-composition claim;
- CPU/CUDA parity, repeated-CUDA determinism, memory fallback, and performance gates pass and are retained under `results/development/`;
- resource manifests prove that only one deconvATAC process ran at a time, every fit used one host thread, preprocessing used at most two workers, and every launch passed the co-tenant load/GPU gate;
- failed runs and optional-baseline gates remain visible in summaries;
- the canonical layout validator passes for configs, work paths, reusable products, references, runnable datasets, registry entries, result scopes, and per-run provenance;
- all result groups, configs, code/input hashes, and environment metadata revalidate; and
- the final report presents negative, null, or positive outcomes under the same predeclared rules.
