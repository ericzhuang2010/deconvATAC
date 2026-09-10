# ShapeMix-ATAC full evaluation report

Status: complete 2026-09-09  
Evaluation version: `shapemix_full_evaluation_v2`  
Canonical result root:
`results/external_validation/shapemix_full_evaluation_v2/`

## 1. Executive summary

The complete ShapeMix evaluation has finished across the original PBMC
benchmark and all four acquired GEO families: GSE129785, GSE194122, GSE205055,
and GSE263333. The evaluation contains 124 runnable dataset units and 372
successful method runs. Each unit received length-aware ShapeMix, the exactly
matched count-only ShapeMix ablation, and NNLS. No evaluation or monitoring
process remains active.

The main scientific conclusion is negative but informative: the frozen
three-bin parent-fragment-length likelihood does **not** provide a reliable
general improvement over count-only ShapeMix.

- The original PBMC benchmark showed clear JSD degradation and no meaningful
  RMSE benefit.
- The ten-donor GSE194122 evaluation showed small, heterogeneous effects whose
  donor-level confidence intervals all included zero.
- GSE129785 length-aware predictions were worse on mean nominal RMSE, JSD, and
  off-target mass for both physical-mixture families.
- The PBMC stress experiment identified a limited low-depth/coarse-label
  niche, but the benefit did not generalize and rare-NK AUPRC was consistently
  worse.
- In eight real-tissue sections, length-aware and count-only maps were highly
  similar. Orthogonal RNA, protein, epigenome, and marker concordance was weak
  and inconsistent, and every real-spatial run raised a reconstruction-mismatch
  warning. These sections have no exact composition truth, so they cannot
  establish accuracy.

For the current fixed-signature, three-bin MAP model, count-only ShapeMix
should remain the default. Fragment length is better treated as a targeted
hypothesis for low-depth or coarsened-label settings that needs independent
confirmation before another production campaign.

## 2. Scope and completed inventory

The final synthesis keeps evidence classes separate. Exact pseudo-spot truth,
nominal physical input ratios, prediction-only preparation comparisons,
backend qualification, and real-spatial proxy evidence are not pooled into one
effect estimate.

| Stage | Dataset family | Evidence class | Runnable units | Successful jobs |
|---|---|---|---:|---:|
| E0 | PBMC protocol-v1 primary | Exact held-out pseudo-spot truth within one donor | 20 | 60 |
| E1 | CUDA qualification v2 | Computational backend qualification | Development arms | Passed |
| E2 | GSE129785 | Nominal physical mixtures and prediction-only cohorts | 16 | 48 |
| E3 | GSE194122 | Exact held-out pseudo-spot truth across ten donors | 40 | 120 |
| E4 | PBMC diagnostic stress v2 | Exact held-out pseudo-spot truth within one donor | 40 | 120 |
| E5 | GSE205055 | Six real-spatial sections without exact composition truth | 6 | 18 |
| E5 | GSE263333 | Two real-spatial sections without exact composition truth | 2 | 6 |
| **Total** | Four GEO families plus PBMC | Five scientific evidence classes | **124** | **372** |

The final cross-family synthesis contains 7 separately labeled evidence
entries, 82 effect rows, and 352 normalized resource rows. The resource table
has 352 rather than 372 rows because the historical PBMC primary resource
summary covers its 40 ShapeMix runs but not its 20 NNLS runs. Those NNLS runs
are present and passed output-manifest validation.

## 3. Data acquisition and preprocessing

### 3.1 Canonical storage policy

All files follow the
[ShapeMix file-organization specification](file_organization.md):

```text
data/work/downloads/<source_id>/             incomplete transfers
data/work/preprocessing/<family>/<job_id>/   restartable intermediates
data/raw/sources/<provider>/<accession>/      immutable provider files
data/processed/shapemix/<family>/             reusable derived caches/audits
data/processed/references/<reference_id>/     standardized references
data/processed/datasets/<dataset_id>/         validated runner inputs
results/<scope>/<campaign_id>/                tracked experiment outputs
```

The large `data/` tree remains ignored. The entire `results/` tree is exposed
to Git. Standardized references live under `data/processed/references/`, not
`data/raw/references/`.

### 3.2 PBMC benchmark

The original human PBMC source contains one healthy donor. The benchmark uses
9,500 labeled cells across 16 immune cell types. Five stratified 70/30 splits
produce 6,644 reference and 2,856 held-out cells per split. Each of the 20
primary datasets contains 1,024 pseudo-spots with approximately ten held-out
cells per spot.

Reference and held-out cells are disjoint. Feature selection uses reference
cells only. Exact pseudo-spot proportions and source-cell membership are
retained, making RMSE and Jensen-Shannon divergence valid accuracy endpoints.

### 3.3 GSE129785: human blood and sorted immune cells

The selected scope contains 30 scATAC fragment samples: 14 physical dilution
mixtures, 9 sorted immune reference populations, 4 PBMC replicates, and 3
fresh/frozen preparation samples.

The acquisition is 40,071,277,814 bytes (37.319 GiB) across 43 frozen
resources. The fragments contain 3,112,258,770 rows and 497,297 retained
barcodes. The standardized reference contains 14,688 author-called cells
across nine immune types and exactly 5,000 reference-only peaks.

The analysis stays in hg19. A coordinate audit selected
`right_cut_offset=-1` and exactly reproduced all 9,528 tested author-matrix
entries. Sixteen runner-ready datasets were created: 14 dilutions, one
four-replicate PBMC cohort, and one three-preparation cohort. Physical ratios
are nominal validation evidence, not exact `truth/proportions.csv`.

### 3.4 GSE194122: human bone-marrow mononuclear cells

The Multiome arm contains 69,249 cells, 10 donors, 13 donor/site batches, 22
author cell labels, 116,490 ATAC peaks, and 13,431 genes. Forty-one source
files total 34,604,880,348 validated bytes. The fragment suite is hard-linked
into reusable processed storage to avoid a second physical copy.

All 69,249 cells map one-to-one through the author's per-barcode metrics bridge
to fragment barcodes. A source-only coordinate audit selected `chromEnd` with
no offset. Ten leave-one-donor-out folds keep every site for a donor on the
same side. Fold-specific features and signatures use only the other nine
donors. Two mixture seeds under two abundance conditions produce 40
exact-truth datasets with 1,024 spots and seven harmonized cell types each.

### 3.5 GSE205055: mouse embryo, mouse brain, and human hippocampus

The complete author-processed family contains 22 samples and 38 supplementary
files. Its immutable raw archive is 8,177,080,320 bytes. Six ATAC-bearing
sections total 29,966 pixels and 481,802,842 ATAC fragments:

| Section | Organ/stage | Resolution or design | ATAC pixels |
|---|---|---|---:|
| GSM6204623 | Mouse brain, P21 | 20 µm | 2,500 |
| GSM6204624 | Mouse whole embryo, E13 | 25 µm | 9,966 |
| GSM6758284 | Mouse brain, P21 replicate | 20 µm | 2,500 |
| GSM6758285 | Mouse brain, P22 replicate | 100-barcode aggregation | 10,000 |
| GSM6801813 | Mouse whole embryo, E13 | 50 µm | 2,500 |
| GSM6206884 | Adult human hippocampus | 50 µm | 2,500 |

Matched RNA, histone/CUT&Tag, coordinates, and other deposited validation
materials were retained but not used to fit ShapeMix.

### 3.6 GSE263333: mouse embryo and EAE brain

The complete series contains 12 samples and 32 supplementary files. Its raw
archive is 3,306,485,760 bytes. Two sections contain spatial ATAC:

| Section | Organ/stage | Resolution | ATAC pixels | ATAC fragments |
|---|---|---:|---:|---:|
| GSM8189706 | Mouse whole embryo, E13 | 50 µm | 2,500 | 108,080,932 |
| GSM8494157 | Five-month EAE mouse brain | 20 µm | 10,000 | 25,588,029 |

RNA, histone, and selected protein measurements were retained for orthogonal
validation. The five-month section is treated as EAE brain based on the
publication and deposited filenames; the conflicting GEO tissue label remains
documented in the source audit.

### 3.7 Real-spatial references

Three independent references were built instead of using a universal atlas:

- GSE216371 E13.5 mouse embryo: 104,194 cells, 12 major types, and 5,000
  reference-only cCREs;
- GSE246791 adult mouse brain: nine broad types and 4,996 positive-support
  intervals; and
- GSE244618 adult human hippocampus: six broad types and 5,000 reference-only
  features.

The EAE section intentionally uses the adult mouse-brain reference with its
disease/stage mismatch recorded. Spatial outcomes did not influence reference
cells, labels, or feature selection.

### 3.8 Three-bin AnnData contract

Every ShapeMix object contains these ordered sparse layers:

```text
fragment_length_lt_100       [0, 100) bp
fragment_length_100_249      [100, 250) bp
fragment_length_ge_250       [250, infinity) bp
```

Each deduplicated fragment contributes two Tn5 cut sites. `readSupport` is not
used as a weight. Each cut is assigned to its containing selected peak. In
every prepared object, `.X` is the exact elementwise sum of the three layers.
The real-spatial audit rechecked this identity across all eight sections.

The AnnData inputs retain the ordered feature axis, cell-type universe,
fragment semantics, feature hash, coordinate convention, spatial coordinates,
and links to the validation evidence required by later analyses.

## 4. Experimental design

### 4.1 Methods

Every unit received:

1. `shapemix_length`: negative-binomial total counts plus a conditional
   multinomial likelihood over the three fragment-length bins;
2. `shapemix_count_only`: the same signatures, prior, initialization,
   optimizer, seeds, and total-count likelihood, with only `use_shape: false`;
3. `nnls`: nonnegative least squares on collapsed peak totals.

The contrast is `shapemix_length - shapemix_count_only`; negative accuracy
effects favor fragment length. Count-only is a controlled ablation, not a
separately tuned method.

### 4.2 Fitting contract

Production ShapeMix uses fixed signatures, a Gamma abundance prior, three
deterministic restarts, Adam, at most 2,000 steps, patience 100, and tolerance
`1e-5`. It reports row-normalized proportions from positive effective
abundance. It has no spatial prior, observed-depth offset, learned background,
or posterior credible intervals.

Each run preserves predictions, abundances, resolved inputs, environment,
restart/convergence/reconstruction diagnostics, and an SHA-256 manifest over
all immutable artifacts.

### 4.3 Evidence-specific endpoints

| Evidence class | Permitted interpretation |
|---|---|
| Exact pseudo-spot truth | RMSE, JSD, per-type error, rare detection, reconstruction, runtime |
| GSE129785 nominal ratios | Descriptive nominal RMSE/JSD, calibration, rare recovery, off-target mass |
| Prediction-only cohorts | Replicate/preparation stability without accuracy claims |
| Real spatial | Map stability, continuity, boundaries, replicate consistency, cross-modality correlation, mismatch warnings |

Spatial pixels are not biological replicates. GSE194122 donors are the only
population-level units. PBMC split and mixture seeds quantify conditional
resampling within one donor.

## 5. GPU qualification and resource management

### 5.1 CUDA qualification

CUDA qualification v2 passed on the NVIDIA GeForce RTX 3080.

| Check | Result |
|---|---:|
| Full-size cached CUDA speedup over CPU | 4.082× |
| CUDA repeat maximum proportion difference | 0 |
| Cached-versus-streamed CUDA maximum difference | 0 |
| Restart-zero CPU/CUDA maximum proportion difference | `9.16081e-5` |
| Restart-zero CPU/CUDA metric difference | `8.03279e-9` |
| Restart-zero CUDA speedup | 4.173× |

The direct three-restart CPU/CUDA comparison failed because the devices chose
different local optima. That comparison remains preserved. A fixed restart-zero
check separated backend arithmetic from non-convex restart selection; no
threshold was relaxed.

Small GSE129785 inputs stayed on paired CPU configs. GSE194122, PBMC stress,
and all real-spatial ShapeMix pairs used deterministic CUDA with full GPU count
caching. Both arms of every pair always used the same backend.

### 5.2 Load containment and monitoring

The initial co-tenant launcher limited CPU threads and refused launches during
competing scheduler activity, excessive load, GPU use, temperature, or memory
pressure. After full-machine use was authorized, execution resumed under the
exclusive profile while retaining the project lock and GPU, temperature, and
memory gates.

Execution resumed only from hash-verified completed runs. Periodic monitoring
tracked processes, CPU, GPU, output growth, current-stage ETA, and total ETA.
Low CPU use during GPU-bound fitting was expected. The final audit found no
remaining driver, worker, or summarizer.

### 5.3 Recorded compute cost

These are sums of per-run wall times, not calendar elapsed time:

| Campaign | Length-aware hours | Count-only hours | NNLS hours recorded |
|---|---:|---:|---:|
| PBMC primary | 19.270 | 11.030 | Not in historical resource table |
| GSE129785 | 0.310 | 0.121 | 0.007 |
| GSE194122 | 7.242 | 3.932 | 0.124 |
| PBMC stress v2 | 6.975 | 3.843 | 0.030 |
| GSE205055 | 13.853 | 7.585 | 0.010 |
| GSE263333 | 9.007 | 4.872 | 0.004 |
| **Recorded total** | **56.657** | **31.382** | **0.174 plus primary NNLS** |

Length-aware fitting was approximately 1.8 times slower in external full-size
campaigns. The largest job was EAE brain: 8.37 hours length-aware versus 4.56
hours count-only. Peak process-tree RSS ranged from about 1.6 GiB in the
historical PBMC summary to about 9.5 GiB for CPU GSE129785.

## 6. Results

### 6.1 PBMC protocol-v1 primary

| Condition | Endpoint | Mean length-minus-count effect | Bootstrap 95% interval | Outer splits favoring length |
|---|---|---:|---:|---:|
| Observed abundance | `rmse_v1` | +0.000025 | [-0.000107, +0.000126] | 2/5 |
| Observed abundance | `jsd_v2` | +0.005767 | [+0.005492, +0.006042] | 0/5 |
| Equal cell type | `rmse_v1` | +0.005230 | [+0.005017, +0.005443] | 0/5 |
| Equal cell type | `jsd_v2` | +0.027054 | [+0.026397, +0.027640] | 0/5 |

The preregistered support rule was not met. Observed-abundance RMSE was nearly
neutral, but JSD worsened in all five outer splits. Both endpoints worsened in
all five equal-cell-type splits. Rare-cell metrics decreased. NNLS was worse
than count-only on both primary endpoints in every outer split.

### 6.2 GSE194122 ten-donor exact-truth evaluation

| Condition | Endpoint | Mean length-minus-count effect | Donor 95% interval | Donors favoring length |
|---|---|---:|---:|---:|
| Observed abundance | `rmse_v1` | +0.000068 | [-0.000958, +0.001094] | 5/10 |
| Observed abundance | `jsd_v2` | +0.000298 | [-0.001828, +0.002424] | 7/10 |
| Equal cell type | `rmse_v1` | +0.000698 | [-0.001521, +0.002916] | 5/10 |
| Equal cell type | `jsd_v2` | +0.001981 | [-0.002729, +0.006692] | 4/10 |

Count-only and length-aware mean RMSE were 0.098945 and 0.099013 under
observed abundance, and 0.111887 and 0.112585 under equal cell type. NNLS mean
RMSE was 0.118936 and 0.130711. Effects were heterogeneous; all four intervals
include zero, so there is no evidence of an average donor-level benefit.

### 6.3 GSE129785 nominal physical mixtures

| Mixture family | Method | Mean nominal RMSE | Mean nominal JSD | Mean off-target mass | Rare levels detected at half nominal |
|---|---|---:|---:|---:|---:|
| CD4-memory/CD8-naive | Length-aware | 0.213103 | 0.154030 | 0.252316 | 5/6 |
| CD4-memory/CD8-naive | Count-only | 0.163002 | 0.114483 | 0.200221 | 3/6 |
| CD4-memory/CD8-naive | NNLS | 0.139709 | 0.095321 | 0.170793 | 2/6 |
| Monocyte/T-cell | Length-aware | 0.065097 | 0.040204 | 0.061927 | 6/6 |
| Monocyte/T-cell | Count-only | 0.049937 | 0.029264 | 0.049225 | 6/6 |
| Monocyte/T-cell | NNLS | 0.045598 | 0.029759 | 0.054060 | 3/6 |

Length-aware RMSE and JSD were worse in all seven CD4-memory/CD8-naive
samples. It improved both endpoints in 4/7 monocyte/T-cell samples, but larger
degradations made its means worse. Extra rare-component detection came with
higher overall error and off-target mass. These remain descriptive nominal,
not exact-truth, comparisons. PBMC replicate and preparation predictions are
retained separately without accuracy claims.

### 6.4 PBMC diagnostic stress v2

Across 40 paired dataset-seed comparisons, length improved RMSE in 14, JSD in
6, and both together in 6. The joint improvements were both seeds at 25%
retained depth, both seeds at `reference_support=all`, and both seeds for
`subtype=broad3`.

The clearest benefit was at 25% depth: mean RMSE effect -0.001486 and JSD
effect -0.003131. At 50% depth, RMSE was neutral and JSD worsened; at 75%, both
means worsened. Length did not help consistently with 1,000/2,500 peaks, five
bins, low reference support, very small spots, or fine CD4 subtypes. Rare-NK
AUPRC was lower than count-only in all eight observed/0.1%/0.5%/1%
comparisons across two seeds.

This suggests a possible low-depth or coarse-label niche, but it does not
overturn the primary and external-average results.

### 6.5 Real-spatial map stability

These sections have no exact composition truth. Map comparisons measure how
much the shape term changes predictions, not which map is accurate.

| Family/section | Informative spots | Median Pearson | Median Spearman | Median mean absolute difference |
|---|---:|---:|---:|---:|
| GSE205055 human hippocampus 50 µm | 2,500 | 0.982163 | 0.976788 | 0.013984 |
| GSE205055 mouse brain P21 primary | 2,500 | 0.955070 | 0.927009 | 0.014307 |
| GSE205055 mouse brain P21 replicate | 2,500 | 0.938412 | 0.944982 | 0.014314 |
| GSE205055 mouse brain P22 100-barcode | 10,000 | 0.957649 | 0.936789 | 0.021333 |
| GSE205055 embryo E13 25 µm | 9,901 | 0.980388 | 0.985274 | 0.006384 |
| GSE205055 embryo E13 50 µm | 2,500 | 0.978130 | 0.980143 | 0.007671 |
| GSE263333 embryo E13 50 µm | 2,500 | 0.979339 | 0.983355 | 0.012992 |
| GSE263333 EAE brain 20 µm | 9,930 | 0.915882 | 0.894848 | 0.014081 |

The arms are highly similar overall. EAE brain changes most, but that does not
establish greater accuracy. Boundary Jaccard medians were 0.678 for GSE205055
and 0.654 for GSE263333. Median Moran's I was 0.258 count-only versus 0.245
length-aware in GSE205055, and 0.280 versus 0.272 in GSE263333. Both ShapeMix
arms were smoother than NNLS without using a spatial prior.

The P21 replicate comparison was weak for all methods: median cell-type
Spearman was -0.010 for NNLS, 0.031 for count-only, and 0.025 for length-aware.

### 6.6 Real-spatial orthogonal validation

Reference marker panels were frozen before spatial outcomes were examined.

For GSE205055, the paired length-minus-count Spearman delta was +0.000177 for
RNA (32/56 improved), -0.005530 for reference-ATAC markers (15/57 improved),
and +0.000481 for epigenome markers (25/45 improved).

For GSE263333, the delta was -0.004657 for RNA (7/21 improved), -0.007082 for
reference-ATAC markers (7/21 improved), and -0.000091 for epigenome evidence
(21/42 improved). Protein evidence covered only three cell types.

The deltas are small and change direction across evidence classes and
families. Absolute correlations are generally low. The evidence does not show
a consistent length-aware advantage.

### 6.7 Reconstruction warnings

All 24 real-spatial method/section combinations crossed the preregistered
warning threshold. Median proxies were about 1.49 for ShapeMix and 0.87 for
NNLS in GSE205055, and 1.47 versus 0.88 in GSE263333. The proxy is not
identified off-reference mass, but universal warnings reinforce caution about
disease, developmental, protocol, reference, and model mismatch.

## 7. Summary-stage issues and resolutions

### 7.1 Exact zero-signal spots

The first spatial summary attempt failed before writing evidence because NNLS
followed `zero_policy: zeros` at 65 GSE205055 embryo 25-µm spots and 70
GSE263333 EAE-brain spots. Every NNLS zero row exactly matched a spatial input
row with zero total count on the registered feature axis and the ShapeMix
NNLS-fallback spot list. No unexplained invalid row existed.

A zero-signal spot has no data-identified composition. The dated amendment
therefore excludes those spots from all spotwise endpoints for all methods and
rebuilds spatial graphs on informative spots. It does not normalize, delete,
or rewrite predictions. Positive-signal rows must remain finite, nonnegative,
and unit-sum within `1e-6`.

Each family summary includes `spot_exclusions.csv` with the spot IDs, input
totals, method row sums, reason, and affected endpoints. The exact-truth metric
contract is unchanged.

### 7.2 Unsupported frozen marker panel

The next fail-closed attempt found that the GSE205055 embryo 50-µm RNA matrix
contains only `Klf1` from the four-marker Erythroid panel, below the two-feature
minimum. No other combination failed the complete marker availability audit.

The panel and threshold were not changed. The endpoint remains as three
explicit unavailable rows, one per method, with present/missing markers and no
correlation. GSE205055 has 474 available and 3 unavailable cross-modality rows;
GSE263333 has 261 available rows and none unavailable.

Both amendments followed summary errors and preceded any written evidence
table or reviewed biological effect. They are recorded in the validation
config, protocol, plan, tests, and summaries.

## 8. Validation and reproducibility

### 8.1 Run-level validation

All 372 run directories passed their SHA-256 output manifests:

| Campaign | Verified runs |
|---|---:|
| PBMC primary | 60 |
| GSE129785 | 48 |
| GSE194122 | 120 |
| PBMC stress v2 | 120 |
| GSE205055 | 18 |
| GSE263333 | 6 |

Across all 24 real-spatial jobs, the final invariant audit confirmed exact
axes; finite, nonnegative outputs; unit-sum informative rows; exact NNLS
zero/input-zero agreement; `.X == sum(layers)`; 16/16 successful ShapeMix
fits with three converged restarts and no nonfinite events; deterministic RTX
3080 full-CUDA execution; exact zero shape likelihood for count-only; and
239,904 complete peak-by-bin residual rows.

### 8.2 Summary and synthesis validation

Both spatial summaries report `complete`, and every declared row count and
hash was recomputed. The final synthesis config, source summaries, outputs,
and row counts also passed hash verification. Canonical file-layout validation
passed for every external, stress, and spatial campaign. New results are
visible to Git rather than ignored.

### 8.3 Tests

- Real-spatial preparation and summary: 34/34 passed.
- Full-evaluation synthesis: 3/3 passed.
- Repository-wide suite: 454 passed, 1 expected skip, 0 failures.
- `git diff --check`: passed.

## 9. Result artifacts

### 9.1 Final synthesis

```text
results/external_validation/shapemix_full_evaluation_v2/
  evidence_summary.yaml
  evidence_table.tsv
  effect_table.tsv
  resource_table.tsv
```

| Output | SHA-256 |
|---|---|
| `evidence_table.tsv` | `62dac1d1d3c898234270f37605bf73fe5eb09b9c520a96ef70ed2828d288b3bc` |
| `effect_table.tsv` | `f1ac083e7a82d64faa56d4c24bf2b5e9a2f44c9a4653f2c4b9ab40584e02a07b` |
| `resource_table.tsv` | `821272e3fa811ad49c1297bf20293d8c43ae410905d924b2716b13d7129580f6` |

### 9.2 Campaign roots

```text
results/primary/shapemix_primary_ablation_protocol_v1_summary/
results/external_validation/shapemix_gse129785_v2/
results/external_validation/shapemix_gse194122_lodo_v1/
results/sensitivity/shapemix_pbmc_stress_v2/
results/real_spatial/shapemix_gse205055_v1/
results/real_spatial/shapemix_gse263333_v1/
```

Spatial campaign summaries include map concordance, continuity, boundary
agreement, cross-modality concordance, replicate consistency, spot exclusions,
reconstruction warnings, resources, and their evidence manifest.

### 9.3 Governing documents

- [Model specification](model_specification.md)
- [Benchmark protocol](benchmark_protocol.md)
- [Full-evaluation protocol](full_evaluation_protocol.md)
- [Implementation plan](implementation_plan.md)
- [File organization](file_organization.md)
- [Concise results summary](results_summary.md)
- [Step 6 primary report](step6_results.md)

## 10. Scientific interpretation

The broader hypothesis is not supported for the current model. The conclusion
now spans one-donor PBMC resampling, ten held-out BMMC donors, physical immune
mixtures, diagnostic perturbations, and eight tissue sections across human
blood, human bone marrow, human hippocampus, mouse brain, diseased mouse brain,
and mouse embryo.

The credible positive signal is the repeatable 25%-depth improvement. Shape
may add information when count evidence is extremely sparse, and coarse-label
benefits may mean length patterns transfer better at broad lineage resolution.
These are diagnostic interpretations, not a basis for production deployment.

Likely limitations of v1 include an overconfident fixed multinomial shape
term, transfer mismatch in fixed signatures, no explicit off-reference
background, weak shape separation among fine subtypes, roughly 80% extra
optimization time, and non-convex restart selection.

## 11. Recommended next work

The completed campaign should remain immutable. Follow-ups need new versioned
configs and result roots.

1. Keep count-only ShapeMix as the operational v1 default.
2. Confirm the 25%-depth result on independent donors.
3. Predeclare a low-depth donor-level experiment with matched count-only and
   NNLS controls and a fixed directional rule.
4. Test a Dirichlet-multinomial or other overdispersed shape likelihood.
5. Test preregistered shape weighting or hierarchical shrinkage toward
   count-only without tuning on these outcomes.
6. Add an identifiable background/off-reference component before stronger
   disease or cross-protocol tissue claims.
7. Improve reference matching for EAE and developmental tissue under new IDs.
8. Seek independent per-spot composition truth or controlled spatial mixtures;
   map agreement and marker correlation cannot rank accuracy.
9. Preserve the negative result, zero-signal audit, unavailable marker
   endpoint, and failed direct multi-restart CPU/CUDA comparison in reports.

## 12. Completion and repository state

The recurring evaluation goal is complete. Recorded usage was 37,441,273
tokens over approximately 7 days and 10 hours of elapsed goal time, including
waiting and monitoring. No ShapeMix process remains active.

The new documentation and result files are present but uncommitted. Committing
or pushing was not requested. The result files are visible to Git and ready
for review.
