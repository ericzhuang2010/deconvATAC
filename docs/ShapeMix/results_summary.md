# ShapeMix full-evaluation results summary

Status: complete 2026-09-09

Canonical synthesis:
`results/external_validation/shapemix_full_evaluation_v2/`

## Main conclusion

The frozen three-bin parent-fragment-length term is **not a reliable overall
improvement** over the matched count-only ShapeMix model. Effects below are
length-aware minus count-only, so negative values favor fragment length.

- In the original one-donor PBMC benchmark, fragment length worsened JSD in
  every outer split and worsened equal-cell-type RMSE in every split.
- In the ten-donor GSE194122 held-out evaluation, average effects were small
  and slightly unfavorable; all donor-level 95% intervals included zero.
- In GSE129785 physical mixtures, length-aware fits had worse mean nominal
  RMSE, JSD, and off-target mass than count-only for both mixture families.
- The PBMC stress campaign found limited niches: both exact-truth endpoints
  improved at 25% retained depth and in the coarse/full-support reference
  views, but JSD improved in only 6/40 dataset-seed comparisons overall and
  rare-NK AUPRC was lower in all eight evaluated settings.
- In eight real spatial sections, length-aware and count-only maps were highly
  similar. Orthogonal RNA, protein, histone, and marker concordance was weak
  and inconsistent, so these data do not establish an accuracy improvement.

The evidence supports retaining the count-only model as the default for this
model version. Fragment length remains a targeted research direction for
low-depth or deliberately coarsened label settings, not a generally validated
replacement.

## Completed evaluation inventory

| Stage | Evidence | Runnable units | Successful jobs |
|---|---|---:|---:|
| PBMC protocol-v1 primary | Exact held-out pseudo-spot truth | 20 | 60 |
| GSE129785 | Nominal physical mixtures and prediction-only preparations | 16 | 48 |
| GSE194122 | Ten-donor held-out pseudo-spot truth | 40 | 120 |
| PBMC diagnostic stress v2 | Exact held-out pseudo-spot truth | 40 | 120 |
| GSE205055 | Six real spatial sections without composition truth | 6 | 18 |
| GSE263333 | Two real spatial sections without composition truth | 2 | 6 |
| **Total** | Five scientific evidence classes | **124** | **372** |

Every runnable unit received `shapemix_length`, `shapemix_count_only`, and
NNLS. All 372 jobs completed successfully and all 372 SHA-256 output manifests
were revalidated after synthesis.

## Exact-truth results

### PBMC protocol-v1 primary

| Condition | Endpoint | Mean effect | Bootstrap 95% interval | Outer splits favoring length |
|---|---|---:|---:|---:|
| Observed abundance | `rmse_v1` | +0.000025 | [-0.000107, +0.000126] | 2/5 |
| Observed abundance | `jsd_v2` | +0.005767 | [+0.005492, +0.006042] | 0/5 |
| Equal cell type | `rmse_v1` | +0.005230 | [+0.005017, +0.005443] | 0/5 |
| Equal cell type | `jsd_v2` | +0.027054 | [+0.026397, +0.027640] | 0/5 |

Neither condition met the preregistered directional-support rule. NNLS was
worse than count-only for both endpoints in every outer split.

### GSE194122 ten-donor held-out evaluation

| Condition | Endpoint | Mean effect | Donor-level 95% interval | Donors favoring length |
|---|---|---:|---:|---:|
| Observed abundance | `rmse_v1` | +0.000068 | [-0.000958, +0.001094] | 5/10 |
| Observed abundance | `jsd_v2` | +0.000298 | [-0.001828, +0.002424] | 7/10 |
| Equal cell type | `rmse_v1` | +0.000698 | [-0.001521, +0.002916] | 5/10 |
| Equal cell type | `jsd_v2` | +0.001981 | [-0.002729, +0.006692] | 4/10 |

Count-only and length-aware mean RMSE were respectively 0.098945 and 0.099013
under observed abundance, and 0.111887 and 0.112585 under equal cell type.
NNLS was worse in both conditions. The donor results therefore show
heterogeneity but no population-level average benefit from fragment length.

### PBMC diagnostic stress v2

Across the 40 paired dataset-seed comparisons, length improved RMSE in 14 and
JSD in 6; both endpoints improved together in 6. Those six were the two seeds
  at 25% retained depth and the two seeds in each of the
coarse/full-support reference views (`reference_support=all` and
`subtype=broad3`). The largest consistent benefit was at 25% depth:
mean RMSE effect -0.001486 and mean JSD effect -0.003131.

The benefit did not generalize. JSD worsened at 50% and 75% retained depth,
with 1,000 or 2,500 peaks, for five bins, at low reference support, and for the
fine CD4-related subtype task. Rare-NK AUPRC was lower for length-aware
ShapeMix than count-only in all eight observed/0.1%/0.5%/1% comparisons.

## GSE129785 nominal physical mixtures

These values compare predictions with nominal input ratios; they are
descriptive and are not exact sorted-cell composition truth.

| Mixture family | Method | Mean nominal RMSE | Mean nominal JSD | Mean off-target mass | Rare levels detected at half nominal |
|---|---|---:|---:|---:|---:|
| CD4-memory/CD8-naive | Length-aware | 0.213103 | 0.154030 | 0.252316 | 5/6 |
| CD4-memory/CD8-naive | Count-only | 0.163002 | 0.114483 | 0.200221 | 3/6 |
| CD4-memory/CD8-naive | NNLS | 0.139709 | 0.095321 | 0.170793 | 2/6 |
| Monocyte/T-cell | Length-aware | 0.065097 | 0.040204 | 0.061927 | 6/6 |
| Monocyte/T-cell | Count-only | 0.049937 | 0.029264 | 0.049225 | 6/6 |
| Monocyte/T-cell | NNLS | 0.045598 | 0.029759 | 0.054060 | 3/6 |

Length-aware RMSE and JSD were worse in all seven CD4-memory/CD8-naive
samples. For monocyte/T-cell mixtures, length improved both metrics in 4/7
samples but had worse means because the degradations were larger.

## Real-spatial results

There is no exact per-spot composition truth for GSE205055 or GSE263333.
Consequently these results measure map stability and orthogonal concordance,
not deconvolution accuracy.

Across the eight sections, the median cell-type map correlation between
length-aware and count-only predictions was high:

- section-level median Pearson range: 0.915882 to 0.982163;
- section-level median Spearman range: 0.894848 to 0.985274; and
- section-level median mean-absolute-difference range: 0.006384 to 0.021333.

The EAE-brain section showed the largest change (median Pearson 0.915882,
Spearman 0.894848), while both embryo families were especially stable.
Boundary-edge Jaccard medians were 0.678 for GSE205055 and 0.654 for GSE263333.
The P21 replicate-pair median Spearman correlations were low for every method:
-0.010 for NNLS, 0.031 for count-only, and 0.025 for length-aware ShapeMix.

Orthogonal concordance did not show consistent gains. In GSE205055, the paired
length-minus-count RNA Spearman delta had median +0.000177 (32/56 cell-type
comparisons improved); the reference-ATAC marker delta had median -0.005530
(15/57 improved). In GSE263333, the RNA delta had median -0.004657 (7/21
improved), the reference-ATAC marker delta -0.007082 (7/21 improved), and the
histone/epigenome delta -0.000091 (21/42 improved). Protein evidence comprised
only three cell types and was not decisive.

Every real-spatial method/section run crossed the preregistered reconstruction
mismatch-warning threshold. These are warning proxies, not identified
off-reference mass, but they reinforce the need for caution about disease,
protocol, and reference mismatch.

The summary explicitly excludes 65 zero-signal spots in the GSE205055 embryo
25-um section and 70 in the GSE263333 EAE-brain section from all spotwise
comparisons for every method. Each excluded spot has zero total input on the
registered feature axis; IDs and per-method row sums are preserved in
`spot_exclusions.csv`. Raw predictions were not changed. One GSE205055 embryo
50-um RNA/Erythroid endpoint is explicitly unavailable because only one of four
frozen markers is present, below the two-feature minimum.

## GPU and execution results

CUDA qualification v2 passed. On the representative full-size input, cached
CUDA was 4.082 times faster than CPU, CUDA repeats were exact, and cached and
streamed CUDA outputs were exact. A deterministic restart-zero comparison
showed a 4.173-times speedup, maximum CPU/CUDA proportion difference
`9.16e-5`, and metric difference `8.03e-9`. Direct multi-restart CPU/CUDA fits
selected different local optima; that failed comparison remains preserved.

Small GSE129785 inputs therefore stayed on paired CPU configs. GSE194122, PBMC
stress, and all real-spatial paired ShapeMix arms used the qualified RTX 3080
CUDA backend with deterministic algorithms and full CUDA count caching. In the
real-spatial campaign, length-aware fits consumed 22.86 cumulative hours versus
12.46 hours for count-only (about 1.83 times longer). The largest single run was
EAE-brain: 8.37 hours length-aware versus 4.56 hours count-only.

## Integrity and artifact map

- Final evidence index: `results/external_validation/shapemix_full_evaluation_v2/evidence_table.tsv`
- Non-pooled effect table: `results/external_validation/shapemix_full_evaluation_v2/effect_table.tsv`
- Normalized resource table: `results/external_validation/shapemix_full_evaluation_v2/resource_table.tsv`
- Final hash manifest: `results/external_validation/shapemix_full_evaluation_v2/evidence_summary.yaml`
- GSE205055 evidence: `results/real_spatial/shapemix_gse205055_v1/shapemix_gse205055_real_spatial_protocol_v1_cuda/`
- GSE263333 evidence: `results/real_spatial/shapemix_gse263333_v1/shapemix_gse263333_real_spatial_protocol_v1_cuda/`

The synthesis contains seven separately labeled campaign entries, 82 effect
rows, and 352 normalized resource rows. Exact pseudo-truth, nominal physical
ratios, backend qualification, and qualitative real-spatial evidence are never
pooled. Canonical file-layout validation passed for all external, stress, and
real-spatial campaigns. The targeted real-spatial and synthesis regression
suites passed 34/34 and 3/3 tests, respectively. The final repository-wide
suite passed 454 tests with one expected skip and no failures.

## Interpretation limits

The primary PBMC and PBMC stress data come from one donor. GSE194122 supports
donor-level inference for its held-out pseudo-spots but does not provide real
spatial tissue truth. GSE129785 ratios are nominal physical inputs, and the
GSE205055/GSE263333 sections have no exact composition truth. The results apply
to the fixed-signature, three-bin MAP model and frozen reference tracks; they
do not rule out other fragment-shape representations, learned signatures,
background models, or stronger matched references.

See [the protocol](full_evaluation_protocol.md),
[the implementation plan](implementation_plan.md), and
[the Step 6 report](step6_results.md) for definitions and provenance.
