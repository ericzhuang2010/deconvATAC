# Composition labels and evaluation cases

## What a composition label means

For a truth-based deconvolution evaluation, the label is a spot-by-cell-type
composition matrix:

```text
rows    = spots or mixture samples
columns = cell types
values  = true cell-type proportions
```

Each row sums to one. This is different from a **single-cell label**, which
only identifies the type of an individual reference cell. Labeled reference
cells let ShapeMix learn cell-type signatures, but they do not by themselves
reveal the cell-type composition of a real spatial spot.

## Which datasets have composition labels?

| Dataset | Spot composition available? | Source of composition | Permitted interpretation |
|---|---|---|---|
| PBMC protocol-v1 primary | Yes, exact pseudo-spot truth | Counts of sampled, held-out labeled PBMC cells | Quantitative accuracy using RMSE, JSD, per-type error, and rare-cell metrics |
| PBMC diagnostic stress-v2 | Yes, exact pseudo-spot truth | Counts of sampled, held-out labeled PBMC cells under each stress condition | Quantitative accuracy within the diagnostic sensitivity design |
| GSE194122 | Yes, exact pseudo-spot truth | Counts of sampled cells from the completely held-out donor | Quantitative accuracy and donor-level effects |
| GSE129785 physical dilutions | Nominal sample-level proportions only | Intended laboratory mixing ratios | Descriptive agreement with nominal ratios, calibration, rare-component recovery, and off-target mass |
| GSE129785 PBMC/preparation cohorts | No | No known mixture ratio or source-cell membership | Replicate and preparation stability only |
| GSE205055 | No | Real spatial tissue with matched proxy modalities | Map stability, spatial behavior, replicate consistency, and orthogonal concordance only |
| GSE263333 | No | Real spatial tissue with matched proxy modalities | Map stability, spatial behavior, and RNA/protein/histone concordance only |

The evidence classes are intentionally kept separate. Exact pseudo-spot truth,
nominal physical inputs, prediction-only cohorts, and real-spatial proxy
evidence must not be pooled into a single accuracy estimate.

## Case 1: exact held-out pseudo-spot truth

### Construction

The PBMC and GSE194122 experiments start with cell-resolved ATAC or Multiome
data in which every retained cell has a fixed cell-type annotation.

```text
labeled single cells
        |
        +--> training/reference cells --> peaks and fixed signatures
        |
        +--> held-out cells -----------> aggregate into pseudo-spots
                                               |
                                               +--> ATAC input for inference
                                               +--> label counts for truth
```

If pseudo-spot `s` contains `N_s` held-out cells and `N_sc` of them have cell
type `c`, its truth is

```text
truth[s, c] = N_sc / N_s.
```

For example, a pseudo-spot containing four CD14 monocytes, one CD4-naive cell,
one CD4-TCM cell, and one Treg has the composition

```text
CD14 Mono = 4/7
CD4 Naive = 1/7
CD4 TCM   = 1/7
Treg      = 1/7
all other declared types = 0
```

During inference, ShapeMix sees only the aggregated pseudo-spot ATAC counts
and fragment-length layers. It does not see the held-out source-cell identities
or their labels. After inference, the saved source-cell labels are used to
construct the truth matrix and score the predictions.

### PBMC protocol-v1 primary benchmark

**Protocol-v1 is the frozen design of the original, primary PBMC benchmark.**
It is a version of the experimental contract, not a cell subset, model arm, or
metric. The contract fixed the data split, pseudo-spot generator, feature
selection, model comparison, metrics, seeds, and interpretation rules before
the primary results were examined.

The PBMC source contains 9,500 retained cells across 16 immune cell types. For
each protocol-v1 outer split, 70% of the cells are used to select peaks and
estimate the reference signatures; the remaining 30% are used to construct
pseudo-spots. Inference is performed on the aggregated pseudo-spots, not on the
individual held-out cells.

The standard protocol-v1 setting is:

| Component | Protocol-v1 setting |
|---|---|
| Source | One PBMC Multiome donor with 9,500 retained labeled cells |
| Cell-type universe | 16 immune cell types |
| Reference/test separation | Stratified 70% reference and 30% held-out cells |
| Outer splits | Five different 70/30 split seeds |
| Mixtures within each split | Two mixture seeds |
| Mixture conditions | Observed cell-type abundance and equal cell-type abundance |
| Pseudo-spots per dataset | 1,024 |
| Mean cells per pseudo-spot | Approximately 10 |
| Features | 5,000 peaks selected using reference cells only |
| Fragment-length representation | Three bins: `[0,100)`, `[100,250)`, and `[250,infinity)` bp |
| Depth | Full observed fragment depth; no thinning |
| Primary comparison | Length-aware ShapeMix versus matched count-only ShapeMix |

Five splits multiplied by two mixture seeds and two mixture conditions produce
20 primary pseudo-spot datasets. Protocol-v1 asks one broad question under
these fixed standard conditions:

> Does adding the three-bin fragment-length term improve cell-type composition
> estimates relative to the otherwise identical count-only model?

The protocol-v1 result remains the primary PBMC result. It was not replaced or
rewritten by the later stress campaign.

### PBMC diagnostic stress-v2

Stress-v2 is a **follow-up diagnostic sensitivity experiment** performed after
the protocol-v1 result was negative. It asks a different question:

> Is there a particular data regime in which fragment length becomes helpful
> or especially harmful?

It keeps the exact-source-cell construction: pseudo-spots are still assembled
from held-out labeled cells, their ATAC data are given to the model without the
labels, and their true compositions are calculated afterward by counting the
sampled source-cell labels. Thus, changing a stress condition changes the
model's information or task difficulty; it does not make the truth unknown.

Stress-v2 uses PBMC split 1103, two new mixture seeds (`307` and `401`), and
1,024 pseudo-spots per dataset. Its standard anchor is:

```text
all 16 cell types
observed-abundance mixtures
mean 10 cells per spot
full fragment depth
5,000 peaks
all available reference cells
three fragment-length bins
```

"One factor at a time" means changing one item while holding the other relevant
settings at their control values. For example, the 25%-depth condition keeps
the same 16-type task, 5,000-peak axis, approximately 10 cells per spot, full
reference, and three length bins, but binomially retains only 25% of the
fragment-derived cut-site counts. The sampled cells and truth proportions stay
unchanged. It is therefore interpretable as a depth change rather than a
simultaneous depth, cell-number, and feature-number change.

The tested factor families are:

| Factor | What is changed | Levels or comparison | Question being tested |
|---|---|---|---|
| Retained depth | Fraction of pseudo-spot fragments retained | 25%, 50%, 75%, or full depth | Does fragment length help more when sequencing is sparse? |
| Cells per spot | Expected number of held-out cells aggregated into a spot | Mean 2, 5, 10, or 20 | Does spot mixture complexity affect the value of length information? |
| Rare NK abundance | Probability assigned to NK cells during mixture generation | 0.1%, 0.5%, 1%, or observed PBMC abundance | Does fragment length improve rare-cell recovery? |
| Peak count | Number of reference-selected accessibility features supplied to the model | 1,000, 2,500, or 5,000 peaks | Does fragment length compensate for a smaller feature set? |
| Cell-type task | Which three cell types must be separated | Broad types: CD14 Mono, CD4 Naive, and CD8 Naive; or related types: CD4 Naive, CD4 TCM, and CD4 TEM | Is fragment length more useful for broad or closely related distinctions? |
| Reference support | Maximum reference cells retained per broad cell type | 50, 100, 250, or all available cells per type | Does fragment length help when the reference is small? |
| Length bins | Resolution of the fragment-length representation | Two, three, or five bins | Is the frozen three-bin representation too coarse or too detailed? |

The depth, cells-per-spot, rare-NK, peak-count, and bin families use the
all-16-type observed-abundance anchor. The cell-type-task family compares two
separate equal-abundance three-type tasks. The reference-support family uses
the broad three-type task and compares 50, 100, or 250 reference cells per type
with all available reference cells. Thus, "one factor" is interpreted within
each predeclared factor family and its matching control.

This is not a complete factorial experiment: it does not, for example, test
every combination of 25% depth, two cells per spot, 1,000 peaks, and five bins.
Such a combination would make it unclear which change caused the result and
would require a much larger experiment.

The `v2` in **stress-v2** does not mean that it supersedes PBMC protocol-v1.
There was an initial stress-v1 execution whose restricted reference variants
contained peaks with zero training-reference support. Stress-v2 repaired that
data-contract problem by using one common, support-safe 5,000-peak axis across
the planned reference variants, replacing 174 unsupported peaks without using
prediction or accuracy results. The valid protocol-v1 primary result remained
unchanged. In short:

```text
protocol-v1 = original primary PBMC benchmark
stress-v2   = support-safe version of the later diagnostic stress campaign
```

### GSE194122

GSE194122 begins with 22 author-provided cell labels that are frozen into seven
broad classes. Each fold holds out one complete donor. The other nine donors
provide the reference, while cells from the held-out donor alone are sampled
into pseudo-spots. Consequently, the spot composition is exact by construction
and the donor, rather than an individual spot, is the population-level analysis
unit.

### Accuracy comparison

For exact-truth datasets, every method produces a predicted spot-by-cell-type
matrix aligned to the truth matrix. Let:

```text
S       = number of spots
C       = number of cell types in the fixed evaluation universe
T[s,c]  = true proportion of cell type c in spot s
P[s,c]  = predicted proportion of cell type c in spot s
```

Both matrices must contain the same spot IDs. Truth is evaluated over the
complete, predeclared cell-type universe, not merely the intersection of truth
and prediction columns. A missing declared prediction type is inserted with
proportion zero so that omission is penalized; an unexpected extra prediction
type is an error. All values must be finite and nonnegative, and every truth
and prediction row must sum to one within `1e-6`. Invalid rows are not silently
renormalized or removed.

#### RMSE (`rmse_v1`)

RMSE is calculated over every entry of the two `S x C` matrices:

```text
rmse_v1 = sqrt(
    (1 / (S * C))
    * sum over spots s and cell types c of (T[s,c] - P[s,c])^2
)
```

The calculation therefore proceeds as follows:

1. Subtract the predicted proportion from the true proportion for every
   spot-cell-type pair.
2. Square every difference, so overprediction and underprediction cannot
   cancel.
3. Average the squared differences over all spots and all declared cell types.
4. Take the square root to return to proportion units.

Every spot-cell-type entry has equal weight. A large error for one type in one
spot contributes more because the difference is squared. A perfect prediction
has RMSE zero; larger values indicate greater entrywise error.

#### Jensen-Shannon divergence (`jsd_v2`)

JSD treats each spot's complete composition vector as a probability
distribution. For each spot, first define its midpoint distribution:

```text
M[s,c] = (T[s,c] + P[s,c]) / 2
```

Then calculate the base-2 Kullback-Leibler divergences from truth and
prediction to that midpoint:

```text
KL2(T[s] || M[s]) = sum over c of T[s,c] * log2(T[s,c] / M[s,c])
KL2(P[s] || M[s]) = sum over c of P[s,c] * log2(P[s,c] / M[s,c])

JSD2[s] = 0.5 * KL2(T[s] || M[s])
          + 0.5 * KL2(P[s] || M[s])
```

Terms with a zero numerator contribute zero. The final reported metric is the
unweighted mean of the spot-level divergences:

```text
jsd_v2 = (1 / S) * sum over spots s of JSD2[s]
```

`jsd_v2` ranges from zero to one because logarithms use base 2. Zero means the
two composition distributions are identical. Unlike an asymmetric KL
divergence, JSD treats truth and prediction symmetrically and remains finite
when a cell type is absent from one distribution. This metric is the squared
base-2 Jensen-Shannon distance, not the repository's historical unsquared
`js_distance_v1`.

#### Small numerical example

For one spot with three declared cell types, suppose

```text
truth      T = [0.50, 0.25, 0.25]
prediction P = [0.40, 0.35, 0.25]
```

For RMSE, the differences are `[0.10, -0.10, 0.00]`, and the squared
differences are `[0.01, 0.01, 0.00]`. Therefore:

```text
RMSE = sqrt((0.01 + 0.01 + 0.00) / 3)
     = 0.08165
```

For JSD, the midpoint is `[0.45, 0.30, 0.25]`. Applying the two base-2 KL
calculations above gives:

```text
JSD = 0.01006
```

The numbers are on different scales and should not be subtracted from each
other. They are co-primary endpoints that describe different aspects of the
same prediction: RMSE measures entrywise magnitude error, whereas JSD measures
how different each complete spot composition distribution is.

#### Comparing length-aware and count-only ShapeMix

RMSE and JSD are calculated independently for the length-aware and count-only
predictions against the exact same truth matrix. For metric `M`, the paired
effect for one matched dataset is

```text
delta = M(length-aware) - M(count-only).
```

Because lower error is better, a negative effect favors the fragment-length
model, zero means a tie, and a positive effect favors count-only ShapeMix. For
example, count-only RMSE `0.052` and length-aware RMSE `0.058` give
`delta = +0.006`, meaning the length-aware prediction is worse by 0.006 RMSE.

For PBMC protocol-v1, the two mixture-seed effects are first averaged within
each outer 70/30 split. The five resulting outer-split effects are the
resampling units used for the reported mean and bootstrap interval. Individual
spots are not treated as independent biological replicates. For GSE194122,
the two mixture seeds are similarly summarized within each held-out donor, and
the ten donors are the population-level units.

The same formulas can be applied to the GSE129785 intended mixture vectors,
but those results are called **nominal RMSE** and **nominal JSD** because the
input ratios are not exact recovered-sample compositions. RMSE and JSD are not
calculated for GSE205055 or GSE263333 because those real-spatial datasets have
no spot-level composition truth.

## Case 2: GSE129785 physical dilution mixtures

The GSE129785 dilution samples are real laboratory mixtures rather than
computer-generated pseudo-spots. Two sorted-cell mixture series were assayed:

1. CD4-memory cells mixed with CD8-naive cells.
2. Monocytes mixed with T cells.

Each series contains the intended component ratios
`0.1:99.9`, `0.5:99.5`, `1:99`, `50:50`, `99:1`, `99.5:0.5`, and
`99.9:0.1`. All retained barcodes from one physical dilution sample are
aggregated into one mixture observation for deconvolution.

These are **nominal** input ratios, not exact recovered-cell proportions.
Sorting impurity, counting and pipetting error, differential survival,
nucleus recovery, capture efficiency, and barcode filtering can make the
post-assay mixture differ from the intended input. The mixed-sample barcodes
also lack experimentally known per-cell component identities. Classifying
those barcodes computationally would create another prediction rather than
independent truth.

The CD4-memory/CD8-naive series maps directly to two columns in the model's
reference: `Memory CD4 T cells` and `Naive CD8 T cells`. Although its ratios
remain nominal rather than exact, the names of both input components and model
outputs have a direct correspondence.

The **nine-type reference** is the labeled single-cell ATAC reference used by
ShapeMix for every GSE129785 sample. It was constructed from nine separately
sorted immune populations. The model learns one accessibility and
fragment-shape signature for each reference type and outputs nine proportions
that sum to one:

| Model output | Role in the monocyte/T-cell dilution |
|---|---|
| Dendritic cells | Off-target |
| Monocytes | Intended component |
| B cells | Off-target |
| Regulatory T cells | Part of the intended broad T-cell component |
| Naive CD4 T cells | Part of the intended broad T-cell component |
| Memory CD4 T cells | Part of the intended broad T-cell component |
| NK cells | Off-target |
| Naive CD8 T cells | Part of the intended broad T-cell component |
| Memory CD8 T cells | Part of the intended broad T-cell component |

Thus, the model always predicts **nine cell-type proportions**. Six of those
outputs are relevant to the two intended mixture components: one monocyte
output plus five T-cell-subtype outputs. The remaining three outputs are
off-target populations.

The monocyte/T-cell experimental record is less specific because it contains
only two intended broad components:

```text
Monocytes
T cells
```

The reported ratio is between the monocyte pool and the **entire broad T-cell
pool**. For example, `99% Monocytes : 1% T cells` means that the intended input
contains 1% T cells in total. It does **not** mean 1% of each T-cell subtype,
and it does not identify one particular T-cell subtype as that 1% component.

The experiment does not report how the broad T-cell pool is divided among
subtypes. In contrast, five of the nine model outputs are separate T-cell
predictions:

```text
Regulatory T cells
Naive CD4 T cells
Memory CD4 T cells
Naive CD8 T cells
Memory CD8 T cells
```

The only subtype-related constraint supplied by the nominal mixture record is
therefore

```text
Treg + Naive CD4 + Memory CD4 + Naive CD8 + Memory CD8
    = nominal total T-cell fraction.
```

For a nominal 1% total-T sample, all of the following subtype predictions have
the same evaluable total:

```text
[1%, 0%, 0%, 0%, 0%]
[0%, 0%, 0%, 1%, 0%]
[0.2%, 0.2%, 0.2%, 0.2%, 0.2%]
```

The nominal evidence cannot determine which of these subtype allocations is
biologically correct. It can only determine whether their sum is near 1%.

Therefore, the evaluation collapses the nine model outputs into three scoring
categories: `Monocytes`, `total T`, and `off-target`. It uses three rather than
two categories so that abundance incorrectly assigned outside the intended
mixture cannot disappear during aggregation:

```text
predicted Monocytes = prediction[Monocytes]

predicted total T = prediction[Regulatory T cells]
                  + prediction[Naive CD4 T cells]
                  + prediction[Memory CD4 T cells]
                  + prediction[Naive CD8 T cells]
                  + prediction[Memory CD8 T cells]

predicted off-target mass = 1 - predicted Monocytes - predicted total T
```

Here, off-target mass consists of predictions assigned to dendritic cells,
B cells, or NK cells, because those populations were not intentional
components of the physical mixture.

For example, a nominal `99% Monocytes : 1% T cells` sample is compared at the
broad level with target vector `[0.99, 0.01, 0.00]` for
`[Monocytes, total T, off-target]`. If a model predicts 0.94 monocytes, a total
of 0.04 across the five T-cell subtypes, and 0.02 off-target mass, the evaluated
prediction is `[0.94, 0.04, 0.02]`.

The experiment can therefore test whether the model recovers **total T-cell
abundance**, but it cannot determine whether the model correctly divided that
T-cell abundance into regulatory, naive, memory, CD4, and CD8 components. A
model could predict the correct total T-cell fraction while assigning all of
it to the wrong T-cell subtype. That is why this series is described as an
exploratory broad monocyte-versus-total-T comparison rather than exact
subtype-level validation.

There is a second limitation: even the broad monocyte-versus-total-T ratio is
the intended laboratory input ratio, not an exact count of the cells recovered
after sequencing and quality control. Differential recovery, sorting impurity,
and barcode filtering may shift the actual assayed ratio. Thus, the analysis
compares the summed prediction with a **nominal broad target**, not exact
subtype truth or exact recovered-sample truth.

The counts can be summarized as follows:

```text
experimental input information: 2 broad components
    Monocytes, T cells

raw ShapeMix prediction: 9 cell types
    1 monocyte + 5 T subtypes + 3 off-target types

nominal scoring vector: 3 categories
    Monocytes, total T, off-target
```

Comparisons with the nominal ratios include descriptive RMSE and JSD,
rare-component recovery, calibration, and predicted mass assigned to cell
types that were not intentionally mixed. These results must not be described
as exact-truth accuracy.

## Case 3: GSE129785 prediction-only cohorts

These are seven real PBMC scATAC-seq samples from GSE129785. They were not
constructed from known numbers of the nine reference cell types, and the
source does not provide a nominal mixture ratio for them. They are therefore
kept separate from the 14 physical dilution samples described above.

### Four unsorted PBMC replicates

The replicate cohort contains these four author-named samples:

| Replicate | GEO accession | Author sample title |
|---|---|---|
| 1 | `GSM3722015` | `PBMC_Rep1` |
| 2 | `GSM3722076` | `PBMC_Rep2` |
| 3 | `GSM3722075` | `PBMC_Rep3` |
| 4 | `GSM3722077` | `PBMC_Rep4` |

Here, **unsorted PBMC** means that the assayed sample is a mixed peripheral
blood mononuclear cell population rather than a population purified to contain
only one of the nine reference types. The word **replicate** refers to the four
separately named scATAC-seq samples. It does not mean that their true cell-type
proportions are known or deliberately equal.

For evaluation, all retained cell barcodes and fragments within one replicate
are summed to form one aggregate deconvolution observation:

```text
many barcodes in PBMC_Rep1 -> one aggregate spot -> one 9-type prediction
many barcodes in PBMC_Rep2 -> one aggregate spot -> one 9-type prediction
many barcodes in PBMC_Rep3 -> one aggregate spot -> one 9-type prediction
many barcodes in PBMC_Rep4 -> one aggregate spot -> one 9-type prediction
```

The analysis therefore has **four sample-level observations**, not one
independent observation per cell barcode. The project record calls them four
independent unsorted-PBMC samples, but it does not use a verified common donor
or known equal starting composition to establish that differences are purely
technical. PBMC composition can genuinely differ between samples, so the four
predicted composition vectors are not expected to be identical.

This limits, but does not eliminate, the value of the cohort. Its useful
questions are:

1. **Does the method operate reliably on several real PBMC samples?** Both
   model arms should converge, return valid nine-type composition vectors, and
   avoid severe reconstruction failures in all four samples.
2. **How much does fragment-length information change a prediction for the
   same sample?** Within each replicate, the length-aware and count-only models
   receive the same aggregated sample. Their paired difference can therefore
   be measured without assuming that Rep1 and Rep2 have the same biology.
3. **Is the model-arm difference consistent across samples?** For example, one
   can ask whether adding fragment length repeatedly increases or decreases the
   same predicted cell type, or whether its effect is erratic and driven by one
   sample.
4. **Is a reported pattern dependent on one PBMC sample?** Repeating the
   analysis across four samples reveals whether a qualitative conclusion is
   broad or is caused by a single outlier.

What the cohort cannot answer is whether either model recovered the correct
composition or which model is more accurate. Raw variation among Rep1--Rep4 is
also not a clean estimate of technical reproducibility: high variation could
be real donor/sample biology, while low variation could occur even if every
prediction is similarly biased. A proper technical-reproducibility experiment
would require matched aliquots of the same input PBMC material processed
repeatedly. An accuracy experiment would additionally require independent
composition measurements or known source-cell labels.

### Fresh/frozen preparation samples

The preparation cohort contains three differently prepared PBMC samples:

| Preparation label | GEO accession | Author sample title | Meaning used in this evaluation |
|---|---|---|---|
| Fresh | `GSM3722040` | `Fresh_pbmc_5k` | PBMC material carrying the author's fresh-preparation label |
| Frozen, sorted | `GSM3722041` | `Frozen_sorted_pbmc_5k` | PBMC material carrying both frozen and sorted labels |
| Frozen, unsorted | `GSM3722042` | `Frozen_unsorted_pbmc_5k` | Frozen PBMC material without the sorted label |

These labels identify the preparation conditions recorded by the source. The
analysis should not infer an undocumented processing order, sorting method, or
exactly matched starting composition from the short sample titles alone. As in
the replicate cohort, all retained barcodes within each sample are aggregated
into one observation. The model therefore produces one nine-type composition
vector for `fresh`, one for `frozen_sorted`, and one for `frozen_unsorted`.

The three predictions can be compared to ask whether the inferred PBMC profile
changes with the recorded preparation condition. For example:

- fresh versus frozen-unsorted is informative about the combined difference
  associated with those two recorded preparation labels;
- frozen-sorted versus frozen-unsorted is informative about the difference
  associated with the recorded sorting status among the frozen samples; and
- length-aware versus count-only predictions can be paired within each of the
  three samples.

These are **preparation-sensitivity comparisons**, not controlled causal
estimates. A prediction shift may be caused by the preparation, differential
cell recovery, a real composition difference between input samples, reference
mismatch, or model instability. Without known starting or recovered cell-type
proportions, those explanations cannot be separated.

### What can and cannot be concluded

The two cohorts can be used to summarize replicate-to-replicate dispersion,
preparation-associated prediction shifts, differences between length-aware and
count-only predictions, convergence, and reconstruction diagnostics. They do
not support RMSE or JSD against composition truth, because neither exact
source-cell membership nor a nominal mixing vector is available. Stable
predictions are evidence of robustness, but are not proof of accuracy;
unstable predictions flag sensitivity, but do not reveal which prediction is
correct.

## Case 4: GSE205055 and GSE263333 real spatial tissue

These datasets contain real spatial ATAC pixels, not mixtures constructed from
known labeled source cells. Independent labeled single-cell references are
available and are used to learn signatures:

- GSE216371 for E13 mouse embryo;
- GSE246791 for adult mouse brain; and
- GSE244618 for adult human hippocampus.

Those reference annotations answer, "What does each cell type look like in the
reference?" They do not answer, "Which cells occupy this spatial pixel?" The
adult healthy mouse-brain reference also has a recorded disease/stage mismatch
when applied to the five-month EAE brain section.

Therefore, GSE205055 and GSE263333 do not receive truth-based RMSE or JSD.
Their comparisons instead include the following.

### Length-aware versus count-only map agreement

For each predicted cell type, the two methods' values across spatial spots are
compared using Pearson correlation, Spearman correlation, and mean absolute
difference. This measures how much fragment length changes a map; it does not
identify which map is correct.

### Orthogonal concordance

Predicted cell-type abundance across spots is compared with independently
measured RNA marker scores, protein markers, histone signals, or
reference-derived accessibility markers. The analysis asks whether the
length-aware or count-only map has stronger spatial concordance with each
proxy.

These signals are not composition truth. Marker expression and chromatin state
can change without a corresponding change in cell abundance, markers may be
shared across cell types, and the orthogonal measurements may themselves be
spot-level mixtures.

### Spatial and replicate behavior

Additional endpoints include spatial continuity, boundary similarity,
replicate-section consistency, and reconstruction-mismatch warnings. They can
detect implausibility or instability, but they cannot establish absolute
deconvolution accuracy.

## Why not replace every evaluation with pseudo-spots?

Pseudo-spots should be generated wherever appropriate labeled, cell-resolved
data exist, but they answer a narrower question than physical or real-spatial
samples.

| Evidence design | Main question | Important limitation |
|---|---|---|
| Exact pseudo-spots | Can the model recover a known mixture under controlled reference/test separation? | Cells are combined computationally after they were assayed |
| Physical dilutions | Does the method tolerate real sorting, mixing, preparation, capture, and recovery? | The intended input ratio may differ from recovered composition |
| Real spatial tissue | Are maps stable and biologically concordant under real anatomy and reference mismatch? | Exact spot composition is unavailable |

A GSE129785 pseudo-spot benchmark could be created by splitting cells from its
nine sorted populations into reference and held-out pools. It would be a useful
additional controlled experiment, but it would not replace the physical
dilutions because both pseudo-spots and signatures would come from closely
related sorted-cell preparations.

Similarly, pseudo-spots could be made from the external references used for
GSE205055 and GSE263333. Such an experiment would test deconvolution on those
reference sources; it would not create truth for the existing real-spatial
pixels or test their tissue, disease, capture, and protocol mismatch.

The strongest evaluation therefore retains all evidence classes while keeping
their interpretations separate:

1. Exact pseudo-spots for controlled quantitative accuracy.
2. Physical mixtures for experimental robustness and rare-component recovery.
3. Real spatial tissue for map stability and orthogonal biological evidence.

## Relevant project records

- [PBMC benchmark protocol](../benchmark_protocol.md)
- [PBMC stress-v2 amendment](../pbmc_stress_v2_amendment.md)
- [Full-evaluation protocol](../full_evaluation_protocol.md)
- [Full-evaluation report](../full_evaluation_report.md)
- [GSE129785 preprocessing record](../datasets/gse129785_preprocessing.md)
- [GSE205055/GSE263333 preprocessing record](../datasets/gse205055_gse263333_preprocessing.md)
- [GSE129785 source configuration](../../../configs/data_sources/shapemix_gse129785.yaml)
- [PBMC stress-v2 dataset configuration](../../../configs/datasets/shapemix_pbmc_stress_v2.yaml)
- [GSE205055 dataset configuration](../../../configs/datasets/shapemix_gse205055_real_spatial_v1.yaml)
- [GSE263333 dataset configuration](../../../configs/datasets/shapemix_gse263333_real_spatial_v1.yaml)
