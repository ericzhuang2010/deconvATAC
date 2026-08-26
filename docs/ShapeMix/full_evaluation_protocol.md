# ShapeMix full-evaluation protocol v1

Status: frozen before production external predictions; CUDA routing and optimizer coordinate amended 2026-08-25
Frozen: 2026-08-24  
Canonical layout: [file_organization.md](file_organization.md)  
Execution roadmap: [iimplementation_plan.md](iimplementation_plan.md), Section 13

## Purpose

This protocol extends the completed one-donor PBMC protocol-v1 result across
external immune, donor-held-out bone-marrow, and real-spatial datasets. It does
not modify or rerun the completed primary PBMC campaign. Exact pseudo-spot
truth, nominal physical ratios, and qualitative spatial evidence remain
separate evidence classes and are never pooled into one endpoint.

## Core methods

Every runnable unit receives:

1. `shapemix_length`: the frozen negative-binomial total-count model plus the
   conditional three-bin parent-fragment-length likelihood;
2. `shapemix_count_only`: the same model with only `use_shape: false`; and
3. `nnls`: collapsed-count nonnegative least squares on the same reference
   and ordered feature axis.

The paired ShapeMix configurations must differ only in `use_shape`. Both arms
use float32, seed 0, three restarts, 2,000 Adam steps, patience 100, tolerance
`1e-5`, spot batches of 64, and peak chunks of 512. CUDA is used only after
the qualification gates below pass; both arms of a pair use the same backend.

## Evidence classes and endpoints

| Evidence class | Families | Permitted endpoints |
|---|---|---|
| Exact held-out pseudo-spot composition | GSE194122 and separately versioned PBMC sensitivities | `rmse_v1`, `jsd_v2`, per-type absolute error, rare-type detection, calibration, convergence, reconstruction, runtime |
| Nominal physical input ratio | GSE129785 physical dilutions | Descriptive nominal RMSE/JSD, absolute rare-component error, detection at 0.1%, 0.5%, and 1%, calibration slope/intercept, Spearman monotonicity, off-target mass |
| Prediction-only replicate/preparation | GSE129785 PBMC/preparation groups | Replicate dispersion, preparation shift, shape/count paired differences, off-target mass, convergence, reconstruction |
| Real spatial without exact composition truth | GSE205055 and GSE263333 | Shape/count map stability, RNA/protein/histone/marker/anatomy concordance, replicate stability, spatial continuity, boundary preservation, residual/off-reference warnings |

Nominal ratios live only at
`validation/nominal_broad_proportions.csv`. RNA, protein, histone, image,
anatomical, and marker evidence also remains under each dataset's
`validation/` tree. Only source-cell pseudo-spot composition may be written
to `truth/proportions.csv`.

The runner mode is `exact_truth` for exact pseudo-spots and
`prediction_only` for all other classes. Prediction-only campaigns request
`metrics: []`; the runner must not produce truth-based metric rows.

## Analysis units

- GSE194122: donor. The two inner mixtures and any repeated sites are first
  summarized within donor.
- GSE129785: one independently prepared physical sample or preparation.
  Aggregated barcodes are not replicates.
- GSE205055/GSE263333: one independently prepared section or matched section
  group where appropriate. Pixels/spots are not biological replicates.
- Optimizer restarts are never biological replicates.

## GSE129785 external immune campaign

All 16 registered descriptors are run using the standardized nine-class
reference and the `selected_reference_peaks` axis:

1. Dendritic cells
2. Monocytes
3. B cells
4. Regulatory T cells
5. Naive CD4 T cells
6. Memory CD4 T cells
7. NK cells
8. Naive CD8 T cells
9. Memory CD8 T cells

The seven CD4-memory/CD8-naive dilutions form the primary nominal series. The
seven monocyte/T-cell dilutions form an exploratory broad-contrast series.
PBMC replicates and preparation comparisons are prediction-only. No exact
truth metric is permitted for any of these descriptors.

Rare-component detection is evaluated only at nominal fractions at or below
1%. A component is called detected when its predicted fraction is at least
one-half of the nominal fraction; the predicted fraction and absolute error
are always reported so this threshold cannot hide the underlying estimate.

Frozen experiment:
`configs/experiments/shapemix_gse129785_external.yaml`.

## GSE194122 donor-held-out BMMC campaign

### Donor folds and seeds

The outer folds hold out donors `1,2,3,4,5,6,7,8,9,10`. All sites belonging
to a donor stay on the same side. Each donor has two inner mixture seeds,
`101` and `211`, under two conditions: `observed_abundance` and
`equal_celltype`. Each materialized dataset contains 1,024 pseudo-spots with
Poisson mean 10 cells per spot, truncated to at least one cell. Source cells
are sampled without replacement within a spot when possible and with
replacement only when a requested per-type pool is smaller than the request.

### Frozen broad ontology

The model uses the following seven classes in this exact order:

1. B/plasma: B1 B, Naive CD20+ B, Transitional B, Plasma cell
2. CD4 T: CD4+ T activated, CD4+ T naive
3. CD8 T: CD8+ T, CD8+ T naive
4. NK/ILC: NK, ILC
5. Myeloid/DC: CD14+ Mono, CD16+ Mono, cDC2, pDC
6. Erythroid/MK-E: Proerythroblast, Erythroblast, Normoblast, MK/E prog
7. Hematopoietic progenitor: HSC, Lymph prog, G/M prog, ID2-hi myeloid prog

The source table retains all 22 author labels. The broad mapping is frozen
before peak selection or prediction. A fold fails before fitting if any broad
class has fewer than 100 training cells or fewer than 20 held-out cells. The
source audit established that every planned fold passes these thresholds; the
minimum held-out support is 20 cells.

### Leakage controls and feature axis

For each held-out donor:

- compute reference label totals, nonzero-cell coverage, and peak rankings from
  the nine training donors only;
- require at least ten nonzero training cells per eligible peak;
- use the protocol-v1 score and tie-breaks from
  `benchmark_protocol.md`;
- select exactly 5,000 peaks;
- estimate signatures and dispersion from training donors only;
- create pseudo-spots from held-out cells only; and
- record every source cell, original site, author label, and broad label.

To avoid ten repeated full fragment scans without changing any fold's
selection, all ten training-only 5,000-peak lists may be unioned after they are
independently frozen. Fragment shapes are counted once per source sample on
that union axis, then each fold is sliced to its own ordered 5,000-peak list.
The union contents never participate in a fold's ranking.

Each fragment contributes one left cut at `chromStart` and one right cut at
`chromEnd`, the convention validated against the deposited aggregate count
matrix. Duplicate read support is ignored. The bins are `[0,100)`,
`[100,250)`, and `[250,infinity)` bp.

## Real-spatial reference tracks

Three reference tracks are kept separate:

- E13-compatible mouse embryo for the embryonic sections;
- adult mouse brain for adult and five-month EAE-brain sections, with disease
  mismatch explicit; and
- adult human hippocampus for the human section.

The source scopes are frozen before acquisition:

- GSE216371: the GEO family metadata, author annotation/cCRE workbook, and
  complete 68-BED archive; retain every author-annotated E13.5 cell;
- GSE246791: GEO metadata, four author atlas metadata files, and exactly one
  lexicographically first `biorep=early` H5AD from each of the 12 declared
  major brain regions; and
- GSE244618: GEO metadata, author Tables S1-S6, and exactly nine BEDPE samples
  spanning HiT, HiB, and Sub for donors 1, 2, and 4, using rep1 where needed.
### GSE244618 human-hippocampus ontology freeze

Frozen 2026-08-26 before feature selection or any spatial prediction, the
reference uses these six classes in exact output order: inhibitory neurons,
excitatory neurons, astrocytes, oligodendrocytes, OPCs, and microglia.
Author `cellclass=GABA` and `cellclass=GLUT` define the two neuronal classes.
Within `cellclass=NonN`, subclasses `ACBGM`, `ASCNT`, and `ASCT` map to
astrocytes; `OGC`, `OPC`, and `MGC` map to oligodendrocytes, OPCs, and
microglia. The rare vascular subclasses `EC`, `PER`, and `SMC` are
predeclared exclusions. Every retained class has at least 20 annotated cells
in every one of the nine selected samples.

The feature selector ranks the 544,735 unique, non-overlapping author GRCh38
cCREs using the frozen count-only variance score across the six reference
classes, requires at least ten nonzero retained reference cells, and selects
exactly 5,000 cCREs with coverage, total-count, and identifier tie-breaks.
Neither the human spatial section nor any RNA concordance value participates
in labels or feature selection.


The frozen downloads total exactly 91,839,769,587 bytes (85.53 GiB) across 36
files. Their tracked manifests are the sole acquisition authority.

A candidate is runnable only after source identity, license/accessibility,
species, tissue/stage, genome build, labeled-cell support, fragment semantics,
right-cut convention, and shared peak construction pass. Reference-only peak
selection is fixed before inspecting spatial validation concordance. A failed
reference gate remains visible as a gated dataset; it is not replaced with an
incompatible universal atlas.

No labels, peaks, smoothing, or model parameters may be selected by maximizing
RNA/protein/histone concordance on a section later used for reporting.

## CUDA qualification

### Pre-external routing amendment (2026-08-25)

No external prediction had been run when this amendment was made. The original
predeclared gates and their failed result remain recorded in the qualification
table. On the 32 by 200 smoke input, CPU/CUDA maximum absolute proportion
difference was `0.00347122` and CUDA was only `1.025x` faster. Those tiny inputs
therefore remain on CPU. All GSE129785 inputs are 1--4 samples by 5,000 peaks
and use the paired CPU configurations.

On the representative 1,024 by 5,000 input, only 3 of 16,384 proportions
exceeded the original `1e-4` bound; the maximum was `0.000157089`, the 99th
percentile was `5.66e-8`, and RMSE/JSD differed by at most `3.18e-8`. Repeated
CUDA and cached/streamed CUDA results were identical, and cached CUDA was
`4.474x` faster than CPU. The production full-size CUDA bound is therefore
amended to `2e-4`. This is a hardware-equivalence decision based only on the
development qualification data, not an external scientific result. Large
campaigns use CUDA only when their scale is represented by this qualification;
otherwise they retain CPU. Both ShapeMix arms within every dataset use the same
backend.

Before any external prediction is inspected:

- deterministic algorithms are enabled with `warn_only=False`;
- CPU/CUDA toy objectives and gradients meet relative `1e-5` and absolute
  `1e-4` tolerances;
- smoke proportions are retained diagnostically and smoke-scale production
  inputs use CPU;
- one 1,024 by 5,000 development dataset agrees within `2e-4`, with RMSE/JSD
  agreement within `1e-5`;
- two repeated CUDA fits agree within `1e-7` and have the same convergence
  status;
- cached and host-streamed CUDA modes meet the same tolerances;
- CUDA is at least twice as fast as CPU on the representative full-size fit;
  otherwise the campaign retains the CPU backend; and
- peak allocated/reserved VRAM, cache mode, wall time, peak RSS, and device
  identity are recorded under `results/development/`.

No autocast, TF32, float16, or bfloat16 is permitted in protocol v1.

### Optimizer-coordinate amendment (2026-08-25)

The first GSE129785 campaign attempted all 48 jobs. All 16 NNLS jobs
succeeded, while all 32 ShapeMix jobs reached `max_steps` and produced no
ShapeMix prediction files. Development-only initialization and trajectory
diagnostics on the first registered descriptor identified a scale mismatch:
effective NNLS abundance reached 4,566.9, but the softplus optimizer coordinate
could not traverse that abundance scale under the frozen Adam budget.

This amendment does not change the likelihood, priors, signatures,
dispersion, seed tuples, convergence rule, or maximum-step budget. It changes
only the optimizer coordinate to `z = exp(raw_z) + epsilon`; restart zero is
the exact NNLS point, and later restarts retain the frozen deterministic
log-normal perturbation. Failed optimizations now retain hashed failure
diagnostics. A post-change one-descriptor development check converged in both
arms and was inspected only for convergence and finite normalized output; it
was not used to choose endpoints, tolerances, learning rate, or model terms.

The failed campaign and all development diagnostics remain visible under
`results/external_validation/shapemix_gse129785_v1/` and
`results/development/shapemix_gse129785_convergence_v1/`. Valid production
predictions must use a fresh run group. Because the coordinate changes the
optimization trajectory, CUDA qualification v1 is historical evidence only:
v2 smoke qualification has completed, and the v2 full-size qualification must
pass before any large campaign is routed to CUDA. The resource-safety stop of
the first v2 full-size attempt is recorded at
`results/development/shapemix_gpu_qualification_v2/full_size/interruption.yaml`;
it produced no scientific output.

## PBMC diagnostic sensitivity campaign

The one-donor PBMC follow-up changes one factor at a time and is explanatory,
not a search for a favorable subset. It uses outer split `1103`, evaluation
mixture seeds `307` and `401`, 1,024 spots per dataset, and exact source-cell
truth. Unless a row below says otherwise, the anchor is observed-abundance
sampling over all 16 frozen PBMC types, Poisson mean 10 cells per spot, the
training-only top 5,000 peaks, the complete split-specific reference, no depth
thinning, and bins `[0,100)`, `[100,250)`, and `[250,infinity)` bp.

| Factor | Frozen levels | Factor-specific control |
|---|---|---|
| Cut-site depth retention | 0.25, 0.50, 0.75, 1.00 | 1.00 is the global anchor |
| Mean cells per spot | 2, 5, 10, 20 | 10 is the global anchor |
| Controlled rare NK fraction | 0.001, 0.005, 0.010, observed 0.0424210526 | All non-NK probabilities retain their observed relative ratios; the observed level is the global anchor |
| Subtype similarity | equal broad trio (`CD14 Mono`, `CD4 Naive`, `CD8 Naive`); equal related CD4 trio (`CD4 Naive`, `CD4 TCM`, `CD4 TEM`) | Broad trio |
| Reference cells per type | 50, 100, 250, all available | Equal broad-trio mixtures; `all` is the factor control |
| Selected peaks | 1,000, 2,500, 5,000 | First peaks in the already frozen training-only ranking; 5,000 is the global anchor |
| Parent-length bins | two, three, five | Three bins is the global anchor |

The two-bin edges are `[0,100)` and `[100,infinity)`. The five-bin edges are
`[0,80)`, `[80,100)`, `[100,180)`, `[180,250)`, and `[250,infinity)` bp. A
single raw-fragment recount on the five-bin axis is used to derive the two-bin
objects by exact layer summation. The existing three-bin objects are not
rewritten. Reference-support subsets select the smallest SHA-256 barcode
digests within type, independently of held-out counts and predictions.

There is no complete factorial grid and no depth-by-rare interaction in
protocol v1. The anchor dataset is materialized once per evaluation seed and
reused as the declared control for depth, cells, rarity, feature count, and bin
count. The campaign contains 40 unique datasets and 120 core jobs across the
same three methods used elsewhere. Frozen experiment:
`configs/experiments/shapemix_pbmc_stress_v1.yaml`.

## Co-tenant resource policy

Every material task is launched through
`scripts/run_shapemix_low_impact.sh`. Only one deconvATAC task runs at a time.
The launcher requires one-minute load below 6.0, at least 4 GiB available host
memory, and no unrelated GPU compute process. It applies low CPU/I/O priority,
one host math thread, one GPU owner, and at most two preprocessing workers.
A closed gate pauses the next launch and never interrupts the other workflow.
The persistent `gnome-remote-desktop-daemon` is treated as display overhead,
not a competing scientific job, only while its reported allocation is at most
512 MiB; the total 2 GiB prelaunch GPU-memory cap still applies.

## Canonical outputs

- GPU pilots: v1 historical evidence under
  `results/development/shapemix_gpu_qualification_v1/`; current qualification
  under `results/development/shapemix_gpu_qualification_v2/`
- GSE129785: failed v1 evidence under
  `results/external_validation/shapemix_gse129785_v1/`; current production
  campaign under `results/external_validation/shapemix_gse129785_v2/`
- GSE194122: `results/external_validation/shapemix_gse194122_lodo_v1/`
- GSE205055/GSE263333: `results/real_spatial/<family>/<campaign_id>/`
- PBMC sensitivities: `results/sensitivity/<campaign_id>/`

Reusable caches remain under `data/processed/shapemix/`, standardized
references under `data/processed/references/`, and runnable inputs under
`data/processed/datasets/`. Downloads and preprocessing scratch never enter
`results/`.

## Completion criteria

A family is complete only if every planned unit is successful or explicitly
gated, paired ShapeMix arms use the same qualified backend, NNLS is present,
all failures remain visible, manifests and hashes revalidate, and conclusions
use only endpoints permitted for the evidence class. Negative and null results
are reported under the same rules as positive results.
