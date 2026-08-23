# ShapeMix-ATAC

ShapeMix-ATAC is a reference-guided method for deconvolving spatial ATAC mixtures. Its narrow primary question is:

> Does the parent-fragment length composition of ATAC cut sites improve cell-type proportion estimates beyond the same cut sites collapsed to peak totals?

The MVP is a non-spatial, fixed-reference-signature, maximum-a-posteriori model. It is intentionally smaller than the broader models discussed in the early tutorials so that the contribution of fragment shape can be tested in an exactly matched ablation.

## Canonical documents

Use these files in this order:

1. [Model specification](model_specification.md) — normative data semantics and statistical model.
2. [Benchmark protocol](benchmark_protocol.md) — normative data split, peak selection, seeds, metrics, and comparison rules.
3. [Current Bayesian tutorial source](tutorials/ShapeMix_ATAC_Bayesian_Model_Tutorial.tex) and [generated PDF](tutorials/ShapeMix_ATAC_Bayesian_Model_Tutorial.pdf) — canonical explanatory derivation of the MVP and labeled extensions.
4. [Implementation plan](implementation_plan.md) — file-level engineering roadmap and acceptance gates.
5. [Results summary](results_summary.md) — concise primary outcome, validation, controls, and interpretation limit.
6. [Step 6 execution report](step6_results.md) — detailed non-normative result, diagnostics, gates, provenance, and limitations.
7. [Research proposal](../research_class/ShapeMix_ATAC_proposal_draft.md) — current motivation, research question, and broader study design.

If two documents disagree about the executable MVP, `model_specification.md` governs the model and `benchmark_protocol.md` governs the experiment. The implementation plan governs sequencing and repository integration. The proposal governs the scientific motivation and was reconciled with the canonical count unit and likelihood on 2026-08-22; the two normative specifications still take precedence over any future proposal drift.

## Frozen MVP contract

- Input rows follow the Cell Ranger ARC 2.0 five-column fragments schema: chromosome, start, end, barcode, and `readSupport`. There is no strand field.
- Each row represents a deduplicated fragment. The primary pipeline ignores `readSupport`, emits its two Tn5 cut sites, and tags both cut sites with the parent fragment's length.
- Parent-fragment length is `end - start` and uses bins `[0,100)`, `[100,250)`, and `[250,∞)` bp.
- Each cut site is assigned independently to the unique non-overlapping Cell Ranger peak that contains it. Unassigned cut sites are not counted.
- Three sparse `AnnData` layers hold cut-site counts by length bin, and `.X` is their exact elementwise sum.
- Reference and held-out mixture cells are disjoint. Peaks and fixed signatures are learned from the reference split only.
- Total peak counts use a negative-binomial likelihood. Conditional bin composition uses a multinomial likelihood.
- The peak-only arm uses the same total-count likelihood, signatures, abundance prior, initialization, optimizer, seeds, and compute budget; only the conditional shape term is disabled.
- No observed spot-depth offset is used. Inferred abundance is an effective reference-cell-equivalent, depth-scaled quantity; normalized abundance is the reported cell-type proportion.
- The primary benchmark has no learned background, spatial smoothing, or uncertainty interval.
- Metrics use the dataset-declared ordered cell-type universe. The primary Jensen–Shannon metric is `jsd_v2`, the mean base-2 divergence, not the repository's historical unsquared distance.

## Document status

The following binary documents are historical, non-normative context. Retain them for provenance, but do not implement directly from them:

- `ShapeMix_high_level.docx`
- `tutorials/ShapeMix_ATAC_Tutorial.docx`
- `tutorials/ShapeMix_ATAC_Tutorial.pdf`
- `tutorials/ShapeMix_ATAC_Bayesian_Model_Tutorial_deprecated.docx`

`tutorials/ShapeMix_ATAC_Bayesian_Model_Tutorial.tex` is the editable source of the detailed Bayesian tutorial and has been reconciled with the frozen MVP; its PDF is the generated explanatory artifact rebuilt from that source on 2026-08-22. The tutorial labels independent-bin likelihoods, latent signatures, a depth offset, background, and a spatial prior as future extensions rather than executable MVP requirements. When the tutorial and normative specification differ, follow the normative specification. Never edit a generated PDF manually.

## Scope boundaries

The MVP does not claim donor-level or cross-dataset generalization: the first PBMC data source contains one donor. Its repeated splits quantify conditional resampling variability within that donor. The following are gated extensions, not reasons to change the primary model in place:

- Dirichlet-multinomial shape overdispersion;
- variational uncertainty;
- learned background or cross-protocol scaling;
- peak-specific dispersion;
- cut-position, TSS, or motif-footprint features;
- spatial smoothing;
- external-donor and real-spatial validation.

Any such extension requires a new configuration, explicit negative controls, and a versioned data schema when the observation axis changes.

## Implementation status

Steps 0 through 6 of the implementation plan were completed in this workspace on 2026-08-23. The scientific and benchmark contracts are frozen; the raw ARC 2.0 fragments/index are checksum-pinned and tabix-validated; the reusable sparse fragment-shape counter and optional AnnData contract are implemented; the held-out PBMC benchmark has been generated; the fixed-signature MAP core is implemented and tested; both nested arms are integrated with the maintained runner; and the strict paired primary ablation has been executed. Official-matrix reconstruction froze `chromEnd` (`right_cut_offset: 0`) as the right Tn5 coordinate.

Step 3 produced five 6,644-reference/2,856-held-out primary splits, one three-type development split, one smoke dataset, and 20 registered primary pseudo-spatial datasets covering the frozen outer seeds, inner seeds, and equal/observed conditions. Reference-only signal audits support retaining the three-bin per-peak representation with planned hierarchical smoothing, while documenting moderate depth-related composition as a technical caveat. These artifacts quantify conditional resampling within one donor; they do not establish cross-donor or cross-dataset generalization.

Step 4 added strict typed configuration; deterministic, reference-only `A`, `omega`, and cross-fitted `phi_ref`; stable negative-binomial and conditional-multinomial likelihoods; chunked diagnostics; and one streamed PyTorch MAP path for both the shape-aware and peak-only arms. The focused synthetic suite exercises all ten model invariants, and both frozen arms converge on the development smoke data through ordered sparse layers. That smoke run validates execution only and is not a scientific comparison.

Step 5 added one lazy `shapemix` adapter for both `use_shape` settings, named `method_runs` with complete config provenance, preflight collision/error checks, standardized outputs, and five compact native diagnostics. Both frozen arms completed through the maintained runner on the 32-spot development dataset and reproduced exactly in a second CPU output root, including proportions, abundances, fitted objectives, and stable native artifacts. Those runs validate integration and determinism only—the historical smoke metrics are not a scientific ShapeMix result.

Step 6 added strict `rmse_v1`/`jsd_v2` evaluation, exact declared-universe and spot validation, nested paired summaries, rare-cell and reconstruction diagnostics, output hash manifests, full provenance, resource reporting, deterministic sharding, and fail-closed resume. All 40 primary fits completed and converged with no failed or unavailable pairs. The effect was shape-aware minus count-only, where negative favors ShapeMix. The directional-support rule was false in both conditions: observed-abundance RMSE was mixed and near zero, while JSD worsened in all five outer splits; both RMSE and JSD worsened in all five equal-cell-type outer splits. See the [Step 6 execution report](step6_results.md) for exact effects and intervals.

Development degeneracy controls passed their exact algebraic checks, and the deterministic permutation diagnostic was frozen. Because the primary direction was not positive, the protocol did not trigger the expensive primary negative-control campaign. NNLS completed 20/20 same-data baseline runs and was worse than count-only on both endpoints in every outer split. Cell2location, RCTD, SpatialDWLS, and all stress/sensitivity datasets remain explicit unactivated gates, not missing result rows.

One pre-launch shape-aware fit on `pbmc_granulocyte_sorted_10k_shapemix_equal_celltype_split_1103_mix_101` was used as a blinded resource pilot: only dimensions, resource demand, convergence, and output shape were inspected. Its retained artifact records 1,268.645 seconds of fit time but not the thread environment. It remained under `/private/tmp`, outside `results/primary` and every strict summary. The completed benchmark still contains one donor only; its five outer splits quantify conditional resampling variability, not biological replication or population generalization. The final repository suite passes with 230 tests and one existing skip.

Because `data/` and primary results are ignored, another clone must fetch the raw files and regenerate the derived artifacts using the [recreation guide](<../recreate_data_directory (important).md>) and [tracked source manifest](../../configs/data_sources/pbmc_granulocyte_sorted_10k_cellranger_arc_2.0.0.yaml), then rerun the frozen experiment and summary configurations.
