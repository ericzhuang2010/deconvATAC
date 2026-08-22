# HSR 2026–27 Proposal — Andy Zhuang

## Project Title

**ShapeMix-ATAC: A Novel Fragment-Shape-Aware Bayesian Approach to Spatial ATAC-seq Deconvolution**

## General Problem and Research Question

Although most cells in an organism contain nearly the same DNA sequence, different cell types use different parts of that DNA. Promoters, enhancers, and transcription-factor binding sites must be accessible before they can help regulate genes. ATAC-seq measures this chromatin accessibility by using the Tn5 transposase to cut and tag exposed DNA [1]. Consequently, ATAC-seq can reveal the regulatory programs that help establish cell identity and state, including information that may not be apparent from RNA abundance alone.

Spatial ATAC-seq is important because the location of a regulatory program within a tissue can be as meaningful as the program itself. Development, immune responses, brain organization, and diseases such as cancer all depend on which cells occupy particular regions and how their regulatory DNA responds to the local environment. Ordinary single-cell ATAC-seq measures individual cells but usually loses their original tissue locations. Spatial ATAC-seq preserves those locations, making it possible to connect open chromatin, cell identity, and tissue structure [2]. Newer spatial multi-omics assays can also measure chromatin accessibility alongside RNA, proteins, or histone modifications, emphasizing the importance of chromatin accessibility as a distinct layer of spatial biology [9].

However, the resolution of a spatial assay does not always equal the resolution of a single cell. One spatial spot or pixel can contain several nuclei. For example, the original spatial-ATAC-seq study reported approximately 1–2 cells in a 10 µm pixel and approximately 25 cells in a 50 µm pixel in mouse embryo tissue [2]. The fragments observed in such a spot are therefore a mixture from multiple cell types.

Deconvolution is important because a mixed signal can hide the biology researchers want to study. A strong accessibility signal from a rare cell type may be diluted by more common neighboring cells, while a change in cell-type composition may be mistaken for a regulatory change within one cell type. Without deconvolution, it can be difficult to determine whether a peak is accessible throughout a tissue region or is highly accessible only in one population of cells. Estimating the cell-type proportions in each spot can reveal rare or spatially restricted populations, improve maps of tissue organization, and make downstream comparisons of healthy and diseased regions more interpretable. In this way, computational deconvolution can extract finer biological information from existing spatial data without requiring a new experiment at every location.

ATAC fragments also contain more information than a total count at each open-chromatin peak. Their lengths and cut positions can reflect nucleosome organization and transcription-factor occupancy [1]. This creates an opportunity for an ATAC-specific deconvolution method to use the physical patterns of fragments rather than treating ATAC data exactly like RNA count data.

Published work has tested spatial ATAC deconvolution mainly by borrowing algorithms developed for spatial RNA data and substituting ATAC peaks for genes. In other words, an RNA method normally receives a spot × gene-expression matrix, while its adapted ATAC version receives a spot × peak-accessibility matrix. A recent benchmark followed this strategy with Cell2location, RCTD, Tangram, SpatialDWLS, and DestVI and showed that several of these RNA-based methods can be transferred successfully to spatial ATAC data [3]. These results establish useful methods and performance baselines for studying the problem.

Simply transferring RNA methods may leave useful information unused because RNA molecules and ATAC fragments are not equivalent measurements. RNA deconvolution models are built around gene-expression counts. ATAC data contain peak counts, but each fragment also has a genomic start position, end position, and length related to the physical organization of chromatin. Existing transferred methods generally collapse all fragments within a peak into one count, removing this extra structure. The benchmark also found that ATAC-based results remained slightly less accurate than RNA-based results, especially for rare cell types [3]. This suggests that borrowing RNA methods is feasible but may not be the final solution. An ATAC-native model could improve deconvolution by using biological information that has no direct equivalent in transcriptome data, especially because quantitative fragment counts preserve useful regulatory information in single-cell ATAC-seq [4].

The existing methods also provide useful comparisons for model design. Cell2location and RCTD use probabilistic count-mixture models, while SpatialDWLS estimates cell-type proportions with a weighted least-squares strategy [5, 6, 8]. ShapeMix-ATAC will retain the useful idea of a reference-guided mixture while changing the data representation to include ATAC fragment length.

**General research question:** Does the parent-fragment length composition of observed, deduplicated Tn5 cut sites improve cell-type deconvolution beyond the total cut-site count in each peak?

## Specific Goal and Hypothesis

**Specific goal:** I will develop and evaluate **ShapeMix-ATAC**, a reference-based computational method that estimates the proportion of each cell type in a mixed spatial ATAC spot using both (1) the total number of deduplicated Tn5 cut sites in each chromatin peak and (2) the distribution of those cut sites among bins defined by their parent fragments' lengths.

**Hypothesis:** If parent-fragment length distributions contain reproducible cell-type information, then ShapeMix-ATAC will estimate known cell-type proportions more accurately than an otherwise matched peak-count-only model. I expect the largest improvement for rare cell types and closely related immune-cell subtypes that have similar total accessibility profiles.

The central test will be an exactly nested ablation experiment. Both arms use the same total-count likelihood, fixed reference signature, abundance prior, parameters, initialization, optimizer, stopping rule, random seed, and compute budget. ShapeMix adds only the conditional likelihood for allocating each peak total among parent-fragment length bins. This design isolates the contribution of parent-fragment length composition.

## Proposed Methods

### 1. Data and Experimental Design

I will begin with a public, annotated single-cell multiome dataset containing both RNA and ATAC measurements from human peripheral blood mononuclear cells (PBMCs). RNA measurements can support cell-type labels, while raw ATAC fragment files provide chromosome, start position, end position, and cell barcode. This dataset contains related T-cell, B-cell, monocyte, dendritic-cell, and natural-killer-cell populations, making it useful for testing fine-subtype and rare-cell recovery.

To avoid evaluating the method on the same cells used to build its reference, I will use deterministic, cell-type-stratified 70/30 reference/test splits. The available PBMC input contains one donor, so donor-level separation and donor-level generalization claims are not possible in the primary experiment. Each outer split seed defines a new 70/30 partition; mixture-generation seeds are nested within that split and never change reference construction. The reference pool will be used only to select peaks and estimate fingerprints and dispersion. Cells from the held-out test pool will be combined into simulated spatial spots, so the exact cell-type proportions in every spot are known. Simulated mixtures are an established strategy for evaluating spatial ATAC deconvolution when experimental ground truth is unavailable [3].

### 2. Building Cell-Type Fragment Fingerprints

For the primary experiment, I will select 5,000 peaks using a deterministic accessibility-variance ranking computed from the reference pool only; held-out cells, shape-bin counts, and simulated mixtures will not influence selection. Experiments with 10,000 and 20,000 peaks will be sensitivity analyses. Highly variable peaks are prioritized because they performed better than simply choosing the most accessible peaks in the deconvATAC benchmark [3].

Each row of the 10x fragments file represents one deduplicated fragment. I will calculate its parent-fragment length as `chromEnd - chromStart` from the adjusted, zero-based half-open interval and assign both Tn5 cut sites from that row to the same one of three bins:

- Short: less than 100 base pairs
- Mononucleosome-like: 100–249 base pairs
- Long: 250 base pairs or longer

Thus, the exact intervals are `[0,100)`, `[100,250)`, and `[250,infinity)` base pairs. Each cut site is assigned to the unique selected, non-overlapping Cell Ranger peak that contains it. The exact right-cut coordinate convention will be frozen only after selected columns of the official peak matrix have been reconstructed. The fragments-file `readSupport` field will be ignored in the primary analysis because it reports read-pair support for an already deduplicated fragment; weighting by it would reintroduce PCR support. It may be studied later as an explicitly separate sensitivity analysis.

For each cell type, I will estimate a total cut-site signature and a conditional parent-length-bin distribution from reference cells only. Starting with three bins will limit sparsity and keep the first model feasible. The primary tensor therefore counts **deduplicated Tn5 cut sites annotated by their parent fragment's length**, not unique fragments. A unique-fragment analysis would require its own versioned dataset and a matched fragment-count baseline; it is not interchangeable with the primary input.

### 3. Creating Pseudo-Spatial Spots

I will construct simulated spatial spots by sampling and summing labeled test cells. The initial benchmark will contain 1,024 spots with an average of approximately 10 cells per spot. I will create both balanced simulations and realistic imbalanced simulations in which some cell types are rare. Additional experiments will vary sequencing depth, cells per spot, cell-type similarity, and rare-cell abundance. Each spot's sampled cells will be recorded to create a ground-truth proportion table.

### 4. Count-Only Baseline and ShapeMix-ATAC Model

Let `Y[s,p,b]` be the cut-site count for spot `s`, peak `p`, and parent-fragment length bin `b`, and let `N[s,p] = sum_b Y[s,p,b]` be the ordinary peak total. Let fixed reference estimates `A[c,p]` and `omega[c,p,b]` denote, respectively, the mean total cut-site rate per reference cell and the conditional length-bin probabilities for cell type `c`. The model infers positive effective abundances `z[s,c]` and defines:

```text
v[s,p,b]   = sum_c z[s,c] * A[c,p] * omega[c,p,b]
mu[s,p]    = sum_b v[s,p,b] = sum_c z[s,c] * A[c,p]
rho[s,p,b] = v[s,p,b] / mu[s,p]
```

The primary same-dataset simulation uses zero background and the factorized likelihood:

```text
N[s,p] ~ NegativeBinomial(
    mean=mu[s,p],
    dispersion=max(sum_c z[s,c], epsilon) * phi_ref
)
Y[s,p,1:B] | N[s,p] ~ Multinomial(N[s,p], rho[s,p,1:B])
```

Here `phi_ref` is one global inverse-dispersion estimated by a deterministic, cell-type-stratified two-fold cross-fit within the training-reference pool. Scaling it by `sum_c z[s,c]` approximates the aggregation of independent reference-depth cells. The count-only arm uses the negative-binomial term and the same prior on `z`; ShapeMix uses those exact terms plus the conditional multinomial. Independent negative-binomial likelihoods for the bins are not the primary model because their sum is generally not the same negative-binomial total and would prevent an exactly nested ablation. They may be evaluated later as a clearly labeled sensitivity model.

The primary model has no observed spot-depth offset: total depth is absorbed into `z`. Therefore `z` is an **effective reference-cell-equivalent, depth-scaled abundance**, not a calibrated nucleus count. The reported cell-type proportions are `pi[s,c] = z[s,c] / sum_c z[s,c]`. Background, cross-protocol scaling, learned signatures, spatial smoothing, and peak-specific dispersion are future extensions rather than components of the primary test.

I will fit the fixed-signature model by maximum a posteriori optimization, inferring only `z`, because this version is easier to identify, implement, and debug. Variational inference and credible intervals are future extensions after the MAP model is validated. The two primary arms will differ only by inclusion of the conditional shape likelihood.

### 5. Evaluation and Statistical Analysis

The primary dependent variable will be the error between the true and estimated cell-type proportions. The fixed primary metrics are RMSE over every spot-by-declared-cell-type entry and `jsd_v2`, the mean per-spot base-2 Jensen--Shannon **divergence**. Operationally, `jsd_v2 = mean(jensenshannon(truth, prediction, base=2)^2)`, with range `[0,1]`; it must not be confused with unsquared Jensen--Shannon distance. The dataset-declared ordered cell-type universe and exact spot set are fixed for every method. Missing declared prediction columns are filled with zero and penalized; extra cell types, non-finite values, negative values, all-zero rows, or mismatched spot sets fail evaluation rather than being dropped or silently normalized. Secondary outcomes will include rare-cell detection F1 score, correlation between true and estimated proportions, and runtime. Credible-interval coverage will be reported only if a later uncertainty model is implemented.

I will repeat the experiment across frozen outer split seeds and nested mixture seeds and report paired differences between models. Inner-mixture effects will first be averaged within each outer split; spots and inner mixtures will not be treated as independent biological replicates. Any resampling interval therefore quantifies conditional variability for this one donor, not donor-level uncertainty or biological generalization. I will also examine results separately for common cell types, rare cell types, and related immune subtypes. A successful result would be a consistent reduction in RMSE or `jsd_v2` for ShapeMix-ATAC without an unacceptable increase in false-positive rare-cell calls. If the two models perform similarly, that will also be informative because it would show that parent-fragment length composition does not add enough signal under the tested conditions.

The matched shape-versus-count ablation will remain the primary test. As a secondary benchmark, I will also compare ShapeMix-ATAC with established count-only implementations available in the project, including nonnegative least squares, Cell2location, RCTD, and SpatialDWLS [3, 5, 6, 8]. This will show whether any gain is meaningful relative to existing methods rather than only relative to a weak baseline.

### 6. Robustness and Optional Real-Data Validation

I will test whether any improvement remains after changing fragment-length cutoffs, peak number, spot depth, and cells per spot. I will inspect residuals and reconstructed fragment-length distributions to identify reference mismatch or technical bias. If the simulated benchmarks support the hypothesis, I will apply the method to a public spatial ATAC or spatial multi-omics dataset and compare its predicted spatial cell-type patterns with known tissue regions and marker peaks. Newer assays such as spatial-Mux-seq demonstrate that chromatin accessibility can be measured together with other molecular layers in the same tissue, which may provide useful cross-modality evidence for qualitative validation [9]. Because exact proportions are not available for most real tissues, this step will be treated as qualitative validation rather than the primary accuracy test.

Spatial smoothing will be considered only after the fragment-length model works independently at each spot. SONAR showed that adaptively borrowing information from similar neighboring spots can improve deconvolution while reducing bias across sharp tissue boundaries [7]. A later ShapeMix-ATAC version could use that principle, but the first experiment will exclude spatial smoothing so it does not confound the test of fragment shape.

## Variables and Controls

| Category | Definition |
| --- | --- |
| **Primary independent variable** | Likelihood: shared total peak-count model alone versus that exact model plus the conditional parent-fragment-length composition term |
| **Secondary independent variables** | Sequencing depth, number of cells per spot, peak set, fragment-length cutoffs, cell-type similarity, and rare-cell abundance |
| **Primary dependent variables** | RMSE and `jsd_v2` between estimated and true cell-type proportions |
| **Secondary dependent variables** | Rare-cell F1 score, proportion correlation, runtime, and uncertainty calibration/coverage |
| **Controlled factors** | Outer reference/test split, nested simulated spots, selected peaks, total-count likelihood, signatures, dispersion, abundance prior, initialization, optimizer/stopping settings, random seeds, compute budget, metric universe, and evaluation code |

## Novelty, Feasibility, and Limitations

### Novelty

Spatial ATAC deconvolution is a new and relatively undeveloped research area compared with spatial RNA deconvolution. Spatial transcriptomics has a much larger collection of established methods, whereas the main published spatial ATAC benchmark evaluated methods transferred from spatial RNA [3]. Those transferred methods demonstrated that deconvolution from ATAC peak counts is possible, but they were not originally designed around the physical properties of ATAC fragments.

An April 2026 survey of the area organizes the existing work into three groups [10]. First, spatial RNA methods—including Cell2location, RCTD, Tangram, SpatialDWLS, and DestVI—have been applied directly to spatial ATAC peak-count matrices [3, 10]. Second, other spatial RNA methods contribute useful ideas such as probabilistic count models, spatial priors, and graph-based smoothing. Third, bulk ATAC deconvolution and adjacent spatial epigenomics methods contribute chromatin-specific feature selection, reference construction, denoising, and accessibility reconstruction. These neighboring methods solve important pieces of the problem, but they do not by themselves provide a fragment-shape-aware model for estimating cell-type proportions in mixed spatial ATAC spots.

The 2025 deconvATAC benchmark stated that dedicated spatial chromatin-accessibility deconvolution methods were not available at the time of that study [3]. The April 2026 survey likewise did not identify a purpose-built, spot-level spatial ATAC deconvolution algorithm in the primary and official sources it reviewed [10]. Because failure to find a method is not proof that none exists, the novelty claim should remain qualified as **“to our knowledge”** and should be checked again before final submission or publication.

The proposed novelty is therefore not the use of a reference-guided mixture model, Bayesian inference, spatial regularization, or ATAC peak counts by themselves. To our knowledge, current spot-level spatial ATAC deconvolution methods do not directly model a spot × peak × fragment-shape tensor. ShapeMix-ATAC's specific contribution is to incorporate cell-type-specific fragment-length fingerprints into the deconvolution likelihood and test whether they add information beyond total peak accessibility. This claim is deliberately narrow and testable: the matched ablation will compare the same model with and without fragment-length bins rather than assume that a more complicated model must perform better.

### Feasibility

The project is feasible within seven months because it is computational, uses public data, and does not require new wet-lab experiments or reagents. The existing deconvATAC workspace already provides annotated PBMC multiome matrices, a pseudo-spot simulator, count-only deconvolution baselines, and evaluation metrics. The remaining implementation work will be processing raw fragments, constructing the shape-count tensor, implementing ShapeMix-ATAC, and conducting controlled benchmarks.

### Limitations and Risk Management

The main risks are sparse per-peak shape counts, technical differences between reference and spatial assays, the one-donor design, and increased computation. I will address these by starting with only three fragment-length bins, selecting well-covered peaks, using shrinkage or peak grouping if necessary, comparing matched and mismatched data conditions, and treating background, cross-protocol scaling, learned signatures, independent-bin likelihoods, motif footprints, spatial smoothing, and uncertainty inference as optional extensions rather than requirements for the first experiment. The deterministic 70/30 partitions prevent cell-level leakage but cannot establish cross-donor generalization.

## Seven-Month Work Plan

1. **Month 1:** Complete the literature review, obtain raw fragment files, finalize data splits, and reproduce count-only baselines.
2. **Month 2:** Implement and validate fragment-length binning and the spot × peak × bin tensor.
3. **Month 3:** Build reference fingerprints and implement the minimal ShapeMix-ATAC model.
4. **Month 4:** Debug the model on toy data and run the first controlled ablation experiment.
5. **Month 5:** Run repeated PBMC benchmarks for rare cells, related subtypes, depth, and spot density.
6. **Month 6:** Perform sensitivity analyses, uncertainty modeling if feasible, and optional real-data validation.
7. **Month 7:** Analyze results, prepare figures, document limitations, and write the final report and presentation.

## Core Five Full-Text PDF References

1. [Buenrostro, J. D. et al. (2013). *Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position*. Nature Methods, 10, 1213–1218.](https://research.stowers.org/cws/CompGenomics/Papers/Buenrostro_2013_NatMeth.pdf)
2. [Deng, Y. et al. (2022). *Spatial profiling of chromatin accessibility in mouse and human tissues*. Nature, 609, 375–383.](https://europepmc.org/articles/PMC9452302?pdf=render)
3. [Ouologuem, S. et al. (2025). *Spatial transcriptomics deconvolution methods generalize well to spatial chromatin accessibility data*. Bioinformatics, 41, i314–i322.](https://europepmc.org/articles/PMC12261446?pdf=render)
4. [Martens, L. D. et al. (2024). *Modeling fragment counts improves single-cell ATAC-seq analysis*. Nature Methods, 21, 28–31.](https://europepmc.org/articles/PMC10776385?pdf=render)
5. [Kleshchevnikov, V. et al. (2020; published in 2022). *Comprehensive mapping of tissue cell architecture via integrated single-cell and spatial transcriptomics* (Cell2location full-text preprint).](https://www.biorxiv.org/content/10.1101/2020.11.15.378125v1.full.pdf)

## Additional Full-Text PDF References

6. [Cable, D. M. et al. (2020; published in 2022). *Robust decomposition of cell type mixtures in spatial transcriptomics* (RCTD full-text preprint).](https://www.biorxiv.org/content/10.1101/2020.05.07.082750v1.full.pdf)
7. [Liu, Z. et al. (2023). *SONAR enables cell type deconvolution with spatially weighted Poisson-Gamma model for spatial transcriptomics*. Nature Communications, 14, 4727.](https://europepmc.org/articles/PMC10406862?pdf=render)
8. [Dong, R. & Yuan, G.-C. (2021). *SpatialDWLS: accurate deconvolution of spatial transcriptomic data*. Genome Biology, 22, 145.](https://europepmc.org/articles/PMC8108367?pdf=render)
9. [Guo, P. et al. (2025). *Multiplexed spatial mapping of chromatin features, transcriptome and proteins in tissues*. Nature Methods, 22, 520–529.](https://europepmc.org/articles/PMC11906265?pdf=render)
10. [*Computational Methods for Deconvolution of Spatial ATAC-seq Data* (internal literature survey, April 26, 2026; not peer reviewed).](<../related papers/Survey of computational Methods for Deconvolution of Spatial ATAC-seq Data.pdf>)
