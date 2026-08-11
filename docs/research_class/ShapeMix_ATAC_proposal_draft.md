# HSR 2026–27 Proposal — Andy Zhuang

## Project Title

**ShapeMix-ATAC: Testing Whether Fragment-Shape Information Improves Cell-Type Deconvolution of Spatial ATAC-seq Data**

## General Problem and Research Question

ATAC-seq measures open regions of DNA by using the Tn5 transposase to cut accessible chromatin. The resulting DNA fragments contain more information than a total count at each open-chromatin peak: their lengths and cut positions can reflect nucleosome organization and transcription-factor occupancy [1]. Spatial ATAC-seq preserves the locations of these measurements in tissue, but one spatial spot or pixel can contain several nuclei. For example, the original spatial-ATAC-seq study reported approximately 1–2 cells in a 10 µm pixel and approximately 25 cells in a 50 µm pixel in mouse embryo tissue [2]. Therefore, the signal measured at one spot can be a mixture of multiple cell types.

Current spatial deconvolution methods estimate the cell types in each spot mainly from a spot-by-feature count matrix. A recent benchmark showed that RNA-based methods such as Cell2location and RCTD can be applied to spatial ATAC peak counts, but ATAC results remained slightly less accurate than RNA results, especially for rare cell types [3]. These methods generally collapse all fragments within a peak into one count. This discards ATAC-specific physical information, even though quantitative fragment counts have been shown to preserve useful regulatory information in single-cell ATAC-seq [4].

The existing methods also provide useful comparisons for model design. Cell2location and RCTD use probabilistic count-mixture models, while SpatialDWLS estimates cell-type proportions with a weighted least-squares strategy [5, 6, 8]. ShapeMix-ATAC will retain the useful idea of a reference-guided mixture while changing the data representation to include ATAC fragment length.

**General research question:** Can cell-type-specific ATAC fragment-length patterns improve the accuracy of cell-type deconvolution compared with using total peak counts alone?

To our knowledge, current spot-level spatial ATAC deconvolution methods do not directly use a spot × peak × fragment-shape representation. This project will test that focused claim rather than assume that a more complicated model must perform better.

## Specific Goal and Hypothesis

**Specific goal:** I will develop and evaluate **ShapeMix-ATAC**, a reference-based computational method that estimates the proportion of each cell type in a mixed spatial ATAC spot using both (1) which chromatin peaks contain fragments and (2) the length distribution of those fragments within each peak.

**Hypothesis:** If fragment-length distributions contain reproducible cell-type information, then ShapeMix-ATAC will estimate known cell-type proportions more accurately than an otherwise matched peak-count-only model. I expect the largest improvement for rare cell types and closely related immune-cell subtypes that have similar total accessibility profiles.

The central test will be an ablation experiment: I will run the same model on the same simulated spots twice, once retaining fragment-length bins and once collapsing those bins into total peak counts. This design isolates the contribution of fragment-shape information.

## Proposed Methods

### 1. Data and Experimental Design

I will begin with a public, annotated single-cell multiome dataset containing both RNA and ATAC measurements from human peripheral blood mononuclear cells (PBMCs). RNA measurements can support cell-type labels, while raw ATAC fragment files provide chromosome, start position, end position, and cell barcode. This dataset contains related T-cell, B-cell, monocyte, dendritic-cell, and natural-killer-cell populations, making it useful for testing fine-subtype and rare-cell recovery.

To avoid evaluating the method on the same cells used to build its reference, I will divide cells into separate reference and test pools, stratified by cell type. A donor-level split will be used when donor labels allow it. The reference pool will be used only to learn cell-type fingerprints. Cells from the test pool will be combined into simulated spatial spots, so the exact cell-type proportions in every spot are known. Simulated mixtures are an established strategy for evaluating spatial ATAC deconvolution when experimental ground truth is unavailable [3].

### 2. Building Cell-Type Fragment Fingerprints

I will select approximately 5,000–20,000 reproducible, variable peaks shared by the reference and test data. Highly variable peaks will be prioritized because they performed better than simply choosing the most accessible peaks in the deconvATAC benchmark [3].

For each fragment, I will calculate its length and initially assign it to one of three bins:

- Short: less than 100 base pairs
- Mononucleosome-like: 100–249 base pairs
- Long: 250 base pairs or longer

The exact cutoffs will be checked in a sensitivity analysis because fragment-size distributions can differ by dataset and protocol. For each cell type, I will then estimate a peak-accessibility and fragment-length fingerprint from the reference cells. Starting with only three bins will limit sparsity and keep the first model feasible. Fragment counts will be retained rather than binarized because count modeling preserves quantitative ATAC information [4].

### 3. Creating Pseudo-Spatial Spots

I will construct simulated spatial spots by sampling and summing labeled test cells. The initial benchmark will contain 1,024 spots with an average of approximately 10 cells per spot. I will create both balanced simulations and realistic imbalanced simulations in which some cell types are rare. Additional experiments will vary sequencing depth, cells per spot, cell-type similarity, and rare-cell abundance. Each spot's sampled cells will be recorded to create a ground-truth proportion table.

### 4. Count-Only Baseline and ShapeMix-ATAC Model

The count-only baseline will use one value for every spot and peak: the total number of fragments in that peak. ShapeMix-ATAC will instead use a three-dimensional count tensor:

> **observed data = spot × peak × fragment-length bin**

For each spot, the model will estimate the nonnegative abundance of each cell type whose reference fingerprint best reconstructs the observed counts. A negative-binomial likelihood will be used to represent noisy sequencing counts, following the general count-mixture strategy used by Bayesian deconvolution methods such as Cell2location [5]. RCTD's Poisson mixture model and correction for differences between reference and spatial technologies provide an additional design comparison [6]. The ShapeMix-ATAC model will include a sequencing-depth offset, a background component for fragments not explained by the reference, and—if cross-protocol bias is detected—a simple reference-to-target scaling term. Estimated abundances will be normalized to produce cell-type proportions that sum to one.

I will first fit a minimal version by maximum a posteriori optimization because it is easier to implement and debug. If time and computational resources permit, I will add variational inference to estimate uncertainty intervals. The same peaks, data splits, simulated spots, likelihood family, optimization budget, and evaluation code will be used for the shape-aware and count-only versions. The only intended difference will be whether the fragment-length bins are retained or collapsed.

### 5. Evaluation and Statistical Analysis

The primary dependent variable will be the error between the true and estimated cell-type proportions. I will measure this using root mean squared error (RMSE) and Jensen–Shannon divergence (JSD), which were also used in the deconvATAC benchmark [3]. Secondary outcomes will include rare-cell detection F1 score, correlation between true and estimated proportions, runtime, and—if uncertainty is implemented—credible-interval coverage.

I will repeat simulations with multiple random seeds and report means, variability, paired differences between models, and bootstrap confidence intervals. I will also examine results separately for common cell types, rare cell types, and related immune subtypes. A successful result would be a consistent reduction in RMSE or JSD for ShapeMix-ATAC without an unacceptable increase in false-positive rare-cell calls. If the two models perform similarly, that will also be informative because it would show that fragment length does not add enough signal under the tested conditions.

The matched shape-versus-count ablation will remain the primary test. As a secondary benchmark, I will also compare ShapeMix-ATAC with established count-only implementations available in the project, including nonnegative least squares, Cell2location, RCTD, and SpatialDWLS [3, 5, 6, 8]. This will show whether any gain is meaningful relative to existing methods rather than only relative to a weak baseline.

### 6. Robustness and Optional Real-Data Validation

I will test whether any improvement remains after changing fragment-length cutoffs, peak number, spot depth, and cells per spot. I will inspect residuals and reconstructed fragment-length distributions to identify reference mismatch or technical bias. If the simulated benchmarks support the hypothesis, I will apply the method to a public spatial ATAC or spatial multi-omics dataset and compare its predicted spatial cell-type patterns with known tissue regions and marker peaks. Newer assays such as spatial-Mux-seq demonstrate that chromatin accessibility can be measured together with other molecular layers in the same tissue, which may provide useful cross-modality evidence for qualitative validation [9]. Because exact proportions are not available for most real tissues, this step will be treated as qualitative validation rather than the primary accuracy test.

Spatial smoothing will be considered only after the fragment-length model works independently at each spot. SONAR showed that adaptively borrowing information from similar neighboring spots can improve deconvolution while reducing bias across sharp tissue boundaries [7]. A later ShapeMix-ATAC version could use that principle, but the first experiment will exclude spatial smoothing so it does not confound the test of fragment shape.

## Variables and Controls

| Category | Definition |
| --- | --- |
| **Primary independent variable** | Input representation: peak-count only versus peak count plus fragment-length bins |
| **Secondary independent variables** | Sequencing depth, number of cells per spot, peak set, fragment-length cutoffs, cell-type similarity, and rare-cell abundance |
| **Primary dependent variables** | RMSE and JSD between estimated and true cell-type proportions |
| **Secondary dependent variables** | Rare-cell F1 score, proportion correlation, runtime, and uncertainty calibration/coverage |
| **Controlled factors** | Reference/test split, simulated spots, selected peaks, model likelihood, optimizer settings, random seeds, and evaluation code |

## Novelty, Feasibility, and Limitations

The proposed novelty is not simply applying an RNA deconvolution method to ATAC peak counts. ShapeMix-ATAC will test whether an ATAC-specific physical measurement—fragment length—adds useful evidence to the deconvolution likelihood. The project is feasible within seven months because it is computational, uses public data, and does not require new wet-lab experiments or reagents. The existing deconvATAC workspace already provides annotated PBMC multiome matrices, a pseudo-spot simulator, count-only deconvolution baselines, and evaluation metrics. The main new work will be processing raw fragments, constructing the shape-count tensor, implementing ShapeMix-ATAC, and conducting controlled benchmarks.

The main risks are sparse per-peak shape counts, technical differences between reference and spatial assays, and increased computation. I will address these by starting with only three fragment-length bins, selecting well-covered peaks, using shrinkage or peak grouping if necessary, comparing matched and mismatched data conditions, and treating motif footprints or spatial smoothing as optional extensions rather than requirements for the first experiment.

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
