# Dataset sizes used in the ShapeMix evaluation

"Dataset size" can mean the number of source cells, experimental samples, or
evaluated spots. These quantities are different and should not be treated as
interchangeable.

| Dataset | Source size | Evaluation size |
|---|---:|---:|
| **10x PBMC 10k** | 9,500 labeled cells from one donor | 20 pseudo-spot datasets x 1,024 spots = **20,480 pseudo-spots** |
| **GSE194122 BMMC** | 69,249 annotated cells from ten donors | 40 pseudo-spot datasets x 1,024 spots = **40,960 pseudo-spots** |
| **GSE129785** | 14 physical mixtures plus nine sorted reference samples | **14 evaluated mixtures**: seven CD4-memory/CD8-naive and seven monocyte/T-cell mixtures |
| **GSE205055** | Six spatial tissue sections | 29,966 input spots; **29,901 evaluated spots** after excluding 65 zero-signal spots |
| **GSE263333** | Two spatial tissue sections | 12,500 input spots; **12,430 evaluated spots** after excluding 70 zero-signal spots |

## Important distinctions

- The 20,480 PBMC and 40,960 BMMC pseudo-spots are simulated observations,
  not independent biological samples.
- PBMC has only **one biological donor**.
- GSE194122 has **ten biological donors**, which supports donor-level
  evaluation.
- Each GSE129785 physical mixture produces one aggregate prediction rather
  than thousands of spatial-spot predictions.
- For GSE205055 and GSE263333, the tissue section is the independent sample;
  spatial spots are nested within each section.

The GSE205055 and GSE263333 input totals include all registered spatial spots.
The evaluated totals exclude spots with zero ATAC signal on the frozen feature
axis from spotwise comparisons. The raw predictions were not changed.

## Why the physical-dilution difference looks much larger

For this comparison, the method labels in
`shapemix_vs_count_only_rmse_jsd.tsv` should be used as written: the values in
the `ShapeMix` row are the ShapeMix results, and the values in the `Count-only`
row are the matched count-only results. Because lower RMSE and JSD are better,
the positive percentages below are reductions in error achieved by ShapeMix:

| Evaluation | RMSE: count-only -> ShapeMix | RMSE reduction | JSD: count-only -> ShapeMix | JSD reduction |
|---|---:|---:|---:|---:|
| PBMC, natural cell frequency | 0.046757 -> 0.046732 | 0.1% | 0.128065 -> 0.122298 | 4.5% |
| PBMC, equal cell frequency | 0.057626 -> 0.052396 | 9.1% | 0.183554 -> 0.156501 | 14.7% |
| GSE194122, natural cell frequency | 0.099013 -> 0.098945 | 0.1% | 0.096842 -> 0.096544 | 0.3% |
| GSE194122, equal cell frequency | 0.112585 -> 0.111887 | 0.6% | 0.124303 -> 0.122322 | 1.6% |
| GSE129785, CD4-memory/CD8-naive mixtures | 0.213103 -> 0.163002 | 23.5% | 0.154030 -> 0.114483 | 25.7% |
| GSE129785, monocyte/T-cell mixtures | 0.065097 -> 0.049937 | 23.3% | 0.040204 -> 0.029264 | 27.2% |

Thus, ShapeMix changes the pseudo-spot errors only modestly in most settings
but reduces both errors by about 23%–27% for the physical dilution mixtures.
Several differences between the evaluations could contribute to this larger
observed improvement:

1. **The truth labels are different in precision.** Each pseudo-spot has exact
   composition truth because the identities of all source cells used to build
   it are known. A physical dilution has a planned input ratio, but sorting,
   handling, library preparation, and recovery can make the molecules measured
   in the final library differ from that nominal ratio.
2. **The physical mixtures are a different measurement setting.** They include
   sorting, mixing, library-preparation, and recovery effects that are absent
   from computer-generated pseudo-spots. These effects may leave ambiguity that
   peak counts alone do not resolve. If fragment-length information supplies
   complementary signal under these conditions, it can produce a larger gain
   over count-only.
3. **The evaluation labels are broader than the prediction labels.** For the
   monocyte/T-cell series, the known target gives only the ratio between
   monocytes and all T cells combined, whereas the reference predicts five
   separate T-cell subtypes plus other possible cell types. The evaluation
   must sum the five T-cell
   predictions and account for off-target predictions. Improvements in broad
   lineage allocation can therefore outweigh mistakes among individual T-cell
   subtypes and produce a substantial change in the aggregate score.
4. **The physical-mixture result is based on very few samples.** Each dilution
   series has only seven mixtures. A few mixtures with large improvements can
   therefore have a strong effect on the reported mean. The pseudo-spot
   averages contain tens of thousands of spot-level observations, although
   those pseudo-spots should not be mistaken for independent biological
   samples.
5. **Relative percentages depend on the baseline.** A similar absolute change
   can produce a larger percentage when the count-only error is small. These
   percentages describe within-dataset method changes and should not be read as
   direct measures of dataset difficulty.

These explanations are plausible interpretations, not proven causes. The
aggregate RMSE and JSD values alone cannot show that physical dilution itself
causes ShapeMix to perform better. Identifying the main reason would require
examining per-mixture errors, off-target predictions, and fragment-length
compatibility between each mixture and its reference.

## Why the percentage change in JSD is larger than in RMSE

A larger percentage change in JSD does **not** mean that its absolute error
changed by more than RMSE. RMSE and JSD use different formulas and scales, so
their raw values—and even their percentage changes—should not be compared as
if they had the same units.

- **RMSE measures entry-by-entry proportion errors.** It squares the difference
  between the predicted and target proportion for every cell type, averages
  those squared differences, and then takes the square root. It mainly answers,
  "On average, how far is each predicted cell-type fraction from its target?"
- **JSD compares each complete composition distribution.** It measures how
  probability mass is distributed across all cell types. It can react strongly
  when mass moves between cell types, when a rare type is missed, or when mass
  is assigned to a cell type that should be absent—even if each individual
  proportion changes by only a modest amount.
- **The percentage uses a metric-specific baseline.** For an error reduction,
  the calculation is
  `(count-only error - ShapeMix error) / count-only error x 100%`.
  Because RMSE and JSD have different starting values and respond differently
  to the same prediction changes, their percentages need not be similar.

For example, a method could move several small amounts of probability from
incorrect or off-target cell types to the correct cell types. The average
entrywise change may be small, producing only a small RMSE change, while the
overall composition distribution becomes noticeably better, producing a
larger relative JSD reduction. The reverse is also possible: redistribution
into wrong or off-target cell types can increase JSD proportionally more than
RMSE.

In these results, the JSD reduction is larger than the RMSE reduction in every
displayed condition. For example, in the PBMC natural-frequency evaluation,
ShapeMix reduces RMSE by about 0.1% and JSD by 4.5%. In the physical
monocyte/T-cell mixtures, it reduces RMSE by 23.3% and JSD by 27.2%. This
suggests that ShapeMix improves how probability is distributed across the full
set of cell types more than it changes the average entrywise proportion error.
It does not mean that JSD is more important, or that a JSD percentage can be
directly compared with an RMSE percentage as an absolute quantity.

## Concise wording for a slide

- **10x PBMC 10k:** 9,500 cells; 20 x 1,024 pseudo-spots
- **GSE194122:** 69,249 cells, ten donors; 40 x 1,024 pseudo-spots
- **GSE129785:** 14 physical mixtures; nine sorted references
- **GSE205055:** six sections; 29,901 evaluated spots
- **GSE263333:** two sections; 12,430 evaluated spots

## Project records

- [PBMC protocol-v1 design](../step6_results.md)
- [Full evaluation inventory and exclusions](../results_summary.md)
- [Current RMSE/JSD comparison table](shapemix_vs_count_only_rmse_jsd.tsv)
- [GSE129785 source configuration](../../../configs/data_sources/shapemix_gse129785.yaml)
- [GSE194122 source configuration](../../../configs/data_sources/shapemix_gse194122.yaml)
