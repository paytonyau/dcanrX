# bicorX

Robust differential co-expression analysis based on biweight
midcorrelation (bicor), extending three previously-published, never
packaged into R before now, methods for compositional and multi-omics
data:

| Method | Source paper | What it does |
|---|---|---|
| **BMHT** | Zheng, Yuan, Sha & Sun (2014), *BMC Bioinformatics* [doi:10.1186/1471-2105-15-S15-S3](https://doi.org/10.1186/1471-2105-15-S15-S3) | Per-gene hub-rewiring score: which genes' overall correlation structure changes most between two conditions |
| **ROS-DET** | Kayano, Takigawa, Shiga, Tsuda & Mamitsuka (2011), *Nucleic Acids Research* [doi:10.1093/nar/gkr130](https://doi.org/10.1093/nar/gkr130) | Per-pair "switching mechanism" detection: which specific gene pairs flip from positively to negatively correlated between conditions |
| **BMKC** | Yuan, Zheng, Xia & Huang (2015), *BioMed Research International* [doi:10.1155/2015/836929](https://doi.org/10.1155/2015/836929) | Module-level detection: which groups of genes form a differentially co-expressed community |

Run `citation("bicorX")` for full citation details, including notes on
where this implementation departs from each paper (see "Honest scope"
below).

## Why bicor instead of Pearson correlation?

Biweight midcorrelation down-weights outliers instead of letting them
dominate the correlation estimate the way Pearson's does. This is a real,
measured advantage, not a marketing claim: under injected sample
corruption, a Pearson-based competitor's precision collapsed from 0.68 to
0.38 while bicorX's held steady or improved (see `NEWS.md`, Phase 5
benchmark).

## Installation

```r
# install.packages("remotes")
remotes::install_github("paytonyau/bicorX")

# To install with vignettes (recommended - they are the main
# documentation; see "Where to look next" below):
remotes::install_github("paytonyau/bicorX", build_vignettes = TRUE)
```

Core dependencies (Rcpp, RcppArmadillo, matrixStats, progress, igraph)
install automatically from CRAN. Optional dependencies for specific input
types (SummarizedExperiment, SingleCellExperiment, Seurat,
MultiAssayExperiment, scuttle) are listed under `Suggests` and are only
needed if you use those input types; several are from Bioconductor, so
install them with `BiocManager::install()` if `install_github` does not
pick them up.

To install from a local source tarball instead:

```r
install.packages("bicorX_0.0.1.tar.gz", repos = NULL, type = "source")
```

**Note on vignettes**: build the package with `R CMD build` (or use
`build_vignettes = TRUE` above) rather than running `R CMD INSTALL`
directly on a cloned source directory — the latter skips vignette
building, leaving `vignette(package = "bicorX")` empty.

**Performance tip**: this package's runtime is dominated by a
correlation-matrix computation that benefits substantially from an
optimized BLAS. Installing OpenBLAS instead of R's default reference BLAS
gives a measured ~1.3-1.4x speedup at zero code changes - see
`DESCRIPTION`'s `SystemRequirements` field for the exact command.

## Quick start

```r
library(bicorX)

# expr_matrix: genes (rows) x samples (columns), any numeric data
# condition: a factor with exactly two levels
res <- run_bmht(expr_matrix, condition, n_permutations = 1000)
head(res$results[order(res$results$P_value), ])

# Per-pair "switching" detection, with a real ~30-50x speedup available
# when you have >=30-50 samples per condition:
res <- run_rosdet(expr_matrix, condition, significance_method = "analytical")

# Module-level detection, with an optional significance test:
res <- run_bmkc(expr_matrix, condition, n_permutations = 1000)
```

See `vignettes/multiomics-workflow.Rmd` and `vignettes/custom-plots.Rmd`
for worked examples, including real (not just synthetic) data.

## What data this supports

| | Supported via |
|---|---|
| Plain matrix / data.frame | direct input |
| Bulk RNA-seq, microarray, any continuous data | `transform = "none"` |
| Compositional data (microbiome 16S/shotgun, ATAC) | `transform = "clr"` / `"rclr"` |
| Single-cell / spatial transcriptomics | `make_pseudobulk()` to aggregate to sample level first |
| Multi-omics / cross-domain (RNA+ATAC, host+microbe) | `run_multiomic()`, or the `layer`/`species` argument directly |
| Bioconductor objects | `SummarizedExperiment`, `SingleCellExperiment`, `Seurat`, `MultiAssayExperiment` all dispatch automatically via `run_bicordcea()` |

**Real limitation, not a bug**: exactly two conditions per comparison,
not multi-group designs.

## Two ways to get a significance value for ROS-DET, and when to use which

| | `significance_method = "permutation"` (default) | `significance_method = "analytical"` |
|---|---|---|
| Assumption | None - exact by construction | Approximate bivariate-normal-derived (Kayano et al.'s own "ECOR" test) |
| Speed | Scales with `n_permutations` | Fixed, small cost per pair regardless of desired resolution |
| Calibration | Always valid | Well-calibrated for **N >= 50 samples/condition**; degrades below that (real, but bounded - not catastrophic) |
| Use when | Small samples, or you want zero distributional assumptions | N >= 30-50 per condition and you want continuous, well-resolved p-values fast |

This isn't a hypothetical tradeoff - on a real published ground-truth
benchmark, switching from permutation to analytical significance testing
took ROS-DET's ranking quality from roughly 3.5x worse than two
established competitors (DGCA, dcanr) to roughly on par with them. See
`NEWS.md`'s "Phase A" entry and `dev/benchmarks/03_sim102_ground_truth.R`.

## What's actually been validated (not just implemented)

This project's development process (documented in full in `NEWS.md`)
put a premium on checking claims before making them, including reporting
several attempted fixes that failed and were rejected rather than
quietly discarded. Headline validated results:

- **Outlier robustness**: real, replicated advantage over Pearson-based
  competitors under sample corruption (Phase 5 benchmark).
- **Competitive ranking quality**: ROS-DET with
  `significance_method = "analytical"` is roughly on par with DGCA/dcanr
  on a real published ground-truth simulation (Phase A) - a fixable
  implementation gap (permutation testing where the source paper
  specifies an analytical test), not an inherent limitation.
- **Compositional correction, on real data**: on a real published lung
  microbiome dataset (33 vs. 33 patients), CLR correction changed a
  specific, spurious naive finding (driven by a per-sample sequencing-
  depth artifact) into a correct non-finding - not just a synthetic-data
  demonstration (Phase E).
- **Multi-omics, on real data**: on real TCGA cancer RNA-seq + miRNA-seq
  data, found a statistically massive (Fisher's exact p = 5e-13), formally
  tested cross-layer pattern - not an eyeballed proportion (see
  `dev/real_data_validation/acc_multiomics_check.R`).
- **BMHT does not show the same competitiveness gap as ROS-DET** did
  before its fix - checked directly against real ground truth, not
  assumed either way.

## Honest scope - what this is not

- **Not a strictly faithful reimplementation of any of the three source
  papers.** Real, documented departures exist (e.g. k-core decomposition
  in place of BMKC's exact k-clique percolation, for scalability). See
  `citation("bicorX")` for a per-method summary and `NEWS.md` for full
  detail.
- **Not validated at genome-wide scale on real data yet** - the real-data
  validations above use hundreds of features, not tens of thousands.
- **BMHT's own permutation-resolution characteristics have not been
  fixed the way ROS-DET's were** - Zheng et al.'s BMHT paper doesn't
  specify an analytical alternative, so this remains an open item if a
  gap is ever found there.
- **Package name has not been checked against the real CRAN/Bioconductor
  namespace** - only against GitHub and a partial Debian mirror, due to
  this development environment's network restrictions. Verify before any
  public release/submission.

## Where to look next

- `NEWS.md` - the full, chronological record of what was built, found,
  fixed, and rejected, including several negative results reported
  honestly rather than discarded.
- `TODO.md` - current open items.
- `vignettes/` - worked examples, including real-data validations.
- `dev/benchmarks/` and `dev/real_data_validation/` - reproducible
  scripts behind every validated claim above.
