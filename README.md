# bicorDCEA

[![R](https://img.shields.io/badge/R-4.2%2B-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-Development-yellowgreen)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Biweight Midcorrelation Differential Co-expression Analysis**

bicorDCEA (Biweight Midcorrelation Differential Co-Expression Analysis) is an R/C++ package designed to identify structural topological vulnerabilities, logic inversions, and rewired causal cliques across two distinct biological states.
---

## 🚀 Key Features

- **Three complementary methods** in one package:
  - **BMKC**: Detects rewired gene **modules** (cliques)
  - **BMHT**: Ranks genes by connectivity **rewiring**
  - **ROS-DET**: Finds specific gene **pairs** that switch correlation sign
- **Native single-cell support** via pseudobulk
- **Microbial-ready**: CLR transformation + **cross-species mode**
- **Robust & Fast**: BiocParallel + S4 methods + comprehensive validation
- **Modern Bioconductor design**

---

## Installation

```r
# Development version (GitHub)
devtools::install_github("yourusername/bicorDCEA")

# After Bioconductor acceptance
BiocManager::install("bicorDCEA")
```

---

## Quick Start

### Bulk RNA-seq

```r
library(bicorDCEA)

# Simulated data
set.seed(42)
expr_matrix <- matrix(rnorm(100 * 60), nrow = 100)
rownames(expr_matrix) <- paste0("Gene", 1:100)
condition <- factor(rep(c("Normal", "Tumor"), each = 30))

# Unified interface
bmht_res <- run_dcanr(expr_matrix, condition, method = "bmht")
rosdet_res <- run_rosdet(expr_matrix, condition, transform = "none")
bmkc_res <- run_bmkc(expr_matrix, condition, T1 = 0.8, T2 = 0.4)

summary(bmht_res)
plot(bmht_res)
```

### Single-Cell RNA-seq

```r
# One-step pseudobulk + analysis
res <- run_dcanr_pseudobulk(sce,
                           condition_col = "treatment",
                           sample_col = "donor",
                           method = "bmht",
                           transform = "clr")   # recommended for scRNA-seq
```

### Microbial / Cross-Species Analysis

```r
# Example with species labels
species_vec <- c(rep("Ecoli", 2500), rep("Pseudomonas", 2200), rep("Staph", 1800))

res <- run_rosdet(expr_matrix, condition,
                  transform = "clr",           # critical for microbes
                  species = species_vec,
                  pair_type_filter = "inter")  # only inter-species pairs

summary(res)
plot(res)
```

---

## New Features (Latest Update)

| Feature                     | Description                                      | Recommended For               |
|---------------------------|--------------------------------------------------|-------------------------------|
| **CLR Transformation** | Compositional data handling                      | Microbial, metatranscriptomics, scRNA-seq |
| **Cross-Species Mode** | Intra- vs Inter-species analysis                 | Bacteria-bacteria interactions |
| **Transform Parameter** | `"none"`, `"clr"`, `"log1p"`                     | All users                     |
| **Species Vector** | Automatic pair/module type annotation            | Multi-species studies         |

---

## Methods Comparison

| Method     | Level       | Best For                               | Key Output               | Cross-Species Support |
|------------|-------------|---------------------------------------|--------------------------|-----------------------|
| **ROS-DET** | Pair        | Specific interaction switches         | Ranked gene pairs        | Yes (intra/inter)    |
| **BMHT** | Gene        | Rewired hub genes                     | Ranked gene list         | Yes                   |
| **BMKC** | Module      | Pathway/module rewiring               | Gene cliques             | Yes                   |

---

## Data Preparation Recommendations

**Bulk RNA-seq**: Use VST or voom normalized data  
**Single-cell**: Use `make_pseudobulk()` or `run_dcanr_pseudobulk()`  
**Microbial data**: Always use `transform = "clr"`

---

## Visualization

```r
plot(rosdet_res)           # Top switching pair
plot(bmht_res)             # Top rewired genes
plot_bmkc_exploration(explorer)  # Threshold optimization
```

---

## Why Choose bicorDCEA?

- **Most robust correlation** (`bicor` vs Pearson)
- **True differential focus** (condition rewiring)
- **Three complementary views** in one package
- **Built for noisy data** (single-cell & microbes)
- **Modern Bioconductor standards** (S4, BiocParallel, SummarizedExperiment)

---

## Citation

If you use **bicorDCEA**, please cite:

> Yau et al. (2026). bicorDCEA: Robust differential co-expression analysis using biweight midcorrelation. R package.

---

## Contributing & Feedback

- Bug reports / feature requests: [Issues](https://github.com/yourusername/bicorDCEA/issues)
- Pull requests are welcome!

---

**License**: MIT  
**Last updated**: April 2026

---

**Ready for robust differential co-expression analysis in bulk, single-cell, and microbial systems.**