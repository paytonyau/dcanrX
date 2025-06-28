# dcanr: Differential Co-expression Analysis in R
### Overview

`dcanr` is an R package that provides a unified, robust implementation of three distinct methods for differential co-expression analysis (DCEA). These algorithms were originally developed for microarray data and are designed to uncover how gene-gene relationships are rewired between two experimental conditions (e.g., normal vs. disease).

A key feature of `dcanr` is that all three implemented methods leverage the **biweight midcorrelation (`bicor`)**, a robust alternative to Pearson correlation that is less sensitive to outliers. This makes the analyses more reliable, especially with noisy biological data.

This package makes these powerful algorithms accessible for use on modern RNA-sequencing data by providing a self-contained set of tools with built-in parallel processing to ensure scalability.

### The `dcanr` Toolbox

The three implemented methods offer complementary views of differential co-expression, allowing researchers to move from systems-level discovery to specific, testable hypotheses.

1.  **BMKC `(run_bmkc)`**: Identifies entire functional **modules** of genes that are co-regulated in one condition but not the other. *Best for a high-level view of pathway dysregulation.*
2.  **BMHT `(run_bmht)`**: Ranks individual **genes** based on their overall change in co-expression connectivity. *Best for identifying key "hub" genes that are significantly rewired.*
3.  **ROS-DET `(run_rosdet)`**: Finds specific gene **pairs** that "switch" their correlation from positive to negative. *Best for pinpointing specific, interpretable interaction changes.*

### Comparison of Methods

| Feature | BMKC (Module-Level) | BMHT (Gene-Level) | ROS-DET (Pair-Level) |
| :--- | :--- | :--- | :--- |
| **Core Similarity Metric**| Biweight Midcorrelation | Biweight Midcorrelation | Biweight Midcorrelation |
| **Primary Goal** | Discover functional modules of genes that lose or gain co-regulation. | Rank individual genes by their overall change in network connectivity. | Identify specific gene pairs that switch correlation from positive to negative. |
| **Final Output** | A set of **gene modules** (cliques). | A ranked list of individual **genes**. | A ranked list of significant **gene pairs**. |
| **Key Strength** | Provides a systems-level view of pathway disruption. | Identifies the most "rewired" hub genes that may be key drivers. | Finds highly specific, interpretable relationships ideal for generating hypotheses. |

### Installation

You can install the development version of `dcanr` from GitHub. First, ensure you have `devtools` installed.

```r
# install.packages("devtools")
devtools::install_github("paytonyau/dcanr")
```


### Important Note on Data Preparation

Correlation-based methods are sensitive to the properties of the input data. For meaningful results, especially with RNA-seq data, please ensure your expression matrix is properly processed:

1.  **Filtering:** Remove genes with low expression/variance across all samples.
2.  **Normalization & Transformation:** Raw counts from RNA-seq **must** be normalized for library size and transformed to stabilize variance. We strongly recommend using either the **Variance Stabilizing Transformation (VST)** or `rlog` from the `DESeq2` package, or the **`voom`** method from the `limma` package.

### Analysis Workflow: A Complete Example

This example demonstrates the intended workflow, from data simulation to running and visualizing the results from all three methods.

#### 1. Load Libraries and Prepare Data

```r
library(dcanr)
library(ggplot2)
library(igraph)

# --- Create a simulated dataset ---
set.seed(42)
n_genes <- 100
n_samples <- 60 # 30 per condition
expr_matrix <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
rownames(expr_matrix) <- paste0("Gene", 1:n_genes)
condition <- factor(c(rep("Normal", 30), rep("Tumor", 30)))

# --- Engineer patterns for each model to find ---
# 1. A switching pair for ROS-DET (Gene5, Gene6)
expr_matrix["Gene6", condition == "Normal"] <- expr_matrix["Gene5", condition == "Normal"] + rnorm(30, 0, 0.5) # Positive
expr_matrix["Gene6", condition == "Tumor"]  <- -expr_matrix["Gene5", condition == "Tumor"] + rnorm(30, 0, 0.5) # Negative

# 2. A co-regulated module for BMKC/BMHT (Genes 10-15) that is lost in Tumor
module_genes <- paste0("Gene", 10:15)
module_seed <- rnorm(30)
for (gene in module_genes) {
  # Strong correlation in Normal
  expr_matrix[gene, condition == "Normal"] <- module_seed + rnorm(30, sd = 0.1)
  # Random noise in Tumor
  expr_matrix[gene, condition == "Tumor"] <- rnorm(30)
}
```

#### 2. Run the Analyses

We will run all three functions using parallel processing to speed up the computations.

```r
# Use 2 workers for parallel processing
N_WORKERS <- 2

# --- Method 1: BMKC ---
# First, explore thresholds to find the optimal parameters
explorer_results <- explore_bmkc_thresholds(
  expr_matrix, condition,
  T1_range = seq(0.7, 0.9, by = 0.1),
  T2_range = seq(0.3, 0.5, by = 0.1),
  workers = N_WORKERS
)

# Visualize the exploration results (optional but recommended)
plot_bmkc_exploration(explorer_results)

# Run the final analysis with chosen thresholds
bmkc_results <- run_bmkc(
  expr_matrix, condition,
  T1 = 0.8, T2 = 0.4, min_clique_size = 5,
  workers = N_WORKERS
)

# --- Method 2: BMHT ---
bmht_results <- run_bmht(
  expr_matrix, condition,
  n_permutations = 500, # Use >=1000 for a real analysis
  workers = N_WORKERS
)

# --- Method 3: ROS-DET ---
rosdet_results <- run_rosdet(
  expr_matrix, condition,
  workers = N_WORKERS
)
```

#### 3. Interpret and Visualize Results

##### **BMKC Module Results**

The BMKC analysis identifies that the module of genes from Gene10 to Gene15 loses its tight co-expression in the tumor condition.

```r
cat("--- BMKC Results ---\n")
print(bmkc_results$modules)
#> [[1]]
#> [1] "Gene10" "Gene11" "Gene12" "Gene13" "Gene14" "Gene15"

# Visualize the largest module's connectivity change
# (Code for this plot is in the package documentation)
```
**Interpretation:** This result suggests that the biological pathway represented by Genes 10-15 is functionally intact in normal tissue but becomes dysregulated and disorganized in tumor tissue.

##### **BMHT Gene Ranking Results**

The BMHT results rank the most "rewired" genes. As expected, the genes from the dysregulated module (10-15) are ranked at the top.

```r
cat("\n--- BMHT Results (Top 6) ---\n")
print(head(bmht_results$results))
#>      Gene  DC_Score P_value   FDR
#> 15 Gene15 0.7719658   0.000 0.000
#> 14 Gene14 0.7679383   0.000 0.000
#> 13 Gene13 0.7621118   0.000 0.000
#> 12 Gene12 0.7588325   0.000 0.000
#> 11 Gene11 0.7554587   0.000 0.000
#> 10 Gene10 0.7495066   0.000 0.000
```
**Interpretation:** This tells us that Genes 10-15 have undergone the most significant changes to their network connectivity, marking them as key players in the disease's regulatory rewiring.

##### **ROS-DET Switching Pair Results**

The ROS-DET analysis pinpoints the specific gene pair that flips its correlation sign.

```r
cat("\n--- ROS-DET Results (Top Result) ---\n")
print(head(rosdet_results$results, 1))
#>   Gene1 Gene2     Score      P_value        r1         r2 c_weight
#> 1 Gene5 Gene6 0.9022416 1.838528e-05 0.6387069 -0.6698188 0.6894086
```
**Interpretation:** This result provides a highly specific and testable hypothesis: the relationship between Gene5 and Gene6 is fundamentally inverted in tumors. This could suggest that a shared regulator is lost or a new one is gained, leading to this dramatic switch in behavior.

### Bug Reports and Contributions

If you encounter a bug or have a suggestion for improvement, please open an issue on the GitHub repository page. We welcome contributions! Please feel free to fork the repository and submit a pull request.

### Citation

If you use the methods implemented in this package for your research, please cite the original publications.

### License

This package is licensed under the MIT License.
