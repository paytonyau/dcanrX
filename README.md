# dcanr: Differential Co-expression Analysis in R

`dcanr` is an R package that provides a unified implementation of three distinct methods for differential co-expression analysis (DCEA). These methods were originally developed for microarray data and are designed to uncover how gene-gene relationships are rewired between two experimental conditions (e.g., normal vs. tumor).

This package makes these powerful algorithms accessible for use on modern RNA-sequencing data by providing a robust, self-contained set of tools. The three implemented methods offer complementary views of differential co-expression:

1.  **ROS-DET**: Identifies specific gene **pairs** that "switch" their correlation from positive to negative.
2.  **BMHT**: Ranks individual **genes** based on their overall change in co-expression connectivity.
3.  **BMKC**: Discovers entire functional **modules** of genes that are co-regulated in one condition but not the other.

## Comparison of Implemented Methods

| Feature | ROS-DET (Kayano et al., 2011) | BMHT (Zheng et al., 2014) | BMKC (Yuan et al., 2015) |
| :--- | :--- | :--- | :--- |
| **Primary Goal** | To identify "switching mechanisms": gene **pairs** that switch from positive to negative correlation. | To identify and rank individual **genes** based on their overall change in co-expression connectivity. | To identify differentially coexpressed gene **modules** or functional communities. |
| **Core Metric** | **Weighted difference of correlations**. The score for a pair is `c * |r1 - r2|`, where `c` penalizes for differences in expression range (range bias). | **Average change in coexpression**. The score (`dc`) for a gene is the average change in its correlation values against a filtered set of neighbors. | **Threshold-based network construction**. A final network is built by keeping edges only if coexpression is high in one condition and low in the other. |
| **Filtering Strategy** | **ECOR hypothesis test**. A statistical test removes pairs whose correlation change is not statistically significant, specifically addressing small sample sizes. | **'Half-thresholding'**. A gene pair is kept if its correlation is high in *at least one* of the two conditions. Significance is assessed via permutation test. | **"Differential Coexpression Threshold"**. A final binary network is created based on hard thresholds for both conditions simultaneously. |
| **Final Output** | A ranked list of significant **gene pairs** with scores and P-values. | A ranked list of individual **genes** with a differential coexpression score, P-value, and FDR. | A set of **gene modules** (cliques), presented as lists of genes. |
| **Key Strength** | Finds specific, highly interpretable relationships between two genes. Robust to range bias and small sample sizes. | Provides a ranked list of the most "rewired" genes in the network. | Identifies groups of functionally related genes, providing a systems-level view of pathway dysregulation. |

## Installation

You can install the development version of `dcanr` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("payton_yau/dcanr")
```

## Important Note on Data Preparation

Correlation-based methods are sensitive to the properties of the input data. For meaningful results, especially with RNA-seq data, please ensure your expression matrix is properly processed:

1.  **Filtering:** Remove genes with low expression/variance across all samples.
2.  **Normalization & Transformation:** Raw counts from RNA-seq must be normalized for library size and transformed to stabilize variance. We strongly recommend using either the **Variance Stabilizing Transformation (VST)** or `rlog` from the `DESeq2` package, or the **`voom`** method from the `limma` package.

## Quick Start Example

Here is a complete example demonstrating how to use the three main functions on a simulated dataset.

```r
# 1. Load the package
library(dcanr)
library(ggplot2)

# 2. Create a simulated dataset
set.seed(42)
n_genes <- 50
n_samples <- 40 # 20 per condition
expr_matrix <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
rownames(expr_matrix) <- paste0("Gene", 1:n_genes)
colnames(expr_matrix) <- paste0("Sample", 1:n_samples)

condition <- factor(c(rep("Normal", 20), rep("Tumor", 20)))

# Engineer patterns for demonstration:
# - A switching pair for ROS-DET (Gene5, Gene6)
expr_matrix["Gene6", 1:20] <- expr_matrix["Gene5", 1:20] + rnorm(20, 0, 0.5) # Positive
expr_matrix["Gene6", 21:40] <- -expr_matrix["Gene5", 21:40] + rnorm(20, 0, 0.5) # Negative

# - A co-regulated module for BMHT/BMKC (Genes 10-15)
module_genes <- paste0("Gene", 10:15)
# Strong correlation in Normal
for (i in 1:(length(module_genes)-1)) {
  expr_matrix[module_genes[i+1], 1:20] <- expr_matrix[module_genes[i], 1:20] + rnorm(20, 0, 0.4)
}
# Weak correlation in Tumor
expr_matrix[module_genes, 21:40] <- rnorm(length(module_genes) * 20)

# 3. Run the Analyses

# Method 1: ROS-DET to find switching pairs
cat("--- Running ROS-DET ---\n")
rosdet_results <- run_rosdet(expr_matrix, condition, significance_level = 0.05)
print(head(rosdet_results))
#>   Gene1 Gene2     Score      P_value        r1         r2  c_weight
#> 1 Gene5 Gene6 0.9634887 2.112814e-05 0.7029671 -0.7303358 0.6720173

# Method 2: BMHT to rank rewired genes
cat("\n--- Running BMHT ---\n")
bmht_results <- run_bmht(expr_matrix, condition, half_threshold = 0.3, n_permutations = 100)
print(head(bmht_results))
#>     Gene    DC_Score P_value       FDR
#> 10 Gene10  0.5593846    0.00 0.0000000
#> 11 Gene11  0.7226503    0.00 0.0000000
#> 12 Gene12  0.8662283    0.00 0.0000000
#> 13 Gene13  0.9839957    0.00 0.0000000
#> 14 Gene14  1.0827250    0.00 0.0000000
#> 15 Gene15  1.1685820    0.00 0.0000000


# Method 3: BMKC to find co-regulated modules
cat("\n--- Running BMKC ---\n")
bmkc_modules <- run_bmkc(expr_matrix, condition, T1 = 0.6, T2 = 0.4, min_clique_size = 4)
print(bmkc_modules)
#> [[1]]
#> [1] "Gene10" "Gene11" "Gene12" "Gene13" "Gene14" "Gene15"

# 4. Visualize a Top Result from ROS-DET
if (nrow(rosdet_results) > 0) {
  top_pair <- rosdet_results[1, ]
  plot_df <- data.frame(
    Gene1_Expr = as.numeric(expr_matrix[top_pair$Gene1, ]),
    Gene2_Expr = as.numeric(expr_matrix[top_pair$Gene2, ]),
    Condition = condition
  )
  
  ggplot(plot_df, aes(x = Gene1_Expr, y = Gene2_Expr, color = Condition)) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste("ROS-DET Top Switching Pair:", top_pair$Gene1, "&", top_pair$Gene2),
      subtitle = paste0("Correlation switches from ", round(top_pair$r1, 2), " to ", round(top_pair$r2, 2)),
      x = paste(top_pair$Gene1, "Expression"),
      y = paste(top_pair$Gene2, "Expression")
    ) +
    theme_bw()
}
```
*Note: The plot image is generated by the code and will not be saved automatically. To include it in your GitHub README, you will need to save the plot and upload it to your repository.*

## Bug Reports and Contributions

If you encounter a bug or have a suggestion for improvement, please open an issue on the GitHub repository page. We welcome contributions! Please feel free to fork the repository and submit a pull request.

## Citation

If you use the methods implemented in this package for your research, please cite the original publications:

* **For ROS-DET**: Kayano, M., Takigawa, I., Shiga, M., Tsuda, K., & Mamitsuka, H. (2011). ROS-DET: robust detector of switching mechanisms in gene expression. *Nucleic acids research*, 39(11), e74.
* **For BMHT**: Zheng, C. H., Yuan, L., Sha, W., & Sun, Z. L. (2014). Gene differential coexpression analysis based on biweight correlation and maximum clique. *BMC bioinformatics*, 15(15), S3.
* **For BMKC**: Yuan, L., Zheng, C. H., Xia, J. F., & Huang, D. S. (2015). Module based differential coexpression analysis method for type 2 diabetes. *BioMed research international*, 2015.

## License

This package is licensed under the MIT License.
