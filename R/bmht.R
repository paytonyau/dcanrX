#' @importFrom WGCNA bicor
#' @importFrom stats p.adjust
NULL

#' Run the BMHT analysis to find differentially coexpressed genes.
#'
#' Implements the BMHT (Biweight Midcorrelation and Half-Thresholding) method.
#' It ranks individual genes based on their overall change in coexpression connectivity.
#'
#' @param expr_matrix A numeric matrix where rows are genes and columns are samples.
#' @param condition A factor or character vector with two levels, matching the columns of expr_matrix.
#' @param half_threshold The correlation threshold for the 'half-thresholding' strategy.
#' @param n_permutations The number of permutations to run for statistical significance testing.
#' @return A data frame with genes ranked by their differential coexpression ('dc') value.
#' @export
run_bmht <- function(expr_matrix, condition, half_threshold = 0.4, n_permutations = 1000) {

  conditions <- unique(condition)
  if (length(conditions) != 2) stop("Condition vector must have exactly two levels.")

  cond1 <- conditions[1]; cond2 <- conditions[2]

  # Calculate bicor matrices for both conditions
  cor_1 <- WGCNA::bicor(t(expr_matrix[, condition == cond1]))
  cor_2 <- WGCNA::bicor(t(expr_matrix[, condition == cond2]))
  diag(cor_1) <- 0; diag(cor_2) <- 0

  # Half-thresholding strategy
  informative_mask <- abs(cor_1) > half_threshold | abs(cor_2) > half_threshold

  # Function to calculate dc scores
  calculate_dc_scores <- function(c1, c2, mask) {
    sapply(1:nrow(c1), function(i) {
      informative_neighbors <- mask[i, ]
      if (sum(informative_neighbors) == 0) return(0)

      diff_sq <- (c1[i, informative_neighbors] - c2[i, informative_neighbors])^2
      # DC formula from Equation (6)
      sqrt(mean(diff_sq))
    })
  }

  # Calculate observed dc scores
  dc_observed <- calculate_dc_scores(cor_1, cor_2, informative_mask)

  # Permutation test for significance
  perm_dc_matrix <- replicate(n_permutations, {
    perm_condition <- sample(condition)
    perm_cor_1 <- WGCNA::bicor(t(expr_matrix[, perm_condition == cond1]))
    perm_cor_2 <- WGCNA::bicor(t(expr_matrix[, perm_condition == cond2]))
    diag(perm_cor_1) <- 0; diag(perm_cor_2) <- 0
    perm_mask <- abs(perm_cor_1) > half_threshold | abs(perm_cor_2) > half_threshold
    calculate_dc_scores(perm_cor_1, perm_cor_2, perm_mask)
  })

  # Calculate p-values based on permutation results
  p_values <- sapply(1:length(dc_observed), function(i) {
    observed_score <- dc_observed[i]
    perm_scores <- perm_dc_matrix[i, ]
    sum(perm_scores >= observed_score) / n_permutations
  })

  # FDR correction
  fdr_values <- stats::p.adjust(p_values, method = "BH")

  results_df <- data.frame(
    Gene = rownames(expr_matrix),
    DC_Score = dc_observed,
    P_value = p_values,
    FDR = fdr_values
  )

  return(results_df[order(-results_df$DC_Score), ])
}
