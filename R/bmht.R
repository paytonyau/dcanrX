#' @importFrom stats p.adjust
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply mclapply
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_minimal reorder
NULL

#' Run the BMHT analysis to find differentially coexpressed genes.
#'
#' Implements the BMHT (Biweight Midcorrelation and Half-Thresholding) method.
#' It ranks individual genes based on their overall change in coexpression connectivity.
#' This improved version uses the robust, built-in 'parallel' package for processing.
#'
#' @param expr_matrix A numeric matrix where rows are genes and columns are samples.
#' @param condition A factor or character vector with two levels, matching the columns of `expr_matrix`.
#' @param half_threshold The correlation threshold for the 'half-thresholding' strategy.
#' @param n_permutations The number of permutations to run for statistical significance testing.
#' @param workers An integer specifying the number of parallel workers to use.
#'   Defaults to 1 (no parallelization).
#' @return An object of class `bmht_result` containing the ranked list of genes.
#'   This object can be used with `summary()` and `plot()`.
#' @export
run_bmht <- function(expr_matrix, condition, half_threshold = 0.4, n_permutations = 1000, workers = 1) {

  # --- Input Validation ---
  conditions <- unique(condition)
  if (length(conditions) != 2) stop("Condition vector must have exactly two levels.")

  cond1 <- conditions[1]; cond2 <- conditions[2]

  # Calculate bicor matrices for both conditions once
  cor_1 <- compute_cor_matrix(expr_matrix[, condition == cond1])
  cor_2 <- compute_cor_matrix(expr_matrix[, condition == cond2])
  diag(cor_1) <- 0; diag(cor_2) <- 0

  # Half-thresholding strategy
  informative_mask <- abs(cor_1) > half_threshold | abs(cor_2) > half_threshold

  # Helper function to calculate dc scores
  calculate_dc_scores <- function(c1, c2, mask) {
    sapply(1:nrow(c1), function(i) {
      informative_neighbors <- mask[i, ]
      if (sum(informative_neighbors) == 0) return(0)

      diff_sq <- (c1[i, informative_neighbors] - c2[i, informative_neighbors])^2
      sqrt(mean(diff_sq))
    })
  }

  # Calculate observed dc scores
  dc_observed <- calculate_dc_scores(cor_1, cor_2, informative_mask)

  # --- Permutation test using the 'parallel' package ---
  # This single function defines the work to be done in each permutation
  run_one_permutation <- function(i) {
    perm_condition <- sample(condition)
    perm_cor_1 <- compute_cor_matrix(expr_matrix[, perm_condition == cond1])
    perm_cor_2 <- compute_cor_matrix(expr_matrix[, perm_condition == cond2])
    diag(perm_cor_1) <- 0; diag(perm_cor_2) <- 0
    perm_mask <- abs(perm_cor_1) > half_threshold | abs(perm_cor_2) > half_threshold
    calculate_dc_scores(perm_cor_1, perm_cor_2, perm_mask)
  }

  # Detect OS and choose the appropriate parallel function
  is_windows <- .Platform$OS.type == "windows"

  perm_results_list <- list()
  if (workers > 1 && !is_windows) {
    # Use mclapply on Mac/Linux (more efficient)
    message(paste("Using", workers, "parallel workers with mclapply."))
    perm_results_list <- parallel::mclapply(1:n_permutations, run_one_permutation, mc.cores = workers)
  } else if (workers > 1 && is_windows) {
    # Use parLapply on Windows
    message(paste("Using", workers, "parallel workers with parLapply (Windows)."))
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl))
    # Export necessary functions and objects to the cluster
    parallel::clusterExport(cl, varlist = c("compute_cor_matrix", "expr_matrix", "condition", "half_threshold", "calculate_dc_scores", "bicor"), envir = environment())
    perm_results_list <- parallel::parLapply(cl, 1:n_permutations, run_one_permutation)
  } else {
    # Fallback to sequential execution
    message("Running permutations sequentially.")
    perm_results_list <- lapply(1:n_permutations, run_one_permutation)
  }

  # Combine list of results into a single matrix
  perm_dc_matrix <- do.call(cbind, perm_results_list)

  # Calculate p-values based on permutation results
  p_values <- sapply(1:length(dc_observed), function(i) {
    sum(perm_dc_matrix[i, ] >= dc_observed[i]) / n_permutations
  })

  # FDR correction
  fdr_values <- stats::p.adjust(p_values, method = "BH")

  # Create the results data frame
  results_df <- data.frame(
    Gene = rownames(expr_matrix),
    DC_Score = dc_observed,
    P_value = p_values,
    FDR = fdr_values
  )

  # Create the custom S3 object
  output <- list(
    results = results_df[order(-results_df$DC_Score), ],
    n_permutations = n_permutations,
    half_threshold = half_threshold
  )

  class(output) <- "bmht_result"
  return(output)
}

#' Print a summary of BMHT results.
#'
#' Provides a concise summary of the results from a `run_bmht` analysis.
#'
#' @param object An object of class `bmht_result`.
#' @param fdr_threshold The FDR cutoff to use for counting significant genes. Defaults to 0.05.
#' @param ... Additional arguments (not used).
#' @return Prints a summary to the console, invisibly returning the object.
#' @export
summary.bmht_result <- function(object, fdr_threshold = 0.05, ...) {
  n_significant <- sum(object$results$FDR < fdr_threshold, na.rm = TRUE)

  cat("--- BMHT Analysis Summary ---\n")
  cat(paste("Number of genes tested:", nrow(object$results), "\n"))
  cat(paste("Number of permutations:", object$n_permutations, "\n"))
  cat(paste("Half-thresholding value:", object$half_threshold, "\n"))
  cat("----------------------------------\n")

  if (n_significant > 0) {
    cat(paste("Found", n_significant, "significant genes at FDR <", fdr_threshold, "\n\n"))
    cat("Top 6 results:\n")
    print(head(object$results))
  } else {
    cat("No significant differentially co-expressed genes found at this FDR level.\n")
  }
  invisible(object)
}

#' Plot the top-ranked genes from a BMHT analysis.
#'
#' Generates a bar plot showing the top differentially co-expressed genes,
#' ranked by their DC_Score.
#'
#' @param x An object of class `bmht_result`.
#' @param n_top The number of top genes to display in the plot. Defaults to 15.
#' @param ... Additional arguments (not used).
#' @return A ggplot2 object.
#' @export
plot.bmht_result <- function(x, n_top = 15, ...) {
  if (nrow(x$results) == 0) {
    message("Cannot plot. No results found in the object.")
    return(invisible(NULL))
  }

  # Filter for genes with a non-zero score for a cleaner plot
  top_genes <- head(x$results[x$results$DC_Score > 0, ], n_top)

  if (nrow(top_genes) == 0) {
    message("No genes with a non-zero DC_Score to plot.")
    return(invisible(NULL))
  }

  # Create the plot
  p <- ggplot2::ggplot(top_genes, ggplot2::aes(x = reorder(Gene, DC_Score), y = DC_Score)) +
    ggplot2::geom_col(fill = "firebrick") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Top", nrow(top_genes), "Differentially Co-expressed Genes (BMHT)"),
      x = "Gene",
      y = "Differential Coexpression (DC) Score"
    ) +
    ggplot2::theme_minimal()

  return(p)
}

#' @export
print.bmht_result <- function(x, ...) {
  summary(x, ...)
}
