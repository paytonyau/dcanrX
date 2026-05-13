#' @title BMHT: Biweight Midcorrelation + Half-Thresholding Method
#' @description Ranks individual genes by their differential co-expression
#' connectivity (rewiring) between two conditions using bicor.
#' Supports CLR transformation for microbiome and cross-species mode.
#'
#' @importFrom stats p.adjust
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam SnowParam bpparam
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_minimal
#' @name bhmt
NULL

#' Run BMHT Analysis
#'
#' @param expr_matrix Numeric matrix (genes × samples)
#' @param condition Factor/character vector with exactly two levels
#' @param half_threshold Minimum absolute correlation to consider an edge informative
#' @param n_permutations Number of permutations for significance testing
#' @param significance_level FDR threshold (default 0.05)
#' @param workers Number of parallel workers
#' @param transform Compositional transformation ("none", "clr", "log1p")
#' @param species Optional character vector indicating species for each gene
#' @param seed Integer seed for exact reproducibility across parallel runs
#' @return An object of class `bmht_result`
#' @export
run_bmht <- function(expr_matrix, condition,
                     half_threshold = 0.4,
                     n_permutations = 1000,
                     significance_level = 0.05,
                     workers = 1,
                     transform = c("none", "clr", "log1p"),
                     species = NULL,
                     seed = 42) {

  # Ensure strict reproducibility for baseline operations
  if (!is.null(seed)) set.seed(seed)

  condition <- .validate_input(expr_matrix, condition)
  transform <- match.arg(transform)

  # Strict validation for species mode to prevent silent failures
  species_mode <- FALSE
  if (!is.null(species)) {
    if (length(species) != nrow(expr_matrix)) {
      stop("Length of `species` must exactly match the number of rows in `expr_matrix`.")
    }
    species_mode <- TRUE
    species <- as.character(species)
  }

  conditions <- levels(condition)

  # Prevent variance calculation crashes from highly unbalanced/small groups
  n_c1 <- sum(condition == conditions[1])
  n_c2 <- sum(condition == conditions[2])
  if (n_c1 < 3 || n_c2 < 3) {
    stop("BMHT requires at least 3 samples per condition to compute valid midvariances.")
  }

  # Apply compositional transformation (Handles Microbiome / Sparsity)
  expr_matrix <- .compositional_transform(expr_matrix, transform = transform, species = species)

  # Compute observed correlation matrices
  message("Computing observed bicor matrices...")
  cor_1 <- compute_cor_matrix(expr_matrix[, condition == conditions[1]], workers = workers)
  cor_2 <- compute_cor_matrix(expr_matrix[, condition == conditions[2]], workers = workers)

  diag(cor_1) <- diag(cor_2) <- 0

  # Safe half-thresholding mask (handling potential NAs properly)
  informative_mask <- (abs(cor_1) > half_threshold) | (abs(cor_2) > half_threshold)
  informative_mask[is.na(informative_mask)] <- FALSE

  # Safe DC score calculation
  calculate_dc_scores <- function(c1, c2, mask) {
    vapply(seq_len(nrow(c1)), function(i) {
      idx <- mask[i, ]
      if (sum(idx, na.rm = TRUE) == 0) return(0)
      diffs <- c1[i, idx] - c2[i, idx]
      # na.rm = TRUE ensures genes don't drop out due to a single NA edge
      sqrt(mean(diffs^2, na.rm = TRUE))
    }, numeric(1))
  }

  dc_observed <- calculate_dc_scores(cor_1, cor_2, informative_mask)

  # ===================================================================
  # Optimized Permutation Test (Windows-Safe Parallel)
  # ===================================================================

  # Extract only the required edges to avoid matrix RAM bloat
  # and prevent statistical bias from re-evaluating the mask under null conditions.
  mask_upper <- informative_mask
  mask_upper[lower.tri(mask_upper, diag = TRUE)] <- FALSE
  pair_indices <- which(mask_upper, arr.ind = TRUE)
  n_pairs <- nrow(pair_indices)
  n_genes <- nrow(expr_matrix)

  if (n_pairs == 0) {
    warning("No edges passed the half-threshold. Returning 0 scores.")
  } else {
    message(sprintf("Running %d permutations on %d informative edges...", n_permutations, n_pairs))
  }

  # Explicitly pass required data to workers to prevent Lexical Scoping crashes on Windows
  run_one_permutation <- function(i, cond, c_levels, em, pairs, n_p, n_g, mask, bicor_fn) {
    tryCatch({
      perm_condition <- sample(cond)
      x_c1 <- em[, perm_condition == c_levels[1]]
      x_c2 <- em[, perm_condition == c_levels[2]]

      p_diff_sq <- numeric(n_p)

      # Fast inner loop: Compute ONLY the specific edges that mattered biologically
      for (p in seq_len(n_p)) {
        g1 <- pairs[p, 1]
        g2 <- pairs[p, 2]

        r1 <- bicor_fn(x_c1[g1, ], x_c1[g2, ])
        r2 <- bicor_fn(x_c2[g1, ], x_c2[g2, ])

        if (is.na(r1) || is.na(r2)) {
          p_diff_sq[p] <- NA_real_
        } else {
          p_diff_sq[p] <- (r1 - r2)^2
        }
      }

      # Aggregate back to a temporary matrix to match the masking logic
      diff_mat <- matrix(NA_real_, nrow = n_g, ncol = n_g)
      diff_mat[pairs] <- p_diff_sq
      diff_mat[pairs[, c(2, 1)]] <- p_diff_sq # Make symmetric

      # Extract permuted DC scores for each gene safely
      p_dc <- vapply(seq_len(n_g), function(g) {
        idx <- mask[g, ]
        if (!any(idx)) return(0)

        vals <- diff_mat[g, idx]
        valid_vals <- vals[!is.na(vals)] # Filter NAs safely to avoid NaN warnings

        if (length(valid_vals) == 0) return(0)
        sqrt(mean(valid_vals))
      }, numeric(1))

      return(p_dc)

    }, error = function(e) {
      return(rep(NA_real_, n_g))
    })
  }

  # Use the Smart Dispatcher
  BPPARAM <- .get_bp_param(workers)

  # Ensure RNGseed is set if using multiple workers for reproducibility
  if(workers > 1){
    BiocParallel::bpRNGseed(BPPARAM) <- seed
  }

  perm_results <- BiocParallel::bplapply(
    seq_len(n_permutations),
    run_one_permutation,
    cond = condition,
    c_levels = conditions,
    em = expr_matrix,
    pairs = pair_indices,
    n_p = n_pairs,
    n_g = n_genes,
    mask = informative_mask,
    bicor_fn = bicor,
    BPPARAM = BPPARAM
  )

  # Column-bind results safely
  perm_dc_matrix <- do.call(cbind, perm_results)

  # Standard pseudo-count added to empirical p-value calculation
  p_values <- vapply(seq_along(dc_observed), function(i) {
    valid_perms <- perm_dc_matrix[i, !is.na(perm_dc_matrix[i, ])]
    if (length(valid_perms) == 0) return(1.0) # Fail safe
    (sum(valid_perms >= dc_observed[i]) + 1) / (length(valid_perms) + 1)
  }, numeric(1))

  # Boost statistical power by only adjusting FDR for tested (non-zero) genes
  fdr_values <- rep(1.0, length(p_values))
  tested_mask <- dc_observed > 0

  if (any(tested_mask)) {
    fdr_values[tested_mask] <- stats::p.adjust(p_values[tested_mask], method = "BH")
  }

  results_df <- data.frame(
    Gene = rownames(expr_matrix),
    DC_Score = dc_observed,
    P_value = p_values,
    FDR = fdr_values,
    stringsAsFactors = FALSE
  )

  if (species_mode) {
    results_df$species <- species
  }

  results_df <- results_df[order(-results_df$DC_Score), ]
  rownames(results_df) <- NULL

  structure(
    list(
      results = results_df,
      n_permutations = n_permutations,
      half_threshold = half_threshold,
      significance_level = significance_level,
      species_mode = species_mode,
      transform = transform,
      data = list(expr_matrix = expr_matrix, condition = condition, species = species)
    ),
    class = "bmht_result"
  )
}

# ===================================================================
# S3 Methods
# ===================================================================

#' Summary for bmht_result
#' @export
summary.bmht_result <- function(object, ...) {
  n_sig <- sum(object$results$FDR < object$significance_level, na.rm = TRUE)

  cat("--- BMHT Analysis Summary ---\n")
  cat(sprintf("Genes tested: %d\n", sum(object$results$DC_Score > 0)))
  cat(sprintf("Permutations performed: %d\n", object$n_permutations))
  cat(sprintf("Half-threshold: %.2f\n", object$half_threshold))
  cat(sprintf("Transform applied: %s\n", toupper(object$transform)))
  if (object$species_mode) cat("Species mode: Enabled\n")
  cat(sprintf("Significant genes (FDR < %.2f): %d\n", object$significance_level, n_sig))
  cat("----------------------------------\n")

  if (n_sig > 0) {
    cat("Top 10 genes by DC Score:\n")
    print(head(object$results, 10))
  } else {
    cat("No significant differentially co-expressed genes found.\n")
  }
  invisible(object)
}

#' Print method
#' @export
print.bmht_result <- function(x, ...) {
  summary(x, ...)
}
