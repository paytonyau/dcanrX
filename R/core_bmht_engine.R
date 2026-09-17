#' @title BMHT: Biweight Midcorrelation + Half-Thresholding Method
#' @description Ranks individual genes by their differential co-expression
#' connectivity (rewiring) between two conditions using bicor. The observed
#' and permutation-null statistics are computed in C++; see bmht_cpp.cpp.
#'
#' @importFrom stats p.adjust mad
#' @importFrom Rcpp sourceCpp
#' @importFrom progress progress_bar
#' @importFrom matrixStats rowMads
#' @importFrom parallel detectCores
#' @useDynLib bicorX, .registration = TRUE
#' @name bhmt
NULL

#' Run BMHT Analysis
#' @details
#' On the one real ground-truth benchmark available for this package
#' (dcanr's `sim102`; paper-fidelity review follow-up), ranking genes by
#' `DC_Score` (the raw effect size) noticeably outperformed ranking by
#' `P_value`/`FDR` for identifying true differential-coexpression genes
#' (AUPRC 0.8895 vs. 0.7931, gene-level ground truth) - the opposite of
#' what happened for ROS-DET's analytical significance method, where
#' p-value-based ranking was dramatically better than the raw score. The
#' likely reason: `P_value` normalizes each gene's `DC_Score` by that
#' gene's own permutation-null variance, which differs across genes
#' (fewer informative edges -> noisier null); on this dataset that
#' normalization reordered genes slightly worse than the raw score did.
#' This is not fully mechanistically isolated (would need a dedicated
#' follow-up to confirm precisely why), but the practical
#' recommendation is clear: prefer `DC_Score` for **ranking**/identifying
#' the best candidate genes, and use `FDR` only for a calibrated
#' significance **cutoff** decision (FDR's calibration itself is
#' unaffected by this - see `tests/testthat/test-bmht-calibration.R`).
#' Unlike ROS-DET (`ecor_pvalue()`), Zheng et al.'s BMHT paper does not
#' specify an analytical alternative to permutation testing, so this
#' ranking behavior is not expected to change without further work.
#' @keywords internal
.run_bmht_matrix <- function(expr_matrix, condition,
                             half_threshold = 0.4,
                             n_permutations = 1000,
                             significance_level = 0.05,
                             workers = 1,
                             transform = c("none", "clr", "rclr", "log1p"),
                             species = NULL,
                             seed = 42) {

  if (!is.null(seed)) set.seed(seed)
  condition <- .validate_input(expr_matrix, condition)
  transform <- match.arg(transform)

  # ===================================================================
  # 0. User Experience Guard Rails & Warning Messages
  # ===================================================================
  phys_cores <- parallel::detectCores(logical = FALSE)
  if (workers > phys_cores) {
    warning(sprintf("Thread Optimization: Requested workers (%d) exceeds physical hardware cores (%d).\n", workers, phys_cores),
            "Resetting workers to match physical cores to optimize cache layout and maximize speed.")
    workers <- max(1, phys_cores)
  }

  feature_mads <- matrixStats::rowMads(expr_matrix, na.rm = TRUE)
  keep <- feature_mads > 0
  if (any(!keep, na.rm = TRUE)) {
    warning("Zero Variance Alert: Some features have a Median Absolute Deviation (MAD) of exactly 0.\n",
            "BMHT has automatically filtered these silent features to prevent matrix arithmetic errors.")
    expr_matrix <- expr_matrix[keep, , drop = FALSE]
    if (!is.null(species)) species <- species[keep]
  }

  group_sizes <- table(condition)
  if (max(group_sizes) / min(group_sizes) > 4) {
    warning("Severe Cohort Imbalance: Group sizes are highly asymmetrical (",
            paste(group_sizes, collapse = " vs "), ").\n",
            "Permutation-derived empirical p-values may exhibit inflated Type I error rates.")
  }

  species_mode <- FALSE
  if (!is.null(species)) {
    if (length(species) != nrow(expr_matrix)) stop("Length of `species` must exactly match rows.")
    species_mode <- TRUE
    species <- as.character(species)
  }

  conditions <- levels(condition)
  if (sum(condition == conditions[1]) < 3 || sum(condition == conditions[2]) < 3) {
    stop("BMHT requires at least 3 samples per condition.")
  }

  expr_matrix <- .compositional_transform(expr_matrix, transform = transform, species = species)
  n_genes <- nrow(expr_matrix)

  # ===================================================================
  # 1. Base Topological Architecture with Scale-Free Normalization
  # ===================================================================
  message("Computing observed bicor matrices with Topological Scaling Factor (MAD Normalization)...")

  # Map observed calculations straight to native C++ layout for strict numeric parity
  cor_1 <- cpp_bmht_observed_bicor(expr_matrix[, condition == conditions[1], drop = FALSE])
  cor_2 <- cpp_bmht_observed_bicor(expr_matrix[, condition == conditions[2], drop = FALSE])
  cor_1[is.na(cor_1)] <- 0
  cor_2[is.na(cor_2)] <- 0

  # Compute observed static background macro-density stabilizers via MAD
  mad_1 <- stats::mad(cor_1[lower.tri(cor_1)], na.rm = TRUE)
  mad_2 <- stats::mad(cor_2[lower.tri(cor_2)], na.rm = TRUE)

  # Floor adjusted to 0.05 to maintain scale separation in noisy sets
  if (is.na(mad_1) || mad_1 < 0.05) mad_1 <- 1.0
  if (is.na(mad_2) || mad_2 < 0.05) mad_2 <- 1.0

  cor_1_scaled <- cor_1 / mad_1
  cor_2_scaled <- cor_2 / mad_2

  informative_mask <- (abs(cor_1) > half_threshold) | (abs(cor_2) > half_threshold)
  informative_mask[is.na(informative_mask)] <- FALSE

  diff_sq_obs <- (cor_1_scaled - cor_2_scaled)^2
  diff_sq_obs[!informative_mask] <- NA
  dc_observed <- sqrt(rowMeans(diff_sq_obs, na.rm = TRUE))
  dc_observed[is.na(dc_observed)] <- 0

  rm(cor_1, cor_2, cor_1_scaled, cor_2_scaled, diff_sq_obs)

  # ===================================================================
  # 2. R-Side Progress Chunking & C++ Execution
  # ===================================================================
  active_genes_idx <- which(rowSums(informative_mask, na.rm = TRUE) > 0)
  n_active <- length(active_genes_idx)

  # Accumulated incrementally, chunk by chunk - never materializes the full
  # n_genes x n_permutations result matrix. That matrix was previously kept
  # around for the whole run but only ever used for a single rowSums
  # reduction at the end (no other consumer) - at the very large
  # permutation counts now often needed for adequate power (Phase 5
  # finding), that was hundreds of MB to GB of R-side memory for something
  # that was purely transient. See tests/testthat/test-bmht-memory.R for
  # the exact-reproduction check this refactor was validated against.
  greater_eq_counts <- integer(n_genes)
  valid_perms_counts <- integer(n_genes)

  if (n_active == 0) {
    warning("No edges passed the half-threshold. Returning 0 scores.")
    # greater_eq_counts/valid_perms_counts both stay all-zero, giving
    # p_values = (0+1)/(0+1) = 1 for every gene - identical to the old
    # behavior (a fully-zero perm_dc_matrix compared against dc_observed,
    # which is also 0 for every gene in this case).
  } else {
    em_sub <- expr_matrix[active_genes_idx, , drop = FALSE]

    rm(informative_mask)
    gc(verbose = FALSE)

    cond_int <- as.integer(condition)

    # Fixed, bounded chunk size regardless of n_permutations - caps peak
    # per-chunk memory (both the shuffled-labels matrix and the C++ output
    # chunk) instead of letting it grow proportionally with however many
    # total permutations the user requests. Progress granularity is a UX
    # nicety, not a correctness requirement, so a fixed size is strictly
    # better: modest runs still get ~10 updates, large runs get more
    # frequent (not less frequent) feedback.
    chunk_size <- min(500, max(10, n_permutations))
    n_chunks <- ceiling(n_permutations / chunk_size)

    dc_observed_active <- dc_observed[active_genes_idx]

    message(sprintf("Deploying C++ OpenMP Engine across %d hardware threads...", workers))
    pb <- progress::progress_bar$new(
      format = "  Permutations Workflow [:bar] :percent | ETA: :eta",
      total = n_chunks, clear = FALSE, width = 65
    )

    for (c in 1:n_chunks) {
      start_p <- ((c - 1) * chunk_size) + 1
      end_p <- min(c * chunk_size, n_permutations)
      curr_b_count <- (end_p - start_p) + 1

      shuffled_chunk <- replicate(curr_b_count, sample(cond_int))

      cpp_chunk_out <- cpp_bmht_permutations(em_sub, shuffled_chunk,
                                             half_threshold, mad_1, mad_2, workers)

      chunk_ge <- rowSums(cpp_chunk_out >= dc_observed_active, na.rm = TRUE)
      chunk_valid <- rowSums(!is.na(cpp_chunk_out))

      greater_eq_counts[active_genes_idx] <- greater_eq_counts[active_genes_idx] + chunk_ge
      valid_perms_counts[active_genes_idx] <- valid_perms_counts[active_genes_idx] + chunk_valid

      pb$tick()
    }
  }

  # ===================================================================
  # 3. Vectorized Significance & FDR
  # ===================================================================
  # greater_eq_counts / valid_perms_counts already accumulated incrementally
  # above - no full perm_dc_matrix to reduce here anymore.
  p_values <- (greater_eq_counts + 1) / (valid_perms_counts + 1)

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
  if (species_mode) results_df$species <- species

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
    class = c("bmht_result", "list")
  )
}
