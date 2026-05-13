#' @title ROS-DET: Robust Switching Differential Co-expression
#' @description Identifies gene pairs that switch correlation sign between conditions.
#' Optimized with precomputed matrices, WGCNA scaling, and cross-species mode.
#'
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam SnowParam bpparam
#' @importFrom stats optimize pchisq p.adjust
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_bw
#' @name rosdet
NULL

# ===================================================================
# Main ROS-DET Function (with Cross-Species Support & Optimization)
# ===================================================================

#' Run ROS-DET Analysis
#'
#' @param expr_matrix Numeric matrix (genes × samples)
#' @param condition Factor/character with exactly two levels
#' @param min_delta Minimum absolute difference in correlation to consider testing
#' @param significance_level FDR threshold (default 0.05)
#' @param workers Number of parallel workers
#' @param transform Compositional transformation ("none", "clr", "log1p")
#' @param species Optional character vector indicating species for each gene
#' @param pair_type_filter Filter results to "intra", "inter", or "all"
#' @param weight_mode "wcor" (original variance penalty) or "unweighted" (for single-cell)
#' @param seed Integer seed for parallel reproducibility
#' @return `rosdet_result` object
#' @export
run_rosdet <- function(expr_matrix, condition,
                       min_delta = 0.4,
                       significance_level = 0.05,
                       workers = 1,
                       transform = c("none", "clr", "log1p"),
                       species = NULL,
                       pair_type_filter = c("all", "intra", "inter"),
                       weight_mode = c("wcor", "unweighted"),
                       seed = 42) {

  if (!is.null(seed)) set.seed(seed)

  condition <- .validate_input(expr_matrix, condition)
  transform <- match.arg(transform)
  pair_type_filter <- match.arg(pair_type_filter)
  weight_mode <- match.arg(weight_mode)

  cond_levels <- levels(condition)

  # Strict Sample Size Guardrail
  n1 <- sum(condition == cond_levels[1])
  n2 <- sum(condition == cond_levels[2])
  if (n1 < 3 || n2 < 3) {
    stop("ROS-DET requires at least 3 samples per condition to compute valid variances and Chi-Square statistics.")
  }

  # Strict Species Validation
  use_species_mode <- FALSE
  if (!is.null(species)) {
    if (length(species) != nrow(expr_matrix)) {
      stop("Length of `species` must exactly match the number of rows in `expr_matrix`.")
    }
    use_species_mode <- TRUE
    species <- as.character(species)
    message(sprintf("Cross-species mode enabled. Filtering mode: %s", toupper(pair_type_filter)))
  } else if (pair_type_filter != "all") {
    warning("pair_type_filter requires 'species' to be provided. Defaulting to 'all'.")
    pair_type_filter <- "all"
  }

  # Apply compositional transformation
  expr_matrix <- .compositional_transform(expr_matrix, transform, species)
  gene_names <- rownames(expr_matrix)
  n_genes <- nrow(expr_matrix)

  expr_1 <- expr_matrix[, condition == cond_levels[1]]
  expr_2 <- expr_matrix[, condition == cond_levels[2]]

  # ===================================================================
  # Step 1: O(N) Precomputations (Massive Speedup)
  # ===================================================================
  message("Precomputing biweight midvariances...")
  var_1 <- apply(expr_1, 1, biweight_midvariance)
  var_2 <- apply(expr_2, 1, biweight_midvariance)

  message("Precomputing observed bicor matrices...")
  cor_1 <- compute_cor_matrix(expr_1, workers = workers)
  cor_2 <- compute_cor_matrix(expr_2, workers = workers)

  # Prevent ECOR division-by-zero crashes by clamping perfect correlations
  cor_1[cor_1 > 0.99] <- 0.99; cor_1[cor_1 < -0.99] <- -0.99
  cor_2[cor_2 > 0.99] <- 0.99; cor_2[cor_2 < -0.99] <- -0.99

  # ===================================================================
  # Step 2: Pre-filter Target Pairs (Independent Filtering)
  # ===================================================================
  message(sprintf("Filtering candidate pairs (min_delta >= %.2f)...", min_delta))

  delta_mat <- abs(cor_1 - cor_2)

  # Safe NA handling: convert missing edges to 0 so they fail the min_delta check safely
  delta_mat[is.na(delta_mat)] <- 0
  delta_mat[lower.tri(delta_mat, diag = TRUE)] <- 0 # Only check upper triangle

  candidate_indices <- which(delta_mat >= min_delta, arr.ind = TRUE)
  n_candidates <- nrow(candidate_indices)

  if (n_candidates == 0) {
    warning("No gene pairs met the min_delta threshold. Try lowering min_delta.")
    return(invisible(NULL))
  }

  message(sprintf("Evaluating %s high-potential pairs with ECOR optimization...", format(n_candidates, big.mark=",")))

  # ===================================================================
  # Step 3: ECOR Parallel Execution (Windows-Safe)
  # ===================================================================

  # Cross-platform dispatcher respecting RNGseed
  if (workers <= 1) {
    BPPARAM <- BiocParallel::SerialParam(progressbar = TRUE, RNGseed = seed)
  } else if (.Platform$OS.type == "windows") {
    BPPARAM <- BiocParallel::SnowParam(workers = workers, progressbar = TRUE, RNGseed = seed)
  } else {
    BPPARAM <- BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE, RNGseed = seed)
  }

  # explicitly pass ALL required variables to prevent scoping errors in Snow workers
  eval_one_pair <- function(i, cand_idx, c1, c2, v1, v2, w_mode, num1, num2) {
    g1_idx <- cand_idx[i, 1]
    g2_idx <- cand_idx[i, 2]

    r1 <- c1[g1_idx, g2_idx]
    r2 <- c2[g1_idx, g2_idx]

    # WCOR Calculation with Single-Cell Override
    if (w_mode == "unweighted") {
      c_weight <- 1.0
    } else {
      max_var1 <- max(v1[g1_idx], v1[g2_idx], na.rm = TRUE)
      max_var2 <- max(v2[g1_idx], v2[g2_idx], na.rm = TRUE)
      c_weight <- if (max_var1 == 0 || max_var2 == 0) 0 else min(max_var1 / max_var2, max_var2 / max_var1)
    }

    score <- c_weight * abs(r1 - r2)

    # ECOR Objective Function
    objective_fn <- function(rho_hat) {
      term1 <- num1 * (r1 - rho_hat) / (1 - rho_hat * r1)
      term2 <- num2 * (r2 - rho_hat) / (1 - rho_hat * r2)
      (term1 + term2)^2
    }

    # Safe optimization
    opt_res <- tryCatch({
      stats::optimize(objective_fn, interval = c(-0.99, 0.99))$minimum
    }, error = function(e) NA_real_)

    if (is.na(opt_res)) return(c(score, NA_real_, r1, r2, c_weight))

    rho_hat <- opt_res

    # Machine-Precision Safety (Epsilon) to prevent log(0) and NaN crashes
    eps <- 1e-12
    num1_safe <- max((1 - r1 * rho_hat)^2, eps)
    den1_safe <- max((1 - r1^2) * (1 - rho_hat^2), eps)

    num2_safe <- max((1 - r2 * rho_hat)^2, eps)
    den2_safe <- max((1 - r2^2) * (1 - rho_hat^2), eps)

    term1_log <- num1 * log(num1_safe / den1_safe)
    term2_log <- num2 * log(num2_safe / den2_safe)

    T_statistic <- sum(c(term1_log, term2_log), na.rm = TRUE)
    p_value <- stats::pchisq(T_statistic, df = 1, lower.tail = FALSE)

    return(c(score, p_value, r1, r2, c_weight))
  }

  # Safely pass all variables directly into the bplapply environment
  res_list <- BiocParallel::bplapply(
    seq_len(n_candidates),
    eval_one_pair,
    cand_idx = candidate_indices,
    c1 = cor_1,
    c2 = cor_2,
    v1 = var_1,
    v2 = var_2,
    w_mode = weight_mode,
    num1 = n1,
    num2 = n2,
    BPPARAM = BPPARAM
  )

  res_mat <- do.call(rbind, res_list)

  # ===================================================================
  # Step 4: Formatting and Independent Filtering (FDR)
  # ===================================================================
  all_results <- data.frame(
    Gene1 = gene_names[candidate_indices[, 1]],
    Gene2 = gene_names[candidate_indices[, 2]],
    Score = res_mat[, 1],
    P_value = res_mat[, 2],
    r1 = res_mat[, 3],
    r2 = res_mat[, 4],
    c_weight = res_mat[, 5],
    stringsAsFactors = FALSE
  )

  # Clean failed optimizations safely
  all_results <- all_results[!is.na(all_results$P_value), ]

  # Apply FDR specifically to the filtered candidate list to maximize statistical power
  all_results$FDR <- stats::p.adjust(all_results$P_value, method = "BH")

  # Species pair type attribution
  if (use_species_mode) {
    sp1 <- species[candidate_indices[, 1]]
    sp2 <- species[candidate_indices[, 2]]
    all_results$pair_type <- ifelse(sp1 == sp2, "intra", "inter")

    if (pair_type_filter != "all") {
      all_results <- all_results[all_results$pair_type == pair_type_filter, ]
    }
  }

  # Rank strictly by ECOR Score
  ranked_results <- all_results[order(-all_results$Score), , drop = FALSE]
  rownames(ranked_results) <- NULL

  structure(
    list(
      results = ranked_results,
      n_tested = n_candidates,
      significance_level = significance_level,
      species_mode = use_species_mode,
      pair_type_filter = pair_type_filter,
      weight_mode = weight_mode,
      data = list(expr_matrix = expr_matrix, condition = condition, species = species)
    ),
    class = "rosdet_result"
  )
}

# ===================================================================
# S3 Methods
# ===================================================================

#' Summary for rosdet_result
#' @export
summary.rosdet_result <- function(object, ...) {
  n_sig <- sum(object$results$FDR < object$significance_level, na.rm = TRUE)

  cat("--- ROS-DET Analysis Summary ---\n")
  cat(sprintf("Gene pairs evaluated (passed min_delta): %s\n", format(object$n_tested, big.mark = ",")))
  cat(sprintf("Weight mode: %s\n", toupper(object$weight_mode)))

  if (object$species_mode) {
    cat(sprintf("Mode: Cross-species enabled | Filter: %s\n", toupper(object$pair_type_filter)))
  }

  cat(sprintf("Significant switching pairs (FDR < %.3f): %s\n",
              object$significance_level, format(n_sig, big.mark = ",")))
  cat("----------------------------------\n")

  if (n_sig > 0) {
    cat("Top 6 results:\n")
    print(head(object$results, 6))
  } else {
    cat("No significant switching gene pairs found.\n")
  }
  invisible(object)
}

#' Print method
#' @export
print.rosdet_result <- function(x, ...) summary(x, ...)
