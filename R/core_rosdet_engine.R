#' @title ROS-DET: Universal Re-wired Differential Co-expression Topology
#' @description Identifies gene pairs undergoing rewiring (both sign-flipping and decoupling)
#' between conditions, using a two-gate design: Gate 1 screens candidate pairs on
#' |delta correlation|, Gate 2 runs a permutation test on the surviving candidates.
#'
#' @importFrom stats p.adjust
#' @importFrom Rcpp sourceCpp
#' @importFrom progress progress_bar
#' @importFrom matrixStats rowMads
#' @importFrom parallel detectCores
#' @useDynLib bicorX, .registration = TRUE
#' @name rosdet
NULL

#' Run ROS-DET Analysis
#'
#' @param weight_mode `"unweighted"` (default as of the Phase 5 benchmark -
#'   see below) uses the raw absolute correlation difference for
#'   `Distance_Score`. `"wcor"` instead weights it by a variance-ratio term
#'   derived from the *observed* data (see Details).
#'
#'   IMPORTANT: this weight is a fixed per-pair constant applied identically
#'   to the observed statistic and to every permutation's null statistic, so
#'   it cancels exactly in the `null_dist >= obs_dist` comparison used to
#'   build p-values. `weight_mode` therefore changes `Distance_Score` (and
#'   the *ranking* used for display/sorting, e.g. "top hits by effect size")
#'   but has **no effect whatsoever on `P_value` or `FDR`** - both are
#'   identical between the two modes for the same seed and permutation
#'   count. Verified empirically in
#'   \code{tests/testthat/test-rosdet-weight-mode.R}.
#'
#'   The default was changed from `"wcor"` to `"unweighted"` after a
#'   benchmark (`dev/phase5_compositional_benchmark.R`) found `"wcor"` can
#'   substantially harm - not just shift - `Distance_Score`-based ranking:
#'   on a synthetic dataset with a genuine 10-feature differential module,
#'   `"wcor"` produced a *lower* mean score for true-signal pairs than for
#'   background pairs (0.144 vs. 0.164), inverting the useful signal
#'   (AUPRC 0.034, at random baseline), while `"unweighted"` on the exact
#'   same data correctly separated them (0.841 vs. 0.335; AUPRC 0.359).
#'
#'   Kayano et al. 2011 validate WCOR specifically against *range bias*
#'   (their Fig. 4b/c: one condition's expression range much smaller than
#'   the other's), not the compositional-closure confound above - so this
#'   was retested directly on a range-bias scenario matching their design
#'   (paper-fidelity review, Phase C). Three findings, carefully separating
#'   variance/range differences from true signal (an easy confound to
#'   introduce by accident - an initial attempt at this retest conflated
#'   them and was discarded before drawing any conclusion from it):
#'   \itemize{
#'     \item `"wcor"` correctly suppresses a *pure* range-bias artifact
#'       (large variance difference between conditions, no real
#'       correlation change) to near-zero score, while `"unweighted"`
#'       does not protect against this at all - matching the paper's
#'       design intent exactly.
#'     \item On *clean* true signal with variance genuinely equalized
#'       between conditions, `"wcor"` and `"unweighted"` perform
#'       comparably (0.843 vs. 0.893 mean score for true vs. background
#'       pairs in the retest) - `"wcor"` is not broadly worse at
#'       detecting real signal in general.
#'     \item The real, narrower failure mode: when true signal is
#'       naturally *accompanied by* a variance/range difference between
#'       conditions (as in the original compositional-confound benchmark,
#'       and plausibly common in real biology where regulatory rewiring
#'       often co-occurs with variance shifts), `"wcor"` suppresses that
#'       signal along with the variance difference, rather than
#'       distinguishing "real signal with incidental variance shift" from
#'       "pure variance artifact with no real signal".
#'   }
#'   `"unweighted"` remains the default given this - it's the safer
#'   general-purpose choice, since signal/variance-shift coupling seems
#'   more common in practice than the narrow scenario `"wcor"` protects
#'   against. `"wcor"` is kept available and is a reasonable choice
#'   specifically when technical/range artifacts (not biological variance
#'   shifts) are the main concern.
#' @param significance_method `"permutation"` (default) uses the exact,
#'   distribution-free permutation test (shuffle condition labels,
#'   recompute bicor, count exceedances) that this engine has always
#'   used. `"analytical"` uses the closed-form-ish chi-squared(1) test
#'   from Kayano et al. 2011 ("ECOR", the significance step of the
#'   original ROS-DET method this engine is based on) - see
#'   `.ecor_pvalue()` for the statistic and its calibration properties.
#'
#'   The tradeoff is speed vs. assumptions, not "better" vs. "worse":
#'   `"analytical"` needs no resampling at all (a fixed, small cost per
#'   pair regardless of desired p-value resolution - practical at the
#'   genome-wide scale the original paper validated it at, ~2.6e10 pairs
#'   across 46 real datasets) but assumes approximate bivariate normality
#'   and is only well-calibrated for **N >= ~50 samples per condition**;
#'   below that it becomes progressively anti-conservative (a warning is
#'   issued below N=30). `"permutation"` makes no distributional
#'   assumption and is exact at any sample size, but its resolution is
#'   capped by `n_permutations` (the source of this package's documented
#'   permutation-floor limitation - see NEWS.md) and its cost scales with
#'   however many permutations that resolution requires.
#' @keywords internal
.run_rosdet_matrix <- function(expr_matrix, condition,
                               min_delta = 0.4,
                               n_permutations = 1000,
                               significance_level = 0.05,
                               workers = 1,
                               transform = c("none", "clr", "rclr", "log1p"),
                               species = NULL,
                               pair_type_filter = c("all", "intra", "inter"),
                               weight_mode = c("unweighted", "wcor"),
                               significance_method = c("permutation", "analytical"),
                               seed = 42,
                               n_bootstraps = NULL) {

  # `n_bootstraps` was the original argument name, but this is a permutation
  # test (group labels are shuffled without replacement via sample(), not
  # resampled with replacement) - kept as a deprecated alias so existing
  # calls keep working.
  if (!is.null(n_bootstraps)) {
    warning("`n_bootstraps` is deprecated; use `n_permutations` instead. ",
            "ROS-DET's null is built by permuting condition labels, not by ",
            "resampling with replacement, so 'permutations' is the accurate name.",
            call. = FALSE)
    n_permutations <- n_bootstraps
  }

  if (!is.null(seed)) set.seed(seed)

  condition <- .validate_input(expr_matrix, condition)
  transform <- match.arg(transform)
  pair_type_filter <- match.arg(pair_type_filter)
  weight_mode <- match.arg(weight_mode)
  significance_method <- match.arg(significance_method)

  cond_levels <- levels(condition)

  # ===================================================================
  # 0. Package-Wide UX Guard Rails & Warning Messages
  # ===================================================================
  phys_cores <- parallel::detectCores(logical = FALSE)
  if (workers > phys_cores) {
    warning(sprintf("Thread Optimization [ROS-DET]: Requested workers (%d) exceeds physical hardware cores (%d).\n", workers, phys_cores),
            "Resetting workers to match physical cores to optimize cache layout and maximize speed.")
    workers <- max(1, phys_cores)
  }

  feature_mads <- matrixStats::rowMads(expr_matrix, na.rm = TRUE)
  keep <- feature_mads > 0
  if (any(!keep, na.rm = TRUE)) {
    warning("Zero Variance Alert [ROS-DET]: Features with a MAD of 0 automatically trimmed\n",
            "to secure matrix algebra stability against division-by-zero errors.")
    expr_matrix <- expr_matrix[keep, , drop = FALSE]
    if (!is.null(species)) species <- species[keep]
  }

  group_sizes <- table(condition)
  if (max(group_sizes) / min(group_sizes) > 4) {
    warning("Severe Cohort Imbalance [ROS-DET]: Highly asymmetrical group layouts (",
            paste(group_sizes, collapse = " vs "), ") detected.\n",
            "Resampling null baselines may exhibit conservative empirical p-value shifts.")
  }

  n1 <- sum(condition == cond_levels[1])
  n2 <- sum(condition == cond_levels[2])
  if (n1 < 5 || n2 < 5) stop("ROS-DET bootstrapping requires at least 5 samples per condition.")

  if (significance_method == "analytical" && (n1 < 30 || n2 < 30)) {
    warning("significance_method = \"analytical\" (ECOR) is only well-calibrated for ",
            "roughly N >= 50 samples per condition; below N ~= 30 it becomes ",
            "meaningfully anti-conservative (see .ecor_pvalue() docs for the ",
            "calibration simulation this is based on). With ", n1, "/", n2,
            " samples per condition here, \"permutation\" is the safer choice.",
            call. = FALSE)
  }

  use_species_mode <- FALSE
  if (!is.null(species)) {
    if (length(species) != nrow(expr_matrix)) stop("Length of `species` must match rows.")
    use_species_mode <- TRUE # FIX: Typo assignment correction to enable parsing validation
    species <- as.character(species)
  }

  expr_matrix <- .compositional_transform(expr_matrix, transform, species)
  gene_names <- rownames(expr_matrix)

  expr_1 <- expr_matrix[, condition == cond_levels[1], drop = FALSE]
  expr_2 <- expr_matrix[, condition == cond_levels[2], drop = FALSE]

  message("Initializing ROS-DET Absolute Distance Workflow...")
  message("Precomputing biweight midvariances...")
  var_1 <- apply(expr_1, 1, biweight_midvariance)
  var_2 <- apply(expr_2, 1, biweight_midvariance)

  message("Precomputing observed cache-aligned bicor matrices...")
  cor_1 <- cpp_rosdet_observed_bicor(expr_1)
  cor_2 <- cpp_rosdet_observed_bicor(expr_2)

  # ===================================================================
  # GATE 1: Fast Absolute Distance Screening Pass
  # ===================================================================
  delta_mat <- abs(cor_1 - cor_2)
  delta_mat[is.na(delta_mat)] <- 0
  delta_mat[lower.tri(delta_mat, diag = TRUE)] <- 0

  candidate_indices <- which(delta_mat >= min_delta, arr.ind = TRUE)

  if (nrow(candidate_indices) > 0 && use_species_mode && pair_type_filter != "all") {
    sp1 <- species[candidate_indices[, 1]]
    sp2 <- species[candidate_indices[, 2]]
    keep_mask <- if(pair_type_filter == "inter") sp1 != sp2 else sp1 == sp2
    candidate_indices <- candidate_indices[keep_mask, , drop = FALSE]
  }

  n_candidates <- nrow(candidate_indices)
  if (n_candidates == 0) {
    warning("No gene pairs cleared your minimum delta screening pass (Gate 1). Summary aborted.")
    return(invisible(NULL))
  }

  # ===================================================================
  # 2. Candidate Linearization & RAM Compaction
  # ===================================================================
  g1_idx <- candidate_indices[, 1]
  g2_idx <- candidate_indices[, 2]

  r1_vec <- cor_1[candidate_indices]
  r2_vec <- cor_2[candidate_indices]

  rm(cor_1, cor_2, delta_mat)
  gc(verbose = FALSE)

  if (weight_mode == "unweighted") {
    c_weight_vec <- rep(1.0, n_candidates)
  } else {
    max_v1 <- pmax(var_1[g1_idx], var_1[g2_idx], na.rm = TRUE)
    max_v2 <- pmax(var_2[g1_idx], var_2[g2_idx], na.rm = TRUE)
    c_weight_vec <- numeric(n_candidates)
    valid_mask <- max_v1 > 0 & max_v2 > 0
    c_weight_vec[valid_mask] <- pmin(max_v1[valid_mask] / (max_v2[valid_mask] + 1e-12),
                                     max_v2[valid_mask] / (max_v1[valid_mask] + 1e-12))
  }

  # NOTE: c_weight_vec is a fixed per-pair constant, applied identically to
  # obs_distance_vec here and to every permutation's null distance inside
  # cpp_rosdet_bootstrap_batch(). It therefore cancels exactly in the
  # null_dist >= obs_dist comparison and has zero effect on P_value/FDR -
  # only on Distance_Score (i.e. ranking/display). See weight_mode docs above.
  obs_distance_vec <- c_weight_vec * abs(r1_vec - r2_vec)

  # ===================================================================
  # GATE 2: Significance Testing (permutation or analytical)
  # ===================================================================
  if (significance_method == "analytical") {
    message(sprintf("Gate 2: Executing analytical ECOR test (Kayano et al. 2011) on %s high-potential pairs...",
                    format(n_candidates, big.mark = ",")))
    p_values <- .ecor_pvalue(r1_vec, r2_vec, n1, n2)

  } else {
    message(sprintf("Gate 2: Executing %s-fold Parallel Permutation Testing on %s high-potential pairs...",
                    n_permutations, format(n_candidates, big.mark=",")))

    cond_int <- as.integer(condition)
    # Fixed, bounded chunk size regardless of n_permutations - same fix as
    # BMHT's engine (see core_bmht_engine.R): caps the shuffled-labels chunk
    # size instead of letting it grow proportionally with n_permutations,
    # and keeps progress-bar updates frequent at very large permutation
    # counts instead of collapsing to a handful of huge chunks.
    chunk_size <- min(500, max(10, n_permutations))
    n_chunks <- ceiling(n_permutations / chunk_size)

    greater_eq_counts <- integer(n_candidates)

    pb <- progress::progress_bar$new(
      format = "  Permutation Progress [:bar] :percent | ETA: :eta",
      total = n_chunks, clear = FALSE, width = 65
    )

    unique_active_genes <- unique(c(g1_idx, g2_idx))
    expr_sub <- expr_matrix[unique_active_genes, , drop = FALSE]

    g1_mapped <- match(g1_idx, unique_active_genes) - 1
    g2_mapped <- match(g2_idx, unique_active_genes) - 1

    for (c in 1:n_chunks) {
      start_b <- ((c - 1) * chunk_size) + 1
      end_b <- min(c * chunk_size, n_permutations)
      curr_b_count <- (end_b - start_b) + 1

      shuffled_labels_chunk <- replicate(curr_b_count, sample(cond_int))

      # Trigger zero-contention parallel batch engine
      cpp_counts_chunk <- cpp_rosdet_bootstrap_batch(
        expr_sub, g1_mapped, g2_mapped, obs_distance_vec,
        c_weight_vec, shuffled_labels_chunk, workers
      )

      greater_eq_counts <- greater_eq_counts + cpp_counts_chunk
      pb$tick()
    }

    p_values <- (greater_eq_counts + 1) / (n_permutations + 1)
  }

  # ===================================================================
  # 5. Final Assembly
  # ===================================================================
  all_results <- data.frame(
    Gene1 = gene_names[g1_idx],
    Gene2 = gene_names[g2_idx],
    Distance_Score = obs_distance_vec,
    P_value = p_values,
    r1 = r1_vec,
    r2 = r2_vec,
    c_weight = c_weight_vec,
    stringsAsFactors = FALSE
  )

  all_results <- all_results[!is.na(all_results$P_value), ]
  all_results$FDR <- stats::p.adjust(all_results$P_value, method = "BH")

  if (use_species_mode) {
    all_results$pair_type <- ifelse(species[g1_idx] == species[g2_idx], "intra", "inter")
  }

  ranked_results <- all_results[order(-all_results$Distance_Score), , drop = FALSE]
  rownames(ranked_results) <- NULL

  structure(list(
    results = ranked_results,
    n_tested = n_candidates,
    n_permutations = if (significance_method == "permutation") n_permutations else NA_integer_,
    significance_method = significance_method,
    significance_level = significance_level,
    species_mode = use_species_mode,
    pair_type_filter = pair_type_filter,
    weight_mode = weight_mode,
    data = list(expr_matrix = expr_matrix, condition = condition, species = species)
  ), class = c("rosdet_result", "list"))
}
