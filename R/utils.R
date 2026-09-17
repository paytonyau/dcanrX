#' @title Utility Functions for bicorX
#' @description Core utilities: biweight midcorrelation (bicor),
#' biweight midvariance, and input validation.
#'
#' @importFrom stats median mad complete.cases cor var sd
#' @name utils
NULL

# ===================================================================
# Input Validation (shared across all methods)
# ===================================================================

#' @keywords internal
.validate_input <- function(expr_matrix, condition) {
  if (!is.matrix(expr_matrix) || !is.numeric(expr_matrix)) {
    stop("`expr_matrix` must be a numeric matrix. If using sparse matrices, please convert via as.matrix().")
  }
  if (is.null(rownames(expr_matrix)) || any(duplicated(rownames(expr_matrix)))) {
    stop("`expr_matrix` must have unique gene names as rownames.")
  }
  if (any(is.na(expr_matrix)) || any(is.infinite(expr_matrix))) {
    stop("`expr_matrix` contains NA, NaN, or Inf values. Please impute or filter these before running.")
  }

  row_vars <- apply(expr_matrix, 1, stats::var)
  n_zero_var <- sum(row_vars == 0, na.rm = TRUE)
  if (n_zero_var > 0) {
    warning(sprintf("Zero Variance Warning: Detected %d zero-variance rows.\n", n_zero_var),
            "Filtering these inert features upstream is highly recommended to protect mathematical stability.")
  }

  if (length(condition) != ncol(expr_matrix)) {
    stop(sprintf("Dimension mismatch: Provided condition vector length (%d) does not match expression matrix columns (%d).",
                 length(condition), ncol(expr_matrix)))
  }
  if (any(is.na(condition))) {
    stop("`condition` vector contains NAs. Please filter out these samples/cells before running.")
  }

  condition <- factor(condition)
  if (length(levels(condition)) != 2) {
    stop("`condition` must contain exactly two unique levels.")
  }

  return(condition)
}

# ===================================================================
# Biweight Midcorrelation (bicor)
# ===================================================================

#' Biweight Midcorrelation
#'
#' Computes the biweight midcorrelation between two numeric vectors, a
#' robust alternative to Pearson correlation that down-weights outlying
#' observations. Falls back to Pearson correlation when the median
#' absolute deviation of either input is zero.
#'
#' @param x,y Numeric vectors of equal length.
#' @return A single numeric value in `[-1, 1]`.
#' @export
#' @export
bicor <- function(x, y) {
  complete_obs <- stats::complete.cases(x, y)
  x <- x[complete_obs]
  y <- y[complete_obs]

  if (length(x) < 2) return(NA_real_)

  median_x <- stats::median(x)
  median_y <- stats::median(y)
  # Raw (unscaled) MAD - NOT stats::mad(), which multiplies by 1.4826 for
  # normal-distribution consistency. The bicor "9" tuning constant already
  # assumes raw MAD (Langfelder & Horvath 2012; WGCNA::bicor uses raw MAD
  # internally). Verified to machine precision against WGCNA::bicor in
  # tests/testthat/test-bicor-parity.R.
  mad_x <- stats::median(abs(x - median_x))
  mad_y <- stats::median(abs(y - median_y))

  if (mad_x == 0 || mad_y == 0) {
    if (stats::sd(x) == 0 || stats::sd(y) == 0) return(0)
    return(stats::cor(x, y, method = "pearson"))
  }

  u <- (x - median_x) / (9 * mad_x)
  v <- (y - median_y) / (9 * mad_y)

  w_x <- ifelse(abs(u) < 1, (1 - u^2)^2, 0)
  w_y <- ifelse(abs(v) < 1, (1 - v^2)^2, 0)

  numerator <- sum((x - median_x) * w_x * (y - median_y) * w_y)
  denom_x_sq <- sum(((x - median_x) * w_x)^2)
  denom_y_sq <- sum(((y - median_y) * w_y)^2)
  denominator <- sqrt(denom_x_sq) * sqrt(denom_y_sq)

  if (denominator == 0) return(0)
  numerator / denominator
}

# ===================================================================
# Biweight Midvariance (internal)
# ===================================================================

#' @keywords internal
biweight_midvariance <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) return(0)

  median_x <- stats::median(x)
  # Raw MAD, not stats::mad() - see bicor() above for why.
  mad_x <- stats::median(abs(x - median_x))

  if (mad_x == 0) return(stats::var(x))

  u <- (x - median_x) / (9 * mad_x)
  valid_u <- abs(u) < 1
  u_valid <- u[valid_u]
  x_valid <- x[valid_u]

  if (length(u_valid) == 0) return(0)

  numerator <- n * sum(((x_valid - median_x)^2) * ((1 - u_valid^2)^4))
  denom_term <- sum((1 - u_valid^2) * (1 - 5 * u_valid^2))

  if (denom_term == 0) return(0)
  numerator / (denom_term^2)
}

# ===================================================================
# Compositional Data Transformation (clr & rclr)
# ===================================================================

#' Apply a compositional transform, optionally per-group
#'
#' @param mat A features x samples matrix.
#' @param transform One of `"none"`, `"clr"`, `"rclr"`, `"log1p"`.
#' @param species Optional grouping vector, one entry per row of `mat`
#'   (aliased as `layer` at the user-facing `run_bicordcea()`/`run_bmht()`/
#'   `run_rosdet()`/`run_bmkc()` level - see there).
#'
#' @details
#' This is the mechanism behind bicorX's cross-domain support (compositional
#' microbiome data, or multiple omics layers combined into one matrix), and is
#' otherwise easy to miss since it's driven entirely by whether `species` is
#' supplied.
#'
#' **Without `species`** (the ordinary case): CLR is computed once, globally,
#' using every row's geometric mean across all samples.
#'
#' **With `species`** (`transform = "clr"` only - `rclr`/`log1p` are always
#' global): rows are split into groups by `species`, and CLR is computed
#' *separately within each group* - each group's log-ratio is centered on
#' its own group mean, not a mean shared across the whole matrix. This
#' matters whenever the matrix mixes domains that shouldn't be normalized
#' against each other:
#'
#' - **Microbiome mode**: `species` holds taxonomic labels (e.g. host vs.
#'   a set of microbial taxa). Host transcript counts and microbial taxon
#'   abundances live on very different scales; a single global CLR would
#'   let the (usually much larger) host compartment dominate the
#'   normalization and distort the microbial compositions.
#' - **Multi-omics mode**: `species`/`layer` holds an omics-layer tag (e.g.
#'   `"RNA"` vs. `"ATAC"`). RNA counts and ATAC peak accessibility have
#'   different sparsity and dynamic range; per-layer CLR keeps one layer's
#'   normalization from bleeding into the other's.
#'
#' Both uses are the identical mechanism - a named partition of the
#' matrix's rows - which is why `run_multiomic()` and the microbiome
#' workflow both just populate this same argument (via `layer`/`species`)
#' rather than needing separate code paths.
#'
#' @keywords internal
.compositional_transform <- function(mat, transform = c("none", "clr", "rclr", "log1p"), species = NULL) {
  transform <- match.arg(transform)

  if (transform == "none") return(mat)

  if (transform == "log1p") {
    if (any(mat < 0, na.rm = TRUE)) stop("Transformation Error: log1p requires non-negative count data.")
    return(log1p(mat))
  }

  if (transform == "rclr") {
    if (any(mat < 0, na.rm = TRUE)) {
      stop("Transformation Error: Robust CLR (rCLR) requires non-negative count matrices.\n",
           "Ensure data does not contain negative values from alternative scaling tools.")
    }

    mat_log <- log(mat)
    mat_log[is.infinite(mat_log)] <- NA

    log_geom_means <- colMeans(mat_log, na.rm = TRUE)
    rclr_mat <- sweep(mat_log, 2, log_geom_means, "-")
    rclr_mat[is.na(rclr_mat)] <- 0
    return(rclr_mat)
  }

  if (transform == "clr") {
    if (any(mat < 0, na.rm = TRUE)) {
      stop("Transformation Error: Standard CLR requires non-negative count matrices.\n",
           "If your matrix features biological log-folds or pre-scaled values, set transform = 'none'.")
    }

    min_nonzero <- min(mat[mat > 0], na.rm = TRUE)
    pseudo_count <- ifelse(is.finite(min_nonzero), min_nonzero / 2, 1e-8)
    mat_pseudo <- mat + pseudo_count

    if (!is.null(species) && length(unique(species)) > 1) {
      species <- as.character(species)
      unique_groups <- unique(species)
      clr_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
      rownames(clr_mat) <- rownames(mat); colnames(clr_mat) <- colnames(mat)

      for (gp in unique_groups) {
        idx <- which(species == gp)
        gp_mat_t <- t(mat_pseudo[idx, , drop = FALSE])
        log_gp_mat <- log(gp_mat_t)
        gp_clr_t <- log_gp_mat - base::rowMeans(log_gp_mat)
        clr_mat[idx, ] <- t(gp_clr_t)
      }
      return(clr_mat)
    } else {
      mat_t <- t(mat_pseudo)
      log_mat <- log(mat_t)
      clr_mat_t <- log_mat - base::rowMeans(log_mat)
      return(t(clr_mat_t))
    }
  }
}

# ===================================================================
# Shared plotting helper: observed bicor for ad hoc visualization use
# ===================================================================

#' Compute an observed bicor matrix for plotting-time recomputation
#'
#' Every plotting function across core_bmht_plots.R, core_bmkc_plots.R, and
#' core_rosdet_plots.R calls `cpp_calc_observed_bicor()` to recompute a
#' correlation matrix for a gene subset at plot time (e.g. for
#' `plot_differential_network()`). That function was never defined or
#' exported anywhere in the package - every one of those ~25 call sites was
#' broken (`could not find function "cpp_calc_observed_bicor"`) until this
#' alias was added.
#'
#' The three engine-specific kernels (`cpp_bmht_observed_bicor`,
#' `cpp_rosdet_observed_bicor`, `cpp_fast_bicor_matrix`) are numerically
#' identical (verified to machine precision in
#' tests/testthat/test-bicor-parity.R), so this simply delegates to one of
#' them rather than requiring a fourth compiled entry point.
#'
#' @keywords internal
#' @export
cpp_calc_observed_bicor <- function(expr) {
  cpp_bmht_observed_bicor(expr)
}

# ===================================================================
# Semi-parametric permutation-null tail approximation
# ===================================================================

#' Extrapolate permutation-null p-values below the raw permutation floor
#'
#' @param obs Numeric vector of observed statistics (one per gene/pair).
#' @param null_matrix Numeric matrix of permutation-null values, one row
#'   per gene/pair (matching `obs`), one column per permutation.
#' @param tail_prob Fraction of each row's null distribution treated as
#'   "the tail" for the GPD fit (default 0.10, i.e. the upper 10\%).
#' @param min_exceedances Minimum number of tail samples required for a
#'   GPD fit to be attempted; below this, falls back to the standard
#'   permutation-counting p-value for that row (default 20 - a common
#'   rule-of-thumb minimum for a stable method-of-moments GPD fit).
#' @return Numeric vector of p-values, one per row of `null_matrix`.
#' @keywords internal
#'
#' @section STATUS: EXPERIMENTAL, NOT USED BY ANY ENGINE, NOT RECOMMENDED:
#' This function is **not called anywhere in the BMHT/ROS-DET/BMKC
#' engines**. It was implemented specifically to address the
#' permutation-resolution gap identified in the Phase 5 benchmark (see
#' `dev/benchmarks/BENCHMARK_REPORT.md`), then tested against a 50,000+
#' replicate self-consistency calibration simulation before any attempt
#' to integrate it - and **failed that validation** under both
#' estimators tried:
#' \itemize{
#'   \item Method-of-moments GPD fit: ~27\% Type I error inflation right at
#'     the extrapolation boundary (a known small-sample bias in that
#'     estimator's shape parameter), though over-conservative deeper in
#'     the tail.
#'   \item Maximum-likelihood GPD fit (the current implementation, tried
#'     specifically to fix the method-of-moments problem): *worse*, not
#'     better - consistently anti-conservative with ratios of observed to
#'     nominal Type I error from 1.2x up to 4x deeper into the tail, and
#'     this persisted even at 10x more permutations (5,000 vs. 500),
#'     ruling out "just needs more permutation samples" as the fix.
#' }
#' This is consistent with a documented limitation of peaks-over-threshold
#' extreme value approximation: light-tailed distributions (Gaussian and
#' similar) converge to their asymptotic GPD limit very slowly, making
#' generic POT/GPD extrapolation unreliable at practically-achievable
#' permutation counts. Fixing this would need either a null-distribution-
#' specific parametric model (rather than a generic nonparametric EVT
#' technique) or a fundamentally different approach - not attempted here
#' given the calibration risk and remaining scope. Left in the codebase,
#' unintegrated, as a documented starting point and a record of what was
#' tried and why it didn't work, so this isn't silently re-attempted the
#' same way later. See `tests/testthat/test-gpd-tail-calibration.R` for
#' the calibration simulation that produced these numbers.
#'
#' @section What the function does (if you want to build on it):
#' Pure permutation-counting p-values (`(count(null >= obs) + 1) /
#' (n_perm + 1)`) have a hard floor of `1 / (n_perm + 1)`. Once corrected
#' for multiple testing across many genes or pairs, this coarse resolution
#' costs real statistical power - Phase 5's benchmark work found bicorX
#' trailing (semi-)analytical competitors (DGCA, dcanr) by roughly 3.5x in
#' ranking quality on a real ground-truth benchmark, traced specifically to
#' this floor rather than to the underlying bicor computation (which
#' Phase 1 confirmed is correct, and which Phase 5 confirmed is fast).
#'
#' This fits a Generalized Pareto Distribution (GPD) to the upper tail of
#' the permutation null (a standard peaks-over-threshold extreme-value
#' technique) and uses it to extrapolate a p-value when the observed
#' statistic exceeds the threshold - i.e. exactly the regime where raw
#' counting saturates at the floor. When the observed statistic does not
#' exceed the threshold, the empirical permutation count is already an
#' accurate estimate and no extrapolation is used. When there isn't enough
#' data in the tail for a stable fit, this falls back to the standard
#' permutation-counting formula rather than risk instability.
.gpd_tail_pvalues <- function(obs, null_matrix, tail_prob = 0.10, min_exceedances = 20) {
  n <- ncol(null_matrix)
  n_rows <- nrow(null_matrix)
  p_values <- numeric(n_rows)

  for (i in seq_len(n_rows)) {
    null_vals <- null_matrix[i, ]
    null_vals <- null_vals[!is.na(null_vals)]
    n_valid <- length(null_vals)
    x <- obs[i]

    if (n_valid == 0) {
      p_values[i] <- 1
      next
    }

    # Standard permutation-counting p-value - always computed, used as the
    # fallback and as the answer whenever obs doesn't reach the tail.
    count_ge <- sum(null_vals >= x)
    p_standard <- (count_ge + 1) / (n_valid + 1)

    threshold <- stats::quantile(null_vals, probs = 1 - tail_prob, names = FALSE, type = 7)
    exceedances <- null_vals[null_vals > threshold] - threshold
    k <- length(exceedances)

    # Only attempt GPD extrapolation when obs is actually in the tail
    # (above the threshold) - otherwise the empirical count is already
    # a good, more direct estimate and extrapolation adds nothing but risk.
    if (x <= threshold || k < min_exceedances) {
      p_values[i] <- p_standard
      next
    }

    mean_exc <- mean(exceedances)
    var_exc <- stats::var(exceedances)

    if (!is.finite(var_exc) || var_exc <= 0) {
      p_values[i] <- p_standard
      next
    }

    # Maximum-likelihood GPD fit. Method-of-moments was tried first and
    # rejected: calibration testing (50,000-replicate self-consistency
    # simulation under a known null) showed a ~27% Type I error inflation
    # right at the extrapolation boundary with method-of-moments, a known
    # small-sample bias issue for that estimator's shape parameter. MLE,
    # initialized from the method-of-moments values, tested clean across
    # the same simulation (see tests/testthat/test-gpd-tail-calibration.R).
    mom_xi <- 0.5 * ((mean_exc^2 / var_exc) - 1)
    mom_sigma <- 0.5 * mean_exc * ((mean_exc^2 / var_exc) + 1)

    if (!is.finite(mom_xi) || !is.finite(mom_sigma) || mom_sigma <= 0) {
      p_values[i] <- p_standard
      next
    }

    neg_log_lik <- function(par) {
      xi_par <- par[1]
      sigma_par <- exp(par[2])  # optimize log(sigma) to keep sigma > 0
      z <- 1 + xi_par * exceedances / sigma_par
      if (any(z <= 0)) return(1e10)
      if (abs(xi_par) < 1e-8) {
        k * log(sigma_par) + sum(exceedances) / sigma_par
      } else {
        k * log(sigma_par) + (1 + 1 / xi_par) * sum(log(z))
      }
    }

    fit <- tryCatch(
      stats::optim(c(mom_xi, log(mom_sigma)), neg_log_lik, method = "Nelder-Mead"),
      error = function(e) NULL
    )

    if (is.null(fit) || fit$convergence != 0) {
      p_values[i] <- p_standard
      next
    }

    xi <- fit$par[1]
    sigma <- exp(fit$par[2])

    if (!is.finite(xi) || !is.finite(sigma) || sigma <= 0) {
      p_values[i] <- p_standard
      next
    }

    exceedance_x <- x - threshold
    z <- 1 + xi * exceedance_x / sigma

    if (xi < 0 && z <= 0) {
      # Outside the GPD's support for a negative shape parameter - the fit
      # says this value is essentially impossible under this tail model,
      # which is a sign the model shouldn't be trusted this far out. Fall
      # back rather than report an unreliable (and non-monotonic) value.
      p_values[i] <- p_standard
      next
    }

    tail_mass <- k / n_valid
    p_gpd <- if (abs(xi) < 1e-8) {
      tail_mass * exp(-exceedance_x / sigma)
    } else {
      tail_mass * z^(-1 / xi)
    }

    # The GPD estimate should never be larger than what plain counting
    # already gives at the floor (it's an extrapolation *below* the floor,
    # not an alternative estimate that could go the other way) - clamp
    # defensively in case of numerical edge cases.
    p_values[i] <- min(p_gpd, p_standard)
  }

  p_values
}

# ===================================================================
# ECOR: analytical significance test for equality of two correlations
# (Kayano, Takigawa, Shiga, Tsuda & Mamitsuka 2011, "ROS-DET: robust
# detector of switching mechanisms in gene expression", Nucleic Acids
# Research 39(11):e74, Equation 4)
# ===================================================================

#' Analytical p-value for the equality of two correlation coefficients
#'
#' Implements the likelihood-ratio chi-squared(1) test from Kayano et al.
#' 2011 (ROS-DET's "ECOR" step), Equation 4. The test statistic requires
#' solving an implicit equation for the pooled correlation estimate under
#' H0 (rho_1 = rho_2 = rho_hat):
#'   N1*(r1 - rho)/(1 - rho*r1) + N2*(r2 - rho)/(1 - rho*r2) = 0
#' then:
#'   T = N1*log[(1-r1*rho_hat)^2 / ((1-r1^2)(1-rho_hat^2))]
#'     + N2*log[(1-r2*rho_hat)^2 / ((1-r2^2)(1-rho_hat^2))]  ~ chisq(1)
#'
#' The paper derives this assuming each condition's data is bivariate
#' normal with `r1`/`r2` as Pearson correlations, then applies the same
#' formula to biweight midcorrelation values as a disclosed approximation
#' (validated empirically in their supplementary material, not claimed to
#' be exact for bicor).
#'
#' @section Calibration (verified before this was wired into any engine):
#' A self-consistency simulation (pure-noise data, bicor computed via this
#' package's own kernel, comparing observed vs. nominal false-positive
#' rates) shows:
#' \itemize{
#'   \item N >= 50 per condition: well-calibrated (~1.0-1.2x nominal at
#'     p < 0.05/0.01/0.001).
#'   \item N ~= 20-30 per condition: mildly anti-conservative (~1.2-1.7x).
#'   \item N < 20 per condition: meaningfully anti-conservative (~1.7-2.2x
#'     at N=10), consistent with the known finite-sample bias of an
#'     asymptotic chi-squared approximation.
#' }
#' This degrades gracefully and predictably with N (unlike the four
#' generic parametric approaches tried and rejected for the permutation-
#' floor problem elsewhere in this package - GPD method-of-moments, GPD
#' MLE, plain normal, skew-normal - which showed unpredictable,
#' sometimes-worsening-in-the-tail miscalibration even at moderate N).
#' This is because ECOR is a theoretically-derived asymptotic result, not
#' a distribution empirically fit to permutation output - see
#' tests/testthat/test-ecor-calibration.R for the validation this claim is
#' based on.
#'
#' @param r1,r2 Observed correlation (or bicor) coefficients for the two
#'   conditions. May be vectors (one call per pair).
#' @param N1,N2 Sample sizes for the two conditions. Recycled against
#'   `r1`/`r2` if scalar.
#' @return Numeric vector of p-values, one per pair. `NA` for a pair where
#'   the root-finder fails to converge (rare; falls back gracefully rather
#'   than erroring the whole batch).
#' @keywords internal
.ecor_pvalue <- function(r1, r2, N1, N2) {
  n <- max(length(r1), length(r2))
  r1 <- rep_len(r1, n); r2 <- rep_len(r2, n)
  N1 <- rep_len(N1, n); N2 <- rep_len(N2, n)
  p_values <- numeric(n)

  pooled_rho_eq <- function(rho, r1, r2, N1, N2) {
    N1 * (r1 - rho) / (1 - rho * r1) + N2 * (r2 - rho) / (1 - rho * r2)
  }

  for (i in seq_len(n)) {
    x1 <- r1[i]; x2 <- r2[i]; n1 <- N1[i]; n2 <- N2[i]

    if (is.na(x1) || is.na(x2)) { p_values[i] <- NA; next }

    # rho_hat lies in (-1, 1); search a bracket around the two observed
    # values first (where it usually is), widening to the full range if
    # that bracket doesn't contain a sign change.
    lo <- max(-0.9999, min(x1, x2) - 0.01)
    hi <- min(0.9999, max(x1, x2) + 0.01)
    f_lo <- pooled_rho_eq(lo, x1, x2, n1, n2)
    f_hi <- pooled_rho_eq(hi, x1, x2, n1, n2)
    if (is.na(f_lo) || is.na(f_hi) || sign(f_lo) == sign(f_hi)) {
      lo <- -0.9999; hi <- 0.9999
      f_lo <- pooled_rho_eq(lo, x1, x2, n1, n2)
      f_hi <- pooled_rho_eq(hi, x1, x2, n1, n2)
      if (is.na(f_lo) || is.na(f_hi) || sign(f_lo) == sign(f_hi)) {
        p_values[i] <- NA
        next
      }
    }

    rho_hat <- tryCatch(
      stats::uniroot(pooled_rho_eq, c(lo, hi), r1 = x1, r2 = x2, N1 = n1, N2 = n2,
                      tol = 1e-10)$root,
      error = function(e) NA
    )
    if (is.na(rho_hat)) { p_values[i] <- NA; next }

    denom1 <- (1 - x1^2) * (1 - rho_hat^2)
    denom2 <- (1 - x2^2) * (1 - rho_hat^2)
    if (denom1 <= 0 || denom2 <= 0) { p_values[i] <- NA; next }

    T_stat <- n1 * log((1 - x1 * rho_hat)^2 / denom1) +
              n2 * log((1 - x2 * rho_hat)^2 / denom2)
    p_values[i] <- 1 - stats::pchisq(T_stat, df = 1)
  }

  p_values
}
