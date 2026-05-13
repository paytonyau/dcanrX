#' @title Utility Functions for bicorDCEA
#' @description Core utilities: biweight midcorrelation (bicor),
#' biweight midvariance, correlation matrix computation, and input validation.
#'
#' @importFrom stats median mad complete.cases cor var sd
#' @importFrom BiocParallel bpparam MulticoreParam SerialParam SnowParam bplapply
#' @name utils
NULL

# ===================================================================
# Smart Parallel Dispatcher (Cross-Platform)
# ===================================================================

#' @keywords internal
.get_bp_param <- function(workers = 1) {
  if (workers <= 1) {
    return(BiocParallel::SerialParam(progressbar = TRUE))
  }

  if (.Platform$OS.type == "windows") {
    # Windows requires Socket/Snow clusters
    return(BiocParallel::SnowParam(workers = workers, progressbar = TRUE))
  } else {
    # Mac/Linux can use faster Forking/Multicore clusters
    return(BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE))
  }
}

# ===================================================================
# Input Validation (shared across all methods)
# ===================================================================

#' @keywords internal
.validate_input <- function(expr_matrix, condition) {
  if (!is.matrix(expr_matrix) || !is.numeric(expr_matrix)) {
    stop("`expr_matrix` must be a standard numeric matrix (genes in rows, samples in columns). If using sparse matrices, please convert via as.matrix().")
  }
  if (is.null(rownames(expr_matrix)) || any(duplicated(rownames(expr_matrix)))) {
    stop("`expr_matrix` must have unique gene names as rownames.")
  }

  # Matrix Integrity Check (Protects C-Engine and CLR math)
  if (any(is.na(expr_matrix)) || any(is.infinite(expr_matrix))) {
    stop("`expr_matrix` contains NA, NaN, or Inf values. Please impute or filter these before running.")
  }

  # Zero-Variance Guardrail (Protects Bulk RNA-seq FDR Power)
  row_vars <- apply(expr_matrix, 1, stats::var)
  n_zero_var <- sum(row_vars == 0, na.rm = TRUE)
  if (n_zero_var > 0) {
    warning(sprintf("Detected %d zero-variance genes. They will safely return 0 correlation, but filtering them upstream is highly recommended to save RAM and boost FDR statistical power.", n_zero_var))
  }

  if (length(condition) != ncol(expr_matrix)) {
    stop("Length of `condition` must equal the number of columns in `expr_matrix`.")
  }

  # Prevent NA poisoning from dirty clinical/metadata
  if (any(is.na(condition))) {
    stop("`condition` vector contains NAs. Please filter out these samples/cells before running.")
  }

  condition <- factor(condition)
  if (length(levels(condition)) != 2) {
    stop("`condition` must contain exactly two unique levels.")
  }
  invisible(condition)
}

# ===================================================================
# Biweight Midcorrelation (bicor)
# ===================================================================

#' Biweight Midcorrelation (bicor) with Pearson Fallback
#'
#' A robust measure of correlation, highly resistant to outliers.
#' Includes a fallback to Pearson correlation for highly sparse vectors
#' (Zero-MAD problem common in single-cell data).
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @return Numeric scalar between -1 and 1 (NA if insufficient data)
#' @export
bicor <- function(x, y) {
  complete_obs <- stats::complete.cases(x, y)
  x <- x[complete_obs]
  y <- y[complete_obs]

  if (length(x) < 2) return(NA_real_)

  median_x <- stats::median(x)
  median_y <- stats::median(y)
  mad_x <- stats::mad(x)
  mad_y <- stats::mad(y)

  # Pearson Fallback for Zero-MAD (Single-cell sparsity protection)
  if (mad_x == 0 || mad_y == 0) {
    # If standard deviation is also zero, return 0 to avoid NA from cor()
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
  mad_x <- stats::mad(x)

  # Standard variance fallback for highly sparse single-cell data
  if (mad_x == 0) {
    return(stats::var(x))
  }

  u <- (x - median_x) / (9 * mad_x)

  # Only consider values where |u| < 1
  valid_u <- abs(u) < 1
  u_valid <- u[valid_u]
  x_valid <- x[valid_u]

  if (length(u_valid) == 0) return(0)

  # Standard Robust Estimation formula (Wilcox Estimator)
  numerator <- n * sum(((x_valid - median_x)^2) * ((1 - u_valid^2)^4))
  denom_term <- sum((1 - u_valid^2) * (1 - 5 * u_valid^2))

  if (denom_term == 0) return(0)

  numerator / (denom_term^2)
}

# ===================================================================
# Fast bicor Correlation Matrix with Parallel & WGCNA Support
# ===================================================================

#' Compute full bicor correlation matrix (internal)
#'
#' @param mat Numeric matrix (genes × samples)
#' @param workers Integer. Number of parallel workers
#' @return Symmetric numeric matrix with gene names
#' @keywords internal
compute_cor_matrix <- function(mat, workers = 1) {
  n_genes <- nrow(mat)

  # Route to WGCNA's C-engine for massive performance gains on big data
  if (requireNamespace("WGCNA", quietly = TRUE)) {
    message(sprintf("Using WGCNA's optimized C-engine for bicor matrix (%d genes)...", n_genes))
    cor_mat <- WGCNA::bicor(
      t(mat),
      maxPOutliers = 0.05,
      nThreads = workers,
      robustX = TRUE,
      pearsonFallback = "individual"
    )

    cor_mat[is.na(cor_mat)] <- 0
    diag(cor_mat) <- 1.0
    return(cor_mat)
  }

  if (n_genes > 5000) {
    warning("WGCNA package is missing. Native R bicor on >5,000 genes will be extremely slow. Please install WGCNA.")
  }

  cor_mat <- matrix(0, nrow = n_genes, ncol = n_genes)
  gene_pairs <- utils::combn(seq_len(n_genes), 2, simplify = FALSE)

  # FIX: Explicitly pass dependencies to the worker nodes to prevent Windows SnowParam crashes
  run_one_pair <- function(p, data_mat, bicor_fn) {
    bicor_fn(data_mat[p[1], ], data_mat[p[2], ])
  }

  BPPARAM <- .get_bp_param(workers)

  # We use bplapply's ability to pass external arguments directly to the workers
  vals <- BiocParallel::bplapply(gene_pairs, run_one_pair,
                                 data_mat = mat,
                                 bicor_fn = bicor,
                                 BPPARAM = BPPARAM)

  cor_mat[lower.tri(cor_mat)] <- unlist(vals)
  cor_mat <- cor_mat + t(cor_mat)
  diag(cor_mat) <- 1.0

  rownames(cor_mat) <- colnames(cor_mat) <- rownames(mat)
  cor_mat
}

# ===================================================================
# Compositional Data Transformation (for microbial / highly sparse data)
# ===================================================================

#' @keywords internal
.compositional_transform <- function(mat, transform = c("none", "clr", "log1p"), species = NULL) {
  transform <- match.arg(transform)

  if (transform == "none") return(mat)

  if (transform == "log1p") {
    # Prevent log of negative numbers
    if (any(mat < 0, na.rm = TRUE)) {
      stop("log1p transformation requires non-negative data. Your matrix contains negative values.")
    }
    return(log1p(mat))
  }

  if (transform == "clr") {
    # CRITICAL QC: CLR is for count/abundance data only
    if (any(mat < 0, na.rm = TRUE)) {
      stop("CLR transformation requires strictly non-negative data (counts or abundances). Your matrix contains negative values (perhaps it is already scaled?).")
    }

    # Data-driven pseudo-count (half the minimum non-zero value)
    min_nonzero <- min(mat[mat > 0], na.rm = TRUE)
    # is.finite checks if min() returned Inf due to all zeros
    pseudo_count <- ifelse(is.finite(min_nonzero), min_nonzero / 2, 1e-8)
    mat_pseudo <- mat + pseudo_count

    # Species-aware CLR for multi-omics (e.g., Host + Microbiome)
    if (!is.null(species) && length(unique(species)) > 1) {
      message("Multi-species detected: Applying CLR transformation independently per species...")
      species <- as.character(species)
      unique_spp <- unique(species)

      clr_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
      rownames(clr_mat) <- rownames(mat)
      colnames(clr_mat) <- colnames(mat)

      for (sp in unique_spp) {
        idx <- which(species == sp)
        # Transpose to (samples x genes)
        sp_mat_t <- t(mat_pseudo[idx, , drop = FALSE])

        # Safe log-space subtraction (prevents exp() underflow)
        log_sp_mat <- log(sp_mat_t)
        sp_clr_t <- log_sp_mat - rowMeans(log_sp_mat)

        clr_mat[idx, ] <- t(sp_clr_t)
      }
      return(clr_mat)

    } else {
      # Standard CLR for single composition
      mat_t <- t(mat_pseudo)

      # Safe log-space subtraction
      log_mat <- log(mat_t)
      clr_mat_t <- log_mat - rowMeans(log_mat)

      return(t(clr_mat_t))
    }
  }
}
