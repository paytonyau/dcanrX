#' Calculate Biweight Midcorrelation (bicor)
#'
#' A robust measure of correlation that is less sensitive to outliers than
#' Pearson correlation. This is an independent implementation based on the
#' formulas described by Kayano et al. (2011), Zheng et al. (2014),
#' and Yuan et al. (2015).
#'
# --- Shared Helper Functions ---

#' Calculate Biweight Midcorrelation (bicor)
#'
#' A robust measure of correlation that is less sensitive to outliers.
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @return The biweight midcorrelation coefficient.
#' @export
bicor <- function(x, y) {
  # ... (code for bicor remains the same)
  complete_obs <- stats::complete.cases(x, y)
  x <- x[complete_obs]; y <- y[complete_obs]
  if(length(x) < 2) return(NA)
  median_x <- stats::median(x); median_y <- stats::median(y)
  mad_x <- stats::mad(x); mad_y <- stats::mad(y)
  if (mad_x == 0 || mad_y == 0) return(0)
  u <- (x - median_x) / (9 * mad_x); v <- (y - median_y) / (9 * mad_y)
  w_x <- ifelse(abs(u) < 1, (1 - u^2)^2, 0)
  w_y <- ifelse(abs(v) < 1, (1 - v^2)^2, 0)
  numerator <- sum((x - median_x) * w_x * (y - median_y) * w_y)
  denom_x_sq <- sum(((x - median_x) * w_x)^2)
  denom_y_sq <- sum(((y - median_y) * w_y)^2)
  denominator <- sqrt(denom_x_sq) * sqrt(denom_y_sq)
  if (denominator == 0) return(0)
  return(numerator / denominator)
}

#' @keywords internal
biweight_midvariance <- function(x) {
  # ... (code for biweight_midvariance remains the same)
  x <- x[!is.na(x)]; median_x <- stats::median(x); mad_x <- stats::mad(x)
  if (mad_x == 0) return(0)
  u <- (x - median_x) / (9 * mad_x); w <- ifelse(abs(u) < 1, (1 - u^2)^2, 0)
  numerator <- sum(w * (x - median_x)^2); denominator <- sum(w)
  if (denominator == 0) return(0)
  return(length(x) * numerator / (denominator^2))
}

# Helper to compute a correlation matrix using the bicor function
# This version uses the stable, built-in 'parallel' package.
compute_cor_matrix <- function(mat, workers = 1) {
  n_genes <- nrow(mat)
  cor_mat <- matrix(0, n_genes, n_genes)
  gene_pairs <- utils::combn(1:n_genes, 2, simplify = FALSE) # Return as a list

  # This function runs on a single gene pair
  run_one_pair <- function(p) {
    return(bicor(mat[p[1], ], mat[p[2], ]))
  }

  # Detect OS and choose the appropriate parallel function
  is_windows <- .Platform$OS.type == "windows"

  vals <- NULL
  if (workers > 1 && !is_windows) {
    # Use mclapply on Mac/Linux
    vals <- parallel::mclapply(gene_pairs, run_one_pair, mc.cores = workers)
  } else if (workers > 1 && is_windows) {
    # Use parLapply on Windows
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, varlist = c("mat", "bicor"), envir = environment())
    vals <- parallel::parLapply(cl, gene_pairs, run_one_pair)
  } else {
    # Fallback to sequential execution
    vals <- lapply(gene_pairs, run_one_pair)
  }

  # Populate the correlation matrix
  cor_mat[lower.tri(cor_mat)] <- unlist(vals)
  cor_mat <- cor_mat + t(cor_mat)
  diag(cor_mat) <- 1
  rownames(cor_mat) <- colnames(cor_mat) <- rownames(mat)
  return(cor_mat)
}
