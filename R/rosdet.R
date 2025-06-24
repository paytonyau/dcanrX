#' @importFrom WGCNA bicor
#' @importFrom stats optimize pchisq
NULL

#' Run the ROS-DET analysis to find switching gene pairs.
#'
#' Implements the methodology from the ROS-DET paper  to identify gene pairs
#' that switch from positive to negative correlation between two conditions.
#' The method is robust to outliers, range bias, and small sample sizes.
#'
#' @param expr_matrix A numeric matrix where rows are genes and columns are samples.
#' @param condition A factor or character vector with two levels, matching the columns of expr_matrix.
#' @param significance_level The alpha level for the ECOR P-value cutoff.
#' @return A data frame with the ranked and filtered gene pairs.
#' @export
run_rosdet <- function(expr_matrix, condition, significance_level = 0.05) {

  conditions <- unique(condition)
  if (length(conditions) != 2) {
    stop("Condition vector must contain exactly two unique conditions.")
  }

  gene_names <- rownames(expr_matrix)
  gene_pairs <- utils::combn(gene_names, 2, simplify = FALSE)

  all_results <- do.call(rbind, lapply(gene_pairs, function(pair) {
    gene1_name <- pair[1]
    gene2_name <- pair[2]

    x_expr <- as.numeric(expr_matrix[gene1_name, ])
    y_expr <- as.numeric(expr_matrix[gene2_name, ])

    cond1 <- conditions[1]
    cond2 <- conditions[2]

    # --- WCOR Step ---
    x1 <- x_expr[condition == cond1]; y1 <- y_expr[condition == cond1]
    x2 <- x_expr[condition == cond2]; y2 <- y_expr[condition == cond2]

    # Biweight Midcorrelation for each condition
    r1 <- WGCNA::bicor(x1, y1)
    r2 <- WGCNA::bicor(x2, y2)

    # Weight `c` to handle range bias, using biweight midvariances
    q_x1 <- biweight_midvariance(x1); q_y1 <- biweight_midvariance(y1)
    q_x2 <- biweight_midvariance(x2); q_y2 <- biweight_midvariance(y2)

    max_var_1 <- max(q_x1, q_y1, na.rm = TRUE)
    max_var_2 <- max(q_x2, q_y2, na.rm = TRUE)

    c_weight <- if (max_var_1 == 0 || max_var_2 == 0) 0 else min(max_var_1 / max_var_2, max_var_2 / max_var_1)

    # Final WCOR score
    score <- c_weight * abs(r1 - r2)

    # --- ECOR Step ---
    n1 <- length(x1); n2 <- length(x2)

    # Objective function to find rho_hat
    objective_fn <- function(rho_hat) {
      term1 <- n1 * (r1 - rho_hat) / (1 - rho_hat * r1)
      term2 <- n2 * (r2 - rho_hat) / (1 - rho_hat * r2)
      return((term1 + term2)^2)
    }

    rho_hat <- stats::optimize(objective_fn, interval = c(-0.99, 0.99))$minimum

    # Chi-square test statistic T
    term1_log <- n1 * log(((1 - r1 * rho_hat)^2) / ((1 - r1^2) * (1 - rho_hat^2)))
    term2_log <- n2 * log(((1 - r2 * rho_hat)^2) / ((1 - r2^2) * (1 - rho_hat^2)))

    T_statistic <- sum(term1_log, term2_log, na.rm = TRUE)

    p_value <- stats::pchisq(T_statistic, df = 1, lower.tail = FALSE)

    data.frame(
      Gene1 = gene1_name, Gene2 = gene2_name,
      Score = score, P_value = p_value,
      r1 = r1, r2 = r2, c_weight = c_weight
    )
  }))

  # Bonferroni correction for multiple testing
  corrected_alpha <- significance_level / nrow(all_results)

  filtered_results <- all_results[all_results$P_value < corrected_alpha, ]
  ranked_results <- filtered_results[order(-filtered_results$Score), ]

  return(ranked_results)
}
