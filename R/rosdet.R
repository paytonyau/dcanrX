#' @importFrom stats optimize pchisq
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_bw
NULL

#' Run the ROS-DET analysis to find switching gene pairs.
#'
#' Implements the methodology from the ROS-DET paper to identify gene pairs
#' that switch from positive to negative correlation between two conditions.
#' This improved version includes parallel processing capabilities.
#'
#' @param expr_matrix A numeric matrix where rows are genes and columns are samples.
#' @param condition A factor or character vector with two levels, matching the columns of `expr_matrix`.
#' @param significance_level The alpha level for the ECOR P-value cutoff.
#' @param workers An integer specifying the number of parallel workers to use.
#'   Defaults to 1 (no parallelization). Set to a higher number to speed up computation.
#' @return An object of class `rosdet_result` containing the ranked and filtered gene pairs.
#'   This object can be used with `summary()` and `plot()`.
#' @export
#' @examples
#' \dontrun{
#' # Load the package and create dummy data
#' library(dcanr)
#' set.seed(42)
#' expr_matrix <- matrix(rnorm(100 * 40), nrow=100)
#' rownames(expr_matrix) <- paste0("Gene", 1:100)
#' condition <- factor(c(rep("Normal", 20), rep("Tumor", 20)))
#'
#' # Run analysis using 2 parallel workers
#' rosdet_results <- run_rosdet(expr_matrix, condition, workers = 2)
#'
#' # Get a high-level summary
#' summary(rosdet_results)
#'
#' # Automatically plot the top result
#' plot(rosdet_results)
#' }
run_rosdet <- function(expr_matrix, condition, significance_level = 0.05, workers = 1) {

  # --- Input Validation ---
  conditions <- unique(condition)
  if (length(conditions) != 2) {
    stop("Condition vector must contain exactly two unique conditions.")
  }

  # --- Setup Parallel Processing ---
  # Set the parallel plan using the 'future' package.
  # This automatically reverts to the previous plan when the function exits.
  if (workers > 1) {
    future::plan(future::multisession, workers = workers)
    message(paste("Using", workers, "parallel workers for ROS-DET analysis."))
  } else {
    future::plan("sequential")
  }
  on.exit(future::plan("sequential"), add = TRUE)

  gene_names <- rownames(expr_matrix)
  gene_pairs <- utils::combn(gene_names, 2, simplify = FALSE)

  # Use future_lapply for parallel execution
  all_results_list <- future.apply::future_lapply(gene_pairs, FUN = function(pair) {
    # This internal code runs for each gene pair, potentially in parallel
    gene1_name <- pair[1]
    gene2_name <- pair[2]

    x_expr <- as.numeric(expr_matrix[gene1_name, ])
    y_expr <- as.numeric(expr_matrix[gene2_name, ])

    cond1 <- conditions[1]
    cond2 <- conditions[2]

    # --- WCOR Step ---
    x1 <- x_expr[condition == cond1]; y1 <- y_expr[condition == cond1]
    x2 <- x_expr[condition == cond2]; y2 <- y_expr[condition == cond2]

    r1 <- bicor(x1, y1)
    r2 <- bicor(x2, y2)

    q_x1 <- biweight_midvariance(x1); q_y1 <- biweight_midvariance(y1)
    q_x2 <- biweight_midvariance(x2); q_y2 <- biweight_midvariance(y2)

    max_var_1 <- max(q_x1, q_y1, na.rm = TRUE)
    max_var_2 <- max(q_x2, q_y2, na.rm = TRUE)

    c_weight <- if (max_var_1 == 0 || max_var_2 == 0) 0 else min(max_var_1 / max_var_2, max_var_2 / max_var_1)

    score <- c_weight * abs(r1 - r2)

    # --- ECOR Step ---
    n1 <- length(x1); n2 <- length(x2)

    objective_fn <- function(rho_hat) {
      term1 <- n1 * (r1 - rho_hat) / (1 - rho_hat * r1)
      term2 <- n2 * (r2 - rho_hat) / (1 - rho_hat * r2)
      return((term1 + term2)^2)
    }

    if (is.na(r1) || is.na(r2)) {
      return(NULL) # Skip pairs where correlation could not be computed
    }

    rho_hat <- stats::optimize(objective_fn, interval = c(-0.99, 0.99))$minimum

    term1_log <- n1 * log(((1 - r1 * rho_hat)^2) / ((1 - r1^2) * (1 - rho_hat^2)))
    term2_log <- n2 * log(((1 - r2 * rho_hat)^2) / ((1 - r2^2) * (1 - rho_hat^2)))

    T_statistic <- sum(term1_log, term2_log, na.rm = TRUE)

    p_value <- stats::pchisq(T_statistic, df = 1, lower.tail = FALSE)

    data.frame(
      Gene1 = gene1_name, Gene2 = gene2_name,
      Score = score, P_value = p_value,
      r1 = r1, r2 = r2, c_weight = c_weight,
      stringsAsFactors = FALSE
    )
  }, future.seed = TRUE) # future.seed = TRUE ensures reproducibility in parallel

  # Combine the list of single-row data frames into one large data frame
  all_results <- do.call(rbind, all_results_list)

  # Bonferroni correction for multiple testing
  num_tests <- nrow(all_results)
  if (num_tests > 0) {
    corrected_alpha <- significance_level / num_tests
    filtered_results <- all_results[all_results$P_value < corrected_alpha, ]
    ranked_results <- if(nrow(filtered_results) > 0) {
      filtered_results[order(-filtered_results$Score), ]
    } else {
      filtered_results
    }
  } else {
    ranked_results <- all_results # Return empty data frame if no pairs were tested
  }

  # Create the custom S3 object
  output <- list(
    results = ranked_results,
    n_tested = num_tests,
    significance_level = significance_level,
    # Store original data needed for plotting
    data = list(
      expr_matrix = expr_matrix,
      condition = condition
    )
  )

  class(output) <- "rosdet_result"
  return(output)
}

#' Print a summary of ROS-DET results.
#'
#' Provides a concise summary of the results from a `run_rosdet` analysis.
#'
#' @param object An object of class `rosdet_result`.
#' @param ... Additional arguments (not used).
#' @return Prints a summary to the console, invisibly returning the object.
#' @export
summary.rosdet_result <- function(object, ...) {
  n_significant <- nrow(object$results)

  cat("--- ROS-DET Analysis Summary ---\n")
  cat(paste("Number of gene pairs tested:", format(object$n_tested, big.mark = ","), "\n"))
  cat(paste("Significance level (alpha):", object$significance_level, "\n"))
  cat("Correction method: Bonferroni\n")
  cat("----------------------------------\n")

  if (n_significant > 0) {
    cat(paste("Found", n_significant, "significant switching gene pairs.\n\n"))
    cat("Top 6 results:\n")
    print(head(object$results))
  } else {
    cat("No significant switching gene pairs were found at this significance level.\n")
  }
  invisible(object)
}

#' Plot the top switching gene pair from a ROS-DET analysis.
#'
#' Generates a scatter plot visualizing the correlation switch for the
#' top-ranked gene pair identified by `run_rosdet`.
#'
#' @param x An object of class `rosdet_result`.
#' @param ... Additional arguments (not used).
#' @return A ggplot2 object.
#' @export
plot.rosdet_result <- function(x, ...) {
  if (nrow(x$results) == 0) {
    message("Cannot plot. No significant results found in the object.")
    return(invisible(NULL))
  }

  # Get the top-ranked pair and its data
  top_pair <- x$results[1, ]
  expr_matrix <- x$data$expr_matrix
  condition <- x$data$condition

  gene1_name <- top_pair$Gene1
  gene2_name <- top_pair$Gene2

  # Create a data frame for plotting
  plot_df <- data.frame(
    Gene1_Expr = as.numeric(expr_matrix[gene1_name, ]),
    Gene2_Expr = as.numeric(expr_matrix[gene2_name, ]),
    Condition = condition
  )

  conditions <- levels(factor(condition))

  # Create the plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Gene1_Expr, y = Gene2_Expr, color = Condition)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
    ggplot2::labs(
      title = paste("ROS-DET Top Switching Pair:", gene1_name, "&", gene2_name),
      subtitle = paste0(
        "Correlation in '", conditions[1], "': ", round(top_pair$r1, 2),
        " | Correlation in '", conditions[2], "': ", round(top_pair$r2, 2)
      ),
      x = paste(gene1_name, "Expression"),
      y = paste(gene2_name, "Expression")
    ) +
    ggplot2::theme_bw()

  return(p)
}

#' @export
print.rosdet_result <- function(x, ...) {
  summary(x, ...)
}
