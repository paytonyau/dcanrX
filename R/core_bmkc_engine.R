#' @title BMKC: Biweight Midcorrelation + k-Clique Method
#' @description Identifies gene modules (cliques) that are strongly co-expressed
#' in one condition but not in the other.
#'
#' Modernized to support cross-species microbiome data, signed vs. unsigned
#' networks, two-sided differential search, and adaptive percentile thresholding
#' for highly sparse single-cell data.
#'
#' @importFrom igraph graph_from_adjacency_matrix max_cliques
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient labs theme_minimal
#' @importFrom stats quantile
#' @name bmkc
NULL

#' Run BMKC Analysis
#'
#' @param expr_matrix Numeric matrix (genes × samples)
#' @param condition Factor/character vector with exactly two levels
#' @param T1 High co-expression threshold for the active condition
#' @param T2 Low co-expression threshold for the inactive condition
#' @param min_clique_size Minimum number of genes in a reported module
#' @param workers Number of parallel workers
#' @param transform Compositional transformation ("none", "clr", "log1p")
#' @param species Optional character vector indicating species for each gene
#' @param network_type "signed" (positive correlation only) or "unsigned" (absolute)
#' @param threshold_method "absolute" (uses T1/T2 directly) or "percentile" (T1/T2 as quantiles)
#' @param .precomputed_cors Internal use only for threshold exploration speedups.
#' @return An object of class `bmkc_result`
#' @export
run_bmkc <- function(expr_matrix, condition,
                     T1 = 0.7, T2 = 0.3,
                     min_clique_size = 4,
                     workers = 1,
                     transform = c("none", "clr", "log1p"),
                     species = NULL,
                     network_type = c("signed", "unsigned"),
                     threshold_method = c("absolute", "percentile"),
                     .precomputed_cors = NULL) {

  condition <- .validate_input(expr_matrix, condition)
  transform <- match.arg(transform)
  network_type <- match.arg(network_type)
  threshold_method <- match.arg(threshold_method)

  # Strict species validation
  species_mode <- FALSE
  if (!is.null(species)) {
    if (length(species) != nrow(expr_matrix)) {
      stop("Length of `species` must exactly match the number of rows in `expr_matrix`.")
    }
    species_mode <- TRUE
    species <- as.character(species)
  }

  # Sparsity Guardrail: Catch raw single-cell/microbiome matrices
  sparsity <- sum(expr_matrix == 0, na.rm = TRUE) / length(expr_matrix)
  if (sparsity > 0.6) {
    warning(sprintf("High sparsity detected (%.1f%% zeros). For BMKC clique detection, it is highly recommended to use `run_dcanr_pseudobulk()` first to restore correlation magnitudes.", sparsity * 100))
  }

  conditions <- levels(condition)

  # Check for precomputed matrices (Massive speedup for threshold exploration)
  if (!is.null(.precomputed_cors)) {
    cor_1 <- .precomputed_cors[[1]]
    cor_2 <- .precomputed_cors[[2]]
  } else {
    expr_matrix <- .compositional_transform(expr_matrix, transform, species)
    message("Computing observed bicor matrices...")
    cor_1 <- compute_cor_matrix(expr_matrix[, condition == conditions[1]], workers = workers)
    cor_2 <- compute_cor_matrix(expr_matrix[, condition == conditions[2]], workers = workers)
  }

  # ===================================================================
  # Threshold Logic (Adaptive vs Absolute)
  # ===================================================================
  if (threshold_method == "percentile") {
    # Auto-correct whole numbers (e.g., 99 to 0.99)
    if (T1 > 1) T1 <- T1 / 100
    if (T2 > 1) T2 <- T2 / 100

    vals1 <- cor_1[upper.tri(cor_1)]
    vals2 <- cor_2[upper.tri(cor_2)]

    if (network_type == "unsigned") {
      T1_c1 <- stats::quantile(abs(vals1), T1, na.rm = TRUE)
      T1_c2 <- stats::quantile(abs(vals2), T1, na.rm = TRUE)
    } else {
      T1_c1 <- stats::quantile(vals1, T1, na.rm = TRUE)
      T1_c2 <- stats::quantile(vals2, T1, na.rm = TRUE)
    }

    # T2 (inactive threshold) is always evaluated on absolute magnitude
    # (a strong negative correlation is still an active biological connection)
    T2_c1 <- stats::quantile(abs(vals1), T2, na.rm = TRUE)
    T2_c2 <- stats::quantile(abs(vals2), T2, na.rm = TRUE)

  } else {
    T1_c1 <- T1_c2 <- T1
    T2_c1 <- T2_c2 <- T2
  }

  # ===================================================================
  # Two-Sided Differential Adjacency Logic
  # ===================================================================
  if (network_type == "unsigned") {
    adj_cond1_active <- (abs(cor_1) >= T1_c1) & (abs(cor_2) <= T2_c2)
    adj_cond2_active <- (abs(cor_2) >= T1_c2) & (abs(cor_1) <= T2_c1)
  } else {
    adj_cond1_active <- (cor_1 >= T1_c1) & (abs(cor_2) <= T2_c2)
    adj_cond2_active <- (cor_2 >= T1_c2) & (abs(cor_1) <= T2_c1)
  }

  # Helper function to safely extract cliques while avoiding NP-Hard freezes
  extract_cliques <- function(adj_logical, label) {
    adj_logical[is.na(adj_logical)] <- FALSE
    diag(adj_logical) <- FALSE

    adj_matrix <- matrix(as.numeric(adj_logical), nrow = nrow(adj_logical))
    rownames(adj_matrix) <- colnames(adj_matrix) <- rownames(expr_matrix)

    n_nodes <- nrow(adj_matrix)
    n_edges <- sum(adj_matrix) / 2
    max_edges <- (n_nodes * (n_nodes - 1)) / 2
    graph_density <- n_edges / max_edges

    if (graph_density > 0.05 && graph_density <= 0.15) {
      warning(sprintf("Graph density for %s is high (%.1f%%). Finding maximal cliques may take time.", label, graph_density * 100))
    }
    if (graph_density > 0.15) {
      warning(sprintf("Graph density for %s exceeds 15%% (%.1f%%). Execution skipped to prevent R from crashing. Raise T1 or use percentile thresholding.", label, graph_density * 100))
      return(list())
    }

    graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
    cliques <- igraph::max_cliques(graph, min = min_clique_size)
    lapply(cliques, function(c) names(c))
  }

  message("Extracting cliques for Condition 1...")
  cliques_c1 <- extract_cliques(adj_cond1_active, conditions[1])

  message("Extracting cliques for Condition 2...")
  cliques_c2 <- extract_cliques(adj_cond2_active, conditions[2])

  # Combine results
  all_modules <- c(cliques_c1, cliques_c2)
  module_conditions <- c(rep(conditions[1], length(cliques_c1)),
                         rep(conditions[2], length(cliques_c2)))

  # Species awareness packaging
  module_info <- NULL
  if (species_mode && length(all_modules) > 0) {
    module_info <- lapply(all_modules, function(genes) {
      sp_in_module <- species[rownames(expr_matrix) %in% genes]
      list(
        genes = genes,
        n_genes = length(genes),
        n_species = length(unique(sp_in_module)),
        species_table = table(sp_in_module)
      )
    })
  }

  structure(
    list(
      modules = all_modules,
      module_conditions = module_conditions,
      module_info = module_info,
      parameters = list(T1 = T1, T2 = T2,
                        min_clique_size = min_clique_size,
                        threshold_method = threshold_method,
                        network_type = network_type),
      n_modules = length(all_modules),
      species_mode = species_mode,
      transform = transform,
      data = list(expr_matrix = expr_matrix, condition = condition, species = species)
    ),
    class = "bmkc_result"
  )
}

# ===================================================================
# Threshold Exploration (Highly Optimized)
# ===================================================================

#' Explore BMKC Thresholds
#'
#' Helps choose optimal T1 and T2 values by calculating correlation matrices
#' exactly once, offering massive speedups over traditional grid searches.
#' @export
explore_bmkc_thresholds <- function(expr_matrix, condition,
                                    T1_range = seq(0.6, 0.9, by = 0.1),
                                    T2_range = seq(0.1, 0.5, by = 0.1),
                                    min_clique_size = 4,
                                    workers = 1,
                                    transform = c("none", "clr", "log1p"),
                                    species = NULL,
                                    network_type = c("signed", "unsigned"),
                                    threshold_method = c("absolute", "percentile")) {

  transform <- match.arg(transform)
  network_type <- match.arg(network_type)
  threshold_method <- match.arg(threshold_method)

  condition <- .validate_input(expr_matrix, condition)
  conditions <- levels(condition)

  message("Precomputing correlation matrices for fast threshold exploration...")
  expr_transformed <- .compositional_transform(expr_matrix, transform, species)
  cor_1 <- compute_cor_matrix(expr_transformed[, condition == conditions[1]], workers = workers)
  cor_2 <- compute_cor_matrix(expr_transformed[, condition == conditions[2]], workers = workers)
  precomp <- list(cor_1, cor_2)

  param_grid <- expand.grid(T1 = T1_range, T2 = T2_range)
  message(sprintf("Testing %d parameter combinations...", nrow(param_grid)))

  param_grid$n_modules <- sapply(seq_len(nrow(param_grid)), function(i) {
    tryCatch({
      # Suppress warnings during exploration grid loop to avoid console spam
      suppressWarnings({
        res <- run_bmkc(
          expr_matrix = expr_matrix,
          condition = condition,
          T1 = param_grid$T1[i],
          T2 = param_grid$T2[i],
          min_clique_size = min_clique_size,
          workers = 1, # Parallelism handled up front, no need here
          transform = transform,
          species = species,
          network_type = network_type,
          threshold_method = threshold_method,
          .precomputed_cors = precomp
        )
      })
      res$n_modules
    }, error = function(e) {
      NA_integer_
    })
  })

  param_grid
}

#' Plot BMKC Threshold Exploration
#' @export
plot_bmkc_exploration <- function(exploration_result) {
  # Drop NAs (combinations that were too dense and aborted)
  plot_df <- exploration_result[!is.na(exploration_result$n_modules), ]

  ggplot2::ggplot(plot_df, ggplot2::aes(x = T1, y = T2, fill = n_modules)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Number of\nModules") +
    ggplot2::labs(
      title = "BMKC Threshold Exploration",
      subtitle = "Number of detected modules across valid T1 and T2 thresholds",
      x = "T1: High correlation threshold (Active Condition)",
      y = "T2: Low correlation threshold (Inactive Condition)"
    ) +
    ggplot2::theme_minimal()
}

# ===================================================================
# S3 Methods
# ===================================================================

#' Summary for bmkc_result
#' @export
summary.bmkc_result <- function(object, ...) {
  cat("--- BMKC Analysis Summary ---\n")
  cat(sprintf("Modules found: %d\n", object$n_modules))
  cat(sprintf("Thresholds: T1 = %.2f | T2 = %.2f (%s method)\n",
              object$parameters$T1, object$parameters$T2, object$parameters$threshold_method))
  cat(sprintf("Network type: %s\n", toupper(object$parameters$network_type)))

  if (object$n_modules > 0) {
    cond_table <- table(object$module_conditions)
    cat("\nModules by Active Condition:\n")
    for (cond in names(cond_table)) {
      cat(sprintf("  - %s: %d modules\n", cond, cond_table[[cond]]))
    }
  }

  if (object$species_mode) {
    cat("\nSpecies mode: Enabled\n")
    multi_sp <- sum(sapply(object$module_info, function(x) x$n_species > 1))
    cat(sprintf("Multi-species modules: %d\n", multi_sp))
  }
  cat("----------------------------------\n")

  if (length(object$modules) > 0) {
    cat("Module sizes (top 10):\n")
    print(head(sapply(object$modules, length), 10))
  }
  invisible(object)
}

#' Print method
#' @export
print.bmkc_result <- function(x, ...) {
  summary(x, ...)
}
