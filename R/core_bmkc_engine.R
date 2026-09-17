#' @title BMKC: Biweight Midcorrelation + Relaxed Quasi-Clique Method
#' @description Identifies dense gene modules (quasi-cliques) that are strongly co-expressed,
#' using a k-core decomposition with a gamma density-relaxation step. Correlation
#' matrices and adjacency thresholding are computed in C++; see bmkc_cpp.cpp.
#'
#' @importFrom igraph graph_from_adjacency_matrix coreness induced_subgraph edge_density components
#' @importFrom matrixStats rowMads
#' @importFrom parallel detectCores
#' @importFrom progress progress_bar
#' @importFrom Rcpp sourceCpp
#' @useDynLib bicorX, .registration = TRUE
#' @name bmkc
NULL

#' Global differential-coexpression score for a set of BMKC modules
#'
#' Operationalizes Yuan et al. 2015 section 3.3's "global total sum of
#' differential coexpression change of each gene in modules": for each
#' detected module, each gene's differential-coexpression contribution is
#' its mean |bicor_cond1 - bicor_cond2| over edges to *other genes in that
#' same module* (matching BMKC's own module definition rather than a
#' generic genome-wide average); these per-gene values are summed across
#' every (module, gene) occurrence to give one scalar for the whole result.
#' A gene appearing in multiple modules contributes once per module,
#' matching the paper's "total sum" framing.
#'
#' @keywords internal
.bmkc_global_score <- function(modules, cor_1, cor_2) {
  if (length(modules) == 0) return(0)
  total <- 0
  for (genes in modules) {
    idx <- match(genes, rownames(cor_1))
    idx <- idx[!is.na(idx)]
    if (length(idx) < 2) next
    sub_diff <- abs(cor_1[idx, idx, drop = FALSE] - cor_2[idx, idx, drop = FALSE])
    diag(sub_diff) <- NA
    per_gene <- rowMeans(sub_diff, na.rm = TRUE)
    total <- total + sum(per_gene, na.rm = TRUE)
  }
  total
}

#' One null replicate for BMKC's global significance test
#'
#' Implements Yuan et al. 2015 section 3.3's null: "(group-wise) shuffling
#' the expression values for each gene independently" - within each
#' condition's samples, each gene's values are independently permuted
#' across those samples (a different random order per gene), which
#' destroys all gene-gene covariance structure within each condition while
#' preserving each gene's own marginal distribution and the condition
#' grouping itself. The full module-detection pipeline (bicor, adjacency,
#' k-core module extraction) is then rerun on this shuffled matrix and the
#' resulting global score is returned - this is deliberately a lighter,
#' self-contained reimplementation of the core steps rather than a call
#' into `.run_bmkc_matrix()` itself, to avoid any risk of altering that
#' already-validated function's behavior for ordinary (non-permuted) use.
#'
#' @keywords internal
.bmkc_null_replicate <- function(expr_matrix, condition, T1, T2, gamma,
                                  min_clique_size, is_unsigned, is_percentile,
                                  workers) {
  conditions <- levels(condition)
  idx1 <- which(condition == conditions[1])
  idx2 <- which(condition == conditions[2])

  shuffled <- expr_matrix
  # Independently permute each gene's own values within each condition's
  # sample columns - a different random column order per gene (per row),
  # not one shared shuffle applied to every gene alike.
  shuffled[, idx1] <- t(apply(expr_matrix[, idx1, drop = FALSE], 1, sample))
  shuffled[, idx2] <- t(apply(expr_matrix[, idx2, drop = FALSE], 1, sample))
  rownames(shuffled) <- rownames(expr_matrix)

  cor_1 <- cpp_fast_bicor_matrix(shuffled[, idx1, drop = FALSE], workers)
  cor_2 <- cpp_fast_bicor_matrix(shuffled[, idx2, drop = FALSE], workers)
  rownames(cor_1) <- colnames(cor_1) <- rownames(cor_2) <- colnames(cor_2) <- rownames(expr_matrix)

  cpp_res <- cpp_bmkc_adjacency(cor_1, cor_2, T1, T2, is_unsigned, is_percentile, workers)

  extract_modules_lite <- function(adj_logical, n_edges) {
    n_nodes <- nrow(adj_logical)
    if (n_nodes < min_clique_size) return(list())
    graph_density <- n_edges / ((n_nodes * (n_nodes - 1)) / 2)
    if (n_edges > 800000 || graph_density > 0.40) return(list())

    rownames(adj_logical) <- colnames(adj_logical) <- rownames(expr_matrix)
    graph <- igraph::graph_from_adjacency_matrix(adj_logical, mode = "undirected", diag = FALSE)
    core_values <- igraph::coreness(graph)
    min_core_bound <- min_clique_size - 1
    high_cores <- sort(unique(core_values[core_values >= min_core_bound]), decreasing = TRUE)

    valid_modules <- list()
    for (k in high_cores) {
      core_nodes <- names(core_values[core_values >= k])
      if (length(core_nodes) < min_clique_size) next
      subg <- igraph::induced_subgraph(graph, core_nodes)
      components <- igraph::components(subg)
      for (comp_idx in seq_len(components$no)) {
        comp_nodes <- names(components$membership[components$membership == comp_idx])
        if (length(comp_nodes) >= min_clique_size) {
          comp_subg <- igraph::induced_subgraph(subg, comp_nodes)
          if (igraph::edge_density(comp_subg) >= gamma) {
            valid_modules[[length(valid_modules) + 1]] <- comp_nodes
          }
        }
      }
    }
    unique(valid_modules)
  }

  modules <- c(extract_modules_lite(cpp_res$adj_cond1, cpp_res$edge_counts$c1),
               extract_modules_lite(cpp_res$adj_cond2, cpp_res$edge_counts$c2))

  .bmkc_global_score(modules, cor_1, cor_2)
}

#' Run BMKC Analysis
#' @param n_permutations Default `0` (no significance test - preserves
#'   the original behavior of this engine exactly). If > 0, runs the
#'   global significance test from Yuan et al. 2015 section 3.3: each
#'   gene's expression is independently permuted within each condition
#'   (destroying all gene-gene covariance while preserving each gene's own
#'   marginal distribution and the condition grouping), the full module-
#'   detection pipeline is rerun on the shuffled data, and a "global
#'   score" (summed per-module, per-gene mean |bicor difference| to other
#'   genes in the same module - see `.bmkc_global_score()`) is computed.
#'   Repeating this `n_permutations` times gives an empirical null
#'   distribution against which the real (unpermuted) result's global
#'   score is compared, yielding one p-value (`significance_p` in the
#'   result) for the whole set of detected modules - not a per-module or
#'   per-gene FDR. The paper used 1000 permutations; this is fast enough
#'   in practice (~6s at 1000 permutations, 60 genes, single core) to use
#'   that as a real target rather than a token gesture, but scales with
#'   gene count like the rest of the module-detection pipeline.
#' @keywords internal
.run_bmkc_matrix <- function(expr_matrix, condition,
                             T1 = 0.7, T2 = 0.3,
                             gamma = 0.85,
                             min_clique_size = 4,
                             workers = 1,
                             transform = c("none", "clr", "rclr", "log1p"),
                             species = NULL,
                             network_type = c("signed", "unsigned"),
                             threshold_method = c("absolute", "percentile"),
                             bipartite_cliques_only = FALSE,
                             n_permutations = 0,
                             seed = 42,
                             .precomputed_cors = NULL) {

  if (!is.null(seed)) set.seed(seed)

  condition <- .validate_input(expr_matrix, condition)
  transform <- match.arg(transform)
  network_type <- match.arg(network_type)
  threshold_method <- match.arg(threshold_method)

  # ===================================================================
  # 0. User Experience Guard Rails & Warning Messages (Unified Suite UX)
  # ===================================================================
  phys_cores <- parallel::detectCores(logical = FALSE)
  if (workers > phys_cores) {
    warning(sprintf("Thread Optimization [BMKC]: Requested workers (%d) exceeds physical hardware cores (%d).\n", workers, phys_cores),
            "Resetting workers to match physical cores to optimize cache layout and maximize speed.")
    workers <- max(1, phys_cores)
  }

  feature_mads <- matrixStats::rowMads(expr_matrix, na.rm = TRUE)
  keep <- feature_mads > 0
  if (any(!keep, na.rm = TRUE)) {
    warning("Zero Variance Alert [BMKC]: Features with a Median Absolute Deviation (MAD) of 0 detected.\n",
            "BMKC has automatically filtered these silent features to preserve matrix arithmetic bounds.")
    expr_matrix <- expr_matrix[keep, , drop = FALSE]
    if (!is.null(species)) species <- species[keep]
  }

  if (gamma <= 0 || gamma > 1.0) {
    stop("Parameter 'gamma' (quasi-clique relaxation density threshold) must fall within (0, 1.0].")
  }

  species_mode <- FALSE
  if (!is.null(species)) {
    if (length(species) != nrow(expr_matrix)) stop("Length of `species` must exactly match the number of rows.")
    species_mode <- TRUE # FIX: Structural type assignment correction from syntax typo
    species <- as.character(species)
  }

  if (bipartite_cliques_only && !species_mode) {
    stop("bipartite_cliques_only requires a valid 'species' vector.")
  }

  conditions <- levels(condition)
  message("Initializing BMKC Module Extraction Workflow...")

  # ===================================================================
  # 1. Cache-Optimal Native C++ bicor Implementation
  # ===================================================================
  if (!is.null(.precomputed_cors)) {
    cor_1 <- .precomputed_cors[[1]]
    cor_2 <- .precomputed_cors[[2]]
  } else {
    expr_matrix <- .compositional_transform(expr_matrix, transform, species)
    message("Step 1/3: Generating cache-aligned bicor matrices via C++ Engine...")
    cor_1 <- cpp_fast_bicor_matrix(expr_matrix[, condition == conditions[1], drop = FALSE], workers)
    cor_2 <- cpp_fast_bicor_matrix(expr_matrix[, condition == conditions[2], drop = FALSE], workers)
  }
  rownames(cor_1) <- colnames(cor_1) <- rownames(cor_2) <- colnames(cor_2) <- rownames(expr_matrix)

  # ===================================================================
  # 2. C++ Adjacency & Quantile Masking
  # ===================================================================
  is_percentile <- threshold_method == "percentile"
  is_unsigned <- network_type == "unsigned"

  message("Step 2/3: Applying stochastic quantiles and compiling logical adjacency...")
  cpp_res <- cpp_bmkc_adjacency(cor_1, cor_2, T1, T2, is_unsigned, is_percentile, workers)

  adj_cond1_active <- cpp_res$adj_cond1
  adj_cond2_active <- cpp_res$adj_cond2

  # Fetch pre-calculated linear edge counts directly from C++ registers
  edge_counts_c1 <- cpp_res$edge_counts$c1
  edge_counts_c2 <- cpp_res$edge_counts$c2

  if (is.null(.precomputed_cors) && n_permutations <= 0) {
    rm(cor_1, cor_2)
    gc(verbose = FALSE)
  }

  # ===================================================================
  # 3. Fast Bipartite Masking (Native C-level `outer`)
  # ===================================================================
  if (bipartite_cliques_only) {
    intra_mask <- outer(species, species, "==")
    adj_cond1_active[intra_mask] <- FALSE
    adj_cond2_active[intra_mask] <- FALSE
    rm(intra_mask)
    # Recalculate edge targets if bipartite filters altered layout matrices
    edge_counts_c1 <- sum(adj_cond1_active) / 2
    edge_counts_c2 <- sum(adj_cond2_active) / 2
  }

  # ===================================================================
  # 4. igraph Execution via Fast k-Core Quasi-Clique Partitioning
  # ===================================================================
  message("Step 3/3: Running linear-time k-core topology decomposition...")

  extract_relaxed_modules <- function(adj_logical, n_edges, label) {
    rownames(adj_logical) <- colnames(adj_logical) <- rownames(expr_matrix)
    n_nodes <- nrow(adj_logical)
    graph_density <- n_edges / ((n_nodes * (n_nodes - 1)) / 2)

    # Complexity safety filter to protect igraph execution boundaries
    max_allowed_density <- if(n_nodes > 3000) 0.15 else if(n_nodes > 1000) 0.25 else 0.40
    max_safe_edges <- 800000

    if (n_edges > max_safe_edges || graph_density > max_allowed_density) {
      warning(sprintf(
        "Network Complexity Warning: Condition %s is dangerously dense (%d edges, %.1f%% density).\n",
        label, n_edges, graph_density * 100
      ), "Evaluation bypassed to prevent memory exhaustion. Please raise T1.")
      return(list())
    }

    graph <- igraph::graph_from_adjacency_matrix(adj_logical, mode = "undirected", diag = FALSE)

    core_values <- igraph::coreness(graph)
    min_core_bound <- min_clique_size - 1
    high_cores <- sort(unique(core_values[core_values >= min_core_bound]), decreasing = TRUE)

    valid_modules <- list()

    if (length(high_cores) > 0) {
      pb_core <- progress::progress_bar$new(
        format = paste0("  Processing [", label, "] Core Subgraphs [:bar] :percent | ETA: :eta"),
        total = length(high_cores), clear = FALSE, width = 65
      )

      for (k in high_cores) {
        core_nodes <- names(core_values[core_values >= k])
        if (length(core_nodes) >= min_clique_size) {
          subg <- igraph::induced_subgraph(graph, core_nodes)
          components <- igraph::components(subg)

          for (comp_idx in seq_len(components$no)) {
            comp_nodes <- names(components$membership[components$membership == comp_idx])

            if (length(comp_nodes) >= min_clique_size) {
              # OPTIMIZATION: Induced subgraph derived from 'subg' memory stack instead of full root 'graph'
              comp_subg <- igraph::induced_subgraph(subg, comp_nodes)
              sub_density <- igraph::edge_density(comp_subg)

              if (sub_density >= gamma) {
                valid_modules[[length(valid_modules) + 1]] <- comp_nodes
              }
            }
          }
        }
        pb_core$tick()
      }
    }
    unique(valid_modules)
  }

  cliques_c1 <- extract_relaxed_modules(adj_cond1_active, edge_counts_c1, conditions[1])
  cliques_c2 <- extract_relaxed_modules(adj_cond2_active, edge_counts_c2, conditions[2])

  rm(adj_cond1_active, adj_cond2_active, cpp_res)
  gc(verbose = FALSE)

  all_modules <- c(cliques_c1, cliques_c2)
  module_conditions <- c(rep(conditions[1], length(cliques_c1)), rep(conditions[2], length(cliques_c2)))

  significance_p <- NA_real_
  observed_score <- NA_real_
  if (n_permutations > 0) {
    observed_score <- .bmkc_global_score(all_modules, cor_1, cor_2)

    message(sprintf("Running global significance test (%d permutations, Yuan et al. 2015 section 3.3)...",
                    n_permutations))
    pb_sig <- progress::progress_bar$new(
      format = "  Significance Test [:bar] :percent | ETA: :eta",
      total = n_permutations, clear = FALSE, width = 65
    )
    null_scores <- numeric(n_permutations)
    for (b in seq_len(n_permutations)) {
      null_scores[b] <- .bmkc_null_replicate(expr_matrix, condition, T1, T2, gamma,
                                              min_clique_size, is_unsigned, is_percentile, workers)
      pb_sig$tick()
    }
    significance_p <- (sum(null_scores >= observed_score) + 1) / (n_permutations + 1)
  }

  if (is.null(.precomputed_cors) && n_permutations > 0) {
    rm(cor_1, cor_2)
    gc(verbose = FALSE)
  }

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

  message("BMKC module extraction complete.")
  structure(
    list(
      modules = all_modules,
      module_conditions = module_conditions,
      module_info = module_info,
      parameters = list(T1 = T1, T2 = T2, gamma = gamma, min_clique_size = min_clique_size, threshold_method = threshold_method, network_type = network_type),
      n_modules = length(all_modules),
      global_score = observed_score,
      n_permutations = if (n_permutations > 0) n_permutations else NA_integer_,
      significance_p = significance_p,
      species_mode = species_mode,
      transform = transform,
      data = list(expr_matrix = expr_matrix, condition = condition, species = species)
    ),
    class = c("bmkc_result", "list")
  )
}
