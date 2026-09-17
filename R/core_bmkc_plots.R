#' @title BMKC Visualization Suite
#' @description Plotting functions for Biweight Midcorrelation K-Clique (BMKC) modules and their dynamics across conditions.
#' @name bmkc_plots
NULL

# Prevent CRAN package check notes regarding unquoted column name structures
utils::globalVariables(c("Condition", "Eigengene", "Var1", "Var2", "Freq", "edge_diff",
                         "weight", "name", "is_driver", "Membership", "DC_Score", "Gene",
                         "FDR", "X", "Y", "Activity"))

# ==============================================================================
# CATEGORY 1: UNIVERSAL DIAGNOSTICS
# ==============================================================================

#' Module Eigengene Boxplot (BMKC)
#' @param bmkc_result A `bmkc_result` object from [run_bmkc()].
#' @param module_index Which detected module to plot. Default `1`.
#' @param ... Currently unused; accepted for S3 compatibility.
#' @return A `ggplot` object.
#' @export
plot_module_eigengene <- function(bmkc_result, module_index = 1, ...) {
  if (module_index > bmkc_result$n_modules) stop("Module index out of bounds.")
  genes <- bmkc_result$modules[[module_index]]
  expr <- as.matrix(bmkc_result$data$expr_matrix[genes, , drop = FALSE])

  row_vars <- apply(expr, 1, stats::var)
  expr <- expr[row_vars > 1e-8, , drop = FALSE]
  if (nrow(expr) < 2) stop("Not enough non-zero variance genes to compute Eigengene.")

  eigengene <- stats::prcomp(t(expr), scale. = TRUE)$x[, 1]

  cond <- bmkc_result$data$condition
  if (mean(eigengene[cond == levels(cond)[1]]) < mean(eigengene[cond == levels(cond)[2]])) {
    eigengene <- -eigengene
  }

  df <- data.frame(Eigengene = eigengene, Condition = cond)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$Condition, y = .data$Eigengene, fill = .data$Condition)) +
    ggplot2::geom_boxplot(alpha = 0.8) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.6) +
    ggplot2::labs(title = paste("Module", module_index, "Relaxed Eigengene Profile"),
                  subtitle = paste("Summarizing", nrow(expr), "active components")) +
    ggplot2::theme_minimal(base_size = 14, base_family = "sans")
}

#' Plot a Module Correlation Heatmap
#'
#' Draws the within-module biweight midcorrelation structure for one
#' detected BMKC module.
#'
#' @param bmkc_result A `bmkc_result` object from [run_bmkc()].
#' @param module_index Which detected module to plot. Default `1`.
#' @param ... Currently unused; accepted for S3 compatibility.
#' @return A `ggplot` object.
#' @export
#' @export
plot_module_heatmap <- function(bmkc_result, module_index = 1, ...) {
  if (module_index > bmkc_result$n_modules) stop("Module index out of bounds.")
  genes <- bmkc_result$modules[[module_index]]
  expr_mat <- as.matrix(bmkc_result$data$expr_matrix[genes, , drop = FALSE])
  cond <- bmkc_result$data$condition
  cond_levels <- levels(cond)

  # PERFORMANCE FIX: Replaced direct base R cor iterations with hyper-fast native C++ functions
  cor1 <- cpp_calc_observed_bicor(expr_mat[, cond == cond_levels[1], drop = FALSE])
  cor2 <- cpp_calc_observed_bicor(expr_mat[, cond == cond_levels[2], drop = FALSE])

  rownames(cor1) <- colnames(cor1) <- rownames(cor2) <- colnames(cor2) <- genes

  df1 <- as.data.frame(as.table(cor1)); df1$Condition <- cond_levels[1]
  df2 <- as.data.frame(as.table(cor2)); df2$Condition <- cond_levels[2]
  df <- rbind(df1, df2)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::facet_wrap(~ Condition) +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "bicor", na.value = "white") +
    ggplot2::theme_minimal(base_family = "sans") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::labs(title = paste("Module", module_index, "Scale-Free Correlation Heatmap"))
}

# ==============================================================================
# CATEGORY 2: BULK RNA-SEQ & CLINICAL ARRAYS
# ==============================================================================

#' Differential Edge Network (BMKC/BMHT) - Direction & Single-Cell Aware
#'
#' Works with both `bmkc_result` and `bmht_result` objects. This used to be
#' defined twice (identically-named, divergent implementations in both
#' core_bmht_plots.R and core_bmkc_plots.R) - whichever file loaded second
#' silently won, and the two copies had drifted apart (one was missing
#' rownames/colnames on the recomputed correlation matrices, which corrupts
#' node labeling). This is now the single canonical definition.
#' @param result A `bmkc_result` or `bmht_result` object.
#' @param module_index Which detected module to plot (BMKC only).
#'   Default `1`.
#' @param n_top Maximum number of nodes to display. Default `20`.
#' @param min_diff Minimum absolute correlation difference for an edge to
#'   be drawn. Default `0.3`.
#' @param interactive If `TRUE`, render an interactive `visNetwork` plot
#'   instead of a static `ggraph` one. Default `FALSE`.
#' @param sc_obj Optional Seurat object, used only for node annotation.
#' @param ... Currently unused; accepted for S3 compatibility.
#' @return A `ggraph`/`ggplot` object, or a `visNetwork` widget when
#'   `interactive = TRUE`.
#' @export
plot_differential_network <- function(result, module_index = 1, n_top = 20, min_diff = 0.3, interactive = FALSE, sc_obj = NULL, ...) {
  driver_genes <- c()

  if (inherits(result, "bmkc_result")) {
    if (module_index > result$n_modules) stop("Module index out of bounds.")
    genes <- result$modules[[module_index]]
    title <- paste("Differential Network - Module", module_index)

    expr_sub <- as.matrix(result$data$expr_matrix[genes, , drop = FALSE])
    cond <- result$data$condition

    # PERFORMANCE FIX: Replaced slow base loops with optimized native C++ solvers
    c1 <- cpp_calc_observed_bicor(expr_sub[, cond == levels(cond)[1], drop = FALSE])
    c2 <- cpp_calc_observed_bicor(expr_sub[, cond == levels(cond)[2], drop = FALSE])

    rownames(c1) <- colnames(c1) <- rownames(c2) <- colnames(c2) <- genes

    hub_scores <- rowSums(abs(c2) - abs(c1), na.rm = TRUE)
    if (length(hub_scores) > 0 && max(hub_scores) > 0) driver_genes <- names(which.max(hub_scores))

  } else if (inherits(result, "bmht_result")) {
    genes <- head(result$results$Gene, n_top)
    title <- paste("Top", n_top, "Rewired Genes Network")
    if (!is.null(result$causal_drivers)) {
      driver_genes <- intersect(result$causal_drivers$Gene, genes)
    }
  } else stop("Unsupported input. Requires bmkc_result or bmht_result.")

  expr_mat <- as.matrix(result$data$expr_matrix[genes, , drop = FALSE])
  cond <- result$data$condition
  cond_levels <- levels(cond)

  cor1 <- cpp_calc_observed_bicor(expr_mat[, cond == cond_levels[1], drop = FALSE])
  cor2 <- cpp_calc_observed_bicor(expr_mat[, cond == cond_levels[2], drop = FALSE])

  rownames(cor1) <- colnames(cor1) <- rownames(cor2) <- colnames(cor2) <- genes

  adj_diff <- cor2 - cor1
  adj_diff[abs(adj_diff) < min_diff | is.na(adj_diff)] <- 0

  is_directed <- length(driver_genes) > 0
  g <- igraph::graph_from_adjacency_matrix(adj_diff, mode = ifelse(is_directed, "directed", "undirected"), weighted = TRUE, diag = FALSE)

  if (length(igraph::E(g)) == 0) stop("No structural edges passed the min_diff filtration limits.")

  if (is_directed) {
    edge_list <- igraph::as_data_frame(g, what = "edges")
    to_swap <- edge_list$to %in% driver_genes & !(edge_list$from %in% driver_genes)
    temp <- edge_list$from[to_swap]
    edge_list$from[to_swap] <- edge_list$to[to_swap]
    edge_list$to[to_swap] <- temp
    g <- igraph::graph_from_data_frame(edge_list, directed = TRUE)
  }

  igraph::E(g)$edge_diff <- igraph::E(g)$weight
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  igraph::V(g)$is_driver <- igraph::V(g)$name %in% driver_genes

  node_colors_var <- "is_driver"
  node_sizes_var <- "is_driver"
  using_sc_aesthetics <- FALSE

  if (!is.null(sc_obj)) {
    valid_sc_genes <- intersect(igraph::V(g)$name, rownames(sc_obj))
    if (length(valid_sc_genes) > 0) {
      avg_expr <- Seurat::AggregateExpression(sc_obj, features = valid_sc_genes, group.by = "ident", slot = "data", verbose = FALSE)[[1]]

      avg_expr_mat <- as.matrix(avg_expr)
      dominant_group <- colnames(avg_expr_mat)[apply(avg_expr_mat, 1, which.max)]

      igraph::V(g)$Specificity <- dominant_group
      igraph::V(g)$MeanExpression <- base::rowMeans(avg_expr_mat)

      node_colors_var <- "Specificity"
      node_sizes_var <- "MeanExpression"
      using_sc_aesthetics <- TRUE
    }
  }

  if (interactive) {
    if (!requireNamespace("visNetwork", quietly = TRUE)) stop("Please install 'visNetwork'.")
    vis_data <- igraph::as_data_frame(g, what = "both")

    node_colors <- ifelse(vis_data$vertices$name %in% driver_genes, "#ff7f00", "#97c2fc")
    nodes <- data.frame(id = vis_data$vertices$name, label = vis_data$vertices$name, title = vis_data$vertices$name, color = node_colors)

    edges <- data.frame(from = vis_data$edges$from, to = vis_data$edges$to,
                        value = vis_data$edges$weight,
                        color = ifelse(vis_data$edges$edge_diff > 0, "#e31a1c", "#1f78b4"))

    if (is_directed) edges$arrows <- "to"

    vn <- visNetwork::visNetwork(nodes, edges, main = title)
    vn <- visNetwork::visIgraphLayout(vn)
    return(visNetwork::visOptions(vn, highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE))
  }

  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Please install 'ggraph' for the static (non-interactive) network plot, ",
         "or set interactive = TRUE to use visNetwork instead.")
  }

  p <- ggraph::ggraph(g, layout = ifelse(length(genes) <= 5, "stress", "fr"))

  if (is_directed) {
    p <- p + ggraph::geom_edge_link(ggplot2::aes(color = .data$edge_diff, width = .data$weight), alpha = 0.75,
                                    arrow = ggplot2::arrow(length = ggplot2::unit(2.5, 'mm')),
                                    end_cap = ggraph::circle(4, 'mm'))
  } else {
    p <- p + ggraph::geom_edge_link(ggplot2::aes(color = .data$edge_diff, width = .data$weight), alpha = 0.75)
  }

  if (using_sc_aesthetics) {
    p <- p + ggraph::geom_node_point(ggplot2::aes(color = .data[[node_colors_var]], size = .data[[node_sizes_var]])) +
      ggplot2::scale_color_brewer(palette = "Set1", name = "Cell Specificity") +
      ggplot2::scale_size_continuous(range = c(5, 10), name = "Expression Vol")
  } else {
    p <- p + ggraph::geom_node_point(ggplot2::aes(color = .data$is_driver, size = .data$is_driver)) +
      ggplot2::scale_color_manual(values = c("TRUE" = "#ff7f00", "FALSE" = "steelblue"), guide = "none") +
      ggplot2::scale_size_manual(values = c("TRUE" = 9, "FALSE" = 6), guide = "none")
  }

  p <- p + ggraph::geom_node_text(ggplot2::aes(label = .data$name), repel = TRUE, size = 4.5,
                                  fontface = ifelse(igraph::V(g)$is_driver, "bold", "plain"),
                                  bg.color = "white", bg.r = 0.15) +
    ggraph::scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "\u0394 Correlation") +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::labs(title = title, subtitle = ifelse(is_directed, "Arrows indicate predicted TF master regulation", ""))

  return(p)
}

# ==============================================================================
# CATEGORY 4: SINGLE-CELL RNA-SEQ (scRNA-seq)
# ==============================================================================

#' Module Eigengene Projected onto a UMAP Embedding
#' @param bmkc_res A `bmkc_result` object from [run_bmkc()].
#' @param sc_obj A normalized Seurat object containing a UMAP (or other)
#'   dimensionality reduction.
#' @param module_index Which detected module to score. Default `1`.
#' @param reduction Name of the dimensionality reduction to plot.
#'   Default `"umap"`.
#' @param pt_size Point size. Default `0.5`.
#' @param split_by Optional metadata column to facet the plot by.
#' @param ... Currently unused; accepted for S3 compatibility.
#' @return A `ggplot` object.
#' @export
plot_umap_eigengene <- function(bmkc_res, sc_obj, module_index = 1, reduction = "umap", pt_size = 0.5, split_by = NULL, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("The 'Seurat' package is required.")

  if (is.null(bmkc_res$modules) || length(bmkc_res$modules) == 0) stop("No modules found.")
  if (module_index > length(bmkc_res$modules)) stop("Module index out of bounds.")
  clique_genes <- bmkc_res$modules[[module_index]]

  valid_genes <- intersect(clique_genes, rownames(sc_obj))
  if (length(valid_genes) < 3) warning("Few genes mapped to the Seurat object. Activity may be noisy.")
  if (length(valid_genes) == 0) stop("No network genes found in Seurat object.")

  score_name <- paste0("BMKC_Mod", module_index, "_")
  sc_obj <- Seurat::AddModuleScore(sc_obj, features = list(valid_genes), name = score_name)
  target_col <- paste0(score_name, "1")
  active_cond <- bmkc_res$module_conditions[module_index]

  Seurat::FeaturePlot(sc_obj,
                      features = target_col,
                      reduction = reduction,
                      pt.size = pt_size,
                      order = TRUE,
                      split.by = split_by) +
    ggplot2::scale_color_gradientn(colors = c("#f0f0f0", "#fddbc7", "#f4a582", "#d6604d", "#b2182b")) +
    ggplot2::ggtitle(sprintf("BMKC Module %d Activation Mapping", module_index)) +
    ggplot2::labs(subtitle = sprintf("Network quasi-clique signature track | Active state target: %s", active_cond),
                  caption = sprintf("Aggregated from %d structurally validated co-expression nodes", length(valid_genes)),
                  color = "Activity Score") +
    ggplot2::theme_minimal(base_family = "sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = 10),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "right")
}

# ==============================================================================
# CORE S3 DISPATCHER
# ==============================================================================

#' Plot Method for BMKC
#'
#' As of the Phase 3 surface reduction, only four plot types are provided
#' directly: `"heatmap"`, `"eigengene"`, `"network"` (default; a differential
#' network shared with BMHT results), and `"umap"` (a module-eigengene
#' overlay, requires a Seurat object via `sc_obj`). For other visualizations
#' (raincloud style, module violin, spatial/tissue overlays, multi-method
#' overlap), see the "Custom plots from bicorX results" vignette.
#' @method plot bmkc_result
#' @param x A `bmkc_result` object from [run_bmkc()].
#' @param type Which plot to draw: `"default"`/`"network"` for the
#'   differential network, `"heatmap"` for the module correlation
#'   heatmap, or `"eigengene"` for the module eigengene plot.
#' @param ... Passed to the underlying plotting function.
#' @return A `ggplot` object (or a `visNetwork` widget for interactive
#'   network plots).
#' @export
plot.bmkc_result = function(x, type = c("default", "network", "heatmap",
                                        "eigengene", "umap"), ...) {
  type <- match.arg(type)

  switch(type,
         "heatmap"   = plot_module_heatmap(x, ...),
         "eigengene" = plot_module_eigengene(x, ...),
         "default"   = plot_differential_network(x, ...),
         "network"   = plot_differential_network(x, ...),
         "umap"      = plot_umap_eigengene(x, ...)
  )
}
