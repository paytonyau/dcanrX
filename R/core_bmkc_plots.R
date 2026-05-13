#' @title BMKC Visualization Suite
#' @description Plotting functions for Biweight Midcorrelation K-Cliques (BMKC) and Module Dynamics.
#' @name bmkc_plots
NULL

# ===================================================================
# ======================  CORE VISUALS  =============================
# ===================================================================

#' Differential Edge Network (BMKC/BMHT) - Direction-Aware
#' @description Plots the physical network of a module or top hubs, highlighting changed edges.
#' @export
plot_differential_network <- function(result, module_index = 1, n_top = 20, min_diff = 0.3, interactive = FALSE) {
  driver_genes <- c()
  if (inherits(result, "bmkc_result")) {
    if (module_index > result$n_modules) stop("Module index out of bounds.")
    genes <- result$modules[[module_index]]
    title <- paste("Differential Network - Module", module_index)

    if (!is.null(result$module_drivers)) {
      drivers <- result$module_drivers[[paste0("Module_", module_index)]]
      if (drivers[1] != "No known TF driver") driver_genes <- drivers
    }
  } else if (inherits(result, "bmht_result")) {
    genes <- head(result$results$Gene, n_top)
    title <- paste("Top", n_top, "Rewired Genes Network")

    if (!is.null(result$causal_drivers)) {
      driver_genes <- intersect(result$causal_drivers$Gene, genes)
    }
  } else stop("Unsupported input. Requires bmkc_result or bmht_result.")

  expr_mat <- result$data$expr_matrix[genes, ]
  cond <- result$data$condition
  cond_levels <- levels(cond)

  cor1 <- compute_cor_matrix(expr_mat[, cond == cond_levels[1]], workers = 1)
  cor2 <- compute_cor_matrix(expr_mat[, cond == cond_levels[2]], workers = 1)

  adj_diff <- cor2 - cor1
  adj_diff[abs(adj_diff) < min_diff] <- 0

  is_directed <- length(driver_genes) > 0
  g <- igraph::graph_from_adjacency_matrix(adj_diff, mode = ifelse(is_directed, "directed", "undirected"), weighted = TRUE, diag = FALSE)

  if (length(igraph::E(g)) == 0) stop("No edges passed the min_diff threshold.")

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

  if (interactive) {
    if (!requireNamespace("visNetwork", quietly = TRUE)) stop("Please install the 'visNetwork' package for interactive plots.")
    vis_data <- igraph::as_data_frame(g, what = "both")

    node_colors <- ifelse(vis_data$vertices$name %in% driver_genes, "#ff7f00", "#97c2fc")
    nodes <- data.frame(id = vis_data$vertices$name, label = vis_data$vertices$name, title = vis_data$vertices$name, color = node_colors)

    edges <- data.frame(from = vis_data$edges$from, to = vis_data$edges$to,
                        value = vis_data$edges$weight,
                        color = ifelse(vis_data$edges$edge_diff > 0, "#e31a1c", "#1f78b4"),
                        arrows = ifelse(is_directed, "to", ""))

    return(visNetwork::visNetwork(nodes, edges, main = title) %>%
             visNetwork::visIgraphLayout() %>%
             visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE))
  }

  p <- ggraph::ggraph(g, layout = "fr")

  if (is_directed) {
    p <- p + ggraph::geom_edge_link(ggplot2::aes(color = edge_diff, width = weight), alpha = 0.75,
                                    arrow = ggplot2::arrow(length = ggplot2::unit(3, 'mm')),
                                    end_cap = ggraph::circle(4, 'mm'))
  } else {
    p <- p + ggraph::geom_edge_link(ggplot2::aes(color = edge_diff, width = weight), alpha = 0.75)
  }

  p + ggraph::geom_node_point(ggplot2::aes(color = is_driver, size = is_driver)) +
    ggplot2::scale_color_manual(values = c("TRUE" = "#ff7f00", "FALSE" = "steelblue"), guide = "none") +
    ggplot2::scale_size_manual(values = c("TRUE" = 9, "FALSE" = 6), guide = "none") +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 4.5, fontface = ifelse(igraph::V(g)$is_driver, "bold", "plain")) +
    ggraph::scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Δ Correlation") +
    ggraph::theme_graph() +
    ggplot2::labs(title = title, subtitle = ifelse(is_directed, "Arrows indicate predicted TF master regulation", ""))
}

#' Module Eigengene Boxplot (BMKC)
#' @description Plots the summary expression (First Principal Component) of the entire module.
#' @export
plot_module_eigengene <- function(bmkc_result, module_index = 1) {
  if (module_index > bmkc_result$n_modules) stop("Module index out of bounds.")
  genes <- bmkc_result$modules[[module_index]]
  expr <- bmkc_result$data$expr_matrix[genes, ]

  row_vars <- apply(expr, 1, stats::var)
  expr <- expr[row_vars > 0, , drop = FALSE]
  if (nrow(expr) < 2) stop("Not enough non-zero variance genes to compute Eigengene.")

  eigengene <- stats::prcomp(t(expr), scale. = TRUE)$x[, 1]
  df <- data.frame(Eigengene = eigengene, Condition = bmkc_result$data$condition)

  ggplot2::ggplot(df, ggplot2::aes(x = Condition, y = Eigengene, fill = Condition)) +
    ggplot2::geom_boxplot(alpha = 0.8) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.6) +
    ggplot2::labs(title = paste("Module", module_index, "Eigengene Profile"),
                  subtitle = paste("Summarizing", nrow(expr), "active genes")) +
    ggplot2::theme_minimal()
}

#' Before-After Module Correlation Heatmap (BMKC)
#' @description Faceted heatmap proving that the module is synchronized in one condition but not the other.
#' @export
plot_module_heatmap <- function(bmkc_result, module_index = 1) {
  if (module_index > bmkc_result$n_modules) stop("Module index out of bounds.")
  genes <- bmkc_result$modules[[module_index]]
  expr_mat <- bmkc_result$data$expr_matrix[genes, ]
  cond <- bmkc_result$data$condition
  cond_levels <- levels(cond)

  cor1 <- compute_cor_matrix(expr_mat[, cond == cond_levels[1]], workers = 1)
  cor2 <- compute_cor_matrix(expr_mat[, cond == cond_levels[2]], workers = 1)

  df1 <- as.data.frame(as.table(cor1)); df1$Condition <- cond_levels[1]
  df2 <- as.data.frame(as.table(cor2)); df2$Condition <- cond_levels[2]
  df <- rbind(df1, df2)

  ggplot2::ggplot(df, ggplot2::aes(x = Var1, y = Var2, fill = Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::facet_wrap(~ Condition) +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "bicor") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::labs(title = paste("Module", module_index, "Correlation Heatmap"))
}

#' Membership vs Rewiring Scatter (Multi-method correlation)
#' @description Compares how tightly a gene belongs to a module vs how heavily it was rewired.
#' @export
plot_membership_vs_rewiring <- function(bmkc_result, bmht_result, module_index = 1) {
  mod_genes <- bmkc_result$modules[[module_index]]
  expr <- bmkc_result$data$expr_matrix[mod_genes, ]

  row_vars <- apply(expr, 1, stats::var)
  expr <- expr[row_vars > 0, , drop = FALSE]
  mod_genes <- rownames(expr)

  eigengene <- stats::prcomp(t(expr), scale. = TRUE)$x[,1]
  membership <- abs(stats::cor(t(expr), eigengene))

  df <- bmht_result$results[bmht_result$results$Gene %in% mod_genes, ]
  df$Membership <- membership[match(df$Gene, mod_genes)]

  ggplot2::ggplot(df, ggplot2::aes(x = Membership, y = DC_Score, label = Gene)) +
    ggplot2::geom_point(ggplot2::aes(color = FDR < bmht_result$significance_level), size = 4) +
    ggrepel::geom_text_repel(size = 3) +
    ggplot2::labs(title = paste("Module", module_index, "- Membership vs Rewiring Strength"),
                  x = "Module Membership (Correlation to Eigengene)", y = "Differential Co-expression Score",
                  color = "Significant Rewiring") +
    ggplot2::theme_minimal()
}

# ===================================================================
# ======================  ADVANCED VISUALS  =========================
# ===================================================================

#' Plot Module Eigengene (Raincloud Style)
#' @description High-impact density ridge + boxplot hybrid.
#' @export
plot_raincloud_eigengene <- function(bmkc_result, module_index = 1) {
  if (module_index > bmkc_result$n_modules) stop("Module index out of bounds.")
  if (!requireNamespace("ggdist", quietly = TRUE)) stop("Please install 'ggdist' for raincloud plots.")

  genes <- bmkc_result$modules[[module_index]]
  expr_mat <- bmkc_result$data$expr_matrix[genes, ]
  cond <- as.factor(bmkc_result$data$condition)

  pca <- stats::prcomp(t(expr_mat), center = TRUE, scale. = TRUE)
  eigengene <- pca$x[, 1]

  if (mean(eigengene[cond == levels(cond)[1]]) < mean(eigengene[cond == levels(cond)[2]])) {
    eigengene <- -eigengene
  }

  df <- data.frame(Condition = cond, Eigengene = eigengene)

  ggplot2::ggplot(df, ggplot2::aes(x = Condition, y = Eigengene, fill = Condition, color = Condition)) +
    ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.2, point_colour = NA, alpha = 0.7) +
    ggplot2::geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_point(size = 1.5, alpha = 0.4, position = ggplot2::position_jitter(width = .05, height = 0)) +
    ggplot2::scale_fill_manual(values = c("#1f78b4", "#e31a1c")) +
    ggplot2::scale_color_manual(values = c("#1f78b4", "#e31a1c")) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(title = paste("Module", module_index, "Raincloud Eigengene"),
                  x = "Condition", y = "First Principal Component (Module Expression)")
}

#' Spatial UMAP/Tissue Eigengene Projection (BMKC Hook)
#' @description Maps module activity onto specific physical tissue coordinates.
#' @export
plot_eigengene_spatial <- function(bmkc_result, spatial_coords, module_index = 1) {
  if (module_index > bmkc_result$n_modules) stop("Module index out of bounds.")
  genes <- bmkc_result$modules[[module_index]]
  expr <- bmkc_result$data$expr_matrix[genes, ]

  if (ncol(expr) != nrow(spatial_coords)) stop("Number of samples in bmkc_result must match rows in spatial_coords.")

  row_vars <- apply(expr, 1, stats::var)
  expr <- expr[row_vars > 0, , drop = FALSE]
  if (nrow(expr) < 2) stop("Not enough non-zero variance genes to compute Eigengene.")

  eigengene <- stats::prcomp(t(expr), scale. = TRUE)$x[, 1]
  eigengene <- scale(eigengene)[,1]

  df <- data.frame(X = spatial_coords[,1], Y = spatial_coords[,2], Eigengene = eigengene)

  ggplot2::ggplot(df, ggplot2::aes(x = X, y = Y, color = Eigengene)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0) +
    ggplot2::labs(title = paste("Spatial Eigengene Projection - Module", module_index),
                  x = "Coordinate 1 (UMAP/Spatial)", y = "Coordinate 2 (UMAP/Spatial)",
                  color = "Module\nActivity") +
    ggplot2::theme_minimal()
}

#' Spatial Density Module Overlay
#' @description 2D Density mapping for spatial transcriptomics.
#' @export
plot_spatial_module <- function(bmkc_result, module_index = 1, spatial_meta) {
  if (module_index > bmkc_result$n_modules) stop("Module index out of bounds.")

  genes <- bmkc_result$modules[[module_index]]
  expr_mat <- bmkc_result$data$expr_matrix[genes, ]

  module_activity <- colMeans(expr_mat)

  plot_df <- data.frame(
    X = spatial_meta$X_coord,
    Y = spatial_meta$Y_coord,
    Activity = module_activity
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = X, y = Y)) +
    ggplot2::geom_point(color = "grey80", size = 1, alpha = 0.5) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ..level.., alpha = ..level..),
                             geom = "polygon", data = plot_df[plot_df$Activity > median(plot_df$Activity), ]) +
    ggplot2::scale_fill_viridis_c(option = "plasma", name = "Module\nActivity") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste("Spatial Projection: Module", module_index),
                  x = "Tissue X Coordinate", y = "Tissue Y Coordinate")
}

#' Plot BMKC Network Module Activity on a Single-Cell Reduction (UMAP/t-SNE)
#'
#' @description Projects the summarized activity score of a specific BMKC network clique
#' onto a Seurat object's dimensionality reduction plot.
#'
#' @param bmkc_res A result object generated by \code{run_bmkc()}.
#' @param sc_obj A Seurat object containing the single-cell data.
#' @param module_index Numeric. Which module/clique to plot. Default is 1 (the top module).
#' @param reduction Character. The dimensionality reduction to plot on. Default is "umap".
#' @param pt_size Numeric. Size of the points on the plot. Default is 1.
#'
#' @return A ggplot2 object representing the FeaturePlot.
#' @import ggplot2
#' @export
plot_umap_eigengene <- function(bmkc_res, sc_obj, module_index = 1, reduction = "umap", pt_size = 1) {

  # 1. Dependency Check
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The 'Seurat' package is required for this plot. Please install it.")
  }

  # 2. Extract results and validate the module
  res <- bmkc_res$results
  if (is.null(res) || nrow(res) == 0) {
    stop("No modules found in the BMKC result object.")
  }

  clique_genes <- res$Gene[res$Module == module_index]
  if (length(clique_genes) == 0) {
    stop(sprintf("Module %d does not exist or has no genes.", module_index))
  }

  # 3. The Safety Net: Force an intersection with the Seurat object
  # This prevents AddModuleScore from crashing if gene names are slightly mismatched
  valid_genes <- intersect(clique_genes, rownames(sc_obj))

  if (length(valid_genes) < 3) {
    warning(sprintf("Only %d genes from Module %d mapped to the Seurat object. The activity score may be noisy.",
                    length(valid_genes), module_index))
  }

  if (length(valid_genes) == 0) {
    stop("None of the network genes could be found in the Seurat object. Check gene naming conventions.")
  }

  # 4. Calculate the Network Activity Score
  # Note: Seurat automatically appends a "1" to whatever name we give it here.
  score_name <- paste0("BMKC_Mod", module_index, "_")
  sc_obj <- Seurat::AddModuleScore(sc_obj,
                                   features = list(valid_genes),
                                   name = score_name)

  # The exact column name Seurat generated in the metadata
  target_col <- paste0(score_name, "1")

  # 5. Generate the Plot
  p <- Seurat::FeaturePlot(sc_obj,
                           features = target_col,
                           reduction = reduction,
                           pt.size = pt_size) +
    ggplot2::scale_color_gradientn(colors = c("lightgrey", "#fc9272", "#de2d26", "#a50f15")) +
    ggplot2::ggtitle(sprintf("BMKC Module %d Activity", module_index)) +
    ggplot2::labs(
      subtitle = sprintf("Network activity projected onto %s", toupper(reduction)),
      caption = sprintf("Calculated from %d mapped genes", length(valid_genes))
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
      legend.position = "right"
    )

  return(p)
}

#' Plot BMKC Network Activity on Physical Tissue (Spatial scRNA)
#'
#' @description Projects a network module's activity directly onto the physical
#' photograph of the tissue slide. Perfect for 'Warzone' analyses.
#'
#' @param bmkc_res A result object from \code{run_bmkc()}.
#' @param spatial_obj A Seurat object with spatial data (contains @images).
#' @param module_index Numeric. The module to project. Default is 1.
#' @param pt_size Numeric. Size of the spots. Default is 1.6.
#' @param alpha Numeric vector. Min/Max transparency for the heat gradient.
#'
#' @return A ggplot2 object from SpatialFeaturePlot.
#' @export
plot_spatial_module <- function(bmkc_res, spatial_obj, module_index = 1, pt_size = 1.6, alpha = c(0.1, 1)) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for spatial plotting.")
  }

  # 1. Extract and validate genes
  res <- bmkc_res$results
  clique_genes <- res$Gene[res$Module == module_index]

  # 2. Safety Intersect
  valid_genes <- intersect(clique_genes, rownames(spatial_obj))

  if (length(valid_genes) == 0) {
    stop("No module genes found in the spatial dataset.")
  }

  # 3. Calculate Activity Score
  score_name <- paste0("Spatial_Mod", module_index, "_")
  spatial_obj <- Seurat::AddModuleScore(spatial_obj,
                                        features = list(valid_genes),
                                        name = score_name)

  target_col <- paste0(score_name, "1")

  # 4. Generate the Spatial Plot
  # Uses 'fill' for spatial spots to create a 'warzone' heat look
  p <- Seurat::SpatialFeaturePlot(spatial_obj,
                                  features = target_col,
                                  pt.size.factor = pt_size,
                                  alpha = alpha) +
    ggplot2::scale_fill_gradientn(colors = c("lightgrey", "#fc9272", "#de2d26", "#a50f15")) +
    ggplot2::ggtitle(sprintf("Spatial Warzone: Module %d activity", module_index)) +
    ggplot2::labs(
      subtitle = "Network activity mapped to physical tissue coordinates",
      fill = "Network Score"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right"
    )

  return(p)
}

# ===================================================================
# ======================  S3 DISPATCHER  ============================
# ===================================================================

#' Plot Method for BMKC
#' @description S3 method to easily route plotting commands based on type.
#' @param x A bmkc_result object.
#' @param type String. Type of plot: "network", "heatmap", "eigengene", "raincloud", "spatial", "density", "overlap", or "umap".
#' @param ... Additional arguments passed to the specific plotting function.
#' @export
plot.bmkc_result <- function(x, type = c("default", "network", "heatmap", "eigengene", "raincloud", "spatial", "density", "overlap", "umap"), ...) {
  type <- match.arg(type)
  switch(type,
         "default"   = plot_differential_network(x, ...),
         "network"   = plot_differential_network(x, ...),
         "heatmap"   = plot_module_heatmap(x, ...),
         "eigengene" = plot_module_eigengene(x, ...),
         "raincloud" = plot_raincloud_eigengene(x, ...),
         "spatial"   = plot_eigengene_spatial(x, ...),
         "density"   = plot_spatial_module(x, ...),
         "overlap"   = plot_membership_vs_rewiring(x, ...),
         "umap"      = plot_umap_eigengene(x, ...), # NEW: UMAP/t-SNE
         "spatial"   = plot_spatial_module(x, ...) # NEW: Tissue slide
  )
}
