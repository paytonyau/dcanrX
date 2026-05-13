#' @title BMHT Visualization Suite
#' @description Plotting functions for Biweight Midcorrelation Hub Topology (BMHT) and Causal Inference.
#' @name bmht_plots
NULL

# ===================================================================
# ======================  CORE VISUALS  =============================
# ===================================================================

#' Rewiring Volcano Plot (BMHT)
#' @description Plots Differential Co-expression (DC) Score vs Significance.
#' @export
plot_rewiring_volcano <- function(bmht_result, n_label = 10) {
  df <- bmht_result$results
  df$negLogFDR <- -log10(df$FDR)
  df$is_sig <- df$FDR < bmht_result$significance_level

  top_labels <- df[order(df$FDR, -df$DC_Score), ]
  top_labels <- head(top_labels[top_labels$is_sig, ], n_label)

  ggplot2::ggplot(df, ggplot2::aes(x = DC_Score, y = negLogFDR)) +
    ggplot2::geom_point(ggplot2::aes(color = is_sig), size = 3, alpha = 0.7) +
    ggrepel::geom_text_repel(data = top_labels, ggplot2::aes(label = Gene),
                             size = 4, box.padding = 0.5, max.overlaps = Inf) +
    ggplot2::scale_color_manual(values = c("TRUE" = "red3", "FALSE" = "grey70")) +
    ggplot2::labs(title = "Rewiring Volcano Plot (BMHT)",
                  x = "Differential Co-expression Score", y = "-log10(FDR)",
                  color = "Significant") +
    ggplot2::theme_minimal(base_size = 14)
}

#' Top Gene Connectivity Bar Chart (BMHT)
#' @description Renders a bar chart of the most heavily rewired genes.
#' @export
plot_top_gene_bar <- function(bmht_result, n_top = 15) {
  top_genes <- head(bmht_result$results, n_top)
  if (!"species" %in% names(top_genes)) top_genes$species <- "All"

  ggplot2::ggplot(top_genes, ggplot2::aes(x = reorder(Gene, DC_Score), y = DC_Score, fill = species)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = paste("Top", n_top, "Most Rewired Genes (BMHT)"),
                  x = "Gene", y = "DC Score") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = ifelse(length(unique(top_genes$species)) > 1, "bottom", "none"))
}

#' Node Degree Topology Distribution (BMHT)
#' @description Plots the Cumulative Distribution Function (ECDF) to prove massive network shifts.
#' @export
plot_degree_distribution <- function(bmht_result, edge_threshold = 0.4, max_genes = 2000) {
  expr_mat <- bmht_result$data$expr_matrix

  if (nrow(expr_mat) > max_genes) {
    message(sprintf("Subsampling %d genes for topology rendering...", max_genes))
    expr_mat <- expr_mat[sample(seq_len(nrow(expr_mat)), max_genes), ]
  }

  cond <- bmht_result$data$condition
  cond_levels <- levels(cond)

  cor1 <- compute_cor_matrix(expr_mat[, cond == cond_levels[1]], workers = 1)
  cor2 <- compute_cor_matrix(expr_mat[, cond == cond_levels[2]], workers = 1)

  deg1 <- rowSums(abs(cor1) > edge_threshold) - 1
  deg2 <- rowSums(abs(cor2) > edge_threshold) - 1

  df <- data.frame(
    Degree = c(deg1, deg2),
    Condition = rep(cond_levels, each = length(deg1))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = Degree, color = Condition)) +
    ggplot2::stat_ecdf(linewidth = 1.2) +
    ggplot2::labs(title = "Network Topology: Degree Distribution",
                  subtitle = paste("Cumulative distribution of node connections (|r| >", edge_threshold, ")"),
                  x = "Node Degree (Connections)", y = "Cumulative Probability") +
    ggplot2::theme_bw(base_size = 14)
}

#' Condition-Specific Hub Networks (Side-by-Side)
#' @description Renders two force-directed networks to visually prove the rewiring of top hubs.
#' @export
plot_condition_hubs <- function(bmht_result, n_top = 12, edge_threshold = 0.5) {
  top_genes <- head(bmht_result$results$Gene, n_top)
  expr_mat <- bmht_result$data$expr_matrix[top_genes, ]
  cond <- bmht_result$data$condition
  cond_levels <- levels(cond)

  cor1 <- compute_cor_matrix(expr_mat[, cond == cond_levels[1]], workers = 1)
  cor2 <- compute_cor_matrix(expr_mat[, cond == cond_levels[2]], workers = 1)

  cor1[abs(cor1) < edge_threshold] <- 0
  cor2[abs(cor2) < edge_threshold] <- 0

  g1 <- igraph::graph_from_adjacency_matrix(abs(cor1), mode = "undirected", weighted = TRUE, diag = FALSE)
  g2 <- igraph::graph_from_adjacency_matrix(abs(cor2), mode = "undirected", weighted = TRUE, diag = FALSE)

  p1 <- ggraph::ggraph(g1, layout = "fr") +
    ggraph::geom_edge_link(ggplot2::aes(width = weight), alpha = 0.5, color = "darkgrey") +
    ggraph::geom_node_point(size = 5, color = "#1f78b4") +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE) +
    ggraph::theme_graph() +
    ggplot2::labs(title = paste("Hub Network -", cond_levels[1]))

  p2 <- ggraph::ggraph(g2, layout = "fr") +
    ggraph::geom_edge_link(ggplot2::aes(width = weight), alpha = 0.5, color = "darkgrey") +
    ggraph::geom_node_point(size = 5, color = "#e31a1c") +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE) +
    ggraph::theme_graph() +
    ggplot2::labs(title = paste("Hub Network -", cond_levels[2]))

  patchwork::wrap_plots(p1, p2, ncol = 2)
}

#' Global Correlation Shift Density Plot
#' @description Maps the density of all pairwise correlations to identify global regulatory strictness.
#' @export
plot_correlation_shift <- function(bmht_result) {
  expr_mat <- bmht_result$data$expr_matrix
  cond <- bmht_result$data$condition
  cond_levels <- levels(cond)

  cor1 <- compute_cor_matrix(expr_mat[, cond == cond_levels[1]], workers = 1)
  cor2 <- compute_cor_matrix(expr_mat[, cond == cond_levels[2]], workers = 1)

  vals1 <- cor1[upper.tri(cor1)]
  vals2 <- cor2[upper.tri(cor2)]

  df <- data.frame(
    Correlation = c(vals1, vals2),
    Condition = rep(cond_levels, each = length(vals1))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = Correlation, fill = Condition)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::scale_fill_manual(values = c("#1f78b4", "#e31a1c")) +
    ggplot2::labs(title = "Global Network Shift",
                  subtitle = "Density distribution of all pairwise correlations",
                  x = "Biweight Midcorrelation (bicor)",
                  y = "Density") +
    ggplot2::theme_minimal(base_size = 14)
}

#' Plot Directed Causal Network
#' @description Plots a directed star-graph highlighting a master regulator and its downstream targets.
#' @export
plot_causal_network <- function(bmht_result, min_edge = 0.4) {
  if (is.null(bmht_result$causal_drivers) || nrow(bmht_result$causal_drivers) == 0) {
    stop("No causal drivers found. Run `identify_causal_drivers()` first.")
  }

  driver_gene <- bmht_result$causal_drivers$Gene[1]
  expr_mat <- bmht_result$data$expr_matrix
  cond <- bmht_result$data$condition

  cor2 <- compute_cor_matrix(expr_mat[, cond == levels(cond)[2]], workers = 1)
  driver_cor <- cor2[driver_gene, ]

  targets <- names(driver_cor[abs(driver_cor) >= min_edge & names(driver_cor) != driver_gene])

  if (length(targets) == 0) stop("The driver gene has no strong connections to plot above the threshold.")

  edges_df <- data.frame(
    from = driver_gene,
    to = targets,
    weight = driver_cor[targets]
  )

  g_directed <- igraph::graph_from_data_frame(edges_df, directed = TRUE)

  ggraph::ggraph(g_directed, layout = "star") +
    ggraph::geom_edge_link(ggplot2::aes(color = weight, width = abs(weight)),
                           arrow = ggplot2::arrow(length = ggplot2::unit(4, 'mm')),
                           end_cap = ggraph::circle(5, 'mm'), alpha = 0.8) +
    ggraph::geom_node_point(size = 7, color = ifelse(igraph::V(g_directed)$name == driver_gene, "#e31a1c", "grey50")) +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
    ggraph::scale_edge_color_gradient2(low = "#1f78b4", mid = "white", high = "#e31a1c", guide = "none") +
    ggraph::scale_edge_width(range = c(0.5, 2), guide = "none") +
    ggraph::theme_graph() +
    ggplot2::labs(title = paste("Predicted Causal Circuit: Driven by", driver_gene),
                  subtitle = "Arrow indicates directed regulatory control")
}

# ===================================================================
# ======================  ADVANCED VISUALS  =========================
# ===================================================================

#' Plot Network Hive Layout (BMHT Hubs)
#' @description A highly structured radial layout for displaying massive hubs without edge overlapping.
#' @export
plot_hive_network <- function(bmht_result, n_top = 30, min_edge = 0.4) {
  genes <- head(bmht_result$results$Gene, n_top)
  expr_mat <- bmht_result$data$expr_matrix[genes, ]
  cond <- bmht_result$data$condition

  cor2 <- compute_cor_matrix(expr_mat[, cond == levels(cond)[2]], workers = 1)
  cor2[abs(cor2) < min_edge] <- 0

  g <- igraph::graph_from_adjacency_matrix(cor2, mode = "undirected", weighted = TRUE, diag = FALSE)

  deg <- igraph::degree(g)
  igraph::V(g)$axis <- as.character(cut(deg, breaks = 3, labels = c("Low", "Medium", "High")))
  igraph::V(g)$degree <- deg

  ggraph::ggraph(g, layout = "linear", circular = TRUE) +
    ggraph::geom_edge_arc(ggplot2::aes(color = weight, width = abs(weight)), alpha = 0.6) +
    ggraph::geom_node_point(ggplot2::aes(color = axis, size = degree)) +
    ggraph::scale_edge_color_gradient2(low = "#1f78b4", mid = "grey90", high = "#e31a1c") +
    ggraph::theme_graph() +
    ggplot2::labs(title = "Hive Hub Architecture (Disease State)",
                  subtitle = "Nodes structured radially by connectivity")
}

#' Plot Eruption (Volcano + Network Hybrid)
#' @description Overlays the physical network wiring directly onto a Volcano plot.
#' @export
plot_eruption <- function(bmht_result, de_results, min_edge = 0.5, fdr_cutoff = 0.05) {
  sig_genes <- bmht_result$results$Gene[bmht_result$results$FDR < fdr_cutoff]
  if (length(sig_genes) == 0) stop("No significant BMHT genes found to plot.")

  expr_mat <- bmht_result$data$expr_matrix[sig_genes, ]
  cond <- bmht_result$data$condition

  cor1 <- compute_cor_matrix(expr_mat[, cond == levels(cond)[1]], workers = 1)
  cor2 <- compute_cor_matrix(expr_mat[, cond == levels(cond)[2]], workers = 1)
  diff_mat <- cor2 - cor1
  diff_mat[abs(diff_mat) < min_edge] <- 0

  g <- igraph::graph_from_adjacency_matrix(diff_mat, mode = "undirected", weighted = TRUE, diag = FALSE)

  node_names <- igraph::V(g)$name
  de_match <- de_results[match(node_names, de_results$Gene), ]

  layout_matrix <- as.matrix(data.frame(
    x = de_match$logFC,
    y = -log10(de_match$pvalue)
  ))

  igraph::E(g)$diff <- igraph::E(g)$weight
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)

  ggraph::ggraph(g, layout = layout_matrix) +
    ggraph::geom_edge_link(ggplot2::aes(color = diff, width = weight), alpha = 0.6) +
    ggraph::geom_node_point(size = 4, color = "black") +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE) +
    ggraph::scale_edge_color_gradient2(low = "blue", mid = "grey", high = "red") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Eruption Plot: Network Wiring mapped onto Differential Expression",
                  x = "Log2 Fold Change (Expression)",
                  y = "-Log10 P-value (Expression)")
}

#' Plot Cinematic Rewiring Animation (BMHT)
#' @description Animates the condition shift using gganimate.
#' @export
plot_network_animation <- function(bmht_result, n_top = 15, edge_threshold = 0.4) {
  if (!requireNamespace("gganimate", quietly = TRUE)) stop("Please install 'gganimate' and 'transformr'.")

  genes <- head(bmht_result$results$Gene, n_top)
  expr_mat <- bmht_result$data$expr_matrix[genes, ]
  cond <- bmht_result$data$condition
  cond_levels <- levels(cond)

  c1 <- compute_cor_matrix(expr_mat[, cond == cond_levels[1]], workers = 1)
  c2 <- compute_cor_matrix(expr_mat[, cond == cond_levels[2]], workers = 1)

  c1[abs(c1) < edge_threshold] <- 0
  c2[abs(c2) < edge_threshold] <- 0

  g_global <- igraph::graph_from_adjacency_matrix(abs(c1) + abs(c2), mode = "undirected", weighted = TRUE, diag = FALSE)
  layout_mat <- igraph::layout_with_fr(g_global)

  build_frame <- function(adj_mat, state_name) {
    g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
    df <- igraph::as_data_frame(g, what = "edges")
    if(nrow(df) == 0) return(data.frame())
    df$State <- state_name
    return(df)
  }

  edges_c1 <- build_frame(c1, cond_levels[1])
  edges_c2 <- build_frame(c2, cond_levels[2])
  edge_df <- rbind(edges_c1, edges_c2)
  edge_df$State <- factor(edge_df$State, levels = cond_levels)

  nodes_df <- data.frame(
    id = igraph::V(g_global)$name,
    x = layout_mat[, 1],
    y = layout_mat[, 2]
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = nodes_df, ggplot2::aes(x = x, y = y), size = 5, color = "steelblue") +
    ggplot2::geom_text(data = nodes_df, ggplot2::aes(x = x, y = y, label = id), vjust = -1) +
    ggplot2::theme_void() +
    ggplot2::labs(title = 'Network State: {closest_state}') +
    gganimate::transition_states(State, transition_length = 2, state_length = 1)

  return(p)
}

#' Plot 3D Network Landscape
#' @description Renders an interactive 3D topography using plotly.
#' @export
plot_3d_network <- function(bmht_result, n_top = 20, edge_threshold = 0.5) {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Please install 'plotly'.")

  genes <- head(bmht_result$results$Gene, n_top)
  expr_mat <- bmht_result$data$expr_matrix[genes, ]
  cond <- bmht_result$data$condition

  cor2 <- compute_cor_matrix(expr_mat[, cond == levels(cond)[2]], workers = 1)
  cor2[abs(cor2) < edge_threshold] <- 0

  g <- igraph::graph_from_adjacency_matrix(abs(cor2), mode = "undirected", weighted = TRUE, diag = FALSE)
  layout_3d <- igraph::layout_with_fr(g, dim = 3)

  nodes <- data.frame(
    name = igraph::V(g)$name,
    x = layout_3d[, 1],
    y = layout_3d[, 2],
    z = igraph::degree(g)
  )

  plotly::plot_ly(nodes, x = ~x, y = ~y, z = ~z, text = ~name, mode = "markers+text",
                  marker = list(size = 8, color = ~z, colorscale = "Viridis"),
                  textposition = "top") %>%
    plotly::layout(title = "3D Hub Landscape (Disease State)",
                   scene = list(xaxis = list(title = "X"),
                                yaxis = list(title = "Y"),
                                zaxis = list(title = "Hub Connectivity (Z)")))
}

#' Plot Taxonomic Network Rewiring (Microbiome Bridge)
#'
#' @description Aggregates node-level network shifts into higher taxonomic ranks
#' (e.g., Phylum, Class) to show community-level rewiring.
#'
#' @param bmht_res A result object from \code{run_bmht()}.
#' @param taxa_df A data.frame mapping Gene/Node IDs to Taxonomic Ranks.
#' Must contain a 'Gene' column to match the BMHT results.
#' @param rank Character. The taxonomic rank to aggregate by (e.g., "Phylum", "Family").
#' @param top_n Numeric. How many top taxonomic groups to display. Default is 10.
#'
#' @return A ggplot2 object.
#' @export
plot_taxa_rewiring <- function(bmht_res, taxa_df, rank = "Phylum", top_n = 10) {

  # 1. Dependency and Input Validation
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")

  res <- bmht_res$results
  if (is.null(res) || nrow(res) == 0) stop("No network results found in the BMHT object.")

  if (!"Gene" %in% colnames(taxa_df)) {
    stop("The taxonomy table (taxa_df) must contain a 'Gene' column to match the network results.")
  }
  if (!rank %in% colnames(taxa_df)) {
    stop(sprintf("Taxonomic rank '%s' not found in the provided taxonomy table.", rank))
  }

  # 2. Merge Network Scores with Taxonomy
  merged_df <- merge(res, taxa_df, by = "Gene", all.x = FALSE)

  if (nrow(merged_df) == 0) {
    stop("Safety Net Triggered: No matching nodes found between the network results and the taxonomy table.")
  }

  # 3. Identify the Score Column (Flexibility for future engine updates)
  score_col <- if ("Delta_Degree" %in% colnames(merged_df)) "Delta_Degree" else "Score"
  if (!score_col %in% colnames(merged_df)) stop("Could not identify the network score column.")

  # 4. Aggregate Shifts by Taxonomic Rank
  # We take the absolute value so negative and positive shifts don't cancel each other out
  merged_df$Abs_Shift <- abs(merged_df[[score_col]])
  agg_df <- aggregate(merged_df$Abs_Shift, by = list(Taxon = merged_df[[rank]]), FUN = sum)
  colnames(agg_df)[2] <- "Total_Rewiring"

  # Sort and grab the top results
  agg_df <- agg_df[order(-agg_df$Total_Rewiring), ]
  agg_df <- head(agg_df, top_n)

  # Lock factor levels so ggplot doesn't automatically alphabetize the bars
  agg_df$Taxon <- factor(agg_df$Taxon, levels = rev(agg_df$Taxon))

  # 5. Generate the Plot
  p <- ggplot2::ggplot(agg_df, ggplot2::aes(x = Taxon, y = Total_Rewiring, fill = Taxon)) +
    ggplot2::geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_d(option = "mako", guide = "none") +
    ggplot2::labs(
      title = paste("Community Network Rewiring by", rank),
      subtitle = "Total absolute connectivity shift per taxonomic group",
      x = rank,
      y = "Cumulative Network Shift Score",
      caption = sprintf("Aggregated from %d mapped nodes", nrow(merged_df))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(face = "italic"),
      axis.text.y = ggplot2::element_text(face = "italic", size = 11) # Taxa are usually italicized
    )

  return(p)
}

# ===================================================================
# ======================  S3 DISPATCHER  ============================
# ===================================================================

#' Plot Method for BMHT Results
#'
#' @description S3 method to easily route plotting commands based on the desired visualization type.
#'
#' @param x A \code{bmht_result} object.
#' @param type Character. The type of plot to generate. Options include:
#' "volcano", "bar", "hubs", "topology", "shift", "causal", "hive",
#' "eruption", "animation", "3d", "taxa", or "microbiome".
#' @param ... Additional arguments passed to the specific plotting function
#' (e.g., \code{taxa_df} for microbiome plots, or \code{n_top} for bar charts).
#'
#' @method plot bmht_result
#' @export
plot.bmht_result <- function(x, type = c("default", "volcano", "bar", "hubs", "topology",
                                         "shift", "causal", "hive", "eruption",
                                         "animation", "3d", "taxa", "microbiome"), ...) {
  # Match argument ensures users can type "micro" and it will auto-complete to "microbiome"
  type <- match.arg(type)

  switch(type,
         "default"    = plot_rewiring_volcano(x, ...),
         "volcano"    = plot_rewiring_volcano(x, ...),
         "bar"        = plot_top_gene_bar(x, ...),
         "hubs"       = plot_condition_hubs(x, ...),
         "topology"   = plot_degree_distribution(x, ...),
         "shift"      = plot_correlation_shift(x, ...),
         "causal"     = plot_causal_network(x, ...),
         "hive"       = plot_hive_network(x, ...),
         "eruption"   = plot_eruption(x, ...),
         "animation"  = plot_network_animation(x, ...),
         "3d"         = plot_3d_network(x, ...),
         "taxa"       = plot_taxa_rewiring(x, ...),
         "microbiome" = plot_taxa_rewiring(x, ...))
}
