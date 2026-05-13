#' @title ROS-DET Visualization Suite
#' @description Plotting functions for Differential Edge Topologies (Correlation Switching).
#' @name rosdet_plots
NULL

# ===================================================================
# ======================  CORE VISUALS  =============================
# ===================================================================

#' Stacked Switching Diagram (ROS-DET "Brick" Plot)
#' @description Shows how many specific relationships completely inverted their logic.
#' @export
plot_switching_stacked <- function(rosdet_result, n_top = 12) {
  if (nrow(rosdet_result$results) == 0) stop("No significant results to plot.")
  top_pairs <- head(rosdet_result$results, n_top)

  links <- data.frame(
    source = ifelse(top_pairs$r1 > 0, "Positive (Cond1)", "Negative (Cond1)"),
    target = ifelse(top_pairs$r2 > 0, "Positive (Cond2)", "Negative (Cond2)"),
    value = 1, # Give each pair a weight of 1 so they stack perfectly
    pair = paste(top_pairs$Gene1, "↔", top_pairs$Gene2)
  )

  ggplot2::ggplot(links, ggplot2::aes(x = source, y = value, fill = target, group = pair)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.5, color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = pair), position = ggplot2::position_stack(vjust = 0.5),
                       size = 3.5, color = "white", fontface = "bold") +
    ggplot2::scale_fill_manual(values = c("Positive (Cond2)" = "#1f78b4", "Negative (Cond2)" = "#e31a1c")) +
    ggplot2::labs(title = "Correlation Switching Breakdown (ROS-DET)",
                  subtitle = paste("Top", n_top, "switching pairs"),
                  x = "Starting State", y = "Number of Pairs") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.major.x = ggplot2::element_blank())
}

#' The Quadrant Scatter (ROS-DET "Flip-Flop" Proof)
#' @description Plots Cond1 vs Cond2 correlation to visualize the complete inversion.
#' @export
plot_rosdet_quadrant <- function(rosdet_result, n_label = 15) {
  df <- rosdet_result$results
  df$is_sig <- df$FDR < rosdet_result$significance_level

  top_labels <- head(df[df$is_sig, ], n_label)
  top_labels$pair_label <- paste(top_labels$Gene1, top_labels$Gene2, sep="\n")

  ggplot2::ggplot(df, ggplot2::aes(x = r1, y = r2)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey80", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(color = is_sig, size = abs(Score)), alpha = 0.7) +
    ggrepel::geom_text_repel(data = top_labels, ggplot2::aes(label = pair_label), size = 3, box.padding = 0.5) +
    ggplot2::scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "grey80")) +
    ggplot2::labs(title = "ROS-DET Quadrant Map",
                  subtitle = "Significant pairs sit in the top-left or bottom-right quadrants",
                  x = "Correlation in Condition 1",
                  y = "Correlation in Condition 2",
                  color = "Significant Switch", size = "Switch Score") +
    ggplot2::theme_bw(base_size = 14)
}

#' Species Bipartite Interaction Network (ROS-DET)
#' @description Creates a bipartite network for Host-Pathogen / Cross-Species modes.
#' @export
plot_species_interaction <- function(rosdet_result) {
  if (!rosdet_result$species_mode) stop("Species mode not enabled.")
  inter <- rosdet_result$results[rosdet_result$results$pair_type == "inter", ]
  if (nrow(inter) == 0) stop("No inter-species pairs found.")

  g <- igraph::graph_from_data_frame(inter[, c("Gene1", "Gene2")], directed = FALSE)

  sp_map <- rosdet_result$data$species
  names(sp_map) <- rownames(rosdet_result$data$expr_matrix)
  igraph::V(g)$species <- sp_map[igraph::V(g)$name]

  # CRITICAL FIX: igraph requires a logical 'type' attribute to render bipartite layouts safely
  unique_species <- unique(igraph::V(g)$species)
  igraph::V(g)$type <- igraph::V(g)$species == unique_species[1]

  ggraph::ggraph(g, layout = "bipartite") +
    ggraph::geom_edge_link(ggplot2::aes(color = inter$Score), alpha = 0.8, width = 1.1) +
    ggraph::geom_node_point(ggplot2::aes(color = species), size = 7) +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 4) +
    ggraph::scale_edge_color_gradient(low = "grey70", high = "darkred", name = "Switch Score") +
    ggraph::theme_graph() +
    ggplot2::labs(title = "Cross-Species Switching Interactions")
}

# ===================================================================
# ======================  ADVANCED VISUALS  =========================
# ===================================================================

#' Plot Differential Edge Volcano
#' @description Plots the change in correlation (Delta r) vs significance for gene pairs.
#' @export
plot_edge_volcano <- function(rosdet_res) {
  df <- rosdet_res$results

  # Highlight the radical switchers (Delta > 0.8 and P < 0.05)
  df$Status <- "Stable"
  df$Status[df$Score > 0.8 & df$FDR < rosdet_res$significance_level] <- "Gained Co-expression"
  df$Status[df$Score < -0.8 & df$FDR < rosdet_res$significance_level] <- "Lost Co-expression (Switched)"

  ggplot2::ggplot(df, ggplot2::aes(x = Score, y = -log10(FDR), color = Status)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::geom_vline(xintercept = c(-0.8, 0.8), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(rosdet_res$significance_level), linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("Stable" = "grey",
                                           "Gained Co-expression" = "#1f78b4",
                                           "Lost Co-expression (Switched)" = "#e31a1c")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Network Edge Volcano Plot",
                  x = "Change in Correlation (\u0394 r)", y = "-log10(FDR)")
}

#' Plot Circular intercellular Diagram (ROS-DET)
#' @description Circular layout perfect for Single-Cell Communication and massive Host-Pathogen interactions.
#' @param rosdet_result A ROS-DET result object.
#' @param n_top Numeric. The maximum number of top switching edges to plot. Default is 30.
#' @export
plot_intercellular_network <- function(rosdet_result, n_top = 30) {
  if (!requireNamespace("circlize", quietly = TRUE)) stop("Please install 'circlize' to use chord diagrams.")

  # 1. Extract results
  res <- rosdet_result$results

  if (is.null(res) || nrow(res) == 0) {
    stop("No significant switching edges found in the ROS-DET result.")
  }

  # 2. Sort to guarantee we plot the most dramatic shifts
  # Checks for column name variations (Score vs Delta_Cor)
  if ("Score" %in% colnames(res)) {
    res <- res[order(-abs(res$Score)), ]
    edge_weight <- abs(res$Score)
  } else if ("Delta_Cor" %in% colnames(res)) {
    res <- res[order(-abs(res$Delta_Cor)), ]
    edge_weight <- abs(res$Delta_Cor)
  } else {
    stop("Could not find 'Score' or 'Delta_Cor' column in results.")
  }

  # 3. Apply the limit
  top_pairs <- head(res, n_top)
  edge_weight <- head(edge_weight, n_top)

  # 4. Calculate the directional shift for the colors
  if (all(c("r1", "r2") %in% colnames(top_pairs))) {
    diff_val <- top_pairs$r2 - top_pairs$r1
  } else {
    diff_val <- top_pairs$Delta_Cor
  }

  # 5. Create adjacency data frame for circlize
  chord_df <- data.frame(
    from = top_pairs$Gene1,
    to = top_pairs$Gene2,
    value = edge_weight, # Thickness of the ribbon
    diff = diff_val      # Color of the ribbon (Red = gained, Blue = lost)
  )

  # Set colors: Red = Positive Shift (Cooperating), Blue = Negative Shift (Competing)
  ribbon_colors <- ifelse(chord_df$diff > 0, "#e31a1c99", "#1f78b499")

  # 6. Clear base R plot parameters
  circlize::circos.clear()
  circlize::circos.par(gap.after = 2)

  # 7. Render
  circlize::chordDiagram(
    chord_df[, 1:3],
    col = ribbon_colors,
    transparency = 0.3,
    annotationTrack = c("name", "grid"),
    preAllocateTracks = list(track.height = 0.1)
  )
  title(main = "Network Rewiring Communication")
}

# ===================================================================
# ======================  S3 DISPATCHER  ============================
# ===================================================================

#' Plot Method for ROS-DET
#' @description S3 method to easily route plotting commands based on type.
#' @param x A rosdet_result object.
#' @param type String. Type of plot: "default", "quadrant", "stacked", "species", "volcano", or "chord".
#' @param ... Additional arguments passed to the specific plotting function.
#' @export
plot.rosdet_result <- function(x, type = c("default", "quadrant", "stacked", "species", "volcano", "chord"), ...) {
  type <- match.arg(type)
  switch(type,
         "default"  = plot_rosdet_quadrant(x, ...),
         "quadrant" = plot_rosdet_quadrant(x, ...),
         "stacked"  = plot_switching_stacked(x, ...),
         "species"  = plot_species_interaction(x, ...),
         "volcano"  = plot_edge_volcano(x, ...),
         "intercellular"    = plot_intercellular_network(x, ...))
}
