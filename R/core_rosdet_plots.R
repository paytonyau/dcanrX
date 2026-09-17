#' @title ROS-DET Visualization Suite
#' @description Plotting functions for differential edge / correlation-switching results from ROS-DET.
#' @name rosdet_plots
NULL

# Prevent CRAN package check notes regarding unquoted column name structures
utils::globalVariables(c("r1", "r2", "is_sig", "Distance_Score", "pair_label",
                         "Delta_r", "FDR", "Status", "source", "value", "target",
                         "pair", "Var1", "Var2", "Freq", "Condition", "shift",
                         "name", "species", "Edge_Type", "abs_weight", "Distance", "Signal"))

# ==============================================================================
# CATEGORY 1: UNIVERSAL DIAGNOSTICS
# ==============================================================================

#' The Quadrant Scatter (ROS-DET "Flip-Flop" Proof)
#' @description Plots Cond1 vs Cond2 correlation to visualize complete edge inversion and decoupling.
#' @param rosdet_result A `rosdet_result` object from [run_rosdet()].
#' @param n_label Number of top pairs to label. Default `15`.
#' @param ... Currently unused; accepted for S3 compatibility.
#' @return A `ggplot` object.
#' @export
plot_rosdet_quadrant <- function(rosdet_result, n_label = 15, ...) {
  df <- rosdet_result$results
  df$is_sig <- df$FDR < rosdet_result$significance_level

  top_labels <- head(df[df$is_sig, ], n_label)
  top_labels$pair_label <- paste(top_labels$Gene1, top_labels$Gene2, sep="\n")

  ggplot2::ggplot(df, ggplot2::aes(x = .data$r1, y = .data$r2)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey80", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(color = .data$is_sig, size = abs(.data$Distance_Score)), alpha = 0.7) +
    ggrepel::geom_text_repel(data = top_labels, ggplot2::aes(label = .data$pair_label), size = 3, box.padding = 0.5) +
    ggplot2::scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "grey80")) +
    ggplot2::labs(title = "ROS-DET Quadrant Map",
                  subtitle = "Significant pairs sit in the off-diagonals (re-wired or decoupled)",
                  x = "Correlation in Condition 1",
                  y = "Correlation in Condition 2",
                  color = "Significant Switch", size = "Distance Score") +
    ggplot2::theme_bw(base_size = 14, base_family = "sans")
}

#' Plot Differential Edge Volcano
#' @description Plots the absolute change in co-expression distance vs significance for gene pairs.
#' @param rosdet_res A `rosdet_result` object from [run_rosdet()].
#' @param ... Currently unused; accepted for S3 compatibility.
#' @return A `ggplot` object.
#' @export
plot_edge_volcano <- function(rosdet_res, ...) {
  df <- rosdet_res$results

  # Prevent -log10(0) from producing Inf and erasing the most significant points
  df$FDR[df$FDR == 0] <- 1e-300

  # IMPLEMENTATION OF THE 4-TIER STABILITY SYSTEM (Continuous Delta Integration)
  df$Status <- "Stable"
  sig_mask <- df$FDR < rosdet_res$significance_level
  df$Status[sig_mask & df$Distance_Score >= 1.2] <- "Severe Disruption"
  df$Status[sig_mask & df$Distance_Score >= 0.8 & df$Distance_Score < 1.2] <- "Moderate Re-wiring"
  df$Status[sig_mask & df$Distance_Score >= 0.4 & df$Distance_Score < 0.8] <- "Weak Decoupling"

  df$Status <- factor(df$Status, levels = c("Stable", "Weak Decoupling", "Moderate Re-wiring", "Severe Disruption"))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$Distance_Score, y = -log10(.data$FDR), color = .data$Status)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::geom_vline(xintercept = c(0.4, 0.8, 1.2), linetype = "dashed", color = "grey60") +
    ggplot2::geom_hline(yintercept = -log10(rosdet_res$significance_level), linetype = "dashed", color = "grey60") +
    ggplot2::scale_color_manual(values = c("Stable" = "grey80",
                                           "Weak Decoupling" = "#74add1",
                                           "Moderate Re-wiring" = "#f46d43",
                                           "Severe Disruption" = "#d73027")) +
    ggplot2::theme_minimal(base_family = "sans") +
    ggplot2::labs(title = "Network Edge Volcano Plot",
                  subtitle = "ROS-DET 4-Tier Topology Disruption Space",
                  x = "Absolute Distance Score (Rewiring Magnitude)", y = "-log10(FDR)",
                  color = "Rewiring Status")
}

# ==============================================================================
# CORE S3 DISPATCHER
# ==============================================================================

#' Plot Method for ROS-DET
#'
#' As of the Phase 3 surface reduction, only two plot types are provided
#' directly: `"quadrant"` (default) and `"volcano"`. For other
#' visualizations (switching heatmaps/networks, species interaction, ego
#' networks, spatial/UMAP overlays), see the "Custom plots from bicorX
#' results" vignette.
#' @method plot rosdet_result
#' @param x A `rosdet_result` object from [run_rosdet()].
#' @param type Which plot to draw: `"default"`/`"quadrant"` for the
#'   correlation quadrant plot, or `"volcano"` for the edge volcano.
#' @param ... Passed to the underlying plotting function.
#' @return A `ggplot` object.
#' @export
plot.rosdet_result = function(x, type = c("default", "quadrant", "volcano"), ...) {
  type <- match.arg(type)

  switch(type,
         "default"  = plot_rosdet_quadrant(x, ...),
         "quadrant" = plot_rosdet_quadrant(x, ...),
         "volcano"  = plot_edge_volcano(x, ...)
  )
}
