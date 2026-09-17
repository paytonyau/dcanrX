#' @title bicorX Utilities and Graphics Theme
#' @description Shared themes, color palettes, and multi-method comparison visualizations.
#' @importFrom utils head write.csv
#' @importFrom stats reorder
#' @importFrom rlang .data
#' @name utils_graphics
NULL

utils::globalVariables(c("Time", "Correlation", "GeneRatio", "Description", "Count"))

# ===================================================================
# ======================  GLOBAL STYLES & THEMES ====================
# ===================================================================

#' Custom bicorX ggplot2 Theme
#' @param base_size Base font size in points. Default `14`.
#' @return A `ggplot2` theme object.
#' @export
theme_bicor <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = base_size - 2, hjust = 0.5, color = "grey30"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "#f0f0f0", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 1)
    )
}

#' Standard Palette for bicorX
#' @export
bicor_palette <- function() {
  c(
    "Condition1" = "#1f78b4",
    "Condition2" = "#e31a1c",
    "Neutral"    = "#969696",
    "Highlight"  = "#ff7f00",
    "Positive"   = "#1f78b4",
    "Negative"   = "#e31a1c"
  )
}

# ===================================================================
# ======================  MULTI-METHOD / OMICS  =====================
# ===================================================================

#' Plot Rolling Trajectory Heatmap
#' @param expr_matrix Numeric matrix, features (rows) by cells (columns).
#' @param geneA,geneB Names of the two features whose rolling correlation
#'   is to be plotted.
#' @param pseudotime_vec Numeric vector of pseudotime values, one per
#'   column of `expr_matrix`.
#' @param window_size Number of cells per sliding window. Default `30`.
#' @return A `ggplot` object.
#' @export
plot_trajectory_heatmap <- function(expr_matrix, geneA, geneB, pseudotime_vec, window_size = 30) {
  order_idx <- order(pseudotime_vec)
  exprA <- expr_matrix[geneA, order_idx]
  exprB <- expr_matrix[geneB, order_idx]
  time_sorted <- pseudotime_vec[order_idx]

  n_cells <- length(exprA)
  if (n_cells < window_size * 2) stop("Not enough cells for the specified window size.")

  rolling_cor <- numeric(n_cells - window_size)
  window_time <- numeric(n_cells - window_size)

  for (i in seq_len(n_cells - window_size)) {
    idx <- i:(i + window_size - 1)
    # PERFORMANCE & COMPILATION FIX: Replaced WGCNA hard-coded hook with package native bicor tracker
    rolling_cor[i] <- bicor(exprA[idx], exprB[idx])
    window_time[i] <- mean(time_sorted[idx])
  }

  plot_df <- data.frame(Time = window_time, Correlation = rolling_cor)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Time, y = .data$Correlation, fill = .data$Correlation)) +
    ggplot2::geom_col(width = diff(range(window_time))/length(window_time)) +
    ggplot2::scale_fill_gradient2(low = bicor_palette()["Condition1"],
                                  mid = "white",
                                  high = bicor_palette()["Condition2"], limits = c(-1, 1)) +
    theme_bicor() +
    ggplot2::labs(title = paste("Dynamic Rewiring:", geneA, "vs", geneB),
                  subtitle = paste("Rolling bicor (Window size =", window_size, "cells)"),
                  x = "Pseudotime / Developmental Trajectory",
                  y = "Correlation Strength")
}
