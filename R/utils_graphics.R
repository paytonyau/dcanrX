#' @title bicorDCEA Utilities and Graphics Theme
#' @description Shared themes, color palettes, and multi-method comparison visualizations.
#' @name utils_graphics
NULL

# ===================================================================
# ======================  GLOBAL STYLES & THEMES ====================
# ===================================================================

#' Custom bicorDCEA ggplot2 Theme
#' @description A standardized, clean theme applied across the package to ensure publication-ready consistency.
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

#' Standard Palette for bicorDCEA
#' @description Standardized hex codes to keep Condition 1 and Condition 2 visually consistent across all plots.
#' @export
bicor_palette <- function() {
  c(
    "Condition1" = "#1f78b4", # Deep Blue (Healthy / Control)
    "Condition2" = "#e31a1c", # Firebrick Red (Disease / Treated)
    "Neutral"    = "#969696", # Grey
    "Highlight"  = "#ff7f00", # Orange (For Master Regulators / Drivers)
    "Positive"   = "#1f78b4",
    "Negative"   = "#e31a1c"
  )
}

# ===================================================================
# ======================  MULTI-METHOD / OMICS  =====================
# ===================================================================

#' Multi-Method UpSet Plot (Intersection of Discoveries)
#' @description Compares the genes identified by ROS-DET, BMHT, and BMKC to find universal network drivers.
#' @export
plot_method_overlap <- function(rosdet_res = NULL, bmht_res = NULL, bmkc_res = NULL) {
  if (!requireNamespace("UpSetR", quietly = TRUE)) stop("Please install the 'UpSetR' package.")

  gene_lists <- list()

  if (!is.null(rosdet_res)) {
    sig_pairs <- rosdet_res$results[rosdet_res$results$FDR < rosdet_res$significance_level, ]
    gene_lists[["ROS-DET\n(Switching)"]] <- unique(c(sig_pairs$Gene1, sig_pairs$Gene2))
  }
  if (!is.null(bmht_res)) {
    sig_hubs <- bmht_res$results[bmht_res$results$FDR < bmht_res$significance_level, ]
    gene_lists[["BMHT\n(Hubs)"]] <- sig_hubs$Gene
  }
  if (!is.null(bmkc_res)) {
    gene_lists[["BMKC\n(Modules)"]] <- unique(unlist(bmkc_res$modules))
  }

  if (length(gene_lists) < 2) stop("Provide at least two method results to compare overlap.")

  # Suppress harmless ggplot2 deprecation warnings originating from the UpSetR package
  suppressWarnings({
    UpSetR::upset(UpSetR::fromList(gene_lists),
                  order.by = "freq",
                  main.bar.color = bicor_palette()["Condition1"],
                  text.scale = 1.5)
  })
}

#' GO Enrichment Dot Plot wrapper
#' @description A clean, standardized dot plot for displaying pathway enrichment results of network hubs.
#' @export
plot_enrichment_dotplot <- function(enrich_res, top_n = 15, title = "Pathway Enrichment") {
  if (!all(c("Description", "GeneRatio", "p.adjust", "Count") %in% names(enrich_res))) {
    stop("enrich_res must be a dataframe containing: Description, GeneRatio, p.adjust, Count.")
  }

  df <- head(enrich_res[order(enrich_res$p.adjust), ], top_n)

  # Convert string fractions ("5/100") to numeric decimals if necessary
  if (is.character(df$GeneRatio)) {
    df$GeneRatio <- sapply(df$GeneRatio, function(x) {
      parts <- as.numeric(strsplit(x, "/")[[1]])
      parts[1] / parts[2]
    })
  }

  ggplot2::ggplot(df, ggplot2::aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
    ggplot2::scale_color_gradient(low = bicor_palette()["Condition2"], high = bicor_palette()["Condition1"]) +
    ggplot2::labs(title = title, x = "Gene Ratio", y = "") +
    theme_bicor()
}

#' Plot Rolling Trajectory Heatmap
#' @description For scRNA-seq trajectory analysis. Plots how the correlation between two genes changes dynamically over pseudotime.
#' @export
plot_trajectory_heatmap <- function(expr_matrix, geneA, geneB, pseudotime_vec, window_size = 30) {
  # Sort cells by pseudotime
  order_idx <- order(pseudotime_vec)
  exprA <- expr_matrix[geneA, order_idx]
  exprB <- expr_matrix[geneB, order_idx]
  time_sorted <- pseudotime_vec[order_idx]

  n_cells <- length(exprA)
  if (n_cells < window_size * 2) stop("Not enough cells for the specified window size.")

  rolling_cor <- numeric(n_cells - window_size)
  window_time <- numeric(n_cells - window_size)

  # Calculate rolling biweight midcorrelation over the moving pseudotime window
  for (i in seq_len(n_cells - window_size)) {
    idx <- i:(i + window_size - 1)
    rolling_cor[i] <- WGCNA::bicor(exprA[idx], exprB[idx])
    window_time[i] <- mean(time_sorted[idx])
  }

  plot_df <- data.frame(Time = window_time, Correlation = rolling_cor)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Time, y = Correlation, fill = Correlation)) +
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
