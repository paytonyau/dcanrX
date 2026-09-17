#' @title BMHT Visualization Suite
#' @description Plotting functions for Biweight Midcorrelation Hub Topology (BMHT) results and causal-driver overlays.
#' @name bmht_plots
NULL

# Fix CRAN package check notes regarding unquoted column name structures
utils::globalVariables(c("DC_Score", "negLogFDR", "is_sig", "Gene", "species"))

# ==============================================================================
# CATEGORY 1: UNIVERSAL DIAGNOSTICS
# ==============================================================================

#' Rewiring Volcano Plot (BMHT)
#' @param bmht_result A `bmht_result` object from [run_bmht()].
#' @param n_label Number of top genes to label. Default `10`.
#' @param ... Currently unused; accepted for S3 compatibility.
#' @return A `ggplot` object.
#' @export
plot_rewiring_volcano <- function(bmht_result, n_label = 10, ...) {
  df <- bmht_result$results
  df$negLogFDR <- -log10(df$FDR)
  df$is_sig <- df$FDR < bmht_result$significance_level

  top_labels <- df[order(df$FDR, -df$DC_Score), ]
  top_labels <- head(top_labels[top_labels$is_sig, ], n_label)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$DC_Score, y = .data$negLogFDR)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$is_sig), size = 3, alpha = 0.7) +
    ggrepel::geom_text_repel(data = top_labels, ggplot2::aes(label = .data$Gene),
                             size = 4, box.padding = 0.5, max.overlaps = Inf) +
    ggplot2::scale_color_manual(values = c("TRUE" = "red3", "FALSE" = "grey70")) +
    ggplot2::labs(title = "Rewiring Volcano Plot (BMHT)",
                  subtitle = "Scale-Free Topological Distance Matrix Profile",
                  x = "Differential Co-expression Score (MAD Normalized)", y = "-log10(FDR)",
                  color = "Significant") +
    ggplot2::theme_minimal(base_size = 14, base_family = "sans")
}

#' Top Gene Connectivity Bar Chart (BMHT)
#' @param bmht_result A `bmht_result` object from [run_bmht()].
#' @param n_top Number of top-ranked genes to display. Default `15`.
#' @param ... Currently unused; accepted for S3 compatibility.
#' @return A `ggplot` object.
#' @export
plot_top_gene_bar <- function(bmht_result, n_top = 15, ...) {
  top_genes <- head(bmht_result$results, n_top)
  if (!"species" %in% names(top_genes)) top_genes$species <- "All"

  ggplot2::ggplot(top_genes, ggplot2::aes(x = reorder(.data$Gene, .data$DC_Score), y = .data$DC_Score, fill = .data$species)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = paste("Top", n_top, "Most Rewired Genes (BMHT)"),
                  subtitle = "Radical drivers ranked by scale-free metric authority",
                  x = "Gene", y = "Normalized DC Score") +
    ggplot2::theme_minimal(base_size = 13, base_family = "sans") +
    ggplot2::theme(legend.position = ifelse(length(unique(top_genes$species)) > 1, "bottom", "none"))
}

# ==============================================================================
# CORE S3 DISPATCHER
# ==============================================================================

#' Plot Method for BMHT Results
#'
#' As of the Phase 3 surface reduction, only three plot types are provided
#' directly: `"volcano"` (default), `"bar"`, and `"topology"` (a differential
#' network shared with BMKC results). For other visualizations (drivers,
#' causal networks, spatial/UMAP overlays, taxa/microbiome breakdowns), see
#' the "Custom plots from bicorX results" vignette, which shows how to
#' build them directly from `x$results` and `x$data` using ggplot2/igraph/
#' Seurat - the same approach the removed convenience functions used
#' internally.
#' @method plot bmht_result
#' @param x A `bmht_result` object from [run_bmht()].
#' @param type Which plot to draw: `"default"`/`"volcano"` for the
#'   rewiring volcano, `"bar"` for the top-gene bar chart, or
#'   `"topology"` for the differential network.
#' @param ... Passed to the underlying plotting function.
#' @return A `ggplot` object.
#' @export
plot.bmht_result <- function(x, type = c("default", "volcano", "bar", "topology"), ...) {
  type <- match.arg(type)

  switch(type,
         "default"  = plot_rewiring_volcano(x, ...),
         "volcano"  = plot_rewiring_volcano(x, ...),
         "bar"      = plot_top_gene_bar(x, ...),
         "topology" = plot_differential_network(x, ...)
  )
}
