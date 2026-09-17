#' @title Pseudobulk Aggregation for Single-Cell & Spatial Transcriptomics
#' @importFrom stats na.omit
#' @description Creates pseudobulk expression profiles by aggregating cells/spots
#' from the same biological sample, condition, or spatial microenvironment.
#' Includes strict QC filters for cell counts and low-expression genes.
#'
#' @name pseudobulk
NULL

#' Create Pseudobulk Matrix from SingleCellExperiment
#' @param sce A `SingleCellExperiment` (or object coercible to one)
#'   containing single-cell or spatial counts.
#' @param sample_col Name of the `colData` column identifying independent
#'   biological replicates (patients, mice, wells) to aggregate within.
#'   This should not be a cluster or cell-type label.
#' @param condition_col Name of the `colData` column giving the two-level
#'   condition for each cell.
#' @param assay Name of the assay to aggregate. Default `"counts"`.
#' @param min_cells Minimum number of cells required for a sample to be
#'   retained. Default `10`.
#' @return A list with `expr_matrix` (features by pseudobulk samples) and
#'   `condition` (one entry per pseudobulk sample).
#' @export
make_pseudobulk <- function(sce,
                            sample_col = "sample",
                            condition_col = "condition",
                            assay = "counts",
                            min_cells = 10) {

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("The 'SingleCellExperiment' package must be installed to use this feature.")
  }

  if (!inherits(sce, "SingleCellExperiment") && !inherits(sce, "SummarizedExperiment")) {
    stop("Input must be a SingleCellExperiment or SummarizedExperiment object.")
  }

  coldata <- as.data.frame(SummarizedExperiment::colData(sce))

  if (!sample_col %in% colnames(coldata)) stop(sprintf("Replicate identifier '%s' absent from metadata.", sample_col))
  if (!condition_col %in% colnames(coldata)) stop(sprintf("Experimental condition '%s' absent from metadata.", condition_col))

  valid_cells <- !is.na(coldata[[sample_col]]) & !is.na(coldata[[condition_col]])
  if (sum(!valid_cells) > 0) {
    sce <- sce[, valid_cells, drop = FALSE]
    coldata <- coldata[valid_cells, , drop = FALSE]
  }

  delim <- "___"
  group_factor <- paste(coldata[[sample_col]], coldata[[condition_col]], sep = delim)

  pseudo_se <- scuttle::summarizeAssayByGroup(
    sce,
    ids = group_factor,
    assay.type = assay,
    statistics = "sum"
  )

  expr_matrix <- SummarizedExperiment::assay(pseudo_se, 1)
  n_cells_vec <- pseudo_se$ncells
  if (is.null(n_cells_vec)) n_cells_vec <- pseudo_se$n

  valid_bulks <- n_cells_vec >= min_cells
  if (sum(!valid_bulks) > 0) {
    message(sprintf("Single-Cell QA Filter: Dropped %d sample pseudobulks with fewer than %d cells.", sum(!valid_bulks), min_cells))
    expr_matrix <- expr_matrix[, valid_bulks, drop = FALSE]
    n_cells_vec <- n_cells_vec[valid_bulks]
  }

  if (ncol(expr_matrix) == 0) stop("Quality Control Stop: All constructed sample pseudobulks failed your min_cells boundary threshold.")

  group_info <- strsplit(colnames(expr_matrix), delim)
  sample_vec <- sapply(group_info, `[`, 1)
  condition_vec <- sapply(group_info, `[`, 2)

  list(
    expr_matrix = expr_matrix,
    condition   = factor(condition_vec),
    sample      = sample_vec,
    n_cells     = n_cells_vec,
    group_names = colnames(expr_matrix)
  )
}

# ===================================================================
# Main Single-Cell / Spatial Wrapper
# ===================================================================

#' Run Differential Co-expression on Single-Cell or Spatial Data
#' @param sce A `SingleCellExperiment` (or coercible) object.
#' @param condition_col Name of the `colData` column giving the two-level
#'   condition. Default `"condition"`.
#' @param sample_col Name of the `colData` column identifying independent
#'   biological replicates to aggregate within. Default `"sample"`.
#' @param method Which method to run: `"rosdet"`, `"bmht"`, or `"bmkc"`.
#' @param assay Name of the assay to aggregate. Default `"counts"`.
#' @param transform Compositional transform to apply: `"none"`, `"clr"`,
#'   `"rclr"`, or `"log1p"`.
#' @param split_by Optional `colData` column to split the analysis by,
#'   running the chosen method separately within each level.
#' @param species Optional grouping vector (one entry per feature) used
#'   for per-group CLR normalization and layer-aware pair filtering.
#' @param min_cells Minimum cells required per pseudobulk sample.
#'   Default `10`.
#' @param min_prop Minimum proportion of samples in which a feature must
#'   be detected to be retained. Default `0.10`.
#' @param ... Further arguments passed to the chosen method's engine.
#' @return The chosen method's result object, or a named list of such
#'   objects when `split_by` is supplied.
#' @export
run_dcanr_pseudobulk <- function(sce,
                                 condition_col = "condition",
                                 sample_col = "sample",
                                 method = c("rosdet", "bmht", "bmkc"),
                                 assay = "counts",
                                 transform = c("none", "clr", "rclr", "log1p"),
                                 split_by = NULL,
                                 species = NULL,
                                 min_cells = 10,
                                 min_prop = 0.10,
                                 ...) {

  method <- match.arg(method)
  transform <- match.arg(transform)

  if (!is.null(split_by)) {
    coldata <- as.data.frame(SummarizedExperiment::colData(sce))
    if (!split_by %in% colnames(coldata)) stop(sprintf("Split-by grouping column '%s' not found.", split_by))

    levels_to_test <- unique(na.omit(coldata[[split_by]]))
    message(sprintf("Stratified Split Ingestion Enabled: Running %s over %d distinct structural subsets of '%s'...",
                    toupper(method), length(levels_to_test), split_by))

    results_list <- list()

    for (lvl in levels_to_test) {
      message(sprintf("\n>>> Processing structural partition block: %s", lvl))
      sce_sub <- sce[, coldata[[split_by]] == lvl, drop = FALSE]

      results_list[[lvl]] <- tryCatch({
        run_dcanr_pseudobulk(
          sce = sce_sub, condition_col = condition_col, sample_col = sample_col,
          method = method, assay = assay, transform = transform, split_by = NULL,
          species = species, min_cells = min_cells, min_prop = min_prop, ...
        )
      }, error = function(e) {
        warning(sprintf("Execution failed for subset component '%s': %s", lvl, e$message))
        return(NULL)
      })
    }
    return(Filter(Negate(is.null), results_list))
  }

  pseudo <- make_pseudobulk(
    sce = sce, sample_col = sample_col, condition_col = condition_col,
    assay = assay, min_cells = min_cells
  )

  cond_table <- table(pseudo$condition)
  if (any(cond_table < 3)) {
    stop(sprintf("Replicate Insufficiency Alert: Found only %s. All analytical configurations require at least 3 high-integrity pseudobulks per cohort layer.",
                 paste(names(cond_table), cond_table, sep=":", collapse=" | ")))
  }

  expressed_prop <- rowMeans(pseudo$expr_matrix > 0)
  keep_genes <- expressed_prop >= min_prop

  if (sum(!keep_genes) > 0) {
    message(sprintf("Expression Sparsity Filter: Removing %d components active in fewer than %.0f%% of pseudobulk cohorts.",
                    sum(!keep_genes), min_prop * 100))
    pseudo$expr_matrix <- pseudo$expr_matrix[keep_genes, , drop = FALSE]
  }

  message(sprintf("Aggregation Final Summary: Unified %d genes across %d structural pseudobulks (aggregating %d physical features) | composition = %s",
                  nrow(pseudo$expr_matrix), length(pseudo$condition), sum(pseudo$n_cells), toupper(transform)))

  species_vec <- species
  if (!is.null(species) && length(species) == 1 && is.character(species)) {
    if (species %in% colnames(SummarizedExperiment::rowData(sce))) {
      species_vec <- SummarizedExperiment::rowData(sce)[[species]][keep_genes]
    }
  } else if (!is.null(species) && length(species) == nrow(sce)) {
    species_vec <- species[keep_genes]
  }

  run_bicordcea(
    object = pseudo$expr_matrix,
    condition = pseudo$condition,
    method = method,
    transform = transform,
    species = species_vec,
    ...
  )
}
