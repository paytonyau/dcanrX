#' @title Pseudobulk Aggregation for Single-Cell & Spatial Transcriptomics
#' @description Creates pseudobulk expression profiles by aggregating cells/spots
#' from the same biological sample, condition, or spatial microenvironment.
#' Includes strict QC filters for cell counts and low-expression genes.
#'
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom SummarizedExperiment colData assay rowData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @name pseudobulk
NULL

#' Create Pseudobulk Matrix from SingleCellExperiment
#'
#' @param sce A `SingleCellExperiment` object.
#' @param sample_col Column name identifying biological replicates.
#' @param condition_col Column name containing the experimental condition.
#' @param assay Which assay to aggregate (default: `"counts"`).
#' @param min_cells Minimum number of cells required to form a valid pseudobulk.
#' @return A list containing the pseudobulked matrix and metadata.
#' @export
make_pseudobulk <- function(sce,
                            sample_col = "sample",
                            condition_col = "condition",
                            assay = "counts",
                            min_cells = 10) {

  # 1. Validation
  if (!inherits(sce, "SingleCellExperiment") && !inherits(sce, "SummarizedExperiment")) {
    stop("Input must be a SingleCellExperiment or SummarizedExperiment object.")
  }

  coldata <- as.data.frame(SummarizedExperiment::colData(sce))

  if (!sample_col %in% colnames(coldata)) stop(sprintf("Column '%s' not found in colData.", sample_col))
  if (!condition_col %in% colnames(coldata)) stop(sprintf("Column '%s' not found in colData.", condition_col))

  # 2. NA Scrubbing (Crucial for Spatial dropouts)
  valid_cells <- !is.na(coldata[[sample_col]]) & !is.na(coldata[[condition_col]])
  if (sum(!valid_cells) > 0) {
    sce <- sce[, valid_cells]
    coldata <- coldata[valid_cells, , drop = FALSE]
  }

  # 3. Create composite grouping factor
  delim <- "___"
  group_factor <- paste(coldata[[sample_col]], coldata[[condition_col]], sep = delim)

  # 4. Perform Aggregation via scuttle
  pseudo_se <- scuttle::summarizeAssayByGroup(
    sce,
    ids = group_factor,
    assay.type = assay,
    statistics = "sum"
  )

  # 5. Extract Data
  expr_matrix <- SummarizedExperiment::assay(pseudo_se, 1)
  n_cells_vec <- pseudo_se$ncells
  if (is.null(n_cells_vec)) n_cells_vec <- pseudo_se$n # fallback

  # 6. QC Filter: Drop pseudobulks with too few cells
  valid_bulks <- n_cells_vec >= min_cells
  if (sum(!valid_bulks) > 0) {
    message(sprintf("Dropped %d pseudobulks formed by fewer than %d cells.", sum(!valid_bulks), min_cells))
    expr_matrix <- expr_matrix[, valid_bulks, drop = FALSE]
    n_cells_vec <- n_cells_vec[valid_bulks]
  }

  if (ncol(expr_matrix) == 0) stop("All pseudobulks failed the `min_cells` filter.")

  # 7. Parse metadata safely
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
#'
#' @param sce SingleCellExperiment object
#' @param condition_col Condition column name
#' @param sample_col Sample column name
#' @param method Analysis method ("rosdet", "bmht", or "bmkc")
#' @param assay Assay to use (default: "counts")
#' @param transform Compositional transformation (default: "clr")
#' @param split_by Optional column name for running independent spatial/cell-type networks.
#' @param species Optional character string or column name for cross-species mode.
#' @param min_cells Minimum cells per pseudobulk (default: 10).
#' @param min_prop Minimum proportion of pseudobulks a gene must be expressed in to be kept (default: 0.10).
#' @param ... Additional arguments passed to the chosen analysis function
#' @export
run_dcanr_pseudobulk <- function(sce,
                                 condition_col = "condition",
                                 sample_col = "sample",
                                 method = c("rosdet", "bmht", "bmkc"),
                                 assay = "counts",
                                 transform = c("none", "clr", "log1p"),
                                 split_by = NULL,
                                 species = NULL,
                                 min_cells = 10,
                                 min_prop = 0.10,
                                 ...) {

  method <- match.arg(method)
  transform <- match.arg(transform)

  # -------------------------------------------------------------------
  # SPATIAL / CELL-TYPE SPLIT ROUTINE
  # -------------------------------------------------------------------
  if (!is.null(split_by)) {
    coldata <- as.data.frame(SummarizedExperiment::colData(sce))
    if (!split_by %in% colnames(coldata)) stop(sprintf("split_by column '%s' not found.", split_by))

    levels_to_test <- unique(na.omit(coldata[[split_by]]))
    message(sprintf("Split mode activated. Running %s across %d levels of '%s'...",
                    toupper(method), length(levels_to_test), split_by))

    results_list <- list()

    for (lvl in levels_to_test) {
      message(sprintf("\n>>> Processing subset: %s", lvl))
      sce_sub <- sce[, coldata[[split_by]] == lvl]

      results_list[[lvl]] <- tryCatch({
        run_dcanr_pseudobulk(
          sce = sce_sub, condition_col = condition_col, sample_col = sample_col,
          method = method, assay = assay, transform = transform, split_by = NULL,
          species = species, min_cells = min_cells, min_prop = min_prop, ...
        )
      }, error = function(e) {
        warning(sprintf("Analysis failed for level '%s': %s", lvl, e$message))
        return(NULL)
      })
    }
    # Clean list of NULLs before returning
    return(Filter(Negate(is.null), results_list))
  }

  # -------------------------------------------------------------------
  # STANDARD ROUTINE (Single Run)
  # -------------------------------------------------------------------

  pseudo <- make_pseudobulk(
    sce = sce, sample_col = sample_col, condition_col = condition_col,
    assay = assay, min_cells = min_cells
  )

  # QC Guardrail: Minimum Samples Check
  cond_table <- table(pseudo$condition)
  if (any(cond_table < 3)) {
    stop(sprintf("Not enough samples after pseudobulking! Found %s. All methods require at least 3 valid pseudobulks per condition. Try lowering 'min_cells' or merging samples.",
                 paste(names(cond_table), cond_table, sep=":", collapse=" | ")))
  }

  # QC Guardrail: Low Expression Filter (Saves RAM, boosts FDR power)
  expressed_prop <- rowMeans(pseudo$expr_matrix > 0)
  keep_genes <- expressed_prop >= min_prop

  if (sum(!keep_genes) > 0) {
    message(sprintf("Filtering %d genes expressed in fewer than %.0f%% of pseudobulks...",
                    sum(!keep_genes), min_prop * 100))
    pseudo$expr_matrix <- pseudo$expr_matrix[keep_genes, , drop = FALSE]
  }

  message(sprintf("Final Matrix: %d genes across %d pseudobulks (representing %d cells) | transform = %s",
                  nrow(pseudo$expr_matrix), length(pseudo$condition), sum(pseudo$n_cells), toupper(transform)))

  # Resolve species logic
  species_vec <- species
  if (!is.null(species) && length(species) == 1 && is.character(species)) {
    if (species %in% colnames(SummarizedExperiment::rowData(sce))) {
      species_vec <- SummarizedExperiment::rowData(sce)[[species]][keep_genes]
    }
  } else if (!is.null(species) && length(species) == nrow(sce)) {
    species_vec <- species[keep_genes]
  }

  # Dispatch to the core generic
  run_dcanr(
    object = pseudo$expr_matrix,
    condition = pseudo$condition,
    method = method,
    transform = transform,
    species = species_vec,
    ...
  )
}
