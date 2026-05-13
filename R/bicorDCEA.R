#' @title Unified Interface for bicorDCEA
#' @description Main entry point for the bicorDCEA package.
#' Provides a consistent interface for both bulk RNA-seq (matrix) and
#' single-cell RNA-seq (via SummarizedExperiment / SingleCellExperiment).
#' Supports CLR transformation and cross-species mode for microbial studies.
#'
#' @name bicorDCEA
#' @importFrom methods setGeneric setMethod
#' @importFrom SummarizedExperiment assay colData
NULL

# ===================================================================
# Generic Definition
# ===================================================================

#' Run Differential Co-expression Analysis with bicor
#'
#' Unified function to run BMKC, BMHT, or ROS-DET methods on either
#' bulk (matrix) or single-cell (SummarizedExperiment) data.
#'
#' @param object Numeric matrix (genes × samples) or a
#'   `SummarizedExperiment` / `SingleCellExperiment` object.
#' @param condition For matrix: factor/character vector with 2 levels.
#'   For SE/SCE: character string naming the column in `colData`.
#' @param method Which method to run: `"rosdet"`, `"bmht"`, or `"bmkc"`.
#' @param assay Character. Name of the assay to use if `object` is a
#'   SummarizedExperiment (default `"logcounts"` or `"counts"`).
#' @param transform Character. Compositional transformation to apply
#'   before analysis. Use `"clr"` for microbial/metatranscriptomic data
#'   or highly sparse single-cell data. Default is `"none"`.
#' @param species Optional character vector indicating species for each gene
#'   (enables cross-species / multi-species analysis).
#' @param ... Additional arguments passed to the specific method
#'   (e.g. `T1`, `T2`, `half_threshold`, `n_permutations`, `workers`,
#'   `pair_type_filter`, etc.)
#' @return Result object from the chosen method
#' @export
setGeneric("run_dcanr", function(object, condition,
                                 method = c("rosdet", "bmht", "bmkc"),
                                 assay = "logcounts",
                                 transform = c("none", "clr", "log1p"),
                                 species = NULL, ...) {
  standardGeneric("run_dcanr")
})

# ===================================================================
# Matrix Method (Bulk RNA-seq)
# ===================================================================

#' @rdname run_dcanr
#' @export
setMethod("run_dcanr", "matrix", function(object, condition,
                                          method = c("rosdet", "bmht", "bmkc"),
                                          assay = "logcounts",
                                          transform = c("none", "clr", "log1p"),
                                          species = NULL, ...) {
  method <- match.arg(method)
  transform <- match.arg(transform)

  switch(method,
         rosdet = run_rosdet(object, condition, transform = transform, species = species, ...),
         bmht   = run_bmht(object, condition, transform = transform, species = species, ...),
         bmkc   = run_bmkc(object, condition, transform = transform, species = species, ...)
  )
})

# ===================================================================
# SummarizedExperiment / SingleCellExperiment Method
# ===================================================================

#' @rdname run_dcanr
#' @export
setMethod("run_dcanr", "SummarizedExperiment", function(object, condition,
                                                        method = c("rosdet", "bmht", "bmkc"),
                                                        assay = "logcounts",
                                                        transform = c("none", "clr", "log1p"),
                                                        species = NULL, ...) {
  method <- match.arg(method)
  transform <- match.arg(transform)

  if (!condition %in% names(colData(object))) {
    stop("Column '", condition, "' not found in colData(object). ",
         "Available columns: ", paste(names(colData(object)), collapse = ", "))
  }

  # Extract data
  expr_mat <- SummarizedExperiment::assay(object, assay)
  cond_vec <- colData(object)[[condition]]

  # Forward to matrix method
  run_dcanr(expr_mat, cond_vec, method = method,
            transform = transform, species = species, ...)
})

# SingleCellExperiment inherits from SummarizedExperiment
#' @rdname run_dcanr
#' @export
setMethod("run_dcanr", "SingleCellExperiment", getMethod("run_dcanr", "SummarizedExperiment"))

# ===================================================================
# Convenience Wrapper Functions (forward all new parameters)
# ===================================================================

#' @rdname run_dcanr
#' @export
run_rosdet <- function(object, condition, transform = "none", species = NULL, ...) {
  run_dcanr(object, condition, method = "rosdet",
            transform = transform, species = species, ...)
}

#' @rdname run_dcanr
#' @export
run_bmht <- function(object, condition, transform = "none", species = NULL, ...) {
  run_dcanr(object, condition, method = "bmht",
            transform = transform, species = species, ...)
}

#' @rdname run_dcanr
#' @export
run_bmkc <- function(object, condition, transform = "none", species = NULL, ...) {
  run_dcanr(object, condition, method = "bmkc",
            transform = transform, species = species, ...)
}

# ===================================================================
# Pseudobulk Convenience Function (also forwards new parameters)
# ===================================================================

#' @rdname run_dcanr
#' @export
run_dcanr_pseudobulk <- function(sce,
                                 condition_col = "condition",
                                 sample_col = "sample",
                                 method = c("rosdet", "bmht", "bmkc"),
                                 assay = "counts",
                                 transform = c("none", "clr", "log1p"),
                                 species = NULL,
                                 ...) {

  method <- match.arg(method)
  transform <- match.arg(transform)

  pseudo <- make_pseudobulk(
    sce,
    sample_col = sample_col,
    condition_col = condition_col,
    assay = assay
  )

  message(sprintf("Created %d pseudobulks (from %d cells) | transform = %s",
                  length(pseudo$condition), sum(pseudo$n_cells), transform))

  run_dcanr(
    object = pseudo$expr_matrix,
    condition = pseudo$condition,
    method = method,
    transform = transform,
    species = species,
    ...
  )
}
