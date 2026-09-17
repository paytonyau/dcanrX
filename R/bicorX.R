#' @title Unified Interface for bicorX
#' @description Main entry point for the bicorX package.
#' Provides a consistent interface for bulk RNA-seq (matrix or data.frame),
#' single-cell transcriptomics (via SummarizedExperiment, SingleCellExperiment,
#' or Seurat), and multi-omic/epigenetic datasets (via MultiAssayExperiment).
#' Supports rCLR transformation and cross-species mode for microbial studies.
#'
#' @name bicorX
NULL

# ===================================================================
# S3 Master Wrapper
# ===================================================================

#' Run Differential Co-expression Analysis with bicor
#'
#' Unified function to run BMKC, BMHT, or ROS-DET methods on any supported
#' input type. Was named `run_dcanr()` in earlier versions; that name
#' collides with the existing Bioconductor package `dcanr`
#' (Bhuva et al.), so this is the new name going forward. `run_dcanr()` is
#' kept as a deprecated alias - it warns once and forwards here.
#'
#' Class-specific extraction is handled by `.extract_inputs()` (see
#' Dispatcher_S3.R) so every method below is a thin wrapper: extract
#' `(expr_matrix, condition, species)`, then hand off to `.run_engine()`,
#' which is the single place that calls into the compiled BMHT/ROS-DET/BMKC
#' engines.
#'
#' @param object Numeric matrix, data.frame, SummarizedExperiment, Seurat, or
#'   MultiAssayExperiment.
#' @param condition Factor/character vector (matrix/data.frame input) or a
#'   metadata column name (SummarizedExperiment/Seurat/MultiAssayExperiment
#'   input).
#' @param ... Additional arguments passed on to the chosen method. Most
#'   importantly `method`, one of `"rosdet"` (default), `"bmht"` or
#'   `"bmkc"`; plus that method's own engine arguments
#'   (e.g. `n_permutations`, `min_delta`, `T1`/`T2`, `seed`). Also accepts
#'   `layer`, a synonym for `species` (see Details) that reads more
#'   naturally for multi-omics use - see `run_multiomic()`.
#' @return Result object from the chosen method.
#' @details
#' `species` (or its synonym `layer`) is a grouping vector, one entry per
#' feature/row, used two ways depending on context: as taxonomic labels for
#' per-group CLR normalization in microbiome data, or as omics-layer tags
#' (e.g. `"RNA"`/`"ATAC"`) for the same per-group CLR *and* to classify
#' pairs as intra- vs. inter-layer for ROS-DET's `pair_type_filter`. Both
#' uses are the same underlying mechanism - see
#' `.compositional_transform()` for the CLR details, and `run_multiomic()`
#' for the multi-omics entry point that populates this automatically.
#' @export
run_bicordcea <- function(object, condition, ...) {
  UseMethod("run_bicordcea")
}

#' @rdname run_bicordcea
#' @export
run_dcanr <- function(object, condition, ...) {
  warning("`run_dcanr()` is deprecated and will be removed in a future ",
          "version; use `run_bicordcea()` instead. The name collided with ",
          "the existing Bioconductor package 'dcanr'.", call. = FALSE)
  run_bicordcea(object = object, condition = condition, ...)
}

# ===================================================================
# Matrix Method (Bulk RNA-seq)
# ===================================================================

#' @export
run_bicordcea.matrix <- function(object, condition,
                                  method = c("rosdet", "bmht", "bmkc"),
                                  transform = c("none", "clr", "rclr", "log1p"),
                                  species = NULL, ...) {
  method <- match.arg(method)
  transform <- match.arg(transform)

  extracted <- .extract_inputs(object, condition, species = species)
  .run_engine(method, extracted$expr_matrix, extracted$condition,
              extracted$species, transform, ...)
}

# ===================================================================
# data.frame Method (previously unsupported despite package docs
# advertising it - see Phase 1 review notes)
# ===================================================================

#' @export
run_bicordcea.data.frame <- function(object, condition,
                                      method = c("rosdet", "bmht", "bmkc"),
                                      transform = c("none", "clr", "rclr", "log1p"),
                                      species = NULL, ...) {
  method <- match.arg(method)
  transform <- match.arg(transform)

  extracted <- .extract_inputs(object, condition, species = species)
  .run_engine(method, extracted$expr_matrix, extracted$condition,
              extracted$species, transform, ...)
}

# ===================================================================
# SummarizedExperiment Method
# ===================================================================

#' @export
run_bicordcea.SummarizedExperiment <- function(object, condition,
                                                method = c("rosdet", "bmht", "bmkc"),
                                                assay_name = "counts",
                                                transform = c("none", "clr", "rclr", "log1p"),
                                                species_col = NULL, ...) {
  method <- match.arg(method)
  transform <- match.arg(transform)

  extracted <- .extract_inputs(object, condition, assay_name = assay_name,
                                species_col = species_col)
  .run_engine(method, extracted$expr_matrix, extracted$condition,
              extracted$species, transform, ...)
}

# ===================================================================
# Seurat Method (Single-Cell RNA-seq)
# ===================================================================

#' @export
run_bicordcea.Seurat <- function(object, condition,
                                  method = c("rosdet", "bmht", "bmkc"),
                                  assay_name = "RNA",
                                  layer_name = "counts",
                                  transform = c("none", "clr", "rclr", "log1p"),
                                  species_col = NULL, ...) {
  method <- match.arg(method)
  transform <- match.arg(transform)

  extracted <- .extract_inputs(object, condition, assay_name = assay_name,
                                layer_name = layer_name, species_col = species_col)
  .run_engine(method, extracted$expr_matrix, extracted$condition,
              extracted$species, transform, ...)
}

# ===================================================================
# MultiAssayExperiment Method (Unified Multi-Omics / Epigenetics)
# ===================================================================

#' @export
run_bicordcea.MultiAssayExperiment <- function(object, condition,
                                                method = c("rosdet", "bmht", "bmkc"),
                                                rna_assay = "RNA",
                                                epi_assay = "ATAC",
                                                assay_name = "counts",
                                                transform = c("none", "clr", "rclr", "log1p"),
                                                ...) {
  method <- match.arg(method)
  transform <- match.arg(transform)

  extracted <- .extract_inputs(object, condition, rna_assay = rna_assay,
                                epi_assay = epi_assay, assay_name = assay_name)
  .run_engine(method, extracted$expr_matrix, extracted$condition,
              extracted$species, transform, ...)
}
