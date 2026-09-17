#' @title Input Extraction and Engine Dispatch for bicorX
#' @description Consolidates object-type-specific input handling (matrix,
#' data.frame, SummarizedExperiment, Seurat, MultiAssayExperiment) into one
#' place, and provides a single internal function that routes to the three
#' compiled engines. This replaces what used to be two parallel, largely
#' duplicated S3 dispatch trees:
#'   - `run_dcanr.*` (methods on `object`) in bicorX.R, and
#'   - `run_bmht.*`/`run_rosdet.*`/`run_bmkc.*` (methods on `expr_matrix`) in
#'     the file formerly called Dispatcher_S3.R,
#' plus separate, near-identical input-munging helpers (`.se_dispatcher`,
#' `.mae_dispatcher`, and inline logic inside `run_dcanr.Seurat`). A matrix
#' call used to travel five hops deep
#' (`run_dcanr()` -> `run_dcanr.matrix()` -> `do.call(run_bmht, args_list)` ->
#' `run_bmht()` -> `run_bmht.default()` -> `.run_bmht_matrix()`) partly
#' because of this duplication. Now: `run_bmht()` -> `run_bicordcea()` ->
#' `run_bicordcea.<class>()` -> `.run_bmht_matrix()`.
#'
#' @name bicorX_Dispatchers
NULL

# ===================================================================
# 1. .extract_inputs(): object-type-specific extraction only.
#    Every method returns list(expr_matrix, condition, species).
# ===================================================================

#' @keywords internal
.extract_inputs <- function(object, condition, ...) {
  UseMethod(".extract_inputs")
}

#' @keywords internal
.extract_inputs.matrix <- function(object, condition, species = NULL, ...) {
  list(expr_matrix = object, condition = condition, species = species)
}

#' @keywords internal
.extract_inputs.data.frame <- function(object, condition, species = NULL, ...) {
  list(expr_matrix = as.matrix(object), condition = condition, species = species)
}

#' @keywords internal
.extract_inputs.SummarizedExperiment <- function(object, condition,
                                                  assay_name = "counts",
                                                  species_col = NULL, ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("The 'SummarizedExperiment' package must be installed to use this feature.")
  }

  if (!assay_name %in% SummarizedExperiment::assayNames(object)) {
    stop(sprintf("Assay '%s' not found. Available assays: %s",
                 assay_name, paste(SummarizedExperiment::assayNames(object), collapse = ", ")))
  }
  mat <- as.matrix(SummarizedExperiment::assay(object, assay_name))

  if (!condition %in% colnames(SummarizedExperiment::colData(object))) {
    stop(sprintf("Condition column '%s' not found in colData.", condition))
  }
  cond_vec <- SummarizedExperiment::colData(object)[[condition]]

  species_vec <- NULL
  if (!is.null(species_col)) {
    if (!species_col %in% colnames(SummarizedExperiment::rowData(object))) {
      stop(sprintf("Species column '%s' not found in rowData.", species_col))
    }
    species_vec <- as.character(SummarizedExperiment::rowData(object)[[species_col]])
  }

  list(expr_matrix = mat, condition = cond_vec, species = species_vec)
}

#' @keywords internal
.extract_inputs.Seurat <- function(object, condition,
                                    assay_name = "RNA", layer_name = "counts",
                                    species_col = NULL, ...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("The 'SeuratObject' package must be installed to natively ingest Seurat data.")
  }

  meta_data <- if (is.list(object) && !isS4(object)) object$meta.data else as.data.frame(object[[]])

  if (!condition %in% colnames(meta_data)) {
    stop(sprintf("Condition column '%s' not found in Seurat metadata.", condition))
  }
  cond_vec <- meta_data[[condition]]

  mat <- if (is.list(object) && !isS4(object)) {
    object$matrix # test-environment fallback (a plain list standing in for a Seurat object)
  } else {
    as.matrix(SeuratObject::GetAssayData(object, assay = assay_name, layer = layer_name))
  }

  species_vec <- NULL
  if (!is.null(species_col)) {
    if (species_col %in% colnames(meta_data)) {
      species_vec <- as.character(meta_data[[species_col]])
    } else if (isS4(object) && species_col %in% colnames(object[[assay_name]]@meta.features)) {
      species_vec <- as.character(object[[assay_name]]@meta.features[[species_col]])
    } else {
      stop(sprintf("Species column '%s' not found in Seurat metadata.", species_col))
    }
  }

  list(expr_matrix = mat, condition = cond_vec, species = species_vec)
}

#' @keywords internal
.extract_inputs.MultiAssayExperiment <- function(object, condition,
                                                  rna_assay = "RNA", epi_assay = "ATAC",
                                                  assay_name = "counts", ...) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("The 'MultiAssayExperiment' package must be installed to use this feature.")
  }

  available_exps <- names(MultiAssayExperiment::experiments(object))

  if (!rna_assay %in% available_exps) {
    stop(sprintf("RNA layer '%s' not found in MultiAssayExperiment assay lists.", rna_assay))
  }
  if (!epi_assay %in% available_exps) {
    stop(sprintf("Epigenetic layer '%s' not found in MultiAssayExperiment assay lists.", epi_assay))
  }

  rna_se <- object[[rna_assay]]
  epi_se <- object[[epi_assay]]

  rna_mat <- as.matrix(SummarizedExperiment::assay(rna_se, assay_name))
  epi_mat <- as.matrix(SummarizedExperiment::assay(epi_se, assay_name))

  rownames(rna_mat) <- paste0("RNA_", rownames(rna_mat))
  rownames(epi_mat) <- paste0("Epi_", rownames(epi_mat))

  if (!condition %in% colnames(MultiAssayExperiment::colData(object))) {
    stop(sprintf("Condition metadata column '%s' not found in colData.", condition))
  }
  cond_vec <- MultiAssayExperiment::colData(object)[[condition]]

  common_samples <- intersect(colnames(rna_mat), colnames(epi_mat))
  if (length(common_samples) < 3) {
    stop("Fewer than 3 samples are shared between the RNA and epigenetic layers.")
  }

  rna_mat <- rna_mat[, common_samples, drop = FALSE]
  epi_mat <- epi_mat[, common_samples, drop = FALSE]

  matched_indices <- match(common_samples, rownames(MultiAssayExperiment::colData(object)))
  cond_vec <- cond_vec[matched_indices]

  combined_matrix <- rbind(rna_mat, epi_mat)
  layer_tags <- c(rep("RNA_Layer", nrow(rna_mat)), rep("Epigenetic_Layer", nrow(epi_mat)))

  list(expr_matrix = combined_matrix, condition = cond_vec, species = layer_tags)
}

# ===================================================================
# 2. .run_engine(): the ONE place that routes to a compiled engine and
#    tags the result class. Used by every run_bicordcea.<class>() method.
# ===================================================================

#' @keywords internal
#'
#' `species` doubles as two different things depending on context: a
#' taxonomic label (microbiome mode, used to compute CLR per taxonomic
#' group) or an omics-layer tag (multi-omics mode, e.g. "RNA" vs "ATAC",
#' used the same way - CLR per layer, and to classify pairs as intra- vs
#' inter-layer for `pair_type_filter`). Both uses are the same underlying
#' mechanism: partition the rows of `expr_matrix` into named groups. `layer`
#' is accepted here as a synonym for `species` so multi-omics callers can
#' use the name that actually describes what they're doing, without a
#' second, duplicated code path. If both are supplied, `layer` wins (with a
#' warning) since it's the more specific name for whichever call site used it.
.run_engine <- function(method, expr_matrix, condition, species, transform, ...) {
  dots <- list(...)
  if (!is.null(dots$layer)) {
    if (!is.null(species)) {
      warning("Both `species` and `layer` were supplied; using `layer`. ",
              "They are the same underlying grouping key - `layer` is just ",
              "the clearer name for multi-omics use.", call. = FALSE)
    }
    species <- dots$layer
    dots$layer <- NULL
  }

  engine_fn <- switch(method,
                       rosdet = .run_rosdet_matrix,
                       bmht   = .run_bmht_matrix,
                       bmkc   = .run_bmkc_matrix,
                       stop("Unknown method: ", method))
  res <- do.call(engine_fn, c(
    list(expr_matrix = expr_matrix, condition = condition,
         species = species, transform = transform),
    dots
  ))
  class(res) <- c(paste0(method, "_result"), "list")
  res
}

# ===================================================================
# 3. Thin, non-generic convenience wrappers.
#    No S3 dispatch happens here anymore - run_bicordcea() does all of it.
#    First-argument name (expr_matrix) is preserved from the original
#    contract so existing calls (including the ones fixed in Phase 1's
#    advanced_wrappers.R) keep working unchanged.
# ===================================================================

#' Run BMHT Analysis
#' @param expr_matrix Numeric matrix (features by samples), data frame, or
#'   a `SummarizedExperiment`, `Seurat` or `MultiAssayExperiment` object.
#' @param condition Two-level factor giving each sample's condition, or
#'   the name of a metadata column when `expr_matrix` is an S4 object.
#' @param ... Further arguments passed to the BMHT engine, e.g.
#'   `n_permutations`, `half_threshold`, `transform`, `species`,
#'   `workers`, `seed`. See Details.
#' @return A `bmht_result` object: a list whose `results` element is a
#'   data frame of per-gene `DC_Score`, `P_value` and `FDR`.
#' @export
run_bmht <- function(expr_matrix, condition, ...) {
  run_bicordcea(object = expr_matrix, condition = condition, method = "bmht", ...)
}

#' Run ROS-DET Analysis
#' @param expr_matrix Numeric matrix (features by samples), data frame, or
#'   a `SummarizedExperiment`, `Seurat` or `MultiAssayExperiment` object.
#' @param condition Two-level factor giving each sample's condition, or
#'   the name of a metadata column when `expr_matrix` is an S4 object.
#' @param ... Further arguments passed to the ROS-DET engine, e.g.
#'   `min_delta`, `n_permutations`, `significance_method`, `weight_mode`,
#'   `transform`, `species`, `pair_type_filter`, `workers`, `seed`.
#' @return A `rosdet_result` object: a list whose `results` element is a
#'   data frame of per-pair `Distance_Score`, `P_value` and `FDR`.
#' @export
run_rosdet <- function(expr_matrix, condition, ...) {
  run_bicordcea(object = expr_matrix, condition = condition, method = "rosdet", ...)
}

#' Run BMKC Analysis
#' @param expr_matrix Numeric matrix (features by samples), data frame, or
#'   a `SummarizedExperiment`, `Seurat` or `MultiAssayExperiment` object.
#' @param condition Two-level factor giving each sample's condition, or
#'   the name of a metadata column when `expr_matrix` is an S4 object.
#' @param ... Further arguments passed to the BMKC engine, e.g. `T1`,
#'   `T2`, `gamma`, `min_clique_size`, `n_permutations`, `transform`,
#'   `species`, `workers`, `seed`.
#' @return A `bmkc_result` object: a list of detected modules with
#'   optional global significance (`significance_p`) when
#'   `n_permutations > 0`.
#' @export
run_bmkc <- function(expr_matrix, condition, ...) {
  run_bicordcea(object = expr_matrix, condition = condition, method = "bmkc", ...)
}
