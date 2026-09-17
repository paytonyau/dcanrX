# ===================================================================
# CAUSAL INFERENCE WRAPPERS (TF-GUIDED DIRECTIONALITY)
# ===================================================================

# Satisfy CRAN namespace checks regarding unquoted column mutations
utils::globalVariables(c("Is_TF", "Causal_Driver", "G1_is_TF", "G2_is_TF", "Direction"))

#' Annotate Causal Drivers (TF-Guided BMHT)
#' @param bmht_result A `bmht_result` object from [run_bmht()].
#' @param tf_list Character vector of transcription factor (or other
#'   regulator) gene names to test for enrichment among top-ranked hubs.
#' @param top_pct Proportion of top-ranked genes (by `DC_Score`) to treat
#'   as hubs. Default `0.05` (top 5\%).
#' @return The `bmht_result`'s results data frame with added `Is_TF` and
#'   `Causal_Driver` logical columns.
#' @export
identify_causal_drivers <- function(bmht_result, tf_list, top_pct = 0.05) {
  res <- bmht_result$results
  res <- res[order(res$DC_Score, decreasing = TRUE), ]

  threshold_idx <- max(1, round(nrow(res) * top_pct))
  top_hubs <- res$Gene[1:threshold_idx]

  res$Is_TF <- res$Gene %in% tf_list
  res$Causal_Driver <- ifelse(res$Gene %in% top_hubs & res$Is_TF, "Yes", "No")

  causal_df <- res[res$Causal_Driver == "Yes", ]

  message(sprintf("BMHT Causal Inference: Identified %d master TFs out of the top %d rewired hubs.",
                  nrow(causal_df), threshold_idx))

  bmht_result$results <- res
  bmht_result$causal_drivers <- causal_df

  return(bmht_result)
}

#' Annotate Direction for ROS-DET Switching Pairs
#' @param rosdet_result A `rosdet_result` object from [run_rosdet()].
#' @param tf_list Character vector of transcription factor (or other
#'   regulator) gene names.
#' @return The results data frame with added `G1_is_TF`, `G2_is_TF` and
#'   `Direction` columns indicating putative regulatory direction.
#' @export
annotate_rosdet_direction <- function(rosdet_result, tf_list) {
  res <- rosdet_result$results

  if (nrow(res) == 0) {
    message("No significant ROS-DET pairs to annotate.")
    return(rosdet_result)
  }

  res$G1_is_TF <- res$Gene1 %in% tf_list
  res$G2_is_TF <- res$Gene2 %in% tf_list

  res$Direction <- "Undirected"
  res$Direction[res$G1_is_TF & !res$G2_is_TF] <- paste(res$Gene1[res$G1_is_TF & !res$G2_is_TF], "->", res$Gene2[res$G1_is_TF & !res$G2_is_TF])
  res$Direction[!res$G1_is_TF & res$G2_is_TF] <- paste(res$Gene2[!res$G1_is_TF & res$G2_is_TF], "->", res$Gene1[!res$G1_is_TF & res$G2_is_TF])
  res$Direction[res$G1_is_TF & res$G2_is_TF]  <- "Co-Regulatory (TF <-> TF)"

  directed_count <- sum(res$Direction != "Undirected")
  message(sprintf("ROS-DET Direction Mapping: Assigned regulatory direction to %d switching pairs.", directed_count))

  rosdet_result$results <- res
  return(rosdet_result)
}

#' Annotate Master Regulators for BMKC Modules
#' @param bmkc_result A `bmkc_result` object from [run_bmkc()].
#' @param tf_list Character vector of transcription factor (or other
#'   regulator) gene names.
#' @return A data frame with one row per detected module, reporting the
#'   regulators found within it.
#' @export
annotate_bmkc_drivers <- function(bmkc_result, tf_list) {
  if (bmkc_result$n_modules == 0) stop("No modules found to annotate.")

  driver_list <- list()

  for (i in seq_len(bmkc_result$n_modules)) {
    module_genes <- bmkc_result$modules[[i]]
    tfs_in_module <- intersect(module_genes, tf_list)

    if (length(tfs_in_module) > 0) {
      driver_list[[paste0("Module_", i)]] <- tfs_in_module
    } else {
      driver_list[[paste0("Module_", i)]] <- "No known TF driver"
    }
  }

  message("\n--- BMKC Network Module Master Regulators ---")
  for (name in names(driver_list)) {
    drivers <- paste(driver_list[[name]], collapse = ", ")
    message(sprintf("%s: %s", name, drivers))
  }

  bmkc_result$module_drivers <- driver_list
  return(bmkc_result)
}

# ===================================================================
# DATA AGGREGATION & EXPANSION WRAPPERS
# ===================================================================

#' Trajectory-Aware ROS-DET (Sliding Window Network Rewiring)
#' @param expr_matrix Numeric matrix, features (rows) by cells/samples
#'   (columns).
#' @param pseudotime Numeric vector of pseudotime values, one per column
#'   of `expr_matrix`.
#' @param window_size Number of cells per sliding pseudotime window.
#'   Default `50`.
#' @param workers Number of OpenMP threads. Default `1`.
#' @return A list of per-window ROS-DET results describing how pairwise
#'   co-expression changes along the trajectory.
#' @export
run_trajectory_rosdet <- function(expr_matrix, pseudotime, window_size = 50, workers = 1) {
  order_idx <- order(pseudotime)
  sorted_expr <- expr_matrix[, order_idx, drop = FALSE]
  sorted_time <- pseudotime[order_idx]

  n_cells <- ncol(sorted_expr)
  if (n_cells < window_size * 2) {
    stop("Not enough cells to isolate distinct temporal windows.")
  }

  window_1_idx <- 1:window_size
  window_2_idx <- (n_cells - window_size + 1):n_cells

  combined_expr <- cbind(sorted_expr[, window_1_idx, drop = FALSE], sorted_expr[, window_2_idx, drop = FALSE])
  cond_vec <- factor(rep(c("Early_State", "Late_State"), each = window_size))

  message(sprintf("Temporal Window Analysis: Comparing Early (Mean: %.2f) vs Late (Mean: %.2f)",
                  mean(sorted_time[window_1_idx]),
                  mean(sorted_time[window_2_idx])))

  # run_rosdet's generic signature is (expr_matrix, condition, ...) - not
  # `object` (that's run_dcanr's first argument). Passing object= here sent
  # expr_matrix missing into UseMethod and failed before dispatch.
  res <- run_rosdet(
    expr_matrix = combined_expr,
    condition = cond_vec,
    min_delta = 0.4,
    transform = "none",
    workers = workers
  )

  return(res)
}

#' Multi-Omic Differential Co-Expression (Cross-Layer)
#'
#' Combines two omics layers (e.g. RNA + ATAC, host + microbe) into a single
#' matrix, tags each feature's originating layer, and runs any of the three
#' bicorX methods restricted to (by default) cross-layer pairs/edges only.
#' This is the first-class entry point for multi-omics analysis - previously
#' this capability existed only as `run_multiomic_rosdet()`, hard-coded to
#' ROS-DET. That function is kept as a thin wrapper around this one for
#' backward compatibility.
#'
#' @param omic_A,omic_B Two feature x sample matrices sharing the same
#'   samples/columns (e.g. RNA counts and ATAC peak accessibility for the
#'   same cells or the same subjects).
#' @param condition Factor/character vector of length `ncol(omic_A)`.
#' @param method Which method to run: `"rosdet"` (default), `"bmht"`, or
#'   `"bmkc"`.
#' @param layer_A_name,layer_B_name Labels used to tag features from each
#'   layer (default `"Layer_A"`/`"Layer_B"`) and to prefix feature names so
#'   they stay unique after combining (e.g. `"RNA_Gene1"`,
#'   `"ATAC_Peak1"`).
#' @param pair_type_filter For `method = "rosdet"`, whether to keep only
#'   cross-layer pairs (`"inter"`, the default and the point of this
#'   function), only within-layer pairs (`"intra"`), or everything
#'   (`"all"`). Ignored for BMHT/BMKC, which don't have a pair-type filter.
#' @param ... Additional arguments passed to the underlying engine (e.g.
#'   `n_permutations`, `min_delta`, `T1`/`T2`).
#' @details
#' Internally, `omic_A` and `omic_B` are row-bound into one matrix and each
#' row is tagged with its originating layer via the `layer` argument (see
#' `run_bicordcea()`). If `transform = "clr"` or `"rclr"` is requested, the
#' compositional transform is computed *separately within each layer* - see
#' `.compositional_transform()` - so a layer with a very different scale or
#' sparsity (e.g. sparse ATAC peaks vs. dense RNA counts) doesn't distort
#' the other layer's normalization.
#' @export
run_multiomic <- function(omic_A, omic_B, condition,
                           method = c("rosdet", "bmht", "bmkc"),
                           layer_A_name = "Layer_A", layer_B_name = "Layer_B",
                           pair_type_filter = c("inter", "intra", "all"), ...) {
  method <- match.arg(method)
  pair_type_filter <- match.arg(pair_type_filter)

  if (ncol(omic_A) != ncol(omic_B)) {
    stop("Omic matrices must have matching sample dimensions.")
  }
  if (length(condition) != ncol(omic_A)) {
    stop("Condition vector length must match the number of samples (columns).")
  }

  layer_tags <- c(rep(layer_A_name, nrow(omic_A)), rep(layer_B_name, nrow(omic_B)))

  rownames(omic_A) <- paste0(layer_A_name, "_", rownames(omic_A))
  rownames(omic_B) <- paste0(layer_B_name, "_", rownames(omic_B))
  combined_matrix <- rbind(omic_A, omic_B)

  message(sprintf("Multi-Omic Matrix Alignment: %d '%s' features + %d '%s' features.",
                  nrow(omic_A), layer_A_name, nrow(omic_B), layer_B_name))

  extra_args <- list(...)
  if (method == "rosdet") {
    extra_args$pair_type_filter <- pair_type_filter
  } else if (!identical(pair_type_filter, "inter")) {
    warning("`pair_type_filter` only applies to method = \"rosdet\"; ignored for \"",
            method, "\".", call. = FALSE)
  }

  do.call(run_bicordcea, c(
    list(object = combined_matrix, condition = condition, method = method,
         layer = layer_tags),
    extra_args
  ))
}

#' Multi-Omic ROS-DET (Cross-Layer Differential Co-Expression)
#'
#' Kept for backward compatibility. Equivalent to
#' `run_multiomic(omic_A, omic_B, condition, method = "rosdet", ...)` - see
#' `run_multiomic()` for a version that also supports BMHT and BMKC.
#' @param omic_A Numeric matrix for the first omics layer, features
#'   (rows) by samples (columns).
#' @param omic_B Numeric matrix for the second omics layer, with the same
#'   samples (columns) in the same order as `omic_A`.
#' @param condition Two-level factor giving the condition of each sample.
#' @param ... Further arguments passed to [run_multiomic()].
#' @return A `rosdet_result` object restricted to cross-layer pairs.
#' @export
run_multiomic_rosdet <- function(omic_A, omic_B, condition, ...) {
  run_multiomic(omic_A, omic_B, condition, method = "rosdet", ...)
}

# ===================================================================
# EXPORT WRAPPERS
# ===================================================================

#' Export ROS-DET Network for Cytoscape
#' @param rosdet_result A `rosdet_result` object from [run_rosdet()].
#' @param output_dir Directory to write the node and edge tables to.
#'   Created if it does not exist. Default `"Cytoscape_Networks/"`.
#' @return Invisibly, the output directory path. Called for its side
#'   effect of writing Cytoscape-compatible node and edge tables.
#' @export
export_cytoscape <- function(rosdet_result, output_dir = "Cytoscape_Networks/") {
  if (!dir.exists(output_dir)) dir.create(output_dir)

  res <- rosdet_result$results
  if (nrow(res) == 0) stop("No significant pairs to export.")

  interaction_vec <- if ("Direction" %in% colnames(res)) res$Direction else rep("Undirected_Switch", nrow(res))

  edges <- data.frame(
    Source = res$Gene1,
    Target = res$Gene2,
    Interaction = interaction_vec,
    Cor_Condition1 = res$r1,
    Cor_Condition2 = res$r2,
    Switch_Score = res$Distance_Score
  )

  unique_genes <- unique(c(res$Gene1, res$Gene2))
  nodes <- data.frame(Node_ID = unique_genes, Node_Type = "Gene")

  if ("G1_is_TF" %in% colnames(res)) {
    tf_genes <- unique(c(res$Gene1[res$G1_is_TF], res$Gene2[res$G2_is_TF]))
    nodes$Is_Transcription_Factor <- ifelse(nodes$Node_ID %in% tf_genes, "Yes", "No")
  }

  write.csv(edges, file.path(output_dir, "Cytoscape_Edges.csv"), row.names = FALSE)
  write.csv(nodes, file.path(output_dir, "Cytoscape_Nodes.csv"), row.names = FALSE)

  message(sprintf("Cytoscape Export Complete: Saved %d Nodes and %d Edges to '%s'",
                  nrow(nodes), nrow(edges), output_dir))
}
