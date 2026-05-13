# ===================================================================
# CAUSAL INFERENCE WRAPPERS (TF-GUIDED DIRECTIONALITY)
# ===================================================================

#' Annotate Causal Drivers (TF-Guided BMHT)
#' @description Identifies putative causal drivers by cross-referencing highly rewired BMHT hubs with known Transcription Factors.
#' @param bmht_result The output object from `run_bmht()`.
#' @param tf_list A character vector of known Transcription Factor gene symbols.
#' @param top_pct Numeric. The top percentage of hubs (by score) to consider as master regulators (default: 0.05 for top 5%).
#' @return An updated `bmht_result` object with a new `$causal_drivers` dataframe.
#' @export
identify_causal_drivers <- function(bmht_result, tf_list, top_pct = 0.05) {
  res <- bmht_result$results

  # FIXED: Swapped 'Score' to 'DC_Score' to match your core engine output
  # Sort by Rewiring Score (descending) to find the biggest hubs
  res <- res[order(res$DC_Score, decreasing = TRUE), ]

  # Identify the threshold for the top X% of hubs
  threshold_idx <- max(1, round(nrow(res) * top_pct))
  top_hubs <- res$Gene[1:threshold_idx]

  # Cross-reference with the TF database
  res$Is_TF <- res$Gene %in% tf_list
  res$Causal_Driver <- ifelse(res$Gene %in% top_hubs & res$Is_TF, "Yes", "No")

  # Extract just the causal drivers for a clean output table
  causal_df <- res[res$Causal_Driver == "Yes", ]

  message(sprintf("BMHT: Identified %d putative causal drivers out of the top %d highly rewired hubs.",
                  nrow(causal_df), threshold_idx))

  # Return the updated object
  bmht_result$results <- res
  bmht_result$causal_drivers <- causal_df

  return(bmht_result)
}

#' Annotate Direction for ROS-DET Switching Pairs
#' @description Assigns directed edges to ROS-DET switching pairs based on prior Transcription Factor knowledge.
#' @param rosdet_result The output object from `run_rosdet()`.
#' @param tf_list A character vector of known Transcription Factor gene symbols.
#' @return An updated `rosdet_result` object with a new `$Direction` column.
#' @export
annotate_rosdet_direction <- function(rosdet_result, tf_list) {
  res <- rosdet_result$results

  if (nrow(res) == 0) {
    message("No significant ROS-DET pairs to annotate.")
    return(rosdet_result)
  }

  # Check TF status for both genes in the pair
  res$G1_is_TF <- res$Gene1 %in% tf_list
  res$G2_is_TF <- res$Gene2 %in% tf_list

  # Assign Directed Arrows natively in R without needing dplyr
  res$Direction <- "Undirected"
  res$Direction[res$G1_is_TF & !res$G2_is_TF] <- paste(res$Gene1[res$G1_is_TF & !res$G2_is_TF], "->", res$Gene2[res$G1_is_TF & !res$G2_is_TF])
  res$Direction[!res$G1_is_TF & res$G2_is_TF] <- paste(res$Gene2[!res$G1_is_TF & res$G2_is_TF], "->", res$Gene1[!res$G1_is_TF & res$G2_is_TF])
  res$Direction[res$G1_is_TF & res$G2_is_TF]  <- "Co-Regulatory (TF <-> TF)"

  # Count the causal switches
  directed_count <- sum(res$Direction != "Undirected")
  message(sprintf("ROS-DET: Assigned regulatory direction to %d switching pairs.", directed_count))

  rosdet_result$results <- res
  return(rosdet_result)
}

#' Annotate Master Regulators for BMKC Modules
#' @description Scans BMKC co-expression modules to identify the embedded Transcription Factors acting as clique master regulators.
#' @param bmkc_result The output object from `run_bmkc()`.
#' @param tf_list A character vector of known Transcription Factor gene symbols.
#' @return An updated `bmkc_result` object with a new `$module_drivers` list.
#' @export
annotate_bmkc_drivers <- function(bmkc_result, tf_list) {
  if (bmkc_result$n_modules == 0) stop("No modules found to annotate.")

  driver_list <- list()

  for (i in seq_len(bmkc_result$n_modules)) {
    module_genes <- bmkc_result$modules[[i]]

    # Find which genes in this specific module are TFs
    tfs_in_module <- intersect(module_genes, tf_list)

    if (length(tfs_in_module) > 0) {
      driver_list[[paste0("Module_", i)]] <- tfs_in_module
    } else {
      driver_list[[paste0("Module_", i)]] <- "No known TF driver"
    }
  }

  # Print a clean summary
  message("\n--- BMKC Module Master Regulators ---")
  for (name in names(driver_list)) {
    drivers <- paste(driver_list[[name]], collapse = ", ")
    message(sprintf("%s: %s", name, drivers))
  }

  # Save to the object
  bmkc_result$module_drivers <- driver_list
  return(bmkc_result)
}


# ===================================================================
# DATA AGGREGATION & EXPANSION WRAPPERS (V1.0 HORIZON)
# ===================================================================

#' Trajectory-Aware ROS-DET (Sliding Window Network Rewiring)
#' @description Identifies network correlation switches across a continuous pseudotime trajectory.
#' @param expr_matrix A dense expression matrix (rows = genes, columns = cells/samples).
#' @param pseudotime A numeric vector representing the timeline or developmental stage of each cell.
#' @param window_size The number of cells to include in a single time window.
#' @param workers Number of cores for parallel processing.
#' @export
run_trajectory_rosdet <- function(expr_matrix, pseudotime, window_size = 50, workers = 1) {

  # 1. Sort the data strictly by time
  order_idx <- order(pseudotime)
  sorted_expr <- expr_matrix[, order_idx]
  sorted_time <- pseudotime[order_idx]

  n_cells <- ncol(sorted_expr)
  if (n_cells < window_size * 2) {
    stop("Not enough cells to create distinct temporal windows.")
  }

  # 2. Define "Early" vs "Late" states based on sliding windows
  window_1_idx <- 1:window_size
  window_2_idx <- (n_cells - window_size + 1):n_cells

  # 3. Create a unified matrix and condition vector for the C-engine
  combined_expr <- cbind(sorted_expr[, window_1_idx], sorted_expr[, window_2_idx])
  cond_vec <- factor(rep(c("Early_State", "Late_State"), each = window_size))

  message(sprintf("Comparing Early (Mean Time: %.2f) vs Late (Mean Time: %.2f)",
                  mean(sorted_time[window_1_idx]),
                  mean(sorted_time[window_2_idx])))

  # 4. Feed the temporal data directly into your existing ROS-DET engine
  res <- run_rosdet(
    expr_matrix = combined_expr,
    condition = cond_vec,
    min_delta = 0.4,
    transform = "none",
    workers = workers
  )

  return(res)
}

#' Multi-Omic ROS-DET (Cross-Layer Differential Co-Expression)
#' @description Identifies correlation switches specifically between two different omic layers (e.g., RNA vs. Protein).
#' @param omic_A A dense matrix for Layer A. Rows = features, Columns = samples.
#' @param omic_B A dense matrix for Layer B. Rows = features, Columns = samples.
#' @param condition A factor vector of length equal to the number of columns, denoting the biological state.
#' @param ... Additional arguments passed to the core `run_rosdet` engine.
#' @export
run_multiomic_rosdet <- function(omic_A, omic_B, condition, ...) {

  if (ncol(omic_A) != ncol(omic_B)) {
    stop("Omic matrices must have the exact same number of samples (columns).")
  }
  if (length(condition) != ncol(omic_A)) {
    stop("Condition vector length must match the number of samples.")
  }

  # Rename rows to track layer origin
  rownames(omic_A) <- paste0("OmicA_", rownames(omic_A))
  rownames(omic_B) <- paste0("OmicB_", rownames(omic_B))

  # Stitch into master matrix
  combined_matrix <- rbind(omic_A, omic_B)

  message(sprintf("Running Multi-Omic ROS-DET: Comparing %d Layer-A features vs %d Layer-B features.",
                  nrow(omic_A), nrow(omic_B)))

  # Pass to core engine utilizing omic_mode (assumes underlying C-code handles cross-layer logic)
  res <- run_rosdet(
    expr_matrix = combined_matrix,
    condition = condition,
    omic_mode = TRUE,
    ...
  )

  return(res)
}


# ===================================================================
# EXPORT WRAPPERS
# ===================================================================

#' Export ROS-DET Network for Cytoscape
#' @description Generates Cytoscape-ready Node and Edge tables from a ROS-DET result object.
#' @param rosdet_result The output object from `run_rosdet()`, preferably annotated with directions.
#' @param output_dir String. The directory to save the CSV files.
#' @export
export_cytoscape <- function(rosdet_result, output_dir = "Cytoscape_Networks/") {
  if (!dir.exists(output_dir)) dir.create(output_dir)

  res <- rosdet_result$results

  if (nrow(res) == 0) stop("No significant pairs to export.")

  # Handle the interaction column safely without truncating length
  interaction_vec <- if ("Direction" %in% colnames(res)) res$Direction else rep("Undirected_Switch", nrow(res))

  # 1. Build the EDGE Table (matching core engine column names: r1, r2, Score)
  edges <- data.frame(
    Source = res$Gene1,
    Target = res$Gene2,
    Interaction = interaction_vec,
    Cor_Condition1 = res$r1,
    Cor_Condition2 = res$r2,
    Switch_Score = res$Score
  )

  # 2. Build the NODE Table
  unique_genes <- unique(c(res$Gene1, res$Gene2))

  nodes <- data.frame(
    Node_ID = unique_genes,
    Node_Type = "Gene"
  )

  # If Causal Wrapper was run, flag TFs
  if ("G1_is_TF" %in% colnames(res)) {
    tf_genes <- unique(c(res$Gene1[res$G1_is_TF], res$Gene2[res$G2_is_TF]))
    nodes$Is_Transcription_Factor <- ifelse(nodes$Node_ID %in% tf_genes, "Yes", "No")
  }

  # 3. Save to disk
  write.csv(edges, file.path(output_dir, "Cytoscape_Edges.csv"), row.names = FALSE)
  write.csv(nodes, file.path(output_dir, "Cytoscape_Nodes.csv"), row.names = FALSE)

  message(sprintf("Successfully exported %d Nodes and %d Edges to '%s'",
                  nrow(nodes), nrow(edges), output_dir))
}
