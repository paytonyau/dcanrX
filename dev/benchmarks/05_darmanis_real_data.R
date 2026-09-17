# Phase 5 benchmark, part 5: bicorX vs DGCA on a real dataset.
#
# Unlike 01-04 (synthetic data) and 03 (a published *simulation* with
# known ground truth), this uses genuine published biological data with
# real technical noise (dropout, etc.) - a complementary check, since it
# has no ground-truth network to compute AUPRC against, but does let us
# compare runtime and result-set overlap on data neither tool was tuned
# for. Pair-key construction here already used a canonical sorted key
# (`paste(sort(x), collapse="|")`) from the start, so this script wasn't
# affected by the key-ordering bug found and fixed in 03.
#
# Data: Darmanis et al. 2015 human brain single-cell RNA-seq (572 genes x
# 158 cells, neuron vs oligodendrocyte), bundled with the DGCA package
# itself (data(darmanis), data(design_mat)) - a genuine published dataset
# with real dropout/technical noise, not synthetic.
#
# DGCA only supports Pearson/Spearman correlation (confirmed by reading
# its source: dCor.R takes corrType %in% c("pearson","spearman") with no
# bicor option) - so this is a real methodological contrast, not just two
# implementations of the same estimator.

suppressMessages({
  library(DGCA)
  library(bicorX)
})

data(darmanis)
data(design_mat)
condition <- factor(ifelse(design_mat[, "neuron"] == 1, "neuron", "oligodendrocyte"))
expr <- as.matrix(darmanis)

results <- list(
  meta = list(
    n_genes = nrow(expr), n_samples = ncol(expr),
    condition_table = table(condition),
    zero_fraction = mean(expr == 0),
    r_version = R.version.string,
    dgca_version = as.character(packageVersion("DGCA")),
    timestamp = Sys.time()
  )
)

n_perms <- 50

# --- DGCA (Pearson, permutation-based FDR) ---
cat("=== Running DGCA ===\n")
results$dgca_time <- system.time({
  results$dgca <- ddcorAll(
    inputMat = expr, design = design_mat,
    compare = c("neuron", "oligodendrocyte"),
    corrType = "pearson", adjust = "perm", nPerms = n_perms,
    verbose = FALSE
  )
})

# --- bicorX ROS-DET (bicor, permutation-based FDR) ---
cat("\n=== Running ROS-DET ===\n")
results$rosdet_time <- system.time({
  results$rosdet <- run_rosdet(expr, condition, min_delta = 0.2,
                                n_permutations = n_perms, workers = 1)
})

# --- bicorX BMHT, for reference (different question - per-gene hub score) ---
cat("\n=== Running BMHT ===\n")
results$bmht_time <- system.time({
  results$bmht <- run_bmht(expr, condition, n_permutations = n_perms, workers = 1)
})

# ===================================================================
# Comparison: merge DGCA and ROS-DET results on the same gene pairs
# ===================================================================
dgca_df <- results$dgca
dgca_df$pair_key <- apply(cbind(dgca_df$Gene1, dgca_df$Gene2), 1,
                           function(x) paste(sort(x), collapse = "|"))

rosdet_df <- results$rosdet$results
rosdet_df$pair_key <- apply(cbind(rosdet_df$Gene1, rosdet_df$Gene2), 1,
                             function(x) paste(sort(x), collapse = "|"))

merged <- merge(dgca_df, rosdet_df, by = "pair_key", suffixes = c("_dgca", "_rosdet"))
cat("\nPairs in common (ROS-DET tested a subset via Gate 1 pre-filtering):", nrow(merged), "\n")

results$merged_n_rows <- nrow(merged)  # kept as a count only - full table was 18MB, not needed
results$rank_correlation <- cor(abs(merged$zScoreDiff), merged$Distance_Score, method = "spearman")
cat("Spearman correlation between |DGCA zScoreDiff| and ROS-DET Distance_Score (shared pairs):",
    round(results$rank_correlation, 3), "\n")

# Overlap of "significant" pairs at a matched nominal threshold (both at
# raw p < 0.05, since exact FDR procedures differ and n_perms=50 limits
# how small any empirical p-value can get on either side).
dgca_sig <- dgca_df$pair_key[dgca_df$pValDiff < 0.05]
rosdet_sig <- rosdet_df$pair_key[rosdet_df$P_value < 0.05]
overlap <- length(intersect(dgca_sig, rosdet_sig))
results$overlap_counts <- list(dgca_sig = length(dgca_sig), rosdet_sig = length(rosdet_sig),
                                overlap = overlap,
                                jaccard = overlap / length(union(dgca_sig, rosdet_sig)))
cat("\nDGCA p<0.05 pairs:", length(dgca_sig), " | ROS-DET p<0.05 pairs:", length(rosdet_sig),
    " | overlap:", overlap, " | Jaccard:", round(results$overlap_counts$jaccard, 3), "\n")

results$bmht_n_significant <- sum(results$bmht$results$FDR < 0.05)
results$bmht_n_genes <- nrow(results$bmht$results)

# Drop the large raw tables before saving - only summary statistics are
# needed for reporting (the full merged table was ~18MB and isn't used
# anywhere downstream).
results$merged <- NULL
results$dgca <- NULL
results$rosdet <- NULL
results$bmht <- NULL

saveRDS(results, "darmanis_summary_results.rds")
cat("\nSaved compact summary to darmanis_summary_results.rds\n")
