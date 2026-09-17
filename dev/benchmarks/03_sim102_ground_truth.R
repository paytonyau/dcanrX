# Phase 5 benchmark, corrected (v2): bicorX (ROS-DET) vs DGCA vs dcanr on
# dcanr's own built-in simulated benchmark (sim102), which ships with a
# known ground-truth differential network. This is the strongest evidence
# available in this pass, since it's an external, published benchmark
# rather than data we constructed ourselves.
#
# UPDATED for the paper-fidelity review (Phase A): added ROS-DET with
# significance_method = "analytical" (ECOR, Kayano et al. 2011 Eq. 4).
# This is the direct test of whether implementing the paper's own
# analytical significance test (instead of the permutation test the
# implementation had been using) closes the Phase 5 gap. Result: AUPRC
# 0.2107 (permutation) -> 0.6871 (ECOR), against DGCA/dcanr's 0.7321 -
# closing the great majority of the previously-reported ~3.5x gap. See
# NEWS.md for the full writeup.
#
# CORRECTS TWO BUGS found in earlier drafts of this script:
#
# 1. `getTrueNetwork()` returns an ASYMMETRIC matrix - 20 entries in the
#    upper triangle and 165 DIFFERENT entries in the lower triangle (not
#    mirrored duplicates; verified via isSymmetric() and triangle-by-
#    triangle counts). An earlier draft checked only the upper triangle,
#    silently treating 165 of 185 true edges as non-edges.
#
# 2. `combn(rownames(emat), 2)` does NOT produce alphabetically-ordered
#    pairs - 34% of pairs (3,751 of 11,026) have Gene1 > Gene2
#    alphabetically, since combn just preserves emat's row order. DGCA's
#    output happens to follow that same "natural" row order (so its key
#    construction was accidentally fine), but an earlier draft's ROS-DET
#    key canonicalized pairs alphabetically (`Gene1 < Gene2`) - a
#    systematic mismatch against the non-alphabetical `pair_key`, which
#    silently turned real ROS-DET scores into NA for roughly a third of
#    all pairs regardless of Gate 1 coverage.
#
# Fixed here by using ONE canonical key (always alphabetically sorted)
# applied identically to the ground truth and to all three tools' scores,
# so no tool's matching depends on incidental traversal order.

suppressMessages({
  library(bicorX)
  library(dcanr)
  library(DGCA)
})

canonical_key <- function(a, b) {
  lo <- pmin(a, b); hi <- pmax(a, b)
  paste(lo, hi)
}

data(sim102)
cond_name <- getConditionNames(sim102)[1]
d1 <- getSimData(sim102, cond.name = cond_name, full = FALSE)
emat <- d1$emat
cond <- factor(d1$condition, labels = c("Cond1", "Cond2"))
truenet <- getTrueNetwork(sim102, cond.name = cond_name)

# --- Corrected ground truth: symmetrize by OR across both triangles, key
#     by the same canonical (alphabetical) convention used for every tool ---
true_edges_mat <- (truenet != 0) | (t(truenet) != 0)
diag(true_edges_mat) <- FALSE

gene_pairs <- t(combn(rownames(emat), 2))
true_label_raw <- true_edges_mat[gene_pairs]
truth_key <- canonical_key(gene_pairs[, 1], gene_pairs[, 2])
# Collapse to one row per canonical pair (combn already gives unique
# unordered pairs, so this is a 1:1 relabelling, not a dedup).
true_label <- setNames(true_label_raw, truth_key)

n_genes <- nrow(emat)
cat(sprintf("Dataset: %d genes, %d samples, %d/%d condition split\n",
            n_genes, ncol(emat), sum(cond == "Cond1"), sum(cond == "Cond2")))
cat(sprintf("Ground truth (corrected): %d true edges / %d candidate pairs (%.2f%%)\n\n",
            sum(true_label), length(true_label), 100 * mean(true_label)))

auprc <- function(scores_by_key) {
  common <- names(true_label)
  scores <- scores_by_key[common]
  labels <- true_label[common]
  valid <- !is.na(scores)
  scores <- scores[valid]; labels <- labels[valid]
  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]
  tp <- cumsum(labels); fp <- cumsum(!labels)
  precision <- tp / (tp + fp); recall <- tp / sum(labels)
  keep <- !duplicated(recall, fromLast = TRUE)
  recall <- recall[keep]; precision <- precision[keep]
  sum(diff(c(0, recall)) * precision)
}

results <- list()

# ============================================================
# 1. bicorX ROS-DET
# ============================================================
t_rosdet <- system.time({
  res_rosdet <- suppressMessages(run_rosdet(emat, cond, min_delta = 0.0,
                                             n_permutations = 200, workers = 1))
})
rosdet_key <- canonical_key(res_rosdet$results$Gene1, res_rosdet$results$Gene2)
score_rosdet_by_key <- tapply(res_rosdet$results$Distance_Score, rosdet_key, max)
n_scored <- sum(!is.na(score_rosdet_by_key[names(true_label)]))
cat(sprintf("ROS-DET: %d / %d candidate pairs scored (missing = excluded by Gate 1 even at min_delta=0)\n",
            n_scored, length(true_label)))
cat(sprintf("  Of the %d true edges, %d received a score.\n",
            sum(true_label), sum(true_label & !is.na(score_rosdet_by_key[names(true_label)]))))
results$rosdet <- list(score = score_rosdet_by_key, time = t_rosdet["elapsed"])

# ============================================================
# 1b. bicorX ROS-DET with significance_method = "analytical" (ECOR,
#     Kayano et al. 2011 Eq. 4) - added when implementing the
#     paper-fidelity update. This dataset (406 samples, well above the
#     N=50 calibration threshold documented in .ecor_pvalue()) is exactly
#     the regime ECOR is meant for. Ranked by -log10(p), not
#     Distance_Score - the whole point of ECOR is p-value resolution
#     beyond what permutation counting can achieve, so Distance_Score
#     (unaffected by which significance method is used) is the wrong
#     ranking to demonstrate that with.
# ============================================================
t_ecor <- system.time({
  res_ecor <- suppressWarnings(suppressMessages(
    run_rosdet(emat, cond, min_delta = 0.0, significance_method = "analytical", workers = 1)
  ))
})
key_ecor <- canonical_key(res_ecor$results$Gene1, res_ecor$results$Gene2)
score_ecor_by_key <- setNames(-log10(pmax(res_ecor$results$P_value, 1e-300)), key_ecor)
results$rosdet_ecor <- list(score = score_ecor_by_key, time = t_ecor["elapsed"])

# ============================================================
# 1c. bicorX BMHT - checking whether BMHT shows its own version of the
#     Phase 5 gap (Zheng et al.'s BMHT paper does not specify an
#     analytical alternative to permutation testing the way Kayano et
#     al.'s ROS-DET does, so this isn't fixed by anything in Phase A).
#     BMHT is a per-GENE statistic, so it needs its own gene-level ground
#     truth (a gene counts as "true" if it participates in >=1 true
#     differential edge) - not directly comparable in absolute AUPRC
#     terms to the pair-level numbers above (very different random
#     baseline: ~28% of genes vs ~1.7% of pairs here).
# ============================================================
true_genes <- unique(c(gene_pairs[true_label, 1], gene_pairs[true_label, 2]))

auprc_genes <- function(scores, true_gene_set, all_genes) {
  labels <- all_genes %in% true_gene_set
  ord <- order(scores, decreasing = TRUE); labels <- labels[ord]
  tp <- cumsum(labels); fp <- cumsum(!labels)
  precision <- tp / (tp + fp); recall <- tp / sum(labels)
  keep <- !duplicated(recall, fromLast = TRUE); recall <- recall[keep]; precision <- precision[keep]
  sum(diff(c(0, recall)) * precision)
}

t_bmht <- system.time({
  res_bmht <- suppressMessages(run_bmht(emat, cond, n_permutations = 200, workers = 1))
})
bmht_raw_score <- setNames(res_bmht$results$DC_Score, res_bmht$results$Gene)[rownames(emat)]
bmht_p_score <- setNames(-log10(pmax(res_bmht$results$P_value, 1e-300)), res_bmht$results$Gene)[rownames(emat)]

gene_random_baseline <- mean(rownames(emat) %in% true_genes)
auprc_bmht_raw <- auprc_genes(bmht_raw_score, true_genes, rownames(emat))
auprc_bmht_p <- auprc_genes(bmht_p_score, true_genes, rownames(emat))

cat(sprintf("\nBMHT gene-level check (random baseline = %.4f):\n", gene_random_baseline))
cat(sprintf("  ranked by DC_Score (raw effect size): AUPRC = %.4f\n", auprc_bmht_raw))
cat(sprintf("  ranked by P_value (permutation-based): AUPRC = %.4f\n", auprc_bmht_p))
cat("  NOTE: unlike ROS-DET, p-value ranking is WORSE than raw-score ranking here -\n")
cat("  see NEWS.md for the finding and recommendation (rank by DC_Score, not P_value/FDR,\n")
cat("  for BMHT specifically; use FDR only for a calibrated significance cutoff).\n")
results$bmht_gene_level <- list(auprc_raw = auprc_bmht_raw, auprc_p = auprc_bmht_p, time = t_bmht["elapsed"])

# ============================================================
# 2. DGCA
# ============================================================
design <- matrix(0, nrow = ncol(emat), ncol = 2,
                  dimnames = list(colnames(emat), c("Cond1", "Cond2")))
design[cond == "Cond1", "Cond1"] <- 1
design[cond == "Cond2", "Cond2"] <- 1

t_dgca <- system.time({
  res_dgca <- suppressMessages(suppressWarnings(
    ddcorAll(inputMat = emat, design = design, compare = c("Cond1", "Cond2"),
             adjust = "none", nPerms = 0, verbose = FALSE)
  ))
})
key_dgca <- canonical_key(res_dgca$Gene1, res_dgca$Gene2)
score_dgca_by_key <- setNames(abs(res_dgca$zScoreDiff), key_dgca)
results$dgca <- list(score = score_dgca_by_key, time = t_dgca["elapsed"])

# Gene-level rollup for DGCA (max |zScoreDiff| over all pairs involving
# that gene) - the natural analog to compare against BMHT's own
# gene-level numbers above, since DGCA itself has no native per-gene
# statistic.
gene_max_z <- sapply(rownames(emat), function(g) {
  rows <- res_dgca$Gene1 == g | res_dgca$Gene2 == g
  if (any(rows)) max(abs(res_dgca$zScoreDiff[rows]), na.rm = TRUE) else 0
})
results$dgca_gene_level <- list(auprc = auprc_genes(gene_max_z, true_genes, rownames(emat)))

# ============================================================
# 3. dcanr (z-score method)
# ============================================================
t_dcanr <- system.time({
  score_mat_dcanr <- suppressWarnings(abs(dcScore(emat, cond, dc.method = "zscore")))
})
score_dcanr_by_key <- setNames(score_mat_dcanr[gene_pairs], truth_key)
results$dcanr <- list(score = score_dcanr_by_key, time = t_dcanr["elapsed"])

cat("\n=== Pair-level AUPRC (ROS-DET, DGCA, dcanr) on corrected ground truth ===\n")
cat("(higher = better; random baseline =", round(mean(true_label), 6), ")\n\n")
for (m in setdiff(names(results), c("bmht_gene_level", "dgca_gene_level"))) {
  cat(sprintf("%-10s AUPRC = %.4f | elapsed = %.2fs\n",
              m, auprc(results[[m]]$score), results[[m]]$time))
}
cat(sprintf("\n%-10s AUPRC(DC_Score) = %.4f | AUPRC(P_value) = %.4f | elapsed = %.2fs  [gene-level, not directly comparable to pair-level numbers above]\n",
            "bmht", results$bmht_gene_level$auprc_raw, results$bmht_gene_level$auprc_p,
            results$bmht_gene_level$time))
cat(sprintf("%-10s AUPRC(max |z| per gene) = %.4f  [gene-level, DGCA has no native per-gene statistic]\n",
            "dgca", results$dgca_gene_level$auprc))

saveRDS(list(true_label = true_label, results = results),
        "sim102_corrected_results.rds")
message("\nSaved to dev/benchmarks/sim102_corrected_results.rds")
