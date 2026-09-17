# Phase 5 benchmark, part 2: does bicorX's compositional (CLR) support
# provide a genuine advantage on data where it should matter - i.e., where
# a dominant feature's condition-dependent shift induces spurious
# correlation changes among unrelated features purely through compositional
# closure (the sum-to-one / fixed-total-reads constraint), which is exactly
# the artifact CLR is designed to remove?
#
# Design: 60 taxa x 24 samples (2 conditions x 12).
#   - Taxa 1-10: a TRUE differential module (biologically correlated in
#     condition A only, independent in condition B) - the true positives.
#   - Taxa 11-59: background, independent of everything and of condition.
#   - Taxon 60: a dominant taxon whose absolute abundance is much higher in
#     condition B than condition A. After compositional closure (fixed
#     total reads per sample, as in real sequencing data), this dominant
#     shift mechanically deflates every other taxon's proportional share
#     more in condition B than A - inducing a spurious, purely
#     closure-driven differential-correlation signal among the 49
#     background taxa (11-59) that has nothing to do with true biology.
#     These are the false positives a compositionally-naive method is
#     exposed to and CLR is supposed to protect against.

suppressMessages({
  library(bicorX)
  library(dcanr)
  library(DGCA)
})

set.seed(11)
n_taxa <- 60
n_samples <- 24
cond <- factor(rep(c("A", "B"), each = 12))

# Raw (pre-closure) absolute abundances, log-normal, all independent.
raw <- matrix(rlnorm(n_taxa * n_samples, meanlog = 3, sdlog = 0.4),
              nrow = n_taxa, ncol = n_samples)
rownames(raw) <- sprintf("Taxon_%02d", seq_len(n_taxa))
colnames(raw) <- sprintf("Sample_%02d", seq_len(n_samples))

# Inject the true module: taxa 1-10 share a latent factor in condition A only.
latent_A <- rnorm(sum(cond == "A"), sd = 1)
for (i in 1:10) {
  raw[i, cond == "A"] <- exp(3 + latent_A * 1.2 + rnorm(sum(cond == "A"), sd = 0.15))
  # condition B stays as independent background noise (already set above)
}

# Inject the dominant taxon: much higher absolute abundance in condition B.
raw[60, cond == "A"] <- rlnorm(sum(cond == "A"), meanlog = 3, sdlog = 0.4)
raw[60, cond == "B"] <- rlnorm(sum(cond == "B"), meanlog = 6.5, sdlog = 0.4)  # ~30x higher

# Compositional closure: fixed total reads per sample (as in sequenced count data).
total_reads <- 100000
closed <- sweep(raw, 2, colSums(raw), "/") * total_reads
closed <- round(closed)

true_module <- sprintf("Taxon_%02d", 1:10)
background <- sprintf("Taxon_%02d", 11:59)  # excludes the dominant taxon itself

gene_pairs_true <- t(combn(true_module, 2))                 # 45 true-positive pairs
gene_pairs_bg   <- t(combn(background, 2))                  # background pairs (null)

auprc <- function(scores, labels) {
  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]
  tp <- cumsum(labels); fp <- cumsum(!labels)
  precision <- tp / (tp + fp); recall <- tp / sum(labels)
  keep <- !duplicated(recall, fromLast = TRUE)
  recall <- recall[keep]; precision <- precision[keep]
  sum(diff(c(0, recall)) * precision)
}

run_eval <- function(label, emat, transform_bicor) {
  cat(sprintf("\n--- %s ---\n", label))

  # bicorX ROS-DET, with (or without) the CLR transform under test
  res_bd <- suppressMessages(run_rosdet(emat, cond, min_delta = 0.0,
                                         n_permutations = 100, workers = 1,
                                         weight_mode = "unweighted",
                                         transform = transform_bicor))
  key_bd <- ifelse(res_bd$results$Gene1 < res_bd$results$Gene2,
                    paste(res_bd$results$Gene1, res_bd$results$Gene2),
                    paste(res_bd$results$Gene2, res_bd$results$Gene1))
  score_bd_by_key <- tapply(res_bd$results$Distance_Score, key_bd, max)

  # DGCA and dcanr always operate on the raw (possibly CLR'd) matrix as
  # given - neither has a built-in compositional transform, so whatever
  # `emat` is handed to them is what they see.
  design_mat <- matrix(0, nrow = ncol(emat), ncol = 2,
                       dimnames = list(colnames(emat), c("A", "B")))
  design_mat[cond == "A", "A"] <- 1
  design_mat[cond == "B", "B"] <- 1
  res_dgca <- suppressWarnings(ddcorAll(inputMat = emat, design = design_mat,
                                        compare = c("A", "B"), adjust = "BH",
                                        nPerm = 0, verbose = FALSE))
  key_dgca <- paste(res_dgca$Gene1, res_dgca$Gene2)
  score_dgca_by_key <- setNames(abs(res_dgca$zScoreDiff), key_dgca)

  score_dcanr_mat <- suppressWarnings(abs(as.matrix(dcScore(emat, as.integer(cond), dc.method = "zscore"))))

  score_pairs <- function(pairs, score_by_key_bd, score_by_key_dgca, score_mat_dcanr) {
    k <- paste(pairs[, 1], pairs[, 2])
    list(
      bicordcea = { s <- score_by_key_bd[k]; s[is.na(s)] <- 0; s },
      dgca      = { s <- score_by_key_dgca[k]; s[is.na(s)] <- 0; s },
      dcanr     = score_mat_dcanr[pairs]
    )
  }

  true_scores <- score_pairs(gene_pairs_true, score_bd_by_key, score_dgca_by_key, score_dcanr_mat)
  bg_scores   <- score_pairs(gene_pairs_bg,   score_bd_by_key, score_dgca_by_key, score_dcanr_mat)

  all_labels <- c(rep(1, nrow(gene_pairs_true)), rep(0, nrow(gene_pairs_bg)))
  for (m in c("bicordcea", "dgca", "dcanr")) {
    all_scores <- c(true_scores[[m]], bg_scores[[m]])
    cat(sprintf("  %-10s AUPRC (true module vs background) = %.4f | mean background score = %.4f\n",
                m, auprc(all_scores, all_labels), mean(bg_scores[[m]])))
  }

  invisible(list(true = true_scores, bg = bg_scores))
}

cat(sprintf("Compositional confound dataset: %d taxa, %d samples, %d/%d split\n",
            n_taxa, n_samples, sum(cond == "A"), sum(cond == "B")))
cat(sprintf("True module: taxa 1-10 (%d pairs). Background: taxa 11-59 (%d pairs), exposed to closure artifact from dominant taxon 60.\n",
            nrow(gene_pairs_true), nrow(gene_pairs_bg)))

res_raw <- run_eval("Raw compositional counts (no CLR anywhere - baseline)", closed, "none")
res_clr <- run_eval("bicorX with CLR transform (DGCA/dcanr still see raw counts)", closed, "clr")

saveRDS(list(raw = res_raw, clr = res_clr), "phase5_compositional_results.rds")
message("\nSaved to dev/phase5_compositional_results.rds")
