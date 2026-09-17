# Phase 5 benchmark, part 1: null calibration comparison.
#
# Question: on data with NO true differential-correlation signal, does
# each tool's FDR-significant call rate stay near the nominal alpha, or
# does it inflate? This directly extends the Phase 1 BMHT calibration
# methodology to a head-to-head comparison against DGCA and dcanr.
#
# Requires DGCA and dcanr installed from source (not on the apt mirror):
#   curl -sL https://github.com/andymckenzie/DGCA/archive/refs/heads/master.tar.gz | tar xz
#   R CMD INSTALL DGCA-master
#   curl -sL https://github.com/DavisLaboratory/dcanr/archive/refs/heads/master.tar.gz | tar xz
#   R CMD INSTALL dcanr-master

suppressMessages({
  library(bicorX)
  library(DGCA)
  library(dcanr)
})

set.seed(2024)
n_genes <- 80
n_samples <- 24
n_reps <- 10  # independent pure-noise datasets

run_one_rep <- function(seed) {
  set.seed(seed)
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
              nrow = n_genes, ncol = n_samples)
  rownames(m) <- sprintf("Gene_%03d", seq_len(n_genes))
  cond <- factor(rep(c("A", "B"), each = n_samples / 2))

  # --- bicorX (BMHT) ---
  bmht <- suppressMessages(run_bmht(m, cond, n_permutations = 200, workers = 1))
  bicor_fpr <- mean(bmht$results$FDR < 0.05)

  # --- DGCA ---
  design <- model.matrix(~0 + cond)
  colnames(design) <- levels(cond)
  dgca_res <- suppressMessages(suppressWarnings(
    ddcorAll(inputMat = m, design = design, compare = levels(cond),
             adjust = "perm", nPerms = 50, nPairs = "all", verbose = FALSE)
  ))
  # Per-gene FPR: a gene counts as a "hit" if any pair involving it is
  # FDR-significant - the closest DGCA analogue to BMHT's per-gene call.
  dgca_sig_genes <- unique(c(dgca_res$Gene1[dgca_res$pValDiff_adj < 0.05],
                             dgca_res$Gene2[dgca_res$pValDiff_adj < 0.05]))
  dgca_fpr <- length(dgca_sig_genes) / n_genes

  # --- dcanr (z-score method) ---
  dc_scores <- dcScore(m, cond, dc.method = "zscore")
  dc_z <- suppressWarnings(dcTest(dc_scores, m, cond, dc.method = "zscore"))
  dc_padj <- dcAdjust(dc_z, f = stats::p.adjust, method = "BH")
  dc_padj[upper.tri(dc_padj, diag = TRUE)] <- NA
  sig_mask <- dc_padj < 0.05
  dcanr_sig_genes <- unique(c(
    rownames(dc_padj)[which(sig_mask, arr.ind = TRUE)[, 1]],
    colnames(dc_padj)[which(sig_mask, arr.ind = TRUE)[, 2]]
  ))
  dcanr_fpr <- length(dcanr_sig_genes) / n_genes

  data.frame(rep = seed, bicorX = bicor_fpr, DGCA = dgca_fpr, dcanr = dcanr_fpr)
}

results <- do.call(rbind, lapply(seq_len(n_reps), run_one_rep))

cat("\n=== Per-gene false-positive rate at FDR < 0.05 on PURE NOISE (no true signal) ===\n")
cat("(nominal target: ~0.05 or below; higher = anti-conservative/miscalibrated)\n\n")
print(round(colMeans(results[, -1]), 4))
cat("\nPer-replicate detail:\n")
print(results)

saveRDS(results, "calibration_comparison_results.rds")
