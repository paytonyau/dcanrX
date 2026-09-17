# Phase 5 benchmark, part 2: outlier robustness.
#
# Question: DGCA (Pearson/Spearman only) and dcanr's z-score method (also
# Pearson-based) have no robust-correlation option. bicorX's core claim
# is that biweight midcorrelation should do better when a real signal is
# present but a handful of samples are corrupted by outliers - a common,
# realistic nuisance (technical artifacts, contamination, mislabeled
# samples) in both bulk and microbiome data. This tests that claim directly
# rather than assuming it.
#
# Design: inject a genuine 15-gene differential-correlation module (as in
# the toy fixture), then corrupt 2 of the 24 samples with extreme outlier
# values specifically in the module genes, and see whether each tool still
# recovers the true signal.

suppressMessages({
  library(bicorX)
  library(DGCA)
  library(dcanr)
})

set.seed(99)
n_genes <- 60
n_samples <- 24
n_module <- 15
true_genes <- sprintf("Gene_%03d", seq_len(n_module))

build_data <- function(n_outlier_samples = 0, outlier_multiplier = 15) {
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.0),
              nrow = n_genes, ncol = n_samples)
  rownames(m) <- sprintf("Gene_%03d", seq_len(n_genes))
  cond <- factor(rep(c("A", "B"), each = n_samples / 2))

  latent <- rnorm(sum(cond == "A"), sd = 2)
  for (g in seq_len(n_module)) {
    m[g, cond == "A"] <- 8 + latent + rnorm(sum(cond == "A"), sd = 0.3)
    m[g, cond == "B"] <- rnorm(sum(cond == "B"), mean = 8, sd = 1.0)
  }

  if (n_outlier_samples > 0) {
    outlier_cols <- sample(which(cond == "A"), n_outlier_samples)
    m[seq_len(n_module), outlier_cols] <- m[seq_len(n_module), outlier_cols] * outlier_multiplier
  }

  list(m = m, cond = cond)
}

evaluate_recovery <- function(sig_genes) {
  if (length(sig_genes) == 0) return(c(precision = NA, recall = 0, n_sig = 0))
  precision <- mean(sig_genes %in% true_genes)
  recall <- mean(true_genes %in% sig_genes)
  c(precision = precision, recall = recall, n_sig = length(sig_genes))
}

run_bicordcea <- function(m, cond) {
  res <- suppressMessages(run_bmht(m, cond, n_permutations = 500, workers = 1))
  evaluate_recovery(res$results$Gene[res$results$FDR < 0.05])
}

run_dgca <- function(m, cond) {
  design <- model.matrix(~0 + cond)
  colnames(design) <- levels(cond)
  res <- suppressMessages(suppressWarnings(
    ddcorAll(inputMat = m, design = design, compare = levels(cond),
             adjust = "perm", nPerms = 50, nPairs = "all", verbose = FALSE)
  ))
  sig <- res[res$pValDiff_adj < 0.05, ]
  evaluate_recovery(unique(c(sig$Gene1, sig$Gene2)))
}

run_dcanr <- function(m, cond) {
  dc_scores <- dcScore(m, cond, dc.method = "zscore")
  dc_z <- suppressWarnings(dcTest(dc_scores, m, cond, dc.method = "zscore"))
  dc_padj <- dcAdjust(dc_z, f = stats::p.adjust, method = "BH")
  dc_padj[upper.tri(dc_padj, diag = TRUE)] <- NA
  sig_mask <- which(dc_padj < 0.05, arr.ind = TRUE)
  sig_genes <- unique(c(rownames(dc_padj)[sig_mask[, 1]], colnames(dc_padj)[sig_mask[, 2]]))
  evaluate_recovery(sig_genes)
}

scenarios <- list(
  "Clean data (0 outlier samples)" = 0,
  "2 of 24 samples corrupted (~8%)" = 2,
  "4 of 24 samples corrupted (~17%)" = 4
)

all_results <- list()
for (label in names(scenarios)) {
  n_out <- scenarios[[label]]
  set.seed(99)
  data <- build_data(n_outlier_samples = n_out)

  row <- data.frame(
    scenario = label,
    tool = c("bicorX", "DGCA", "dcanr"),
    rbind(
      run_bicordcea(data$m, data$cond),
      run_dgca(data$m, data$cond),
      run_dcanr(data$m, data$cond)
    )
  )
  all_results[[label]] <- row
  cat("\n---", label, "---\n")
  print(row, row.names = FALSE)
}

final <- do.call(rbind, all_results)
rownames(final) <- NULL
saveRDS(final, "outlier_robustness_results.rds")

cat("\n\n=== Summary: recall (fraction of true 15-gene module recovered) ===\n")
print(reshape(final[, c("scenario", "tool", "recall")],
              idvar = "scenario", timevar = "tool", direction = "wide"))
