test_that("BMHT is not anti-conservative on pure-noise data (no false positives at FDR<0.05)", {
  # Regression test for the fixed-mask bug: an earlier version derived
  # informative_mask from the observed data and reused it, unchanged, for
  # every permutation. That leaked the true group structure into the null
  # and produced ~80% "significant" genes (FDR < 0.05) on pure noise with no
  # true differential signal at all. The mask must be recomputed inside each
  # permutation from that permutation's own correlation matrices.
  set.seed(7)
  n_genes <- 60
  n_samples <- 24
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
              nrow = n_genes, ncol = n_samples)
  cond <- factor(rep(c("A", "B"), each = n_samples / 2))
  rownames(m) <- sprintf("Gene_%03d", seq_len(n_genes))

  res <- suppressMessages(run_bmht(m, cond, n_permutations = 200, workers = 1))

  # Zero true signal -> essentially no genes should survive FDR < 0.05.
  # Allow a very small tolerance for permutation-floor noise rather than
  # requiring exactly zero, since this is a stochastic test.
  n_sig <- sum(res$results$FDR < 0.05)
  expect_lte(n_sig, ceiling(0.10 * n_genes))  # was 48/60 (80%) before the fix
})

test_that("BMHT still detects a real injected differential-correlation module", {
  # Companion test to the calibration check above: confirms the fix didn't
  # destroy power. Uses the Phase 0 toy fixture (20 genes with real signal
  # out of 200).
  toy_path <- system.file("extdata", "toy_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "toy_dataset.RData not found - run inst/extdata/make_toy_data.R")
  load(toy_path)

  res <- suppressMessages(run_bmht(toy_expr, toy_condition, n_permutations = 500, workers = 1))
  true_genes <- sprintf("Gene_%03d", 1:20)

  sig_genes <- res$results$Gene[res$results$FDR < 0.05]
  expect_gt(length(sig_genes), 0)

  # Nearly all significant hits should be true signal genes (allow a little
  # slack for the stochastic permutation test).
  precision <- mean(sig_genes %in% true_genes)
  expect_gte(precision, 0.9)

  # And most of the true signal should be recovered.
  recall <- mean(true_genes %in% sig_genes)
  expect_gte(recall, 0.5)
})
