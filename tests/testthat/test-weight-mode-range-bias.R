test_that("weight_mode='wcor' correctly suppresses a pure range-bias artifact (Kayano et al. 2011 Fig 4b/c scenario)", {
  # Paper-fidelity review, Phase C: Kayano et al. validate WCOR
  # specifically against range bias (one condition's expression range
  # much smaller than the other's), not the compositional-closure
  # confound the original weight_mode default change was based on. This
  # locks in the range-bias-specific finding: wcor correctly suppresses
  # a PURE range-bias artifact (no real signal, just a variance
  # difference between conditions) to near-zero, while unweighted does
  # not protect against this at all - matching the paper's design intent.
  set.seed(6)
  n_samples <- 24
  cond <- factor(rep(c("A", "B"), each = 12))
  n_genes <- 20
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.0), nrow = n_genes, ncol = n_samples)
  rownames(m) <- sprintf("Gene_%02d", seq_len(n_genes))

  # Pure range-bias artifact: no true correlation signal, but condition B
  # has a much smaller range than condition A.
  for (g in 1:15) {
    m[g, cond == "A"] <- 8 + rnorm(12, sd = 1.0)
    m[g, cond == "B"] <- 8 + rnorm(12, sd = 0.02)
  }
  artifact_genes <- sprintf("Gene_%02d", 1:15)

  res_wcor <- suppressMessages(run_rosdet(m, cond, min_delta = 0.2, weight_mode = "wcor",
                                           n_permutations = 200, seed = 1))
  res_unw <- suppressMessages(run_rosdet(m, cond, min_delta = 0.2, weight_mode = "unweighted",
                                          n_permutations = 200, seed = 1))

  score_for_group <- function(res, genes) {
    r <- res$results
    mask <- (r$Gene1 %in% genes) & (r$Gene2 %in% genes)
    mean(r$Distance_Score[mask])
  }

  wcor_artifact_score <- score_for_group(res_wcor, artifact_genes)
  unw_artifact_score <- score_for_group(res_unw, artifact_genes)

  # wcor should suppress the pure artifact close to zero; unweighted
  # should not (it has no mechanism to detect range bias at all).
  expect_lt(wcor_artifact_score, 0.05)
  expect_gt(unw_artifact_score, wcor_artifact_score)
})

test_that("weight_mode='wcor' and 'unweighted' perform comparably on clean signal with no variance confound", {
  # The other half of the Phase C characterization: once variance is
  # genuinely equalized between conditions (no confound at all), wcor is
  # NOT broadly worse than unweighted at detecting real signal - the
  # earlier finding that motivated changing the default was specific to
  # signal *coupled with* a variance/range shift, not a general wcor
  # weakness.
  set.seed(6)
  n_samples <- 24
  cond <- factor(rep(c("A", "B"), each = 12))
  n_genes <- 20
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.0), nrow = n_genes, ncol = n_samples)
  rownames(m) <- sprintf("Gene_%02d", seq_len(n_genes))

  latent <- rnorm(12, sd = 1)
  for (g in 1:10) {
    raw_a <- latent + rnorm(12, sd = 0.3)
    raw_b <- rnorm(12)
    # Rescale both to unit variance - removes variance as a confound entirely.
    m[g, cond == "A"] <- 8 + scale(raw_a)[, 1]
    m[g, cond == "B"] <- 8 + scale(raw_b)[, 1]
  }
  clean_signal_genes <- sprintf("Gene_%02d", 1:10)

  # Sanity check the confound was actually removed.
  va <- apply(m[1:10, cond == "A"], 1, var)
  vb <- apply(m[1:10, cond == "B"], 1, var)
  expect_equal(mean(va), 1, tolerance = 1e-6)
  expect_equal(mean(vb), 1, tolerance = 1e-6)

  res_wcor <- suppressMessages(run_rosdet(m, cond, min_delta = 0.2, weight_mode = "wcor",
                                           n_permutations = 200, seed = 1))
  res_unw <- suppressMessages(run_rosdet(m, cond, min_delta = 0.2, weight_mode = "unweighted",
                                          n_permutations = 200, seed = 1))

  score_for_group <- function(res, genes, in_group) {
    r <- res$results
    mask <- (r$Gene1 %in% genes) & (r$Gene2 %in% genes)
    mean(r$Distance_Score[if (in_group) mask else !mask])
  }

  # Both modes should clearly separate clean signal from background -
  # wcor should not be dramatically worse than unweighted here.
  wcor_signal <- score_for_group(res_wcor, clean_signal_genes, TRUE)
  wcor_bg <- score_for_group(res_wcor, clean_signal_genes, FALSE)
  unw_signal <- score_for_group(res_unw, clean_signal_genes, TRUE)
  unw_bg <- score_for_group(res_unw, clean_signal_genes, FALSE)

  expect_gt(wcor_signal, wcor_bg)  # wcor correctly elevates true signal here
  expect_gt(unw_signal, unw_bg)    # unweighted does too
  # Neither mode should be dramatically worse than the other at this task
  # once the variance confound is removed.
  expect_lt(abs(wcor_signal - unw_signal) / max(wcor_signal, unw_signal), 0.5)
})
