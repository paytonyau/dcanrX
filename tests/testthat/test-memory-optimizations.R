test_that("BMHT's incremental count accumulation gives identical results to the old full-matrix approach", {
  # Regression test for a memory refactor: the R engine previously
  # retained a full n_genes x n_permutations matrix of every permutation's
  # DC score for the whole run, only ever used for a single rowSums
  # reduction at the very end. At the very large permutation counts often
  # needed for adequate power (Phase 5 finding), that was hundreds of MB
  # to GB of transient R-side memory. Replaced with incremental per-chunk
  # accumulation of just the two counts actually needed
  # (greater_eq_counts, valid_perms_counts).
  #
  # This checks against the exact reference values established across
  # many earlier validations in this project's history (Phase 1's
  # calibration fix, Phase 3/4/5's repeated re-checks) - if this refactor
  # changed anything numerically, this is where it would show up.
  toy_path <- system.file("extdata", "toy_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "toy_dataset.RData not found")
  load(toy_path)

  res <- suppressMessages(run_bmht(toy_expr, toy_condition, n_permutations = 500,
                                    workers = 1, seed = 42))

  expect_equal(sum(res$results$FDR < 0.05), 20)

  gene7 <- res$results[res$results$Gene == "Gene_007", ]
  expect_equal(gene7$DC_Score, 2.549791, tolerance = 1e-5)
  expect_equal(gene7$P_value, 0.003992016, tolerance = 1e-8)
  expect_equal(gene7$FDR, 0.03992016, tolerance = 1e-6)
})

test_that("BMHT peak memory no longer scales with n_permutations", {
  # Directly verifies the refactor's actual purpose: memory usage should
  # now be dominated by fixed-size objects (chunk_size-bounded), not by
  # n_permutations. Compares peak memory at two very different permutation
  # counts and confirms it does NOT grow proportionally.
  set.seed(1)
  n_genes <- 100
  n_samples <- 24
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5), nrow = n_genes, ncol = n_samples)
  rownames(m) <- sprintf("Gene_%03d", seq_len(n_genes))
  cond <- factor(rep(c("A", "B"), each = 12))

  measure_peak_alloc <- function(n_perm) {
    gc(reset = TRUE, full = TRUE)
    invisible(suppressMessages(run_bmht(m, cond, n_permutations = n_perm, workers = 1)))
    g <- gc(full = TRUE)
    sum(g[, 6])  # "max used" column, in Mb, across Ncells/Vcells rows
  }

  peak_small <- measure_peak_alloc(200)
  peak_large <- measure_peak_alloc(20000)

  # If memory scaled with n_permutations (the old behavior: a
  # n_genes x n_permutations matrix retained for the whole run), a 100x
  # increase in permutations would show a roughly proportional increase in
  # peak memory. With the fix, peak memory should be roughly flat -
  # allow a generous margin (5x) since gc() reporting is coarse and other
  # allocations (chunk buffers, progress bar, etc.) still exist.
  ratio <- peak_large / peak_small
  expect_lt(ratio, 5)
})

test_that("BMKC's single-extraction quantile refactor gives identical results to the old double-extraction approach", {
  # Regression test for a memory/compute refactor: get_quantile() (now
  # split into extract_and_subsample() + quantile_from_vals()) previously
  # re-extracted the full upper-triangle vector from scratch for both the
  # T1 and T2 threshold calls on the same matrix - duplicating an O(n^2)
  # extraction (and, at n_genes=20,000, a ~1.6GB intermediate vector) for
  # no reason, since the extraction doesn't depend on which threshold is
  # being computed. This checks against the exact module recovered in
  # Phase 3's BMKC validation.
  set.seed(42)
  n_genes <- 60
  n_samples <- 24
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.0), nrow = n_genes, ncol = n_samples)
  rownames(m) <- sprintf("Gene_%02d", seq_len(n_genes))
  cond <- factor(rep(c("A", "B"), each = 12))
  latent <- rnorm(12, sd = 3)
  for (g in 1:15) {
    m[g, cond == "A"] <- 8 + latent + rnorm(12, sd = 0.2)
    m[g, cond == "B"] <- 8 + rnorm(12, sd = 1.0)
  }

  res <- suppressMessages(run_bmkc(m, cond, workers = 1, T2 = 0.3, gamma = 0.5))

  expect_equal(res$n_modules, 1)
  expect_identical(sort(res$modules[[1]]), sprintf("Gene_%02d", 1:15))
})

test_that("BMKC percentile threshold mode still runs correctly after the quantile refactor", {
  # The percentile code path is the only one that exercises
  # extract_and_subsample()/quantile_from_vals() at all (absolute
  # threshold mode never calls it) - a smoke test to make sure the
  # refactor didn't break that path structurally.
  set.seed(3)
  n_genes <- 50
  n_samples <- 24
  m <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(m) <- sprintf("Gene_%03d", seq_len(n_genes))
  cond <- factor(rep(c("A", "B"), each = 12))

  expect_no_error(
    res <- run_bmkc(m, cond, workers = 1, threshold_method = "percentile",
                     T1 = 80, T2 = 20, seed = 7)
  )
  expect_true(is.list(res$parameters))
})

test_that("ROS-DET's bounded chunk_size (same fix as BMHT) gives identical results to before", {
  # Same fixed-chunk-size fix as BMHT's engine, applied here for
  # consistency: chunk_size no longer grows proportionally with
  # n_permutations (previously max(10, ceiling(n_permutations/10)), now
  # min(500, max(10, n_permutations))). The impact is smaller here than
  # for BMHT (this only bounds the shuffled-labels chunk, not a
  # n_genes-scale output matrix - ROS-DET's C++ kernel already reduces to
  # counts internally and never returns a full null-value matrix to R),
  # but it keeps progress-bar granularity from collapsing at very large
  # permutation counts and closes the same class of issue for
  # consistency. Checked against the exact reference values established
  # across this project's history.
  toy_path <- system.file("extdata", "toy_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "toy_dataset.RData not found")
  load(toy_path)

  res <- suppressMessages(run_rosdet(toy_expr, toy_condition, min_delta = 0.3,
                                      n_permutations = 500, workers = 1, seed = 42))

  expect_equal(res$n_tested, 9518)
  expect_equal(sum(res$results$FDR < 0.05), 0)

  top <- res$results[order(res$results$P_value), ][1, ]
  expect_equal(top$Gene1, "Gene_018")
  expect_equal(top$Gene2, "Gene_020")
  expect_equal(top$Distance_Score, 1.690249, tolerance = 1e-5)
})
