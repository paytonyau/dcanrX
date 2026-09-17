test_that(".ecor_pvalue() is well-calibrated at N>=50 on both Pearson and bicor data", {
  # Validates the analytical significance method added in Phase A of the
  # paper-fidelity update (Kayano et al. 2011, Eq. 4). Unlike the four
  # generic parametric approaches tried and rejected for the permutation-
  # floor problem elsewhere (GPD-MoM, GPD-MLE, plain normal, skew-normal -
  # see NEWS.md and test-gpd-tail-calibration.R), this is a theoretically
  # derived asymptotic result, not an empirical curve fit, and it shows a
  # materially different (much better) calibration profile.
  # ~3000 reps x 2 scenarios - a couple seconds, runs fine as a normal test.

  set.seed(1)
  n_reps <- 3000
  N <- 50
  true_rho <- 0.3

  # --- True Pearson data (the paper's own exact assumptions) ---
  Sigma <- matrix(c(1, true_rho, true_rho, 1), 2)
  L <- chol(Sigma)
  p_pearson <- numeric(n_reps)
  for (i in seq_len(n_reps)) {
    d1 <- matrix(rnorm(N * 2), N, 2) %*% L
    d2 <- matrix(rnorm(N * 2), N, 2) %*% L
    p_pearson[i] <- bicorX:::.ecor_pvalue(cor(d1[, 1], d1[, 2]), cor(d2[, 1], d2[, 2]), N, N)
  }
  expect_lt(mean(p_pearson < 0.05, na.rm = TRUE) / 0.05, 1.5)
  expect_lt(mean(p_pearson < 0.01, na.rm = TRUE) / 0.01, 1.5)

  # --- bicor on pure noise (the paper's disclosed approximation, and
  #     what ROS-DET actually needs) ---
  p_bicor <- numeric(n_reps)
  for (i in seq_len(n_reps)) {
    m <- matrix(rnorm(2 * 2 * N), nrow = 2)
    d1 <- m[, 1:N]; d2 <- m[, (N + 1):(2 * N)]
    r1 <- bicorX::bicor(d1[1, ], d1[2, ])
    r2 <- bicorX::bicor(d2[1, ], d2[2, ])
    p_bicor[i] <- bicorX:::.ecor_pvalue(r1, r2, N, N)
  }
  # Much stricter bound than the Pearson case above - this is the
  # calibration result that actually matters for ROS-DET, and it should be
  # close to nominal, not just "not terrible". Bounds account for Monte
  # Carlo noise at this rep count (expected count ~30 at p<0.01, so
  # ~18% relative sampling SD is normal).
  expect_lt(mean(p_bicor < 0.05, na.rm = TRUE) / 0.05, 1.3)
  expect_lt(mean(p_bicor < 0.01, na.rm = TRUE) / 0.01, 1.6)
})

test_that(".ecor_pvalue() degrades gracefully (not catastrophically) at small N", {
  # Documents the known, predictable small-sample bias rather than
  # asserting it away - this is the basis for the N<30 warning in
  # .run_rosdet_matrix().
  set.seed(2)
  n_reps <- 2000

  check_inflation <- function(N) {
    p <- numeric(n_reps)
    for (i in seq_len(n_reps)) {
      m <- matrix(rnorm(2 * 2 * N), nrow = 2)
      d1 <- m[, 1:N]; d2 <- m[, (N + 1):(2 * N)]
      p[i] <- bicorX:::.ecor_pvalue(bicorX::bicor(d1[1, ], d1[2, ]),
                                        bicorX::bicor(d2[1, ], d2[2, ]), N, N)
    }
    mean(p < 0.05, na.rm = TRUE) / 0.05
  }

  # At N=10, real inflation is expected (this is *why* the warning exists)
  # - but it should be bounded, not unbounded/wild like the rejected
  # generic-parametric attempts were.
  ratio_n10 <- check_inflation(10)
  expect_gt(ratio_n10, 1.1)   # inflation is real...
  expect_lt(ratio_n10, 3.0)   # ...but bounded, not catastrophic
})

test_that("significance_method = \"analytical\" warns below N=30 and not at/above it", {
  set.seed(3)
  n_genes <- 20
  m_small <- matrix(rnorm(n_genes * 20), nrow = n_genes, ncol = 20)
  rownames(m_small) <- sprintf("Gene_%02d", seq_len(n_genes))
  cond_small <- factor(rep(c("A", "B"), each = 10))

  expect_warning(
    run_rosdet(m_small, cond_small, min_delta = 0.2, significance_method = "analytical"),
    "N >= 50"
  )

  m_big <- matrix(rnorm(n_genes * 100), nrow = n_genes, ncol = 100)
  rownames(m_big) <- sprintf("Gene_%02d", seq_len(n_genes))
  cond_big <- factor(rep(c("A", "B"), each = 50))

  expect_no_warning(
    run_rosdet(m_big, cond_big, min_delta = 0.2, significance_method = "analytical")
  )
})

test_that("significance_method = \"analytical\" produces a valid rosdet_result with no permutation machinery invoked", {
  toy_path <- system.file("extdata", "multiomics_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "multiomics_dataset.RData not found")
  load(toy_path)

  res <- suppressWarnings(suppressMessages(run_rosdet(
    multiomics_rna, multiomics_condition, min_delta = 0.3, significance_method = "analytical"
  )))

  expect_s3_class(res, "rosdet_result")
  expect_equal(res$significance_method, "analytical")
  expect_true(is.na(res$n_permutations))
  expect_true(all(c("P_value", "FDR") %in% names(res$results)))
  expect_true(all(res$results$P_value >= 0 & res$results$P_value <= 1, na.rm = TRUE))

  # The known true signal (genes 1-5) should dominate the top ranks.
  top10 <- head(res$results[order(res$results$P_value), ], 10)
  true_genes <- sprintf("Gene_%02d", 1:5)
  frac_true <- mean(top10$Gene1 %in% true_genes | top10$Gene2 %in% true_genes)
  expect_gte(frac_true, 0.7)
})

test_that("significance_method = \"permutation\" (default) is completely unaffected by adding analytical mode", {
  # Regression guard: adding the analytical branch must not change the
  # default path's behavior at all.
  toy_path <- system.file("extdata", "toy_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "toy_dataset.RData not found")
  load(toy_path)

  res <- suppressMessages(run_rosdet(toy_expr, toy_condition, min_delta = 0.3,
                                      n_permutations = 500, workers = 1, seed = 42))
  expect_equal(res$n_tested, 9518)
  expect_equal(sum(res$results$FDR < 0.05), 0)
  expect_equal(res$significance_method, "permutation")
})
