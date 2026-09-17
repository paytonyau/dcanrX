test_that("cpp_bmht_permutations (workspace-based) matches an independent gold-standard computation", {
  # Regression test for the buffer-reuse performance refactor: the hot
  # permutation loop was restructured to reuse per-thread workspace buffers
  # (BicorWorkspace/fast_bicor_ws in bmht_cpp.cpp) instead of allocating
  # fresh matrices every permutation, since A/B testing showed this saves
  # ~35-40% of runtime at 600+ genes (negligible at the 50-300 gene scale
  # used elsewhere in this test suite, which is why this needs its own
  # dedicated check rather than relying on existing tests to catch a subtle
  # buffer-reuse bug, e.g. stale data leaking between iterations).
  #
  # This independently reimplements what the C++ kernel should compute,
  # using only the untouched, separately-validated cpp_bmht_observed_bicor
  # path (which does NOT share any workspace/buffer-reuse code with the
  # permutation kernel), and checks they agree to machine precision.
  set.seed(5)
  n_genes <- 40
  n_samples <- 20
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
              nrow = n_genes, ncol = n_samples)
  cond_int <- as.integer(rep(1:2, each = 10))

  n_perm <- 15
  set.seed(99)
  shuffled <- replicate(n_perm, sample(cond_int))
  half_threshold <- 0.4
  mad_1 <- 1.0
  mad_2 <- 1.0

  out_new <- bicorX:::cpp_bmht_permutations(m, shuffled, half_threshold, mad_1, mad_2, 1L)

  gold <- matrix(0, n_genes, n_perm)
  for (p in seq_len(n_perm)) {
    labels <- shuffled[, p]
    idx1 <- which(labels == 1)
    idx2 <- which(labels == 2)
    r1 <- bicorX:::cpp_bmht_observed_bicor(m[, idx1])
    r2 <- bicorX:::cpp_bmht_observed_bicor(m[, idx2])
    mask <- (abs(r1) > half_threshold) | (abs(r2) > half_threshold)
    r1s <- r1 / mad_1
    r2s <- r2 / mad_2
    diff_sq <- (r1s - r2s)^2
    for (g in seq_len(n_genes)) {
      active <- mask[g, ]
      gold[g, p] <- if (sum(active) > 0) sqrt(mean(diff_sq[g, active])) else 0
    }
  }

  expect_equal(out_new, gold, tolerance = 1e-10)
})

test_that("cpp_bmht_permutations gives identical results across repeated calls (no stale-buffer leakage)", {
  # A buffer-reuse bug would most plausibly show up as results depending on
  # call history (e.g. a workspace not being fully re-zeroed between
  # permutations) - this checks that two independent calls with identical
  # inputs give bit-identical output.
  set.seed(3)
  n_genes <- 30
  n_samples <- 16
  m <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  cond_int <- as.integer(rep(1:2, each = 8))
  shuffled <- replicate(20, sample(cond_int))

  out1 <- bicorX:::cpp_bmht_permutations(m, shuffled, 0.4, 1.0, 1.0, 1L)
  out2 <- bicorX:::cpp_bmht_permutations(m, shuffled, 0.4, 1.0, 1.0, 1L)

  expect_identical(out1, out2)
})

test_that("cpp_rosdet_bootstrap_batch (workspace-based) matches an independent gold-standard computation", {
  # Same buffer-reuse refactor and rationale as the BMHT test above,
  # applied to rosdet_cpp.cpp's fast_local_bicor_ws/BicorWorkspace.
  set.seed(6)
  n_genes <- 20
  n_samples <- 20
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
              nrow = n_genes, ncol = n_samples)
  cond_int <- as.integer(rep(1:2, each = 10))

  pairs <- combn(0:(n_genes - 1), 2)
  g1 <- pairs[1, ]
  g2 <- pairs[2, ]
  n_pairs <- length(g1)

  obs_r1 <- bicorX:::cpp_rosdet_observed_bicor(m[, cond_int == 1])
  obs_r2 <- bicorX:::cpp_rosdet_observed_bicor(m[, cond_int == 2])
  obs_dist <- abs(obs_r1[cbind(g1 + 1, g2 + 1)] - obs_r2[cbind(g1 + 1, g2 + 1)])
  weights <- rep(1.0, n_pairs)

  n_perm <- 12
  set.seed(101)
  shuffled <- replicate(n_perm, sample(cond_int))

  out_new <- bicorX:::cpp_rosdet_bootstrap_batch(m, g1, g2, obs_dist, weights, shuffled, 1L)

  gold_counts <- integer(n_pairs)
  for (p in seq_len(n_perm)) {
    labels <- shuffled[, p]
    idx1 <- which(labels == 1)
    idx2 <- which(labels == 2)
    r1n <- bicorX:::cpp_rosdet_observed_bicor(m[, idx1])
    r2n <- bicorX:::cpp_rosdet_observed_bicor(m[, idx2])
    null_dist <- weights * abs(r1n[cbind(g1 + 1, g2 + 1)] - r2n[cbind(g1 + 1, g2 + 1)])
    gold_counts <- gold_counts + as.integer(null_dist >= obs_dist)
  }

  expect_identical(as.integer(out_new), gold_counts)
})

test_that("cpp_rosdet_bootstrap_batch gives identical results across repeated calls", {
  set.seed(8)
  n_genes <- 15
  n_samples <- 16
  m <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  cond_int <- as.integer(rep(1:2, each = 8))
  pairs <- combn(0:(n_genes - 1), 2)
  g1 <- pairs[1, ]; g2 <- pairs[2, ]
  obs_dist <- rep(0.1, length(g1))
  weights <- rep(1.0, length(g1))
  shuffled <- replicate(10, sample(cond_int))

  out1 <- bicorX:::cpp_rosdet_bootstrap_batch(m, g1, g2, obs_dist, weights, shuffled, 1L)
  out2 <- bicorX:::cpp_rosdet_bootstrap_batch(m, g1, g2, obs_dist, weights, shuffled, 1L)

  expect_identical(out1, out2)
})
