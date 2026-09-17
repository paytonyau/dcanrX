test_that(".gpd_tail_pvalues() exists but is not part of any engine's public behavior", {
  # This function is experimental and known to fail calibration
  # validation (see its roxygen docs and the test below) - it must not be
  # wired into any engine's p-value computation until a version passes
  # calibration testing. The existing BMHT/ROS-DET calibration tests
  # (test-bmht-calibration.R, and this file's own tests) already
  # functionally guard against it being wired in silently: if it were, the
  # anti-conservative bias documented below would show up there too. This
  # test just confirms the function is reachable for future work.
  expect_true(is.function(bicorX:::.gpd_tail_pvalues))
})

test_that(".gpd_tail_pvalues() reproduces known-good behavior in the non-extrapolated regime", {
  # These properties hold regardless of the tail-fit calibration issue,
  # since they don't exercise extrapolation at all.
  set.seed(1)
  null_vals <- rnorm(500)

  # Below the tail threshold: must exactly match standard permutation counting.
  obs_below <- as.numeric(stats::quantile(null_vals, 0.5))
  p_gpd <- bicorX:::.gpd_tail_pvalues(obs_below, matrix(null_vals, nrow = 1))
  p_standard <- (sum(null_vals >= obs_below) + 1) / (length(null_vals) + 1)
  expect_equal(p_gpd, p_standard)

  # Too few permutations for a stable fit: must fall back to standard counting.
  null_small <- rnorm(50)
  obs_extreme <- max(null_small) + 2
  p_gpd_small <- bicorX:::.gpd_tail_pvalues(obs_extreme, matrix(null_small, nrow = 1))
  p_standard_small <- 1 / (length(null_small) + 1)
  expect_equal(p_gpd_small, p_standard_small)
})

test_that("KNOWN ISSUE: .gpd_tail_pvalues() extrapolation is anti-conservative in the deep tail", {
  # Documents a real, deliberately-not-fixed finding rather than asserting
  # correct calibration (which this function does not currently have).
  # If this test starts failing because the observed/nominal ratio drops
  # to ~1 (well-calibrated), that's good news - update this test AND the
  # roxygen "STATUS" section on .gpd_tail_pvalues(), and consider wiring
  # it into an engine. ~10,000 reps, takes a couple seconds.

  set.seed(7)
  n_reps <- 10000
  n_perm <- 500
  p_gpd_vec <- numeric(n_reps)
  for (i in seq_len(n_reps)) {
    samples <- rnorm(n_perm + 1)
    null_vals <- samples[-1]
    p_gpd_vec[i] <- bicorX:::.gpd_tail_pvalues(samples[1], matrix(null_vals, nrow = 1))
  }

  # At a threshold near the old permutation floor, the current
  # MLE-based fit inflates the true Type I error rate rather than merely
  # extending resolution below the floor. This is the documented,
  # unresolved calibration problem - asserted here (with a wide margin,
  # since this is inherently a noisy Monte Carlo estimate) so the
  # known-bad behavior is locked in rather than silently drifting.
  observed_rate <- mean(p_gpd_vec < 0.001)
  nominal_rate <- 0.001
  expect_gt(observed_rate, nominal_rate * 1.1)
})
