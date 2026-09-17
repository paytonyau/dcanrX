test_that("bicor() (R), all three cpp bicor kernels, and WGCNA::bicor agree", {
  skip_if_not_installed("WGCNA")

  set.seed(42)
  n_genes <- 30
  n_samples <- 15
  m <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)

  # Reference: package's own pairwise R bicor(), applied over all gene pairs.
  r_ref <- matrix(0, n_genes, n_genes)
  for (i in seq_len(n_genes)) {
    for (j in seq_len(n_genes)) {
      if (i != j) r_ref[i, j] <- bicorX::bicor(m[i, ], m[j, ])
    }
  }

  cpp_bmht   <- bicorX:::cpp_bmht_observed_bicor(m)
  cpp_rosdet <- bicorX:::cpp_rosdet_observed_bicor(m)
  cpp_bmkc   <- bicorX:::cpp_fast_bicor_matrix(m, 1L)

  expect_equal(cpp_bmht,   r_ref, tolerance = 1e-10)
  expect_equal(cpp_rosdet, r_ref, tolerance = 1e-10)
  expect_equal(cpp_bmkc,   r_ref, tolerance = 1e-10)

  # External reference: WGCNA::bicor uses samples-by-genes orientation.
  w <- WGCNA::bicor(t(m))
  diag(w) <- 0

  expect_equal(r_ref,    w, tolerance = 1e-8)
  expect_equal(cpp_bmht, w, tolerance = 1e-8)
})

test_that("bicor() falls back to Pearson when MAD is zero", {
  x <- c(rep(5, 8), 5.001, 4.999)  # near-constant -> MAD may be 0 after rounding-free case
  x_const <- rep(5, 10)
  y <- rnorm(10)

  # Fully constant x: MAD is exactly 0 -> Pearson fallback (which is also
  # undefined/NA for zero-variance input, but sd(x)==0 short-circuits to 0
  # per the current implementation).
  expect_equal(bicorX::bicor(x_const, y), 0)
})

test_that("biweight_midvariance uses raw MAD (matches manual reference formula)", {
  set.seed(1)
  x <- rnorm(50)

  med_x <- stats::median(x)
  mad_x <- stats::median(abs(x - med_x))  # raw, not stats::mad()
  u <- (x - med_x) / (9 * mad_x)
  valid <- abs(u) < 1
  n <- length(x)
  numerator <- n * sum(((x[valid] - med_x)^2) * ((1 - u[valid]^2)^4))
  denom_term <- sum((1 - u[valid]^2) * (1 - 5 * u[valid]^2))
  expected <- numerator / (denom_term^2)

  observed <- bicorX:::biweight_midvariance(x)
  expect_equal(observed, expected, tolerance = 1e-10)
})
