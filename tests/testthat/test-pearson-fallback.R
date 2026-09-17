test_that("C++ bicor kernels fall back to Pearson for zero-MAD genes, matching R bicor()", {
  # Regression test: a gene with mad_x == 0 (>=50% of values tied at the
  # median, but not necessarily constant) got an all-zero X_tilde column in
  # all three C++ kernels, forcing every correlation involving that gene to
  # exactly 0 - a real divergence from R's bicor()/WGCNA::bicor(), both of
  # which fall back to standard Pearson correlation for that gene
  # (pearsonFallback = "individual"). The zero-MAD prefilter each engine
  # applies upstream only catches genes constant across ALL samples; a gene
  # constant within just one condition (exactly what a differential test
  # should care about) reaches these kernels with mad_x == 0 and was
  # previously silently zeroed out instead of falling back correctly.
  set.seed(11)
  n_genes <- 10
  n_samples <- 15
  m <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)

  # Gene 3: >=50% tied at the median (mad == 0) but has nonzero variance ->
  # should fall back to Pearson, NOT be forced to 0.
  m[3, ] <- c(rep(5, 10), rnorm(5))
  # Gene 7: fully constant (sd == 0 too) -> R's bicor() returns exactly 0
  # for this case; the fallback must not produce NaN/Inf.
  m[7, ] <- rep(2.0, n_samples)

  ref <- matrix(0, n_genes, n_genes)
  for (i in seq_len(n_genes)) {
    for (j in seq_len(n_genes)) {
      if (i != j) ref[i, j] <- bicorX::bicor(m[i, ], m[j, ])
    }
  }

  cpp_bmht   <- bicorX:::cpp_bmht_observed_bicor(m)
  cpp_rosdet <- bicorX:::cpp_rosdet_observed_bicor(m)
  cpp_bmkc   <- bicorX:::cpp_fast_bicor_matrix(m, 1L)

  expect_equal(cpp_bmht,   ref, tolerance = 1e-10)
  expect_equal(cpp_rosdet, ref, tolerance = 1e-10)
  expect_equal(cpp_bmkc,   ref, tolerance = 1e-10)

  # Gene 3's row should NOT be all zero (that would mean the fallback isn't
  # firing and it's still being forced to 0).
  expect_false(all(cpp_bmht[3, -3] == 0))

  # Gene 7 (truly constant) should be exactly 0 everywhere, and finite.
  expect_true(all(cpp_bmht[7, -7] == 0))
  expect_true(all(is.finite(cpp_bmht[7, ])))
})
