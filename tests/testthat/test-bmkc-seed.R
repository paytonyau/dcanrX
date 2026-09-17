test_that("run_bmkc is reproducible under a fixed seed and does not warn about RNG", {
  # Regression test: an earlier version called arma_rng::set_seed_random()
  # inside get_quantile(), which reseeds from system entropy on every call -
  # this ignored R's set.seed() entirely and made percentile-threshold runs
  # non-reproducible. A later fix attempt called arma_rng::set_seed()
  # explicitly from C++, which is also wrong for this RcppArmadillo build
  # (its RNG delegates to R's own generator) and triggers Rcpp's own
  # "RNG seed has to be set at the R level" warning. The correct fix needs
  # no C++-side seeding at all: R's set.seed(), called once at the top of
  # .run_bmkc_matrix, is sufficient and is the only mechanism now used.
  toy_path <- system.file("extdata", "toy_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "toy_dataset.RData not found")
  load(toy_path)

  expect_no_warning({
    res1 <- run_bmkc(toy_expr, toy_condition, threshold_method = "percentile",
                      T1 = 80, T2 = 20, seed = 123, workers = 1)
  })
  res2 <- run_bmkc(toy_expr, toy_condition, threshold_method = "percentile",
                    T1 = 80, T2 = 20, seed = 123, workers = 1)

  expect_identical(res1$modules, res2$modules)
  expect_identical(res1$parameters, res2$parameters)
})
