test_that("weight_mode default is 'unweighted', and 'wcor' can invert ranking on real signal", {
  # Phase 5 benchmark finding: 'wcor' weighting can produce a LOWER mean
  # Distance_Score for true differential pairs than for background/null
  # pairs, actively inverting the useful ranking signal - not just
  # rescaling it. This is why the default changed from 'wcor' to
  # 'unweighted'. Reproduces a bounded version of the compositional
  # benchmark scenario (dev/phase5_compositional_benchmark.R) that
  # triggered the finding: a true differential module where signal genes
  # also undergo a genuine variance shift between conditions (which
  # 'wcor' penalizes) alongside unrelated background pairs.
  expect_equal(eval(formals(bicorX:::.run_rosdet_matrix)$weight_mode)[1], "unweighted")

  set.seed(11)
  n_features <- 30
  n_samples <- 24
  cond <- factor(rep(c("A", "B"), each = 12))
  m <- matrix(rlnorm(n_features * n_samples, meanlog = 3, sdlog = 0.4),
              nrow = n_features, ncol = n_samples)
  rownames(m) <- sprintf("Feature_%02d", seq_len(n_features))

  # True module: correlated + lower-variance in condition A; independent +
  # higher-variance (natural background level) in condition B - a genuine
  # variance shift alongside the correlation shift, exactly what 'wcor'
  # penalizes.
  latent_A <- rnorm(12, sd = 1)
  for (i in 1:8) {
    m[i, cond == "A"] <- exp(3 + latent_A * 1.2 + rnorm(12, sd = 0.15))
  }

  true_module <- sprintf("Feature_%02d", 1:8)
  background <- sprintf("Feature_%02d", 9:30)
  true_pairs <- t(combn(true_module, 2))
  bg_pairs <- t(combn(background, 2))

  get_scores <- function(weight_mode) {
    res <- suppressMessages(run_rosdet(m, cond, min_delta = 0.0, n_permutations = 30,
                                        weight_mode = weight_mode, workers = 1))
    key <- ifelse(res$results$Gene1 < res$results$Gene2,
                  paste(res$results$Gene1, res$results$Gene2),
                  paste(res$results$Gene2, res$results$Gene1))
    score_by_key <- tapply(res$results$Distance_Score, key, max)
    list(
      true_mean = mean(score_by_key[paste(true_pairs[, 1], true_pairs[, 2])], na.rm = TRUE),
      bg_mean   = mean(score_by_key[paste(bg_pairs[, 1], bg_pairs[, 2])], na.rm = TRUE)
    )
  }

  unweighted <- get_scores("unweighted")
  wcor <- get_scores("wcor")

  # Unweighted correctly separates true signal from background.
  expect_gt(unweighted$true_mean, unweighted$bg_mean)

  # 'wcor' at minimum substantially closes or inverts that gap on this
  # kind of data - this is the point of the test, not a coincidence.
  unweighted_gap <- unweighted$true_mean - unweighted$bg_mean
  wcor_gap <- wcor$true_mean - wcor$bg_mean
  expect_lt(wcor_gap, unweighted_gap)
})

test_that("weight_mode affects Distance_Score ranking but not P_value/FDR", {
  # This is documented, intentional behavior (see .run_rosdet_matrix's
  # @param weight_mode docs) - not a bug, but a surprising one, so it's
  # locked in here to prevent silent drift in either direction:
  #   - if this test starts failing because p-values now DO differ, that's
  #     a real behavior change and should be called out in NEWS.md.
  #   - if Distance_Score stops differing, weight_mode has become a no-op
  #     entirely and should probably be removed.
  toy_path <- system.file("extdata", "toy_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "toy_dataset.RData not found")
  load(toy_path)

  res_wcor <- suppressMessages(run_rosdet(
    toy_expr, toy_condition, min_delta = 0.3, n_permutations = 50,
    weight_mode = "wcor", workers = 1, seed = 1
  ))
  res_unw <- suppressMessages(run_rosdet(
    toy_expr, toy_condition, min_delta = 0.3, n_permutations = 50,
    weight_mode = "unweighted", workers = 1, seed = 1
  ))

  key <- function(df) paste(df$Gene1, df$Gene2)
  a <- res_wcor$results[order(key(res_wcor$results)), ]
  b <- res_unw$results[order(key(res_unw$results)), ]

  expect_identical(a$P_value, b$P_value)
  expect_identical(a$FDR, b$FDR)
  expect_false(identical(a$Distance_Score, b$Distance_Score))
})
