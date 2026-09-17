test_that("run_multiomic() works across all three methods", {
  set.seed(5)
  rna <- matrix(rnorm(20 * 24, mean = 8, sd = 1.5), nrow = 20, ncol = 24,
                dimnames = list(sprintf("Gene%02d", 1:20), NULL))
  atac <- matrix(rnorm(15 * 24, mean = 8, sd = 1.5), nrow = 15, ncol = 24,
                 dimnames = list(sprintf("Peak%02d", 1:15), NULL))
  cond <- factor(rep(c("A", "B"), each = 12))

  res_rosdet <- suppressMessages(run_multiomic(
    rna, atac, cond, method = "rosdet", min_delta = 0.3, n_permutations = 20, workers = 1
  ))
  expect_s3_class(res_rosdet, "rosdet_result")

  res_bmht <- suppressMessages(run_multiomic(
    rna, atac, cond, method = "bmht", n_permutations = 20, workers = 1
  ))
  expect_s3_class(res_bmht, "bmht_result")

  res_bmkc <- suppressMessages(run_multiomic(
    rna, atac, cond, method = "bmkc", workers = 1
  ))
  expect_s3_class(res_bmkc, "bmkc_result")

  # pair_type_filter = "inter" (the default) should restrict ROS-DET output
  # to cross-layer pairs only.
  is_layer_a <- function(x) grepl("^Layer_A", x)
  expect_true(all(is_layer_a(res_rosdet$results$Gene1) != is_layer_a(res_rosdet$results$Gene2)))

  # pair_type_filter should warn (not error) if supplied for a non-ROS-DET method.
  expect_warning(
    run_multiomic(rna, atac, cond, method = "bmht", pair_type_filter = "intra",
                  n_permutations = 20, workers = 1),
    "pair_type_filter"
  )
})

test_that("run_multiomic_rosdet() (backward-compat alias) matches run_multiomic(method = 'rosdet')", {
  set.seed(5)
  rna <- matrix(rnorm(15 * 24, mean = 8, sd = 1.5), nrow = 15, ncol = 24,
                dimnames = list(sprintf("Gene%02d", 1:15), NULL))
  atac <- matrix(rnorm(10 * 24, mean = 8, sd = 1.5), nrow = 10, ncol = 24,
                 dimnames = list(sprintf("Peak%02d", 1:10), NULL))
  cond <- factor(rep(c("A", "B"), each = 12))

  res_new <- suppressMessages(run_multiomic(
    rna, atac, cond, min_delta = 0.3, n_permutations = 20, workers = 1
  ))
  res_old <- suppressMessages(run_multiomic_rosdet(
    rna, atac, cond, min_delta = 0.3, n_permutations = 20, workers = 1
  ))
  expect_identical(res_new$results, res_old$results)
})

test_that("`layer` is a working synonym for `species`, and takes precedence if both given", {
  set.seed(5)
  rna <- matrix(rnorm(15 * 24, mean = 8, sd = 1.5), nrow = 15, ncol = 24,
                dimnames = list(sprintf("RNA_%02d", 1:15), NULL))
  atac <- matrix(rnorm(10 * 24, mean = 8, sd = 1.5), nrow = 10, ncol = 24,
                 dimnames = list(sprintf("ATAC_%02d", 1:10), NULL))
  cond <- factor(rep(c("A", "B"), each = 12))
  combined <- rbind(rna, atac)
  layer_tags <- c(rep("RNA", 15), rep("ATAC", 10))

  res_layer <- suppressMessages(run_bmht(combined, cond, layer = layer_tags,
                                          n_permutations = 20, workers = 1))
  res_species <- suppressMessages(run_bmht(combined, cond, species = layer_tags,
                                            n_permutations = 20, workers = 1))
  expect_identical(res_layer$results, res_species$results)

  expect_warning(
    res_both <- run_bmht(combined, cond, species = rep("wrong", 25), layer = layer_tags,
                          n_permutations = 20, workers = 1),
    "layer"
  )
  expect_identical(res_layer$results, res_both$results)
})

test_that("the multi-omics fixture demonstrates real cross-layer detection power", {
  fixture_path <- system.file("extdata", "multiomics_dataset.RData", package = "bicorX")
  skip_if(fixture_path == "", "multiomics_dataset.RData not found")
  load(fixture_path)

  res <- suppressMessages(run_multiomic(
    multiomics_rna, multiomics_atac, multiomics_condition,
    min_delta = 0.3, n_permutations = 2000, workers = 1
  ))

  true_rna <- paste0("Layer_A_", sprintf("Gene_%02d", 1:5))
  true_atac <- paste0("Layer_B_", sprintf("Peak_%02d", 1:5))

  sig <- res$results[res$results$FDR < 0.05, ]
  expect_gt(nrow(sig), 0)

  is_true_pair <- (sig$Gene1 %in% true_rna & sig$Gene2 %in% true_atac) |
                  (sig$Gene2 %in% true_rna & sig$Gene1 %in% true_atac)
  # Most significant hits should be genuine cross-layer signal, not noise.
  expect_gte(mean(is_true_pair), 0.5)
})
