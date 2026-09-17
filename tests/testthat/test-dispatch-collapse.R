test_that("run_bicordcea() dispatches correctly across all supported input types", {
  # Regression test for the Phase 2 dispatch collapse: run_dcanr.* and
  # run_bmht.*/run_rosdet.*/run_bmkc.* used to be two parallel, largely
  # duplicated S3 trees. Now run_bicordcea() is the only generic; the
  # matrix/data.frame/SummarizedExperiment/Seurat/MultiAssayExperiment
  # methods all funnel through .extract_inputs() + .run_engine(), and
  # run_bmht()/run_rosdet()/run_bmkc() are thin non-generic wrappers around
  # run_bicordcea(method = ...).
  toy_path <- system.file("extdata", "toy_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "toy_dataset.RData not found")
  load(toy_path)

  # --- matrix, via run_bicordcea() directly and via run_bmht() ---
  res_direct <- suppressMessages(run_bicordcea(
    toy_expr, toy_condition, method = "bmht", n_permutations = 30, workers = 1
  ))
  res_wrapper <- suppressMessages(run_bmht(
    toy_expr, toy_condition, n_permutations = 30, workers = 1
  ))
  expect_s3_class(res_direct, "bmht_result")
  expect_identical(res_direct$results, res_wrapper$results)

  # --- expr_matrix= named argument still works (unchanged contract) ---
  res_named <- suppressMessages(run_bmht(
    expr_matrix = toy_expr, condition = toy_condition, n_permutations = 30, workers = 1
  ))
  expect_identical(res_direct$results, res_named$results)

  # --- data.frame input (previously unsupported entirely) ---
  res_df <- suppressMessages(run_bicordcea(
    as.data.frame(toy_expr), toy_condition, method = "bmht",
    n_permutations = 30, workers = 1
  ))
  expect_s3_class(res_df, "bmht_result")

  # --- run_dcanr() deprecated alias: warns, but produces identical output ---
  expect_warning(
    res_deprecated <- run_dcanr(
      toy_expr, toy_condition, method = "bmht", n_permutations = 30, workers = 1
    ),
    "deprecated"
  )
  expect_identical(res_direct$results, res_deprecated$results)
})

test_that("run_bicordcea() dispatches SummarizedExperiment input correctly", {
  skip_if_not_installed("SummarizedExperiment")
  toy_path <- system.file("extdata", "toy_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "toy_dataset.RData not found")
  load(toy_path)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = toy_expr),
    colData = S4Vectors::DataFrame(group = toy_condition)
  )

  res <- suppressMessages(run_bicordcea(
    se, condition = "group", method = "bmht", assay_name = "counts",
    n_permutations = 30, workers = 1
  ))
  res_wrapper <- suppressMessages(run_bmht(
    se, condition = "group", assay_name = "counts", n_permutations = 30, workers = 1
  ))

  expect_s3_class(res, "bmht_result")
  expect_identical(res$results, res_wrapper$results)
})

test_that("run_bicordcea() dispatches MultiAssayExperiment input correctly", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  set.seed(5)
  rna <- matrix(rnorm(15 * 24, mean = 8, sd = 1.5), nrow = 15, ncol = 24,
                dimnames = list(sprintf("Gene_%02d", 1:15), sprintf("S%02d", 1:24)))
  atac <- matrix(rnorm(10 * 24, mean = 8, sd = 1.5), nrow = 10, ncol = 24,
                 dimnames = list(sprintf("Peak_%02d", 1:10), sprintf("S%02d", 1:24)))
  cond <- factor(rep(c("A", "B"), each = 12))

  rna_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = rna))
  atac_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = atac))
  coldata <- S4Vectors::DataFrame(group = cond, row.names = sprintf("S%02d", 1:24))

  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(RNA = rna_se, ATAC = atac_se),
    colData = coldata
  )

  res <- suppressMessages(run_bicordcea(
    mae, condition = "group", method = "rosdet",
    rna_assay = "RNA", epi_assay = "ATAC",
    min_delta = 0.3, n_permutations = 20, workers = 1
  ))
  res_wrapper <- suppressMessages(run_rosdet(
    mae, condition = "group", rna_assay = "RNA", epi_assay = "ATAC",
    min_delta = 0.3, n_permutations = 20, workers = 1
  ))

  expect_s3_class(res, "rosdet_result")
  expect_gt(res$n_tested, 0)
  expect_identical(res$results, res_wrapper$results)
})
