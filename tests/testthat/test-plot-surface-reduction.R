test_that("the Phase 3 core-8 plotting functions exist and are exported", {
  core_plots <- c(
    "plot_rewiring_volcano", "plot_top_gene_bar",
    "plot_module_heatmap", "plot_module_eigengene", "plot_umap_eigengene",
    "plot_differential_network",
    "plot_rosdet_quadrant", "plot_edge_volcano"
  )
  for (fn in core_plots) {
    expect_true(exists(fn, where = asNamespace("bicorX"), mode = "function"),
                info = fn)
    expect_true(fn %in% getNamespaceExports("bicorX"), info = fn)
  }
})

test_that("removed plotting functions are actually gone (not just unexported)", {
  # Spot-check a representative few of the ~30 functions removed in the
  # Phase 3 surface reduction - regression guard against silently
  # reintroducing them without a deliberate decision.
  removed <- c(
    "plot_3d_network", "plot_hive_network", "plot_intercellular_network",
    "plot_eruption", "plot_hidden_drivers", "plot_causal_network",
    "plot_ego_switching", "plot_raincloud_eigengene", "plot_taxa_rewiring",
    "plot_method_overlap", "plot_enrichment_dotplot",
    "plot_switching_network", "plot_spatial_coords"
  )
  for (fn in removed) {
    expect_false(exists(fn, where = asNamespace("bicorX"), mode = "function"),
                 info = fn)
  }
})

test_that("plot.<result> dispatchers only offer type= options for surviving functions", {
  expect_setequal(eval(formals(bicorX:::plot.bmht_result)$type),
                   c("default", "volcano", "bar", "topology"))
  expect_setequal(eval(formals(bicorX:::plot.bmkc_result)$type),
                   c("default", "network", "heatmap", "eigengene", "umap"))
  expect_setequal(eval(formals(bicorX:::plot.rosdet_result)$type),
                   c("default", "quadrant", "volcano"))
})

test_that("plot.<result> dispatchers render successfully end-to-end for every surviving type", {
  toy_path <- system.file("extdata", "toy_dataset.RData", package = "bicorX")
  skip_if(toy_path == "", "toy_dataset.RData not found")
  load(toy_path)

  res_bmht <- suppressMessages(run_bmht(toy_expr, toy_condition, n_permutations = 20, workers = 1))
  res_rosdet <- suppressMessages(run_rosdet(toy_expr, toy_condition, min_delta = 0.3,
                                             n_permutations = 20, workers = 1))

  for (t in c("default", "volcano", "bar", "topology")) {
    expect_no_error(plot(res_bmht, type = t))
  }
  for (t in c("default", "quadrant", "volcano")) {
    expect_no_error(plot(res_rosdet, type = t))
  }
})
