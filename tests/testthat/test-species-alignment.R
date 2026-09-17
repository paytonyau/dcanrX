test_that("species stays aligned with expr_matrix after the zero-MAD filter (BMHT/BMKC/ROS-DET)", {
  # Regression test: species was validated/used *after* zero-MAD rows were
  # dropped from expr_matrix, but was never subset by the same mask - a
  # length check could pass while species[i] no longer corresponded to
  # rownames(expr_matrix)[i]. BMKC's outer(species, species, "==") bipartite
  # mask and ROS-DET's species[g1_idx]/species[g2_idx] would then be silently
  # misaligned.
  set.seed(1)
  n_genes <- 40
  n_samples <- 20
  m <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
              nrow = n_genes, ncol = n_samples)
  rownames(m) <- sprintf("Gene_%02d", seq_len(n_genes))
  zero_mad_rows <- c(5, 20, 35)
  m[zero_mad_rows, ] <- 7.0  # constant rows -> MAD == 0

  cond <- factor(rep(c("A", "B"), each = n_samples / 2))
  species <- rep(c("host", "microbe"), length.out = n_genes)
  species[zero_mad_rows] <- "SHOULD_BE_DROPPED"

  check_alignment <- function(result) {
    expect_equal(nrow(result$data$expr_matrix), length(result$data$species))
    expect_false(any(result$data$species == "SHOULD_BE_DROPPED"))
    expect_false(any(rownames(result$data$expr_matrix) %in%
                        sprintf("Gene_%02d", zero_mad_rows)))
  }

  res_bmht <- suppressWarnings(suppressMessages(
    run_bmht(m, cond, species = species, n_permutations = 20, workers = 1)
  ))
  check_alignment(res_bmht)

  res_bmkc <- suppressWarnings(suppressMessages(
    run_bmkc(m, cond, species = species, bipartite_cliques_only = FALSE, workers = 1)
  ))
  check_alignment(res_bmkc)

  res_rosdet <- suppressWarnings(suppressMessages(
    run_rosdet(m, cond, species = species, min_delta = 0.3, n_permutations = 20, workers = 1)
  ))
  check_alignment(res_rosdet)
})
