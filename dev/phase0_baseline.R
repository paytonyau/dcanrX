# Phase 0 baseline run.
# Purpose: capture "before" timings and output shapes on a fixed toy dataset
# so that later refactors (Phase 2+) can be checked for behavioural drift.
# Re-run this after any change to the engines and diff against
# dev/phase0_baseline_results.rds.

library(bicorX)

load(system.file("extdata", "toy_dataset.RData", package = "bicorX"))

set.seed(123)

baseline <- list(
  session_info = utils::sessionInfo(),
  r_version    = R.version.string,
  timestamp    = Sys.time()
)

cat("=== BMHT ===\n")
baseline$bmht_time <- system.time({
  baseline$bmht <- run_bmht(toy_expr, toy_condition, n_permutations = 500, workers = 1)
})

cat("\n=== BMKC ===\n")
baseline$bmkc_time <- system.time({
  baseline$bmkc <- run_bmkc(toy_expr, toy_condition, workers = 1)
})

cat("\n=== ROS-DET ===\n")
baseline$rosdet_time <- system.time({
  baseline$rosdet <- run_rosdet(toy_expr, toy_condition, min_delta = 0.3,
                                 n_permutations = 500, workers = 1)
})

# --- Summary ---
cat("\n=== Baseline summary ===\n")
cat(sprintf("BMHT:   %.3fs elapsed | %d genes scored | %d significant (FDR < %.2f)\n",
            baseline$bmht_time["elapsed"],
            nrow(baseline$bmht$results),
            sum(baseline$bmht$results$FDR < baseline$bmht$significance_level, na.rm = TRUE),
            baseline$bmht$significance_level))

cat(sprintf("BMKC:   %.3fs elapsed | %d modules found\n",
            baseline$bmkc_time["elapsed"], baseline$bmkc$n_modules))

cat(sprintf("ROS-DET: %.3fs elapsed | %d pairs tested | %d significant (FDR < %.2f)\n",
            baseline$rosdet_time["elapsed"],
            baseline$rosdet$n_tested,
            sum(baseline$rosdet$results$FDR < baseline$rosdet$significance_level, na.rm = TRUE),
            baseline$rosdet$significance_level))

# Run this script with the working directory set to dev/
saveRDS(baseline, file = "phase0_baseline_results.rds", compress = TRUE)
message("Saved baseline to: ", normalizePath("phase0_baseline_results.rds"))
