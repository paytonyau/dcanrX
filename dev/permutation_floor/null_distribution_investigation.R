# Investigation: does BMHT's DC_Score null distribution have a
# recognizable parametric shape that could support a distribution-specific
# p-value model (Tier 1, next step (a) from TODO.md, following the failed
# generic-GPD attempt documented in NEWS.md / test-gpd-tail-calibration.R)?
#
# Requires: sn (skew-normal/skew-t distributions) - not otherwise a
# package dependency, install via:
#   apt-get install -y r-cran-sn   (or install.packages("sn") with CRAN access)
#
# SUMMARY OF FINDINGS (see NEWS.md for the full writeup):
# 1. DC_Score's null distribution is NOT chi-squared/gamma-like as
#    initially hypothesized (it's an average, not a sum, of many
#    per-edge squared terms, so CLT pulls it toward Gaussian-like shape).
#    It is instead close to Gaussian for most genes, with a real but
#    variable amount of right skew for others (skewness ranging from
#    ~-0.01 to ~1.34 across genes in the toy dataset).
# 2. Skewness does NOT correlate cleanly with the number of "informative"
#    edges contributing to a gene's score (correlation ~0.15, weak and in
#    the opposite direction from the "fewer edges -> more chi-like"
#    hypothesis) - no simple design-derived parameter explains it.
# 3. A plain normal fit to the bulk of the permutation null looks
#    deceptively good in a quick quantile-by-eye comparison (differences
#    of 0.01-0.09 in absolute quantile position) but is BADLY
#    anti-conservative once evaluated properly for p-value calibration:
#    14.5x too many p-values below 0.001 in a self-consistency simulation
#    using the actual empirical null shape from real data (not a
#    synthetic proxy this time). Small-looking quantile errors translate
#    into large tail-probability errors.
# 4. A skew-normal fit (3-parameter: location, scale, shape) is a real,
#    substantial improvement over plain normal (3.5x vs. 14.5x inflation
#    at p<0.001) but still clearly fails calibration - meaning genes near
#    that threshold would get roughly 3.5x too many false positives.
# 5. This is now the FOURTH distinct parametric/semi-parametric approach
#    tested (GPD method-of-moments, GPD MLE, plain normal, skew-normal),
#    and all four show the same systematic pattern: underestimating how
#    thin the true tail actually is, to different degrees. Skew-normal is
#    the best of the four, but "best of four failures" is not a basis for
#    shipping this - none of the four is being integrated into any
#    engine.
#
# CONCLUSION: distribution-specific parametric modeling does not appear
# to close the permutation-floor gap safely with the approaches tried so
# far. TODO.md's recommendation (c) - defaulting to a much higher
# n_permutations, now that per-permutation cost and memory are both
# reduced (see the Performance/RAM NEWS.md sections) - remains the
# safest path. A genuinely different idea (e.g. an analytically-derived
# asymptotic distribution for bicor-based DC-scores, rather than an
# empirically-fit one) is a bigger undertaking not attempted here.

suppressMessages({
  library(bicorX)
  library(sn)
})

load(system.file("extdata", "toy_dataset.RData", package = "bicorX"))
cond_int <- as.integer(toy_condition)

cor_1 <- bicorX:::cpp_bmht_observed_bicor(toy_expr[, cond_int == 1])
cor_2 <- bicorX:::cpp_bmht_observed_bicor(toy_expr[, cond_int == 2])
mad_1 <- stats::mad(cor_1[lower.tri(cor_1)], constant = 1)
mad_2 <- stats::mad(cor_2[lower.tri(cor_2)], constant = 1)

set.seed(1)
n_perm <- 20000
shuffled <- replicate(n_perm, sample(cond_int))
null_mat <- bicorX:::cpp_bmht_permutations(toy_expr, shuffled, 0.4, mad_1, mad_2, 1L)

# --- Step 1: shape survey across genes ---
skews <- apply(null_mat, 1, function(x) mean((x - mean(x))^3) / sd(x)^3)
cat("Skewness range across genes:", paste(round(range(skews), 3), collapse = " to "), "\n")

informative_mask <- (abs(cor_1) > 0.4) | (abs(cor_2) > 0.4)
n_edges <- rowSums(informative_mask, na.rm = TRUE)
cat("Correlation(skewness, n_informative_edges):", round(cor(n_edges, skews), 3), "\n\n")

# --- Step 2: calibration self-consistency check on the most-skewed gene ---
worst_gene <- which.max(skews)
pool <- null_mat[worst_gene, ]
n_pool <- length(pool)

set.seed(3)
n_reps <- 2000
n_null_for_fit <- 2000
p_sn <- numeric(n_reps)
p_normal <- numeric(n_reps)
p_standard <- numeric(n_reps)

for (i in seq_len(n_reps)) {
  idx <- sample(n_pool, n_null_for_fit + 1)
  obs <- pool[idx[1]]
  ns <- pool[idx[-1]]

  mu <- mean(ns); sigma <- sd(ns)
  p_normal[i] <- 1 - pnorm(obs, mu, sigma)

  fit <- sn::sn.mple(y = ns, opt.method = "nlminb")
  dp <- sn::cp2dp(fit$cp, family = "SN")
  p_sn[i] <- 1 - sn::psn(obs, xi = dp["xi"], omega = dp["omega"], alpha = dp["alpha"])

  p_standard[i] <- (sum(ns >= obs) + 1) / (length(ns) + 1)
}

cat("=== Calibration check on gene", worst_gene, "(skewness =", round(skews[worst_gene], 3), ") ===\n")
for (t in c(0.05, 0.01, 0.001)) {
  cat(sprintf("p<%.3f -- standard: %.2fx | normal: %.2fx | skew-normal: %.2fx  (ratio to nominal)\n",
              t, mean(p_standard < t) / t, mean(p_normal < t) / t, mean(p_sn < t) / t))
}
