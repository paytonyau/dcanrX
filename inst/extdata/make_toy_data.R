# Generates the toy fixture used for Phase 0 baseline checks and later tests.
# 200 genes x 24 samples, 2 conditions x 12 replicates, with:
#   - 20 genes carrying a real differential-correlation signal (Module_A)
#   - the rest background noise genes
# Deterministic under set.seed(42).

set.seed(42)
n_genes   <- 200
n_samples <- 24
cond      <- factor(rep(c("A", "B"), each = n_samples / 2))

base_expr <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
                     nrow = n_genes, ncol = n_samples)

# Inject a differentially co-expressed module: genes 1-20 correlated in A,
# uncorrelated in B.
latent_A <- rnorm(sum(cond == "A"))
for (g in 1:20) {
  base_expr[g, cond == "A"] <- latent_A * 2 + rnorm(sum(cond == "A"), sd = 0.5) + 8
  base_expr[g, cond == "B"] <- rnorm(sum(cond == "B"), mean = 8, sd = 1.5)
}

rownames(base_expr) <- sprintf("Gene_%03d", seq_len(n_genes))
colnames(base_expr) <- sprintf("Sample_%02d", seq_len(n_samples))

toy_expr <- base_expr
toy_condition <- cond

# Run this script with the working directory set to inst/extdata/
# (or pass an explicit `out_path` variable before sourcing).
out_path <- if (exists("out_path")) out_path else "toy_dataset.RData"
save(toy_expr, toy_condition, file = out_path)
message("Saved toy dataset to: ", normalizePath(out_path))
