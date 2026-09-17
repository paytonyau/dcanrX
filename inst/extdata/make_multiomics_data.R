# Generates a synthetic RNA + ATAC multi-omics fixture with a genuine
# cross-layer differential signal, for use with run_multiomic() /
# run_bicordcea(..., layer = ...).
#
# Design: 5 RNA genes and 5 ATAC peaks share a common latent factor in
# condition A only (as if a regulatory element's accessibility drives a
# gene's expression, but only under condition A) - uncorrelated in
# condition B. This is exactly the kind of signal pair_type_filter =
# "inter" is meant to surface: a cross-layer relationship that appears in
# one condition and not the other.
#
# Deterministic under set.seed(7).

set.seed(7)

n_genes  <- 30
n_peaks  <- 20
n_samples <- 24
cond <- factor(rep(c("A", "B"), each = n_samples / 2))

rna  <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
               nrow = n_genes, ncol = n_samples)
atac <- matrix(rnorm(n_peaks * n_samples, mean = 5, sd = 1.0),
               nrow = n_peaks, ncol = n_samples)

# Inject the cross-layer module: RNA genes 1-5 and ATAC peaks 1-5 share a
# latent factor in condition A only.
latent_A <- rnorm(sum(cond == "A"))
for (g in 1:5) {
  rna[g, cond == "A"]  <- 8 + latent_A * 2 + rnorm(sum(cond == "A"), sd = 0.4)
  rna[g, cond == "B"]  <- rnorm(sum(cond == "B"), mean = 8, sd = 1.5)
}
for (p in 1:5) {
  atac[p, cond == "A"] <- 5 + latent_A * 1.5 + rnorm(sum(cond == "A"), sd = 0.4)
  atac[p, cond == "B"] <- rnorm(sum(cond == "B"), mean = 5, sd = 1.0)
}

rownames(rna)  <- sprintf("Gene_%02d", seq_len(n_genes))
rownames(atac) <- sprintf("Peak_%02d", seq_len(n_peaks))
colnames(rna)  <- colnames(atac) <- sprintf("Sample_%02d", seq_len(n_samples))

multiomics_rna  <- rna
multiomics_atac <- atac
multiomics_condition <- cond

# Run this script with the working directory set to inst/extdata/
out_path <- if (exists("out_path")) out_path else "multiomics_dataset.RData"
save(multiomics_rna, multiomics_atac, multiomics_condition, file = out_path)
message("Saved multi-omics dataset to: ", normalizePath(out_path))
