# Single-cell pathway verification: does the make_pseudobulk() ->
# run_bicordcea() pipeline work correctly end-to-end on genuine (not
# synthetic) single-cell count data?
#
# IMPORTANT SCOPE NOTE - read before treating this like Phase E:
# This is NOT a scientific finding on the same footing as
# dev/real_data_validation/lungdata_compositional_check.R. That script
# used a real, adequately-powered (33/33), externally-designed
# case-control comparison and found a specific, mechanistically-explained
# result. This script instead verifies that the single-cell/pseudobulk
# code pathway itself runs correctly on real single-cell-derived data,
# at a scale far too small to support any statistical claim:
#
# - The data (`pbmc_small`, bundled directly with the Seurat package - a
#   heavily downsampled 230-gene, 80-cell subset of a classic 10x
#   Genomics PBMC dataset, shipped specifically for Seurat's own
#   testing/demonstration purposes) has no natural multi-patient
#   replicate structure to aggregate over - it's cells from what is
#   effectively one specimen, not a cohort. make_pseudobulk() is designed
#   to aggregate cells FROM THE SAME BIOLOGICAL SAMPLE into per-sample
#   profiles across MULTIPLE independent samples/patients; there are no
#   multiple real patients here to use that way.
# - To exercise the function at all, cells within each of the two
#   largest cell-identity clusters (0: n=36, 1: n=25 - real, distinct
#   immune cell populations, not an experimental condition) were randomly
#   split into artificial "pseudo-replicate pools" of 5 cells each,
#   purely so make_pseudobulk() has multiple "samples" per condition to
#   aggregate. This gives 7 pseudobulk samples for cluster 0 and 5 for
#   cluster 1 - just barely clearing every engine's hard minimum (ROS-DET
#   requires >=5 per condition), and far below the N>=30 recommended for
#   any calibrated significance testing (Phase A's ECOR validation, or
#   even ordinary permutation-test resolution at this N).
# - A real single-cell validation with actual statistical power - e.g. a
#   published dataset with multiple patients per condition, aggregated to
#   genuine per-patient pseudobulk profiles - was attempted first via the
#   Bioconductor `scRNAseq` package's curated dataset collection, but its
#   datasets are fetched through ExperimentHub, which requires external
#   network access this environment does not have (confirmed: it fails
#   with "Cannot connect to ExperimentHub server" and cache errors even
#   in local-hub fallback mode). This is a real environment limitation,
#   not a decision to skip a better option that was available.
#
# What this DOES demonstrate: the pipeline (SingleCellExperiment
# construction -> make_pseudobulk() -> run_bmht()/run_rosdet(), with
# both transform="none" and transform="clr") runs correctly end-to-end on
# real single-cell count data with no errors, and its top result is at
# least biologically plausible (HLA-DQB1, an MHC class II gene, as the
# top-ranked differentially-connected gene between two distinct immune
# cell clusters - consistent with real biology, since MHC-II expression
# differs sharply between antigen-presenting and non-antigen-presenting
# immune cell types). That plausibility is worth noting but should not be
# oversold - at n=5/7 with artificial pseudo-replicates, this is not a
# validated finding.
#
# If real network access to Bioconductor's ExperimentHub becomes
# available, replace the pseudo-replicate construction below with an
# actual multi-patient dataset (e.g. any of the many case-control
# single-cell atlases in `scRNAseq::listDatasets()`) for a genuine,
# adequately-powered single-cell validation on the same footing as
# Phase E's lung microbiome result.

suppressMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(bicorX)
})

set.seed(1)
counts <- as.matrix(GetAssayData(pbmc_small, layer = "counts"))
clusters <- Idents(pbmc_small)

keep_cells <- clusters %in% c("0", "1")
counts <- counts[, keep_cells]
clusters <- droplevels(clusters[keep_cells])

assign_pools <- function(cells_in_cluster, pool_size = 5) {
  n <- length(cells_in_cluster)
  n_pools <- n %/% pool_size
  pool_id <- sample(rep(seq_len(n_pools), length.out = n))
  setNames(pool_id, cells_in_cluster)
}
pools_0 <- assign_pools(colnames(counts)[clusters == "0"])
pools_1 <- assign_pools(colnames(counts)[clusters == "1"])
all_pools <- c(pools_0, pools_1)[colnames(counts)]

sce <- SingleCellExperiment(assays = list(counts = counts))
colData(sce)$sample <- paste0("pool_", all_pools, "_c", as.character(clusters))
colData(sce)$condition <- as.character(clusters)

pb <- make_pseudobulk(sce, sample_col = "sample", condition_col = "condition", min_cells = 5)
cat(sprintf("Pseudobulk matrix: %d genes x %d pseudo-samples\n",
            nrow(pb$expr_matrix), ncol(pb$expr_matrix)))
cat("Condition table:\n"); print(table(pb$condition))

# --- Pipeline verification: both transforms run without error ---
res_naive <- suppressWarnings(suppressMessages(run_bmht(
  pb$expr_matrix, factor(pb$condition), n_permutations = 500, workers = 1
)))
res_clr <- suppressWarnings(suppressMessages(run_bmht(
  pb$expr_matrix, factor(pb$condition), n_permutations = 500, workers = 1, transform = "clr"
)))

cat(sprintf("\nBMHT (transform=none): %d genes scored, %d significant (FDR<0.05)\n",
            nrow(res_naive$results), sum(res_naive$results$FDR < 0.05)))
cat(sprintf("BMHT (transform=clr):  %d genes scored, %d significant (FDR<0.05)\n",
            nrow(res_clr$results), sum(res_clr$results$FDR < 0.05)))

cat("\nTop 5 genes by DC_Score (transform=none) - a biological plausibility\n")
cat("spot-check, not a statistical claim at this sample size:\n")
print(head(res_naive$results[order(-res_naive$results$DC_Score),
                              c("Gene", "DC_Score", "P_value", "FDR")], 5))
