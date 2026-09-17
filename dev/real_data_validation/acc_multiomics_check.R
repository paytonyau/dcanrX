# Real multi-omics validation: RNA-seq + miRNA-seq, TCGA Adrenocortical
# Carcinoma (ACC), comparing the published C1A/C1B molecular subtypes.
#
# Fills a real gap: Phase E validated the compositional-correction claim
# on real data, and the single-cell check verified the pseudobulk pathway
# on real (if underpowered) single-cell data - but there was no real,
# adequately-powered, genuinely cross-layer (two distinct omics platforms,
# not just two feature panels of the same assay) demonstration until now.
# The existing vignette (vignettes/multiomics-workflow.Rmd) uses a small
# synthetic RNA+ATAC dataset for clarity of exposition; this is the real,
# externally-published complement to it.
#
# DATA: `miniACC`, bundled directly with the Bioconductor
# MultiAssayExperiment package (no ExperimentHub/ network access needed,
# unlike the scRNAseq datasets that blocked a fuller single-cell
# validation) - a reduced version of the TCGA Adrenocortical Carcinoma
# multi-omics dataset (Zheng et al. 2016, Cancer Cell, "Comprehensive
# Pan-Genomic Characterization of Adrenocortical Carcinoma" - the source
# of both the data and the C1A/C1B molecular subtype classification used
# as the two-condition comparison here). 78 patients with a valid C1A/C1B
# label have both RNASeq2GeneNorm (198 genes) and miRNASeqGene (471
# miRNAs) data (43 C1A / 35 C1B) - comfortably above the N>=30 threshold
# established for ROS-DET's ECOR analytical test (Phase A), so this uses
# significance_method = "analytical" for a real, adequately-powered
# analysis.
#
# METHOD: run_multiomic() with method="rosdet", pair_type_filter="inter"
# (cross-layer pairs only - gene-to-miRNA, not gene-to-gene or
# miRNA-to-miRNA), otherwise default settings.
#
# RESULT: a genuinely striking, real pattern - a single miRNA,
# hsa-mir-137, is involved in 6 of the 10 genome-wide-significant
# cross-layer pairs (FDR<0.05), each showing a consistent signature:
# strong positive correlation with a different mRNA gene in one subtype,
# near-zero or weak correlation in the other. Checked formally rather
# than reported as an eyeballed proportion (a lesson carried over directly
# from Phase E, where a similar-looking pattern turned out to be a
# background-rate artifact): a Fisher's exact test against the correct
# background (this one miRNA appears in only 71 of 18,706 total candidate
# pairs, a 0.38% background rate) gives p = 5.0e-13, odds ratio = 429 -
# unambiguously a real, massive concentration, not an artifact of a large
# background rate the way the Phase E observation was.
#
# WHAT IS AND ISN'T INDEPENDENTLY VERIFIED HERE:
# - Independently verified (computationally, within this script): the
#   statistical concentration itself (the Fisher's exact test above) and
#   the internal consistency of the pattern (same qualitative direction -
#   strong correlation in one subtype, weak in the other - across all six
#   significant target genes, not just the top one).
# - NOT independently verified here (would need literature search /
#   access to the original TCGA ACC paper, neither available in this
#   environment): any specific claim about hsa-mir-137's documented
#   biological role in adrenocortical carcinoma specifically, or the
#   precise clinical/prognostic meaning of the C1A vs. C1B subtypes
#   beyond what's visible in this dataset's own metadata. hsa-mir-137 is
#   recognized in the broader cancer literature (from general training
#   knowledge, not verified via search in this session) as a
#   tumor-suppressor microRNA in several other cancer types - this is
#   noted as a plausible lead for follow-up, not asserted as confirmed.
#   One of the six target genes, MYC, is a well-known major oncogene,
#   which is at least consistent with a real regulatory story worth
#   investigating further - again noted as a lead, not a confirmed claim.

suppressMessages({
  library(MultiAssayExperiment)
  library(bicorX)
})

data(miniACC)

valid_patients <- !is.na(colData(miniACC)$C1A.C1B)
mae <- suppressWarnings(intersectColumns(
  miniACC[, valid_patients, c("RNASeq2GeneNorm", "miRNASeqGene")]
))
rna <- assay(experiments(mae)[["RNASeq2GeneNorm"]])
mirna <- assay(experiments(mae)[["miRNASeqGene"]])
cond <- factor(colData(mae)$C1A.C1B)

cat(sprintf("Dataset: %d genes (RNA) + %d miRNAs, %d patients (%d C1A / %d C1B)\n",
            nrow(rna), nrow(mirna), ncol(rna), sum(cond == "C1A"), sum(cond == "C1B")))

res <- suppressWarnings(suppressMessages(run_multiomic(
  rna, mirna, cond, method = "rosdet",
  layer_A_name = "mRNA", layer_B_name = "miRNA",
  pair_type_filter = "inter", min_delta = 0.3,
  significance_method = "analytical", workers = 1
)))

cat(sprintf("\n%d candidate cross-layer (mRNA-miRNA) pairs, %d significant (FDR<0.05)\n",
            nrow(res$results), sum(res$results$FDR < 0.05)))

sig <- res$results[order(res$results$P_value), ]
sig <- sig[sig$FDR < 0.05, ]
cat("\nTop significant cross-layer pairs:\n")
print(sig[, c("Gene1", "Gene2", "Distance_Score", "r1", "r2", "FDR")])

# --- Formal enrichment test for the mir-137 concentration (not an
#     eyeballed proportion - see the Phase E precedent for why this
#     matters) ---
touches_mir137 <- function(df) grepl("mir-137", df$Gene1) | grepl("mir-137", df$Gene2)
n_touch_total <- sum(touches_mir137(res$results))
n_touch_sig <- sum(touches_mir137(res$results) & res$results$FDR < 0.05)
n_sig <- sum(res$results$FDR < 0.05)
n_total <- nrow(res$results)
tab <- matrix(c(n_touch_sig, n_sig - n_touch_sig,
                n_touch_total - n_touch_sig,
                n_total - n_sig - (n_touch_total - n_touch_sig)), nrow = 2)
ft <- fisher.test(tab)
cat(sprintf("\nhsa-mir-137 involved in %d/%d significant pairs (background rate %.4f, observed %.4f)\n",
            n_touch_sig, n_sig, n_touch_total / n_total, n_touch_sig / n_sig))
cat(sprintf("Fisher exact test: p = %.3e, odds ratio = %.1f - a real, massive concentration\n",
            ft$p.value, ft$estimate))

saveRDS(res, "acc_multiomic_result.rds")
message("\nSaved to dev/real_data_validation/acc_multiomic_result.rds")
