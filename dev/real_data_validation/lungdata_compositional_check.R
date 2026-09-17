# Phase E: real-data validation of the compositional-correction claim.
#
# Goal (per the update plan): find a published dataset with a genuine
# compositional structure, run a compositionally-naive analysis (no
# transform) alongside bicorX's CLR correction, and find one specific,
# checkable result where CLR changes the answer in a way independently
# supported by other evidence - not just "the numbers changed."
#
# Requires: phyloseq (used only to explore candidate datasets during
# development - not used in the final analysis below) and metagenomeSeq
# (for lungData, the dataset actually used). Install via:
#   apt-get install -y r-bioc-phyloseq r-bioc-metagenomeseq
#
# DATA: Charlson et al. lung/airway microbiome study (Smoker vs
# NonSmoker), bundled directly with the metagenomeSeq package as
# `lungData` - a real, published, peer-reviewed dataset, not something
# constructed for this test. 78 samples total; 66 have a valid
# SmokingStatus label (33 Smoker / 33 NonSmoker) after excluding 12
# sterile/blank control samples. This sample size comfortably clears the
# N>=30-per-condition threshold established for ROS-DET's ECOR analytical
# test (Phase A), so this analysis uses significance_method = "analytical"
# rather than permutation testing - a real, adequately-powered comparison,
# not a token qualitative check.
#
# METHOD: OTUs aggregated to genus level (51,891 raw OTUs is far too
# sparse/zero-inflated for correlation-based analysis at n=66; genus-level
# aggregation is both more robust and more biologically interpretable).
# Samples are filtered to the 66 with a valid SmokingStatus label FIRST,
# then genera are filtered to those present in >=15 of those 66 samples
# (168 genera survive) - prevalence must be computed on the samples
# actually used downstream, not on a larger set that includes samples
# later dropped anyway.
#
# RESULT SUMMARY (see NEWS.md for the full writeup):
# - naive: 872 candidate pairs, 20 significant (FDR<0.05)
# - clr:   1054 candidate pairs, 123 significant (FDR<0.05)
# - Only 2 pairs are significant under BOTH methods - the two transforms
#   produce substantially different result sets, as expected for a real
#   compositional dataset.
#
# ONE well-investigated, mechanistically clear case (naive-significant,
# NOT significant under CLR): Acidithiobacillus_thiooxidans vs
# Bacteroides. Naive: r(Smoker)=0.76, r(NonSmoker)=-0.15, FDR=5.5e-06.
# CLR: r(Smoker)=0.28, r(NonSmoker)=-0.31 (same qualitative direction of
# difference, much smaller magnitude - no longer significant).
# Investigated why: Acidithiobacillus's raw counts correlate strongly with
# each sample's own total library size (r=0.88) - a classic sign that its
# apparent abundance is driven by per-sample sequencing depth rather than
# true relative abundance. Two taxa that both scale with a sample's own
# depth will show spurious raw-count correlation with each other
# regardless of real biology, and this effect can differ across
# conditions even without a significant *average* library-size difference
# between groups (confirmed: t-test on library size vs. condition,
# p=0.51 - not itself a group-level confound) simply from within-group
# variance in per-sample depth. CLR's per-sample geometric-mean centering
# directly removes this per-sample scale effect, which is exactly why the
# correlation shrinks under CLR rather than disappearing outright.
#
# A SECOND claim was investigated and DID NOT SURVIVE PROPER TESTING -
# reported honestly rather than discarded quietly. Initial observation:
# 65 of 121 CLR-only-significant pairs (54%) involve Streptococcus or
# Veillonella, the two genera independently identified (via simple
# relative-abundance comparison, before any differential-correlation
# analysis) as showing the largest real compositional shift between
# Smoker and NonSmoker. This looked like "CLR recovers more signal
# involving the biologically relevant taxa" - but a proper enrichment
# test (Fisher's exact, against the correct background rate) shows this
# is NOT an enrichment at all: these two genera already make up 32.7% of
# all 168 genera in the analysis, which mechanically translates to a
# 52.8% background rate of ANY candidate pair touching one of them
# (since a pair touches the group if EITHER of its two genes matches).
# Observed rate among CLR-only hits (53.7%) is statistically
# indistinguishable from this background (p=0.85, odds ratio=1.04). This
# was an artifact of an uncorrected proportion, not a real finding - a
# smaller, Bacteroides-specific follow-up check (2/45 naive-significant
# vs 7/45 CLR-significant pairs) was also examined but the counts are too
# small to support a confident claim either way, and was not pursued
# further as a headline result.

suppressMessages({
  library(bicorX)
  library(metagenomeSeq)
})

data(lungData)

# --- Build the genus-level analysis matrix ---
taxa_str <- as.character(fData(lungData)$taxa)
parts <- strsplit(taxa_str, ";")
genus <- sapply(parts, function(p) if (length(p) >= 6) trimws(p[6]) else NA)

raw_counts <- MRcounts(lungData)
cond_raw <- pData(lungData)$SmokingStatus

valid_taxa <- !is.na(genus) & genus != ""
agg <- rowsum(raw_counts[valid_taxa, ], group = genus[valid_taxa])

keep_samples <- !is.na(cond_raw)
counts <- agg[, keep_samples]
cond <- droplevels(cond_raw[keep_samples])

prevalence <- rowSums(counts > 0)
counts <- counts[prevalence >= 15, ]
rownames(counts) <- make.unique(gsub("[^A-Za-z0-9_]", "_", rownames(counts)))

cat(sprintf("Dataset: %d genera x %d samples (%d Smoker / %d NonSmoker)\n",
            nrow(counts), ncol(counts), sum(cond == "Smoker"), sum(cond == "NonSmoker")))

# --- Run naive vs CLR-corrected ROS-DET, using the validated analytical
#     (ECOR) significance test since N=33/33 clears its calibration
#     threshold from Phase A ---
res_naive <- suppressWarnings(suppressMessages(run_rosdet(
  counts, cond, min_delta = 0.2, transform = "none",
  significance_method = "analytical", workers = 1
)))
res_clr <- suppressWarnings(suppressMessages(run_rosdet(
  counts, cond, min_delta = 0.2, transform = "clr",
  significance_method = "analytical", workers = 1
)))

cat(sprintf("\nnaive: %d candidate pairs, %d significant (FDR<0.05)\n",
            nrow(res_naive$results), sum(res_naive$results$FDR < 0.05)))
cat(sprintf("clr:   %d candidate pairs, %d significant (FDR<0.05)\n",
            nrow(res_clr$results), sum(res_clr$results$FDR < 0.05)))

# --- The well-investigated case: Acidithiobacillus_thiooxidans vs Bacteroides ---
key <- function(df) paste(pmin(df$Gene1, df$Gene2), pmax(df$Gene1, df$Gene2))
naive_sig <- res_naive$results[res_naive$results$FDR < 0.05, ]
clr_sig <- res_clr$results[res_clr$results$FDR < 0.05, ]
cat("\nOverlap between naive-significant and clr-significant pairs:",
    length(intersect(key(naive_sig), key(clr_sig))), "\n")

target <- naive_sig[grepl("Acidithiobacillus", naive_sig$Gene1) & grepl("Bacteroides", naive_sig$Gene2) |
                     grepl("Bacteroides", naive_sig$Gene1) & grepl("Acidithiobacillus", naive_sig$Gene2), ]
cat("\nAcidithiobacillus_thiooxidans vs Bacteroides (naive):\n")
print(target[, c("Gene1", "Gene2", "Distance_Score", "r1", "r2", "FDR")])

in_clr <- key(clr_sig) %in% key(target)
cat("Same pair significant under CLR?", any(in_clr), "\n")

# --- The claim that did NOT survive testing, reproduced honestly ---
clr_only <- clr_sig[!(key(clr_sig) %in% key(naive_sig)), ]
touches_sv <- function(df) grepl("Streptococcus|Veillonella", df$Gene1) | grepl("Streptococcus|Veillonella", df$Gene2)
bg_rate <- mean(touches_sv(res_clr$results))
obs_rate <- mean(touches_sv(clr_only))
n_touch_total <- sum(touches_sv(res_clr$results))
n_touch_clronly <- sum(touches_sv(clr_only))
tab <- matrix(c(n_touch_clronly, nrow(clr_only) - n_touch_clronly,
                n_touch_total - n_touch_clronly,
                nrow(res_clr$results) - nrow(clr_only) - (n_touch_total - n_touch_clronly)), nrow = 2)
ft <- fisher.test(tab)
cat(sprintf("\n[Checked and REJECTED] Streptococcus/Veillonella 'enrichment' among CLR-only hits:\n"))
cat(sprintf("  background rate = %.3f | observed rate = %.3f | Fisher p = %.3f (not significant)\n",
            bg_rate, obs_rate, ft$p.value))

saveRDS(list(res_naive = res_naive, res_clr = res_clr, counts = counts, cond = cond),
        "lungdata_validation_results.rds")
message("\nSaved to dev/real_data_validation/lungdata_validation_results.rds")
