# Phase 5 Benchmark Report: bicorX vs. DGCA vs. dcanr

**Purpose**: per the pickup roadmap, this is the decision gate - "if this
benchmark produces nothing, stop before further polish." It produced a
real, specific, actionable finding, but not the clean win the roadmap
was hoping for. Reported here exactly as found, including two bugs
discovered in the benchmark harness itself along the way (both disclosed
and fixed before drawing conclusions) and the parts of the evidence that
argue against the package as currently implemented.

**Tools compared**: `DGCA` (Zhang & McKenzie, CRAN) and `dcanr` (Bhuva et
al., Bioconductor) - neither was available via this sandbox's apt mirror;
both were installed from their GitHub source (`andymckenzie/DGCA`,
`DavisLaboratory/dcanr`).

**Environment**: single core, R 4.3.3. Five separate experiments, scripts
`01`-`05` in this directory, described in the order they inform the
verdict rather than the order they were run.

---

## Finding 1: Null calibration - no difference found

`01_calibration_comparison.R`. On pure-noise data (80 genes, 24 samples,
zero true signal, 10 replicates), per-gene false-positive rate at
FDR < 0.05:

| Tool | Mean FPR |
|---|---|
| bicorX (BMHT) | 0.0000 |
| DGCA | 0.0025 |
| dcanr | 0.0050 |

All three control Type I error well below nominal. This mainly confirms
Phase 1's BMHT calibration fix is holding and that bicorX isn't
worse-calibrated than either competitor - it doesn't distinguish them.

## Finding 2: Precision under sample corruption - a real, demonstrated advantage

`02_outlier_robustness.R`. 15 genes given a genuine shared-latent-factor
correlation structure in one condition only; 0, 2, or 4 of 24 samples
then corrupted with 15x-magnitude outliers restricted to those 15 genes
(the worst case for a non-robust correlation measure - corruption lands
exactly on the signal genes).

| Scenario | bicorX precision | DGCA precision | dcanr precision |
|---|---|---|---|
| Clean | 0.00 (no calls) | 0.68 | 0.48 |
| 2/24 corrupted (~8%) | **0.92** | 1.00 | 0.45 |
| 4/24 corrupted (~17%) | **1.00** | **0.38** | 0.47 |

DGCA's precision collapses under corruption (0.68 -> 0.38); at 17%
corruption, over 60% of its calls are false positives. bicorX's
precision holds or improves under the same corruption. This is the
clearest, most mechanistically-explained result in this report: biweight
midcorrelation's outlier-downweighting does what it's supposed to do,
specifically when outliers are present. See
`precision_under_corruption.png`.

## Finding 3: Recall/power on the same synthetic data - a real weakness, with a specific cause

Same experiment, recall (of the 15 true genes):

| Scenario | bicorX recall | DGCA recall | dcanr recall |
|---|---|---|---|
| Clean | **0.00** | 1.00 | 1.00 |
| 2/24 corrupted | 0.80 | 1.00 | 1.00 |
| 4/24 corrupted | 0.87 | 1.00 | 1.00 |

On clean data with an unambiguous signal, bicorX found **zero**
significant genes at a typical permutation count (500), while both
competitors found everything. Traced to a specific cause: the
BH-adjusted FDR floor at `n_permutations=500` across 60 genes was 0.0599
- a hair's-width miss from the permutation-count floor
(`1/(n_permutations+1)`), not absence of signal. Recall improves with
more permutations but plateaus well below 1.0:

| n_permutations | Genes recovered (of 15) |
|---|---|
| 500 | 0 |
| 2,000 | 6 |
| 5,000 | 7 |
| 10,000 | 8 |

**Root cause**: DGCA and dcanr use (semi-)analytical null distributions
with continuous p-value resolution; bicorX's BMHT/ROS-DET use pure
permutation counting, which needs substantially more permutations - and
runtime - to reach comparable resolution once corrected across many
tests.

## Finding 4: A compositional confound - CLR helps, but doesn't close the gap

`04_compositional_confound.R`. 60 taxa, 24 samples: taxa 1-10 form a true
differential module; a 60th "dominant" taxon shifts sharply between
conditions, which after compositional closure (fixed total reads per
sample) mechanically induces spurious differential correlation among the
49 unrelated background taxa - exactly the artifact CLR exists to
correct. AUPRC (true module vs. background, ranked by each tool's raw
score - DGCA/dcanr have no CLR option and always see the same closed
counts regardless of what bicorX is given):

| | bicorX (no CLR) | bicorX (CLR) | DGCA | dcanr |
|---|---|---|---|---|
| AUPRC | 0.359 | **0.550** | 0.807 | 0.807 |

CLR gives bicorX a real, meaningful improvement on its own terms
(0.359 -> 0.550) - the transform is doing what it's supposed to do. But
even with CLR, bicorX's separation of true signal from the
closure-induced confound remains well behind DGCA/dcanr's performance on
the same *uncorrected* data. This is consistent with Finding 3: it looks
like a general ranking/power gap that shows up here too, not a
compositional-specific failure - CLR helps, but isn't enough on its own
to close a gap that exists for other reasons.

## Finding 5: A real published benchmark with known ground truth - the most important, least comfortable result

`03_sim102_ground_truth.R`, using dcanr's own bundled simulated benchmark
(`sim102`) with a genuine ground-truth differential network - the
strongest evidence in this report, since it's external and not built by
us to make a point either way.

**Two bugs were found and fixed in this benchmark's own harness before
trusting its output**, both disclosed here because they impeach earlier
draft numbers that should not be trusted (a pre-fix version had ROS-DET
scoring *below random* - that number was wrong, and is not reported
further):

1. `getTrueNetwork()` returns an **asymmetric** matrix - 20 entries in the
   upper triangle and 165 *different* entries in the lower triangle, not
   mirrored duplicates. An earlier draft checked only the upper triangle,
   undercounting true edges 9x (20 vs. the correct 185).
2. `combn(rownames(emat), 2)` does not produce alphabetically-ordered
   pairs (34% of pairs have `Gene1 > Gene2`). An earlier draft's ROS-DET
   key was canonicalized alphabetically while the ground-truth key
   wasn't, silently turning real ROS-DET scores into `NA` for roughly a
   third of all pairs - unrelated to Gate 1 coverage, which is actually
   complete (`min_delta=0` correctly scores all 11,026 candidate pairs
   once the key bug is fixed).

**Corrected result** (149 genes, 406 samples, 185 true edges / 11,026
candidate pairs):

| Tool | AUPRC | Wall time |
|---|---|---|
| bicorX (ROS-DET) | **0.211** | 1.28s |
| DGCA | **0.732** | 0.03s |
| dcanr | **0.732** | 0.00s |

Random baseline AUPRC is 0.017, so ROS-DET is well above chance (~12x) -
but DGCA and dcanr are both roughly 3.5x better than ROS-DET, and around
40x faster. This is the single most important number in this report: on
external, real ground truth, bicorX's ranking quality trails both
competitors substantially, consistent with Findings 3-4 rather than
contradicting them.

## Finding 6: Real biological data, no ground truth - runtime advantage, and genuine disagreement between methods

`05_darmanis_real_data.R`, using the Darmanis et al. 2015 human brain
single-cell dataset (572 genes, 158 cells, neuron vs. oligodendrocyte)
bundled with DGCA itself - real dropout and technical noise, no synthetic
ground truth available, so this checks concordance and speed rather than
accuracy.

| | Value |
|---|---|
| DGCA runtime | 16.6s |
| ROS-DET runtime | **1.3s** |
| BMHT runtime | **1.4s** |
| Spearman rank correlation (DGCA \|zScoreDiff\| vs. ROS-DET Distance_Score, 71,072 shared pairs) | 0.19 |
| "Significant" pairs at raw p<0.05: DGCA / ROS-DET | 22,011 / 21,790 |
| Overlap (Jaccard) | 0.16 |

bicorX is **~12x faster** than DGCA on real data at this scale - a
genuine, substantial advantage the synthetic experiments above didn't
fully surface (DGCA's permutation loop in pure R is slow; bicorX's C++
kernels are fast, consistent with Finding 3's per-permutation-cost
observation). The two methods call similarly-sized "significant" sets
(~22k each) but mostly *different* pairs (Jaccard 0.16, weak rank
correlation) - they are not interchangeable, and which one is "right" on
real data with no ground truth available can't be settled by this
experiment alone.

---

## Verdict

Not a clean "ship it," not a clean "abandon it," and the picture is more
serious than an earlier draft of this report suggested before Findings
4-6 (from work already in progress that hadn't been fully reconciled
yet) were incorporated.

**What's validated**: the core idea - robust correlation as a defense
against outlier-driven false positives - is real and replicated
(Finding 2). bicorX's compiled engine is also genuinely fast per unit
of work (Findings 3 and 6).

**What's a real, current weakness, not a footnote**: on the one
benchmark with actual external ground truth (Finding 5), bicorX
trails both competitors by roughly 3.5x in ranking quality and 40x in
speed. That gap shows up consistently across the compositional-confound
test (Finding 4) and the clean-data recall test (Finding 3) too, all
pointing at the same root cause: bicorX's pure permutation-testing
machinery has a coarse resolution floor that (semi-)analytical
competitors don't share, and this is costing real detection power, not
just requiring "a few more permutations."

**Recommended next step, if this package continues**: fix the
permutation-resolution problem before investing further in packaging
polish. A semi-parametric approximation to the permutation null
(extrapolating tail p-values below the raw permutation floor, in the
spirit of what DGCA already does) is a well-scoped, specific piece of
engineering - not a redesign of the underlying method, which Finding 2
suggests is worth keeping.

## Limitations of this benchmark

- Findings 1-4 are synthetic, with a ground truth we constructed
  ourselves - directional evidence, not definitive. Findings 5-6 use
  external data (a published simulation with ground truth, and a
  published real dataset without one) and are the more trustworthy
  results in this report for that reason.
- BMHT (a per-gene hub-detection statistic) was used for Findings 2-3
  rather than ROS-DET (a per-pair test, more directly analogous to what
  DGCA/dcanr compute) for reasons of time. Finding 5 does use ROS-DET and
  shows the same qualitative pattern, which is reassuring but doesn't
  fully substitute for redoing 2-3 with ROS-DET directly.
- Only one signal shape was tested in Findings 2-4 (a symmetric,
  densely-correlated gene block). BMHT's hub-centric design may be
  structurally advantaged or disadvantaged by other signal shapes (e.g. a
  true hub gene with many weakly-changed edges) not explored here.
- Two bugs were found in this benchmark's own harness (Finding 5) and
  fixed before reporting - a reminder that benchmark code deserves the
  same scrutiny as the code being benchmarked, and a caution that further
  undiscovered harness issues remain possible.
