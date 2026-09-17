# bicorX — Updated Plan (Post Paper-Fidelity Review)

**STATUS UPDATE**: All phases A-F described below are now **complete**.
This document is kept as the historical record of the plan as it was
written before execution - see NEWS.md and TODO.md for what actually
happened at each phase, including places where execution diverged from
or refined the plan (e.g. Phase C's range-bias retest needed two
attempts before a properly-controlled result; Phase E used a different
real dataset than originally described here, since the specific
compositional-confound dataset used wasn't identified until execution;
"Phase F" concluded with the package renamed to bicorX, which is why
this file uses that name throughout despite being written when the
package was still called bicorDCEA). Read this file for the original
reasoning and sequencing logic; read NEWS.md for what was actually found
and done.

**Supersedes**: the original `bicorX_pickup_roadmap.md` phases 5-6 sequencing.
Phases 0-4 (re-entry, correctness, structural collapse, surface reduction,
multi-omics) are complete and unaffected by this update. This document
reflects what changed once the three source papers were read against the
implementation.

**The one-sentence change**: the Phase 5 benchmark gap (bicorX trailing
DGCA/dcanr by ~3.5x in ranking quality) was treated as an open statistical
research problem. It is now understood to be, for ROS-DET at least, a
**fixable implementation gap** — the source paper specifies a fast
analytical significance test; the implementation uses a slow permutation
test instead. Four attempts to patch the permutation approach (GPD x2,
normal, skew-normal) all failed calibration testing. None of that was
wasted - it rules out an entire family of fixes - but the actual fix was
in the paper the whole time.

---

## Phase A — ROS-DET: implement ECOR (highest priority, do first)

**Why first**: this is the one item that could change the Phase 5
benchmark verdict directly, and it's grounded in a cited, disclosed
method rather than another empirical curve-fit attempt.

| # | Task | Notes |
|---|---|---|
| A1 | Implement Kayano et al. 2011 Eq. 4 (chi-squared(1) likelihood-ratio test for equality of two correlation coefficients) as an alternative `significance_method` for ROS-DET | Requires numerically solving for the pooled MLE `ρ̂` under H0 (implicit equation, not closed-form) |
| A2 | Apply the same statistic to bicor values (as the paper does), explicitly flagged as an approximation - the paper itself validates this empirically rather than claiming it's exact | Document this caveat prominently, not just in code comments |
| A3 | Offer both modes side by side: `"permutation"` (current, exact, slow) and `"analytical"` (new, approximate, fast) - not a replacement | Preserves the permutation option for cases where the parametric approximation is doubtful (e.g. very non-normal data) |
| A4 | Calibration-test the new analytical mode with the same rigor as everything else this project has done - self-consistency simulation, not just an eyeball check | This is non-negotiable given this project's track record: every "looks fine" parametric shortcut so far has failed under real testing |
| A5 | Benchmark the analytical mode's speed at the paper's stated real-data scale (46 datasets, ~2.6×10¹⁰ pairs) or a scaled-down equivalent feasible in this environment | Confirms the speed claim, not just the statistic |
| A6 | Re-run `dev/benchmarks/01-05` with ECOR available | This is the number that actually answers "does this close the Phase 5 gap" |

**Exit criterion**: ECOR passes calibration testing and a rerun of the
sim102 ground-truth benchmark shows ROS-DET's AUPRC/speed materially
closer to DGCA/dcanr than the current 0.21 vs 0.73 gap.

**If it fails calibration**: document why (same discipline as the GPD
writeup), keep permutation-only, and the Phase 5 verdict stands unchanged.

---

## Phase B — BMKC: close the significance-test gap

**Why this order**: independent of Phase A, real regardless of paper
fidelity, and currently the more clear-cut defensibility problem (a
module-detection method with literally no notion of "could this be
chance" is a harder sell than a slow-but-correct significance test).

| # | Task | Notes |
|---|---|---|
| B1 | Implement Yuan et al. 2015 §3.3's global permutation test as a first pass: shuffle each gene's expression independently across samples (destroys all covariance, not just condition labels), rerun full module-detection pipeline, repeat ~1000x, compare aggregate score | Matches the paper; gives one p-value for the whole module set, not per-module |
| B2 | Evaluate whether a per-module (not just global) significance measure would be more useful in practice, and whether it's worth designing one beyond what the paper does | This is where "modernize, don't just replicate" is legitimate - the paper's coarse global test may not be what users actually want |
| B3 | Calibration-test whichever version ships | Same discipline as everywhere else |

**Exit criterion**: `run_bmkc()` returns some defensible measure of
significance alongside modules, not modules alone.

---

## Phase C — Re-examine `weight_mode` on the scenario it was actually built for

**Why this matters**: the current default (`"unweighted"`) was set after
finding `"wcor"` inverts rankings under a *compositional-closure* confound
— a scenario the paper never claims to address. The paper validates WCOR
specifically against *range bias* (Fig. 4b/c), and shows it clearly
outperforming unweighted bicor there.

| # | Task | Notes |
|---|---|---|
| C1 | Build a range-bias test scenario matching the paper's own synthetic design (one condition's expression range much smaller than the other's) | Distinct from the compositional-confound scenario already tested |
| C2 | Compare `"wcor"` vs `"unweighted"` on that scenario specifically | If `"wcor"` wins here as the paper predicts, this confirms it's scenario-dependent, not simply "worse" |
| C3 | Document both scenarios and their outcomes clearly, and consider whether the *default* should depend on detectable data characteristics rather than being a single fixed choice | Avoids quietly picking a "wrong for some real use cases" default a second time |

**Exit criterion**: a documented, evidence-based statement of when to use
`"wcor"` vs `"unweighted"`, not a single default presented as universally better.

---

## Phase D — Documentation and citation (do alongside A-C, finalize after)

| # | Task | Notes |
|---|---|---|
| D1 | Cite all three papers (Zheng et al. 2014, Kayano et al. 2011, Yuan et al. 2015) in `DESCRIPTION`/`CITATION`/vignettes | Currently missing entirely |
| D2 | Publish an honest per-method fidelity table in the package docs (condensed version of the tables from this conversation) | Faithful / modernized / gap, per finding, per method - not a blanket "implements X" claim |
| D3 | Document BMHT's missing maximum-clique stage explicitly as current scope (ranking only), pointing to BMKC for module-level analysis | Lower priority than A-C but cheap and removes a real ambiguity |
| D4 | Document the k-core (BMKC) vs k-clique-percolation (paper) distinction explicitly - different mathematical objects, different guarantees | Rename any code/docs that currently imply the paper's exact structure |
| D5 | Document BMHT's undocumented `mad_1`/`mad_2` normalization as a package addition beyond the published method | Small but should not be silently present |

---

## Phase E — Real-data validation (sequenced after A, not before)

Unchanged in substance from earlier discussion, but now deliberately
sequenced *after* Phase A rather than done in parallel, since ECOR may
change what "the compositional/multi-omics advantage" even looks like for
ROS-DET specifically.

| # | Task |
|---|---|
| E1 | Find a published microbiome/host-microbe dataset with a documented compositional-confound case |
| E2 | Run a compositionally-naive tool alongside bicorX (with whichever significance method - permutation or ECOR - proves best from Phase A) |
| E3 | Show one specific, checkable result where CLR correction changes the answer in a way independently supported by other evidence |

---

## Phase F — Naming / positioning (last, not first)

Deliberately deferred until A-C have real outcomes. A name chosen now
would be guessing; a name chosen after Phase A's benchmark rerun and
Phase B's significance-test fix would be informed by what the package
actually, verifiably does. Candidate framing once ready: *"modernized
implementations of three bicor-based differential coexpression methods
(cite all three), extended for compositional and multi-omics data, with
explicit stated departures from the original methods where scale or the
new data regime required them."*

---

## Summary priority order

1. **Phase A** (ROS-DET/ECOR) — highest expected impact, directly
   testable against the existing Phase 5 benchmark harness.
2. **Phase B** (BMKC significance test) — real gap, independent of A.
3. **Phase C** (`weight_mode` range-bias retest) — cheap, resolves an
   open question from the last default change.
4. **Phase D** (documentation/citation) — do continuously alongside A-C,
   finalize once A-C land.
5. **Phase E** (real-data validation) — after A, since A may change what's
   being validated.
6. **Phase F** (naming) — last, informed by A-E rather than guessed at.

## What does NOT change

- Everything already shipped and validated (Phase 1 correctness fixes,
  Phase 2 dispatch collapse, Phase 3 plot-surface reduction, Phase 4
  multi-omics support, the BLAS/buffer-reuse/RAM performance work) stays
  as is - none of it is paper-fidelity-related and none of it is called
  into question by this review.
- The four rejected permutation-floor fixes (GPD-MoM, GPD-MLE, normal,
  skew-normal) remain rejected regardless of Phase A's outcome for
  BMHT specifically, since BMHT's own paper does not specify an
  analytical alternative - only ROS-DET's does. BMHT's permutation-floor
  problem, if the benchmark still shows one after Phase A, remains open
  and would need Phase 1 (Post-Phase-5)'s original recommendation: either
  a BMHT-specific analytically-derived null (harder, since edges aren't
  independent - noted previously) or simply higher `n_permutations`,
  now cheaper thanks to the performance work already shipped.
