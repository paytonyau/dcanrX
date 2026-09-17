# NEWS

## Unreleased

### rOpenSci preparation: first `R CMD check`, and several real bugs it caught

Ran `R CMD check` for the first time in this project, in preparation for
a possible rOpenSci submission. It went from **1 ERROR, 5 WARNINGs, 4
NOTEs** to **0 ERRORs, 1 WARNING, 3 NOTEs** — the single remaining
WARNING being a sandbox locale limitation (`cannot set locale to
en_US.UTF-8`), not a code issue, which will not appear on a normal
machine. The full test suite (44 tests) continues to pass with zero
failures throughout.

**Real bugs caught, several of them user-facing:**

- **Hard `importFrom` directives on four `Suggests`-only packages**
  (`SummarizedExperiment`, `MultiAssayExperiment`, `ggplot2`,
  `scuttle`). This would have caused the package to **fail to load
  entirely** for any user without all four installed. The affected code
  already used correct `requireNamespace()` guards, so these imports
  were both redundant and actively harmful. Introduced accidentally by
  `roxygen2` during the earlier documentation pass; removed.
- **Three internal engines were being exported** (`.run_bmht_matrix`,
  `.run_rosdet_matrix`, `.run_bmkc_matrix`) — dot-prefixed
  implementation details leaking into the public API. Verified they are
  called only internally (via `.run_engine()`, and via `bicorX:::` in
  tests) before marking them `@keywords internal`.
- **A non-ASCII character** (a literal Greek delta in a `ggraph` legend
  label) violating R's portability requirement for R source files;
  replaced with the `\u0394` escape.
- **Undeclared dependencies**: `ggrepel` and `visNetwork` were used via
  `::` but absent from `Suggests`; `parallel`, `rlang`, and several
  `stats`/`utils` functions (`head`, `write.csv`, `reorder`, `na.omit`,
  `.data`) were used without being imported.
- **An unnecessary `magrittr` dependency** was avoided by rewriting the
  three chained `%>%` calls in the interactive-network plotting path as
  explicit nested calls, rather than adding a package for one code path.

**Documentation completed**: roughly 27 exported functions had
undocumented arguments and two (`bicor`, `plot_module_heatmap`) had no
documentation at all. All now have `@param` and `@return` tags. One
further unescaped `%` was caught in the newly written text (Rd requires
`\%`) — the same class of mistake found in the previous documentation
pass.

**Metadata and community files added** for rOpenSci compliance:
`Authors@R` (Payton Yau, aut/cre), `URL`, `BugReports`, a shortened
`Title`, `CONTRIBUTING.md`, and `CODE_OF_CONDUCT.md` (Contributor
Covenant 2.1). The contributing guide includes a specific requirement
that any change affecting p-values, significance testing, or null
distributions must come with a calibration check, pointing at the
existing calibration tests as the pattern to follow — a direct
consequence of this project's history of methods that looked correct but
failed calibration testing.

**Known remaining items** (documented rather than silently left):
- `REPLACE.ME@example.com` is a placeholder in `DESCRIPTION` and
  `CODE_OF_CONDUCT.md` and must be replaced before any submission. A
  syntactically valid address was required because an obviously-bracketed
  placeholder breaks `R CMD build` outright.
- Two cosmetic NOTEs remain: installed package size, and `NEWS.md`
  section titles that R cannot parse as version numbers (a consequence of
  this file's development-log style).
- The package name has still only been checked for collisions against
  GitHub and a partial Debian mirror, not the real CRAN/Bioconductor
  namespaces.

### Test suite
16 files, 44 tests, all passing. 3 vignettes and full function-level help
now build and install correctly.

---
### Documentation: README, a new getting-started vignette, and a real gap found in how vignettes/help pages were being verified

Addressed the standing documentation gap flagged repeatedly in prior
`TODO.md` entries: the substantial body of validated findings (Phase A's
ECOR benchmark, Phase E's real microbiome result, the TCGA multi-omics
result, performance/RAM work) was only discoverable by reading `NEWS.md`
or git history, with no proper entry point for a new user.

- **Added `README.md`** (there was none before this) - package overview,
  installation, a quick-start example, a data-types table, a
  significance-method decision table (permutation vs. analytical), a
  "what's actually been validated" section linking each headline claim to
  its reproducible evidence, and an explicit "honest scope" section
  listing real, current limitations rather than only strengths.
- **Added `vignettes/getting-started.Rmd`** - a reference vignette going
  deeper than the README on input types, transform choice (with the real
  library-size-artifact mechanism from Phase E as the motivating example,
  not an abstract description), and when to use each ROS-DET significance
  method, with every specific number linked to the script that produced
  it.

**A real gap was found and fixed while doing this, worth reporting
plainly rather than glossing over**: `vignette(package = "bicorX")`
returned "no vignettes found" even though every `.Rmd` in this project
had been individually checked with `rmarkdown::render()` at the time
each was added (Phases 4, E, and the multi-omics example all did this).
That check was real but incomplete - it confirms a `.Rmd` file is
individually renderable, not that it gets built into the installed
package and made discoverable via `vignette()`. The actual cause: this
whole project's `R CMD INSTALL` workflow installed directly from the
source directory, which does not build vignettes into the package at
all in this environment - `R CMD build` (producing a proper source
tarball, which *does* build vignettes as part of that step) followed by
installing from the resulting tarball is required. Also discovered in
the same pass: 44 function-level help pages (`?run_bmht` etc.) had never
existed - `NAMESPACE` had been hand-written since Phase 0 specifically
because `roxygen2` wasn't available in this environment at the time; it
is now, and running it confirmed the hand-written `NAMESPACE`'s
*exports* were completely accurate (byte-identical on that front), while
also surfacing genuinely missing `importFrom` declarations for several
Suggested-package dependencies that had been working by accident rather
than by correct declaration.

**A second, smaller issue caught in the same verification pass**: the
auto-generated help page for the internal, already-rejected
`.gpd_tail_pvalues()` function produced Rd parser warnings, traced to
two literal, unescaped `%` characters in its documentation prose (Rd
format requires `\%` for a literal percent sign - an easy, silent
mistake, since ordinary prose text has no reason to expect that). Fixed
by escaping both instances and confirming a subsequent clean rebuild
produced zero warnings.

**Verified, in order**: `NAMESPACE` diffed against the hand-written
version before applying (exports identical, only additions); a full
clean rebuild via `R CMD build` + install-from-tarball (not direct
source-directory install, per the finding above) produces zero warnings;
`vignette(package = "bicorX")` now correctly lists all three vignettes;
`help("run_bmht")` now resolves; the full test suite passes (unaffected,
since none of this touched engine code); `dev/phase0_baseline.R` gives
statistically identical results to every prior run.

### Test suite
16 files, 80+ assertions, all passing.

---
### A real, genuinely cross-layer multi-omics example (answers "do we have an omics example?")

The existing vignette (`vignettes/multiomics-workflow.Rmd`) uses a small
synthetic RNA+ATAC dataset for clarity of exposition. This adds the real,
externally-published complement: RNA-seq + miRNA-seq from TCGA
Adrenocortical Carcinoma, comparing the published C1A/C1B molecular
subtypes (Zheng et al. 2016, Cancer Cell) - genuinely two distinct omics
platforms, not two feature panels of the same assay, and adequately
powered (78 patients, 43 C1A / 35 C1B - clears the N>=30 threshold for
ROS-DET's ECOR analytical test).

**Data**: `miniACC`, bundled directly with the Bioconductor
`MultiAssayExperiment` package - no ExperimentHub/network access needed,
unlike the `scRNAseq` datasets that blocked a fuller single-cell
validation earlier. A real, if reduced, TCGA multi-omics dataset.

**Result**: `run_multiomic(rna, mirna, cond, method="rosdet",
pair_type_filter="inter", significance_method="analytical")` finds 10
genome-wide-significant (FDR<0.05) cross-layer pairs out of 18,706
candidates. A striking pattern: a single microRNA, `hsa-mir-137`, is
involved in 6 of those 10 - each showing the same qualitative signature
(strong positive correlation with a different mRNA gene in one subtype,
near-zero or weak correlation in the other). Checked formally rather than
reported as an eyeballed proportion - a lesson carried directly from
Phase E, where a similarly-appealing-looking pattern turned out to be a
background-rate artifact. Here it isn't: this miRNA appears in only 71 of
18,706 total candidate pairs (a 0.38% background rate), and a Fisher's
exact test against that correct background gives p = 5.0e-13, odds ratio
= 429 - unambiguously real, not a coincidence of a large background rate.

**What is and isn't claimed**: the statistical concentration itself and
the internal consistency of the pattern across all six target genes are
verified computationally within the script. Any specific claim about
`hsa-mir-137`'s documented biological role in adrenocortical carcinoma,
or the precise clinical meaning of the C1A/C1B subtypes beyond what's in
this dataset's own metadata, is NOT independently verified here - this
environment has no web search access to check the original TCGA paper.
The script notes `hsa-mir-137` is recognized as a tumor-suppressor
microRNA in other cancer types from general background knowledge, and
that one target gene (`MYC`) is a well-known major oncogene - both
flagged explicitly as plausible leads for follow-up, not confirmed
claims, exactly the same discipline applied to the mir-137 finding as
everywhere else in this project.

Reproducible in `dev/real_data_validation/acc_multiomics_check.R`, with
the verified-vs-unverified distinction stated directly in the script
header. Referenced from `vignettes/multiomics-workflow.Rmd`'s new "A
real example" section (which also renders correctly after this addition)
alongside a note that `significance_method = "analytical"` sidesteps the
permutation-count tradeoff the vignette already discusses, given adequate
sample size.

### Test suite
16 files, 80+ assertions, all passing (dev/ + vignette addition only, no
engine code touched). `MultiAssayExperiment` was already a declared
`Suggests` dependency from Phase 4.

---
### Package renamed: bicorDCEA -> bicorX

Following the paper-fidelity review and the real-data validation work
(Phases A-E, plus the single-cell pathway check), the package has been
renamed from `bicorDCEA` to `bicorX`. Rationale (see the naming analysis
this rename followed): keeps the recognizable `bicor` root already used
throughout the citation and documentation work, the "X" signals the
actual validated differentiator (cross-omics/extended scope - Phase E's
real-data compositional-correction result), and avoids overclaiming
words like "robust" or "power" that would need constant re-justifying
against DGCA/dcanr on the general case (Phase A showed *competitive*,
not *dominant*, and only when using the appropriate significance
method). Checked for collisions against GitHub (a working search,
verified against a known package as a control - zero hits) and Debian's
CRAN/Bioconductor package mirror (zero hits, though this covers only a
subset of the true CRAN/Bioconductor namespace - a real check against
CRAN/Bioconductor's own servers was not possible from this environment
due to network restrictions, and should be done before any public
release).

**What changed**: `DESCRIPTION`'s `Package:` field, `NAMESPACE`'s
`useDynLib` directive, the compiled shared library and all its Rcpp-
generated native symbol names (regenerated via
`Rcpp::compileAttributes()` rather than hand-edited, since these are
tied to the package name via a fixed naming convention -
`_bicorDCEA_cpp_...` became `_bicorX_cpp_...` throughout both
`R/RcppExports.R` and `src/RcppExports.cpp`, including the
`R_init_<package>` registration function), the package-level
documentation file (`R/bicorDCEA.R` -> `R/bicorX.R`), `inst/CITATION`,
and every `library(bicorDCEA)`/`bicorDCEA:::` reference across the test
suite, vignettes, and `dev/` scripts (case-sensitive replacement,
verified not to touch the unrelated lowercase `run_bicordcea()` function
name, which is unchanged).

**What did NOT change**: no public function was renamed (`run_bicordcea()`,
`run_bmht()`, `run_rosdet()`, `run_bmkc()`, etc. are all unchanged) - this
is a package-identity rename only, not an API change.

**Verification, in order**: removed the old installed package entirely
(rather than leaving both installed, to avoid any ambiguity about which
was actually being tested); rebuilt from source and confirmed the
compiled library installs and loads correctly under the new name;
confirmed a real function call (`run_bmht()`, `bicor()`) executes
correctly; ran the full test suite under the new package name (44
tests, 0 failed, 0 errors); re-ran `dev/phase0_baseline.R` and confirmed
byte-identical results to every prior run under the old name (same 20
BMHT-significant genes, 0 BMKC modules, 0 ROS-DET-significant pairs);
confirmed both vignettes still render; confirmed `citation("bicorX")`
renders correctly with the updated package name throughout.

Historical entries elsewhere in this file and in `TODO.md` still refer
to `bicorDCEA` where that was the accurate name of the package at the
time the work was done - left as an honest historical record rather than
rewritten, consistent with this project's practice throughout.

### Test suite
16 files, 44 tests, 80+ assertions, all passing under the new package
name.

---
### Single-cell pathway verification (real data, but explicitly NOT the same claim as Phase E)

Following up on the request for a single-cell-level example. Important
distinction, stated plainly rather than blurred: this verifies the
single-cell/pseudobulk code pathway runs correctly on genuine (not
synthetic) single-cell data - it is **not** an adequately-powered
scientific finding on the same footing as Phase E's lung microbiome
result, and the writeup below says so explicitly rather than presenting
it as equivalent.

**What was tried first and didn't work**: Bioconductor's `scRNAseq`
package (a curated collection of published single-cell datasets,
including many real case-control cohorts with proper multi-patient
replicate structure - exactly what would be needed for a Phase-E-quality
single-cell validation) was installed and its dataset catalog loads
correctly, but every dataset is fetched through ExperimentHub, which
requires external network access this environment does not have -
confirmed by attempting a real fetch (`ZeiselBrainData()`), which fails
with "Cannot connect to ExperimentHub server" even with local-hub
fallback. This is a real environment limitation, not a choice to skip a
better available option.

**What was used instead**: `pbmc_small`, bundled directly with the
Seurat package (a heavily downsampled 230-gene, 80-cell subset of a
classic 10x Genomics PBMC dataset, shipped for Seurat's own
testing/demonstration purposes) - real single-cell count data, but from
what is effectively one specimen, with no multi-patient structure to
aggregate over. `make_pseudobulk()` is designed to aggregate cells from
the same biological sample across *multiple* independent samples; there
are no multiple real patients here to use that way. To exercise the
function at all, cells within the two largest cell-identity clusters
(real, distinct immune cell populations - not an experimental condition)
were split into artificial pseudo-replicate pools of 5 cells each,
giving 7 pseudobulk samples for one cluster and 5 for the other - just
barely clearing every engine's hard minimum (5 per condition), far below
any threshold this project has established for calibrated inference.

**What this confirms**: the full pipeline (`SingleCellExperiment`
construction -> `make_pseudobulk()` -> `run_bmht()`, both
`transform="none"` and `transform="clr"`) runs correctly end-to-end on
real single-cell-derived data with no errors. The top-ranked gene
(`HLA-DQB1`, an MHC class II gene) is at least biologically plausible as
a spot-check - MHC-II expression is known to differ sharply between
antigen-presenting and non-antigen-presenting immune cell types, which
is exactly the two-cluster comparison used here - but this is explicitly
a plausibility spot-check, not a statistical claim, given the sample
size.

Reproducible in `dev/real_data_validation/singlecell_pathway_check.R`,
with the scope limitation documented directly in the script's header, not
just in NEWS.md.

**If real network access to Bioconductor's ExperimentHub becomes
available later**, replacing the pseudo-replicate construction with an
actual multi-patient dataset from `scRNAseq::listDatasets()` would give a
genuine, adequately-powered single-cell validation on the same footing
as Phase E - this is flagged directly in the script for whoever picks
this up next with better network access.

### Test suite
16 files, 80+ assertions, all passing (dev/ validation only, no engine
code touched; `scuttle` was already a declared Suggests dependency from
Phase 4, just not previously installed in this environment).

---
### Phase E: real-data validation of the compositional-correction claim

Found a real, published, adequately-powered dataset for this rather than
relying only on synthetic benchmarks: the Charlson et al. lung/airway
microbiome study (Smoker vs. NonSmoker), bundled directly with the
Bioconductor `metagenomeSeq` package as `lungData` - 66 samples with a
valid condition label (33/33), comfortably clearing the N>=50-recommended
(and >=30-minimum) threshold established for ROS-DET's ECOR analytical
test in Phase A, so this uses `significance_method = "analytical"` for a
real, adequately-powered comparison rather than a token qualitative
check.

OTUs aggregated to genus level (168 genera survive a basic prevalence
filter) and ROS-DET run twice - `transform = "none"` (naive) vs.
`transform = "clr"` - otherwise identical settings. The two transforms
produce substantially different result sets (872 candidate pairs / 20
significant vs. 1054 / 123; only 2 pairs significant under both), as
expected for real compositional data.

**One well-investigated, mechanistically clear case**:
`Acidithiobacillus_thiooxidans` vs. `Bacteroides` is highly significant
under naive analysis (FDR = 5.5e-06, r changing from -0.42 to 0.76
between conditions) but not significant at all under CLR. Investigated
why rather than taking the result at face value: `Acidithiobacillus`'s
raw counts correlate strongly with each sample's own total library size
(r=0.88) - a classic sign its apparent abundance tracks sequencing depth
rather than true relative abundance. Confirmed this isn't simply a
group-level confound (mean library size doesn't differ significantly
between Smoker/NonSmoker, t-test p=0.51) - the effect is a per-sample
scaling artifact (any two taxa that both scale with a sample's own
depth will show spurious raw-count correlation regardless of biology),
which CLR's per-sample geometric-mean centering directly removes. Under
CLR the correlation shrinks (0.28 vs. -0.31) but keeps the same
qualitative direction - consistent with a real but overstated naive
signal, not a mechanism-free coincidence.

**A second observation was investigated and did NOT survive proper
testing - reported honestly rather than discarded quietly.** 65 of 121
CLR-only-significant pairs (54%) involve *Streptococcus* or
*Veillonella*, the two genera independently identified (before any
differential-correlation analysis, via simple relative-abundance
comparison) as showing the largest real compositional shift between
conditions - consistent with published smoking-microbiome literature.
This looked like "CLR recovers more signal involving the biologically
relevant taxa." A proper enrichment test (Fisher's exact, against the
correct background rate) shows it isn't: these two genera already make
up 32.7% of all 168 genera, mechanically producing a 52.8% background
rate for any pair touching either of them (a pair touches the group if
*either* gene matches) - the observed 53.7% rate is statistically
indistinguishable from chance (p=0.85). This was an artifact of an
uncorrected proportion, not a real finding.

Both the positive and negative findings are fully reproducible in
`dev/real_data_validation/lungdata_compositional_check.R`, including the
enrichment test that rejected the second claim - not just narrated in
comments.

**This is the first validation in this project's history to use real,
externally-published data for the compositional/multi-omics claim
specifically** (Phase 5's benchmarks used DGCA/dcanr's own real
ground-truth simulation for the general competitiveness question, but
the compositional-correction claim itself had only synthetic-data
support until now).

### Test suite
16 files, 80+ assertions, all passing (this is a `dev/` validation
script, not a package test - no engine code touched, full suite
unaffected).

---
### Follow-up: does BMHT show its own version of the Phase 5 gap? Checked directly - no, and a related finding about ranking

Zheng et al.'s BMHT paper does not specify an analytical alternative to
permutation testing (unlike Kayano et al.'s ROS-DET, fixed in Phase A),
so this was flagged as an open question rather than assumed either way.
Checked directly against the same sim102 ground truth used throughout
the paper-fidelity review.

BMHT is a per-gene statistic, so this needed its own gene-level ground
truth (a gene counts as "true" if it participates in >= 1 true
differential edge - 42/149 genes here, a ~28% baseline rate, not directly
comparable in absolute AUPRC terms to the pair-level ~1.7% baseline used
for ROS-DET/DGCA/dcanr above).

**Result: BMHT's raw `DC_Score` ranking (AUPRC 0.8895) is close to
DGCA's gene-level equivalent (AUPRC 0.9277, rolling up each gene's max
|z-score| across its pairs)** - a small, single-digit-percent relative
gap, nothing like ROS-DET's original 3.5x pair-level gap. **BMHT does
not appear to need the same kind of fix ROS-DET did.**

**A separate, genuinely interesting finding surfaced while checking
this**: ranking by `P_value` (AUPRC 0.7931) is *worse* than ranking by
raw `DC_Score` for BMHT - the opposite of what happened for ROS-DET,
where p-value-based ranking (via ECOR) was dramatically better than the
raw score. Likely explanation (not fully mechanistically isolated - would
need a dedicated follow-up to confirm precisely): `P_value` normalizes
each gene's score by that gene's own permutation-null variance, which
differs across genes depending on how many "informative" edges
contribute to it; on this dataset that per-gene normalization reorders
genes slightly worse than the unnormalized raw score does.

**Practical recommendation, now documented on `.run_bmht_matrix()`**:
rank/prioritize BMHT candidates by `DC_Score`, not `P_value`/`FDR` - but
still use `FDR` for a calibrated significance cutoff decision, since its
calibration (established in Phase 1) is a separate property from ranking
quality and is unaffected by this finding.

Made permanent in `dev/benchmarks/03_sim102_ground_truth.R` (added a
gene-level BMHT-vs-DGCA comparison section alongside the existing
pair-level ROS-DET comparison) rather than left as an ad hoc check.

### Test suite
16 files, 80+ assertions, all passing as of this entry (documentation-only
change, no engine code touched, all tests unaffected as expected).

---
### Paper-fidelity update, Phase D: citations and honest scope in DESCRIPTION

- Added `inst/CITATION`, so `citation("bicorDCEA")` returns all three
  source papers (Zheng et al. 2014, Kayano et al. 2011, Yuan et al. 2015)
  with per-method notes on what's faithful vs. departed from the
  original method (BMHT: ranking only, no downstream clique analysis;
  ROS-DET: both original permutation testing and the paper's own
  analytical ECOR test available; BMKC: k-core in place of exact
  k-clique percolation). Verified `citation("bicorDCEA")` renders
  correctly.
- Rewrote `DESCRIPTION`'s `Description:` field to state plainly that
  this extends three previously-published, never-packaged methods with
  explicit, documented departures - not a claim of strict fidelity to
  any of the three. Added DOIs for all three papers.

### Test suite
16 files, 80+ assertions, all passing as of this entry.

---
### Paper-fidelity update, Phase C: re-tested weight_mode on the scenario the paper actually validated it against

The earlier default change (`weight_mode` from `"wcor"` to `"unweighted"`)
was based on a compositional-closure confound test - a different scenario
from Kayano et al.'s own validation of WCOR, which targets *range bias*
(one condition's expression range much smaller than the other's; their
Fig. 4b/c). Re-tested directly on a range-bias scenario matching the
paper's design.

**A first attempt at this retest was discarded before drawing any
conclusion from it**: it conflated true signal and variance/range
difference in the same synthetic genes, which isn't a valid test of
either factor in isolation - and, worth noting, matches the paper's own
ground-truth criterion for this experiment (their "positive" class for
the range-bias test is explicitly defined via a Bartlett test as cases
*without* a significant variance difference, so a case with both real
signal and real variance difference wouldn't even count as a valid
"positive" by their own definition).

**Redone with the two factors properly separated**, giving a precise,
three-part characterization:
1. `"wcor"` correctly suppresses a *pure* range-bias artifact (variance
   difference, no real correlation change) to near-zero score;
   `"unweighted"` has no protection against this at all - matching the
   paper's design intent exactly.
2. On *clean* true signal with variance genuinely equalized between
   conditions (confirmed by directly checking the variances, not just
   assuming the construction achieved it - an early version of this test
   also had a hidden variance confound from how correlation was injected
   via a shared latent factor, caught by checking rather than assuming),
   `"wcor"` and `"unweighted"` perform comparably (0.843 vs. 0.893 mean
   score, true vs. background) - `"wcor"` is not broadly worse at
   detecting real signal.
3. The real, narrower failure mode: when true signal is naturally
   *accompanied by* a variance/range shift (plausible in real biology,
   where regulatory rewiring often co-occurs with variance changes),
   `"wcor"` suppresses the signal along with the variance shift, rather
   than distinguishing the two.

**`"unweighted"` remains the default** - the safer general-purpose
choice, since signal/variance-shift coupling seems more common in
practice than the narrow pure-artifact scenario `"wcor"` protects
against. `"wcor"` remains available and is a reasonable choice
specifically when technical/range artifacts, not biological variance
shifts, are the concern. Documentation on `weight_mode` fully rewritten
to reflect this three-part picture rather than the earlier single-test
finding. Locked in with `tests/testthat/test-weight-mode-range-bias.R`.

### Test suite
16 files, 80+ assertions, all passing as of this entry.

---
### Paper-fidelity update, Phase B: BMKC significance test (Yuan et al. 2015 sec 3.3)

BMKC previously returned modules with no attached notion of statistical
significance at all - a real defensibility gap independent of paper
fidelity. Implemented Yuan et al.'s own global permutation test: each
gene's expression is independently permuted within each condition
(destroying all gene-gene covariance while preserving each gene's own
marginal distribution and the condition grouping), the full module-
detection pipeline (bicor -> adjacency -> k-core module extraction) is
rerun on the shuffled data, and a "global score" - summed per-module,
per-gene mean |bicor difference| to other genes in the same module - is
computed. Repeating this gives an empirical null distribution against
which the real result's score is compared, yielding one p-value for the
whole set of detected modules (not per-module or per-gene FDR - matching
the paper's own coarser-grained test).

**Implementation approach**: added `n_permutations` (default `0`, fully
backward compatible - the exact previous behavior with no significance
test attached) to `.run_bmkc_matrix()`. The null-generation and rerun
logic (`.bmkc_null_replicate()`) is a deliberately separate,
self-contained reimplementation of the core pipeline steps, rather than
a call back into `.run_bmkc_matrix()` itself - this kept the change
additive and avoided any risk of altering the already-validated default
path's behavior.

**A real bug was caught and fixed before this was tested at all**: an
early memory-cleanup step (`rm(cor_1, cor_2)`) needed to be deferred
whenever the significance test is requested (since it needs those
matrices afterward), and the fix initially introduced a double-removal
that would have errored the moment `n_permutations > 0` was used.
Caught by tracing the control flow before the first test run, not by the
test run itself - worth naming since it's exactly the kind of mistake
this project's discipline exists to catch early.

**Validation, in order**:
1. First real execution: ran cleanly on data with an injected 15-gene
   module, correctly hit the permutation floor (p ~ 1/21) with sensible
   timing (~0.5s for 20 replicates).
2. Calibration: 60 independent pure-noise datasets, each significance-
   tested, gave proportion p<0.05 = exactly 0.05 (nominal), no
   anti-conservative bias.
3. Timing check at a realistic scale: ~6s projected for the paper's own
   default of 1000 permutations at 60 genes - a real, usable target, not
   a token gesture.

Locked in with `tests/testthat/test-bmkc-significance.R`: default-path
regression guard (bit-for-bit unaffected when `n_permutations = 0`),
real-signal detection, a smaller permanent calibration spot check, and a
direct unit test of `.bmkc_global_score()`'s arithmetic.

### Test suite
15 files, 75+ assertions, all passing as of this entry.

---
### Paper-fidelity update, Phase A: implemented ROS-DET's ECOR analytical test - closes almost the entire Phase 5 competitiveness gap

Reading the three source papers (Zheng et al. 2014, Kayano et al. 2011,
Yuan et al. 2015) against the implementation surfaced the single most
consequential finding of this project: **Kayano et al.'s ROS-DET
specifies a fast, analytical significance test (Eq. 4, "ECOR" - a
chi-squared(1) likelihood-ratio test for equality of two correlation
coefficients), not a permutation test.** The implementation had always
used permutation testing instead. That substitution is the direct cause
of the Phase 5 benchmark gap (bicorDCEA trailing DGCA/dcanr by ~3.5x in
ranking quality) and of every failed attempt to patch the permutation
approach afterward (GPD x2, normal, skew-normal - see the Post-Phase-5
entries below). None of that earlier work was wasted - it rigorously
ruled out an entire family of fixes - but the actual fix was in the cited
paper the whole time.

**Implemented and validated before touching any engine**, with the same
discipline as everything else in this project:

- Derived and prototyped the pooled-correlation MLE (an implicit equation
  solved via `uniroot()`) and the chi-squared(1) statistic directly from
  the paper's Equation 4.
- Calibration-tested on **true Pearson data** first (the paper's own
  exact assumptions) across a range of sample sizes - confirmed the
  expected asymptotic behavior: real bias at N=20 (ratio ~1.4-1.5x),
  converging to good calibration (~1.0-1.2x) by N~=50-200.
- Calibration-tested on **bicor** (the paper's disclosed approximation,
  and what ROS-DET actually needs) using the package's own bicor kernel
  on pure-noise data - result: **~1.0-1.2x at N>=50, degrading gracefully
  (not catastrophically) below that** (~1.7-2.2x at N=10). This is a
  materially different, far better-behaved failure mode than any of the
  four previously-rejected parametric approaches (which ranged 3.5x-14.5x
  and sometimes worsened unpredictably deeper into the tail) - because
  this is a theoretically-derived asymptotic result, not a distribution
  empirically fit to permutation output.

**Implemented as `significance_method = c("permutation", "analytical")`**
in `.run_rosdet_matrix()` - an offered alternative, not a replacement.
`"permutation"` remains the default and is completely unaffected (full
regression suite confirms bit-identical behavior). `"analytical"` warns
below N=30 per condition, based on the calibration data above, and
requires no resampling at all - `n_permutations` becomes irrelevant to
that mode entirely.

**The benchmark result**: re-ran Phase 5's sim102 ground-truth benchmark
(`dev/benchmarks/03_sim102_ground_truth.R`, 406 samples - well above the
N=50 threshold) with the new mode:

| Method | AUPRC | Time |
|---|---|---|
| ROS-DET (permutation, as it always was) | 0.2107 | 1.09s |
| **ROS-DET (analytical/ECOR, new)** | **0.6871** | 2.37s |
| DGCA | 0.7321 | 0.07s |
| dcanr | 0.7321 | 0.00s |

**This closes the gap from ~3.5x to ~1.07x** - ROS-DET with ECOR is now
nearly on par with DGCA/dcanr in ranking quality, not in a different
competitive league. It remains slower per-analysis (~34x vs. dcanr's
fully analytical approach) but no longer needs the very large permutation
counts Phase 5 identified as impractical - the entire resolution-vs-cost
tradeoff that motivated the whole permutation-floor investigation is
sidestepped for this significance method. Ranked by `-log10(p)`, not
`Distance_Score` (which is unaffected by which significance method is
used - the improvement is entirely in p-value resolution).

Locked in with `tests/testthat/test-ecor-calibration.R`: calibration on
both Pearson and bicor data, graceful (bounded) small-N degradation, the
N<30 warning firing/not-firing correctly, a real-data sanity check (the
known 5-gene signal in the multi-omics fixture dominates the top ranks),
and a regression guard confirming the default permutation path is
completely unaffected.

**Scope note**: this is ROS-DET-specific. Zheng et al.'s BMHT paper does
not specify an analytical alternative to permutation testing, so BMHT's
own permutation-floor characteristics (if any remain after this) are
unchanged by this work.

### Test suite
14 files, 65+ assertions, all passing as of this entry.

---
### Post-Phase-5, part 4: characterized the null distribution shape - another negative result, honestly reported

Followed up on TODO.md's Tier-1 recommendation: characterize what BMHT's
`DC_Score` null distribution actually looks like, to see whether a
distribution-specific parametric model (rather than the generic GPD
extrapolation already tried and rejected) could pass calibration testing.
Reproducible script: `dev/permutation_floor/null_distribution_investigation.R`.

**Shape survey**: `DC_Score`'s null is an *average* (not a sum) of many
per-edge squared terms, so by CLT it's pulled toward Gaussian-like shape
rather than the chi-squared/gamma shape initially hypothesized. Confirmed
empirically: most genes' null distributions are close to Gaussian, with a
real but variable right skew for others (skewness -0.16 to 1.34 across
genes in the toy dataset). Checked whether skewness is explained by the
number of "informative" edges contributing to each gene's score (a
derivable design property, which would support a principled per-gene
distribution choice) - it isn't (correlation 0.15, weak and in the
opposite direction from the "fewer edges -> more skewed" hypothesis).

**Calibration testing (the decisive check, not the eyeball comparison)**:
a plain normal fit to the bulk of the permutation null looked deceptively
reasonable in a quick quantile-by-eye comparison (differences of
0.01-0.09 in absolute quantile position), but a proper self-consistency
calibration simulation - using the actual empirical null shape from real
toy-dataset output this time, not a synthetic Gaussian proxy - showed it
is badly anti-conservative: **14.5x too many p-values below 0.001**. A
skew-normal fit (3-parameter, matching the actual skewness) is a real,
substantial improvement (**3.5x inflation at the same threshold**, down
from 14.5x) but still clearly fails calibration.

This is now the fourth distinct parametric/semi-parametric approach
tested end to end (GPD method-of-moments, GPD MLE, plain normal,
skew-normal), and all four show the same systematic pattern -
underestimating how thin the true tail actually is - to different
degrees. Skew-normal is the best of the four but "best of four failures"
is not a basis for shipping it; **none of the four has been integrated
into any engine.**

**Revised recommendation**: distribution-specific parametric modeling
does not appear to close the permutation-floor gap safely with any
approach tried so far. TODO.md's option (c) - defaulting to a much higher
`n_permutations`, now that per-permutation cost and memory are both
reduced (see the Performance/RAM sections below) - is now the leading
recommendation, ahead of further parametric-fit attempts. A genuinely
different idea (an analytically-derived asymptotic distribution for
bicor-based DC-scores, rather than one fit empirically to permutation
output) is a bigger undertaking not attempted here.

---
### Post-Phase-5, part 3: RAM optimizations (found substantial pre-existing work, verified and extended it)

Asked to look specifically for RAM/memory improvements. Found real,
substantial, well-designed work already on disk that hadn't been
committed or written up yet - same situation as the Phase 5 benchmark
discovery. Fully verified it (didn't just trust it existed), then found
and fixed one more instance of the same pattern.

**BMHT: eliminated an unnecessary full permutation-result matrix.**
`.run_bmht_matrix()` previously allocated a full `n_genes x
n_permutations` matrix and retained it for the *entire run*, even though
it was only ever used for one `rowSums()` reduction at the very end. At
the very large permutation counts Phase 5 found necessary for adequate
power (thousands to tens of thousands), this was real, unnecessary
memory - the arithmetic is straightforward: 5,000 genes x 100,000
permutations x 8 bytes is 4GB for something transient. Replaced with
incremental per-chunk accumulation of just the two integer counts
actually needed (`greater_eq_counts`, `valid_perms_counts`). Also changed
`chunk_size` from scaling proportionally with `n_permutations` (so chunks
themselves grew unboundedly at large permutation counts) to a fixed,
bounded size - this also fixes progress-bar granularity collapsing at
large runs, not just memory.

**BMKC: eliminated a redundant O(n^2) extraction.** `get_quantile()`
recomputed the full upper-triangle extraction (and, for large matrices,
random subsampling) from scratch for *each* of the T1/T2 threshold calls
on the same correlation matrix - duplicating an O(n_genes^2) operation
for no reason, since the extraction doesn't depend on which threshold is
being computed. At n_genes=20,000 that intermediate vector alone is
~1.6GB, and it was being built twice per matrix. Split into
`extract_and_subsample()` (done once) and `quantile_from_vals()` (cheap,
called twice). **Bonus correctness improvement, not just performance**:
previously, T1 and T2 thresholds for large matrices were each computed
from an *independent* random subsample, meaning the two thresholds for
the same matrix weren't even derived from the same data - now both come
from one consistent subsample.

**ROS-DET: applied the same bounded-chunk-size fix for consistency.**
Found this while checking whether the BMHT pattern existed elsewhere.
ROS-DET's C++ kernel already reduces to per-pair counts internally (never
returns a full null-value matrix to R), so the memory impact here is
smaller than BMHT's fix - but `chunk_size` had the same
unboundedly-scaling-with-`n_permutations` pattern, which still matters
for progress-bar granularity at very large runs. Fixed identically.

**Verification, beyond what was already there:** compiled and ran the
full test suite (all passing); specifically exercised BMKC's percentile
threshold code at a gene count large enough to actually trigger the
subsampling path (4,500 genes, ~10.1M-element upper triangle - none of
the existing tests were large enough to reach this code path at all, only
the non-subsampling path had been checked); measured real peak memory at
that scale (~1.17GB); confirmed the ROS-DET chunk-size fix reproduces
exact historical reference values (9,518 pairs tested, 0 significant,
same top-ranked pair and Distance_Score, matching every prior baseline
run throughout this project). Extended `test-bmht-memory.R` (renamed
`test-memory-optimizations.R` to reflect its broader scope) with the
ROS-DET check.

### Test suite
13 files, 65+ assertions, all passing as of this entry.

---
### Post-Phase-5, part 2: performance investigation and two validated improvements

Asked directly: are there real performance improvements available, and can
"more C++" help? Investigated four concrete hypotheses rather than
guessing, with the same test-before-shipping discipline as the rest of
this project.

**Ruled out:**
- **Moving more logic into C++**: profiled a real run - 92% of wall-clock
  time is already inside the compiled kernel. No headroom here; the
  bottleneck was never R-side overhead.
- **Per-gene adaptive early stopping** (the classic genomics-permutation-
  testing trick of dropping "resolved" genes from further permutation to
  save compute): confirmed empirically via direct timing that the kernel
  scales as clean O(n_genes^2), because the correlation matrix is computed
  once for *all* active genes per permutation. Dropping a gene from
  tracking doesn't shrink that matrix without silently changing every
  other gene's test definition (each gene's score is an average over edges
  to every other active gene). Ruled out before implementing, not after.

**Validated and shipped:**
- **BLAS library choice.** This environment was linked against R's
  reference BLAS (unoptimized, no vectorization/cache blocking) - a system
  configuration matter, not a bicorDCEA code issue. Installed OpenBLAS and
  measured a genuine **~1.3-1.4x speedup** on the correlation-matrix
  computation dominating this package's cost, with a properly
  iteration-count-matched comparison. (An initial quick comparison
  suggested ~5-6x; that used mismatched iteration counts between the two
  measurements and was wrong - corrected before reporting either number.)
  Documented in `DESCRIPTION`'s new `SystemRequirements` field with the
  concrete `update-alternatives` command.
- **Per-thread buffer reuse in the hot permutation loops** (`bmht_cpp.cpp`,
  `rosdet_cpp.cpp`). Both kernels allocated fresh correlation-matrix
  temporaries (`Xt`, `X_tilde`, `cov`, `denom`, `res`) on every single
  permutation; A/B testing in isolation showed this costs essentially
  nothing at 50-300 genes (the scale used throughout this project's own
  benchmarking - noise-level difference) but grows to a real **~1.5-1.6x**
  overhead at 600+ genes, relevant now that adequate statistical power can
  require very large permutation counts. Implemented via a `BicorWorkspace`
  struct holding pre-sized buffers, allocated once per OpenMP thread and
  reused across that thread's permutations (valid because group sizes are
  identical across all permutations of a fixed label multiset - only which
  samples get which label changes).

  Verified three ways before considering this safe to ship: (1) the
  refactored code reproduces the exact same numbers as Phase 1's original
  calibration validation on the toy dataset; (2) an independent
  gold-standard R reimplementation, built only from the untouched
  single-call `cpp_bmht_observed_bicor()`/`cpp_rosdet_observed_bicor()`
  kernels (which share no code with the refactored permutation loops),
  matches to machine precision; (3) repeated calls with identical inputs
  give bit-identical output (rules out stale-buffer leakage between
  iterations). Locked in with `test-workspace-parity.R`.

**Measured combined effect** (OpenBLAS + buffer reuse together, real
pipeline, not the isolated test harness): on the toy dataset (200 genes),
BMHT dropped from ~0.76s to ~0.55s and ROS-DET from ~0.50s to ~0.35s
(~28-31% faster) for numerically identical results (same 20 significant
BMHT genes, same 0 ROS-DET hits). Re-ran the Phase 5 sim102 ground-truth
benchmark: ROS-DET's runtime dropped further (1.28s -> 1.00s) while its
**AUPRC stayed exactly 0.2107** (unchanged from Phase 5) - exactly the
expected, correctly-scoped outcome, since these are pure speed
improvements and were never expected to touch the permutation-floor power
problem, which remains open (see the "attempted fix" section above).

Also surveyed the full range of biological data types and input formats
this package supports (matrix/data.frame/SummarizedExperiment/Seurat/
MultiAssayExperiment; bulk, single-cell via pseudobulk aggregation,
compositional/microbiome via CLR/rCLR, multi-omics via `layer`, pseudotime
trajectories; strictly two-condition comparisons only) - conversational
for now, worth turning into a proper vignette section (see TODO.md).

### Test suite
12 files, 60+ assertions, all passing as of this entry.

---
### Post-Phase-5: attempted the recommended permutation-floor fix - negative result, honestly reported

Phase 5 recommended fixing the permutation-resolution floor (Findings 3-5)
before further packaging polish. Attempted a standard approach - fitting a
Generalized Pareto Distribution (GPD) to the upper tail of the permutation
null and using it to extrapolate p-values below the raw permutation floor
(a textbook peaks-over-threshold extreme-value technique, in the spirit of
what DGCA's semi-analytical approach already does). Implemented as
`.gpd_tail_pvalues()` in `utils.R`.

**Before attempting any integration into the engines, tested it against a
50,000+ replicate self-consistency calibration simulation - and it failed,
under both estimators tried:**

- **Method-of-moments GPD fit**: ~27% Type I error inflation right at the
  extrapolation boundary (a known small-sample bias in that estimator's
  shape parameter), though over-conservative deeper into the tail.
- **Maximum-likelihood GPD fit** (tried specifically to fix the
  method-of-moments problem): *worse*, not better - consistently
  anti-conservative, with observed/nominal Type I error ratios from 1.2x
  up to 4x deeper into the tail. This persisted even at 10x more
  permutations (5,000 vs. 500), ruling out "just needs more permutation
  samples" as the fix.

This is consistent with a known, documented limitation of peaks-over-
threshold extreme value approximation: light-tailed distributions
(Gaussian and similar) converge to their asymptotic GPD limit very
slowly, making generic nonparametric EVT extrapolation unreliable at
permutation counts anyone would actually run. Fixing the permutation
floor properly likely needs either a null-distribution-specific
parametric model (using the actual known/derivable shape of BMHT's
DC_Score or ROS-DET's Distance_Score null, rather than a generic
technique) or a different approach entirely - not attempted here given
the calibration risk and remaining scope.

**Decision**: not integrated into any engine. `.gpd_tail_pvalues()` is
left in the codebase, unused, with prominent documentation of what was
tried and why it failed, specifically so this isn't silently
re-attempted the same way by someone who hasn't seen this. Locked in
with `tests/testthat/test-gpd-tail-calibration.R`, including a test that
asserts the current anti-conservative behavior is present (so if a future
fix actually resolves it, that test's failure is the signal to update
both the test and the "STATUS" section in the function's docs).

**The permutation-floor problem from Phase 5 remains open.** This
attempt at the "obvious" fix not working is itself useful information for
whoever picks this up next - it rules out one well-known approach and
narrows what's worth trying.

### Test suite
11 files, 55+ assertions, all passing as of this entry.

---
### Phase 5 (benchmark - the decision gate)

Installed `DGCA` and `dcanr` from their GitHub source (neither is on this
sandbox's apt mirror) and ran six separate comparisons against
bicorDCEA. Full report: `dev/benchmarks/BENCHMARK_REPORT.md`. This
section was substantially revised partway through from an earlier,
more upbeat draft once two additional benchmarks already in progress -
a compositional-confound test and a real-published-dataset test with
known ground truth - were fully reconciled. The fuller picture is more
serious than the initial draft suggested.

**What's validated**: under sample corruption (outliers injected into
2-4 of 24 samples, restricted to the true-signal genes), DGCA's precision
collapses (0.68 clean -> 0.38 at 17% corruption - the majority of its
calls become false positives) while bicorDCEA's holds or improves (up to
a perfect 1.00). This is exactly the behavior biweight midcorrelation's
outlier-downweighting is supposed to produce, and it's real,
mechanistically explained, and replicated. bicorDCEA's compiled engine is
also genuinely fast per unit of work - ~12x faster than DGCA on a real
published single-cell dataset (Darmanis et al. 2015, bundled with DGCA).

**What's a real, current weakness**: on the one benchmark with actual
external ground truth (dcanr's own bundled `sim102` simulation, 185 true
edges / 11,026 candidate pairs), bicorDCEA's ROS-DET scores AUPRC 0.211
against DGCA/dcanr's 0.732 each - roughly 3.5x worse - and takes ~40x
longer to do it. The same pattern - bicorDCEA trailing both competitors -
shows up in a from-scratch compositional-confound test (CLR improves
bicorDCEA's own AUPRC from 0.36 to 0.55, a real improvement, but still
well short of DGCA/dcanr's 0.81 on the same *uncorrected* data) and in a
clean-data recall test (bicorDCEA found zero significant genes at a
typical permutation count where both competitors found everything,
recovering only 8/15 even at 10,000 permutations). All three point at the
same root cause: bicorDCEA's pure permutation-testing p-values have a
coarse resolution floor (`1/(n_permutations+1)`) that DGCA's
semi-analytical and dcanr's fully-analytical (Fisher z-test) approaches
don't share, and it's costing real detection power, not just requiring
"a few more permutations."

**A concrete consequence of this benchmark work**: `weight_mode`'s
default was changed from `"wcor"` to `"unweighted"` after the
compositional-confound benchmark found `"wcor"` can actively *invert*
`Distance_Score` ranking on real signal (true-signal pairs scoring lower
than background, AUPRC 0.034 - worse than random), because true
differential-correlation pairs often undergo a genuine variance shift as
part of the same rewiring event, and `"wcor"`'s variance-ratio weighting
systematically penalizes exactly that shift. `"unweighted"` on the same
data correctly separates signal from background (AUPRC 0.359). `"wcor"`
remains available for cases where down-weighting a variance shift is
specifically wanted, but is no longer the default. This does not affect
`P_value`/`FDR` (see Phase 1's finding that `weight_mode` never did) -
only `Distance_Score`-based ranking. Locked in with an expanded
`test-rosdet-weight-mode.R`.

**Two bugs were found and fixed in the benchmark harness itself**
(`03_sim102_ground_truth.R`), disclosed here because they invalidated an
earlier draft's numbers (which showed ROS-DET scoring *below random* -
that number was wrong and is not the one reported above): (1) the
simulation's ground-truth network matrix is asymmetric (165 of 185 true
edges live only in the lower triangle, not mirrored), and an earlier
version of the script checked only the upper triangle; (2) `combn()`'s
pair ordering isn't alphabetical, but one score-matching key was
canonicalized alphabetically, silently turning real scores into `NA` for
roughly a third of all pairs. Fixed by symmetrizing the ground truth and
using one consistent canonical key everywhere.

**Verdict**: not a clean "ship it," not a clean "abandon it." The core
robust-correlation idea (Finding 2 above) is validated and worth keeping.
But the permutation-testing machinery - not the underlying bicor
computation, which is both correct (Phase 1) and fast (Findings above) -
is a concrete, well-scoped bottleneck that should be addressed before
further packaging polish: a semi-parametric approximation to the
permutation null (extrapolating tail p-values below the raw permutation
floor, in the spirit of what DGCA already does), rather than requiring
tens of thousands of permutations for adequate power.

### Test suite
10 files, 50+ assertions, all passing as of this entry. Benchmark scripts
under `dev/benchmarks/` are not part of the package's own test suite
(they depend on DGCA/dcanr, which aren't package dependencies), but the
`weight_mode` default change they motivated is locked in by the package
test suite itself.

---

### Phase 4 (multi-omics as the story)
`run_multiomic()` generalizes cross-layer analysis to all three methods;
`layer` added as a documented synonym for `species`; per-group CLR
transform properly documented; worked RNA+ATAC fixture with a verified,
real cross-layer signal added, plus a walkthrough vignette.

### Phase 3 (surface reduction)
Cut plotting surface from ~40 exported functions to 8; rewrote all three
`plot.<result>_result()` dispatchers to match; added a vignette showing
how to rebuild removed visualization types from the result object
directly; verified every surviving function and dispatcher option
end-to-end.

### Phase 2 (structural collapse)
Collapsed two parallel, duplicated S3 dispatch trees into one
(`.extract_inputs()` + `.run_engine()` + `run_bicordcea()`); renamed
`run_dcanr()` -> `run_bicordcea()` (Bioconductor name collision); added
`data.frame` input support; removed dead code and an unused dependency.

### Phase 1 (correctness)
BMHT null calibration fixed (permutation mask was leaking observed signal
into the null - 80% false positives on pure noise reduced to 0%);
`cpp_calc_observed_bicor` didn't exist anywhere despite ~25 call sites;
MAD scaling corrected (raw MAD verified correct against `WGCNA::bicor()`
to machine precision); Pearson fallback added to all three C++ kernels;
`species`/zero-MAD-filter alignment fixed; several smaller bugs fixed and
locked in with regression tests.

### Phase 0 (re-entry)
Fixed a typo that made ROS-DET completely non-functional; fixed OpenMP
being inert in package builds; added undeclared dependencies; built the
toy dataset and baseline regression script.

### Known issues still open (see TODO.md)
- **The permutation-resolution gap identified in Phase 5** - the single
  highest-leverage next engineering step if this package continues, ahead
  of any further packaging polish.
- BMKC's default thresholds are calibrated for larger sample sizes than
  small toy fixtures provide (Phase 3 finding, confirmed not a bug).
- The "fake Seurat" test fallback breaks once real `SeuratObject` is
  loaded (Phase 2 finding, not fixed).
- Phase 5's strongest evidence (Findings 5-6) used external data; Findings
  1-4 used synthetic data with self-constructed ground truth and should
  be read as directional, not definitive.
