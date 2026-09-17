# Contributing to bicorX

Thanks for your interest in contributing.

## Reporting bugs

Please open an issue at
<https://github.com/paytonyau/bicorX/issues> including:

- a minimal reproducible example (the smallest data and code that shows
  the problem),
- the output of `sessionInfo()`,
- what you expected to happen and what happened instead.

## Suggesting features

Open an issue describing the use case before writing code, so we can
discuss whether and how it fits the package's scope.

## Pull requests

1. Fork the repository and create a branch from `main`.
2. Make your change, following the existing code style.
3. **Add or update tests.** This package has a test suite under
   `tests/testthat/`; changes to statistical behaviour in particular are
   expected to come with a test that would fail without the change.
4. Run `R CMD check` locally and make sure it passes without new
   warnings.
5. Update `NEWS.md` describing your change.
6. Open the pull request, referencing any related issue.

## A note on statistical changes

Several parts of this package's history involve statistical methods that
looked correct but failed calibration testing (see `NEWS.md` for the full
record, including approaches that were tried and rejected). If you are
changing anything that affects p-values, significance testing, or null
distributions, please include a calibration check — for example, a
simulation under a known null showing that observed Type I error rates
match nominal ones. Several existing tests
(`test-ecor-calibration.R`, `test-gpd-tail-calibration.R`,
`test-bmht-calibration.R`) show the pattern used elsewhere in the
package.

## Code of Conduct

Please note that this project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating you
agree to abide by its terms.
