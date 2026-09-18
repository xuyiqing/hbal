## Submission

hbal 1.3.0, an update to hbal 1.2.15, which is currently on CRAN.

## Reason for this submission

The check page for hbal 1.2.15 shows a WARNING on
r-devel-linux-x86_64-fedora-gcc: three of att()'s documented examples emit
"Setting row names on a tibble is deprecated.". The cause is 'estimatr' 2.0.0,
published on CRAN on 2026-09-16, whose tidy() method now returns a tibble;
att() then set row names on it. This release fixes that: att() selects the
treatment row and the seven display columns by name and returns a plain
data.frame, so no row names are ever set on a tibble. It was checked here under
both 'estimatr' 2.0.0 and 'estimatr' 1.0.6, with no warnings from the examples
or the test suite under either.

## Test environments

* local: macOS 26.6 (aarch64), R 4.6.1, checked twice: once with 'estimatr'
  2.0.0 (the current CRAN version) and once with 'estimatr' 1.0.6
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  windows-latest (R-release), macOS-latest (R-release)
* win-builder (R-devel): 0 errors | 0 warnings | 0 notes
* macbuilder (R-release, macOS arm64): 0 errors | 0 warnings | 0 notes

## R CMD check results

0 errors | 0 warnings | 0 notes

The examples produce no warnings under either 'estimatr' version.

## Reverse dependencies

None. tools::package_dependencies("hbal", reverse = TRUE) returns character(0),
so no package on CRAN is affected by the change described below.

## Bundled data

The package ships two datasets, both in data/hbal.RData. contenderJudges is the
circuit-court judges data of Black and Owens, deposited at Harvard Dataverse as
doi:10.7910/DVN/25302. lalonde is the LaLonde (1986) / Dehejia and Wahba (1999)
National Supported Work data as distributed in the replication archive of Xu and
Yang (2022), deposited at Harvard Dataverse as doi:10.7910/DVN/QI2WP9. Both
deposits are released under the Creative Commons CC0 1.0 Universal Public Domain
Dedication, which permits redistribution and commercial use. Each dataset's help
page now documents its source and license in a \source section; see ?lalonde and
?contenderJudges.

## Notes for the reviewer

* att()'s default estimator has changed in this release. The default is now
  method = "abw", a cross-fitted, Neyman-orthogonal augmented balancing weights
  estimator. The previous default is unchanged and remains available as
  att(x, method = "lm_robust"). "lm_lin" and "elnet" are unchanged except under
  dr = FALSE, where they now return the weighted difference in means instead of
  silently ignoring the argument; see item 15 of the 1.3.0 section of NEWS.md.
  The change is user visible and is the first entry of the 1.3.0 section of
  NEWS.md.
* hbal() no longer sets a random seed by default. Its `seed` argument has
  defaulted to the hard-coded value 94035 since version 1.1.1, so every call ran
  set.seed(94035) and overwrote the user's .Random.seed. The default is now
  NULL, and a default call leaves the random number generator untouched. A call
  that passes a seed explicitly behaves exactly as before: hbal(..., cv = TRUE,
  seed = 94035) reproduces the cross-validated results of version 1.2.15
  exactly. Results without cv = TRUE are unchanged either way, since the
  non-cross-validated path draws no random numbers.
* hbal() renamed the covariates in X to X1..Xk internally, and the generated
  names could collide with columns of `data` that were not being renamed, or be
  assigned in the wrong order. Either way the function balanced covariates other
  than the ones requested, without a warning. Both are fixed. Results change only
  for calls that hit one of those two conditions; every other call is bit-identical,
  which is checked by a numeric fixture in the test suite.
* This release adds one exported function, balanceData(), which reshapes an
  hbal object's balance table for plotting. It adds no dependency.
