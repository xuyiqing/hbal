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
* TODO: win-builder, R-devel
* TODO: macbuilder

## R CMD check results

0 errors | 0 warnings | 0 notes

The examples produce no warnings under either 'estimatr' version.

## Reverse dependencies

None. tools::package_dependencies("hbal", reverse = TRUE) returns character(0),
so no package on CRAN is affected by the change described below.

## Notes for the reviewer

* att()'s default estimator has changed in this release. The default is now
  method = "abw", a cross-fitted, Neyman-orthogonal augmented balancing weights
  estimator. The previous default is unchanged and remains available as
  att(x, method = "lm_robust"); "lm_lin" and "elnet" are unchanged as well.
  The change is user visible and is the first entry of the 1.3.0 section of
  NEWS.md.
* hbal() takes a user-facing `seed` argument that it passes to set.seed(), so
  that its optional cross-validation is reproducible. The default is a fixed
  number and `seed = NULL` switches the call off entirely. The argument and its
  default are unchanged from version 1.2.15.
