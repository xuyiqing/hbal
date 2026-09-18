# hbal 1.3.0
1. **`att()`'s default estimator has changed.** The default is now a cross-fitted, Neyman-orthogonal augmented balancing weights estimator (`method = "abw"`); the previous default is unchanged and available via `att(out, method = "lm_robust")`.
2. The `"abw"` outcome model is a ridge regression on the standardized columns of `mat`; its penalty is chosen once from the control sample by generalized cross-validation on a fixed grid with a cap on effective degrees of freedom (deterministic, no random numbers).
3. `att()` gains `seed` and `nfolds` arguments for `method = "abw"` (ignored by the other methods); `method = "abw"` does not accept `...` arguments such as `se_type` or `clusters`, which still work with `"lm_robust"` and `"lm_lin"`.
4. `att()` now errors on an unrecognized `method` value instead of failing with an opaque internal error.
5. `att()` with `method = "abw"` now warns when `hbalobject$converged` is 0 (the entropy-balancing weights did not converge), since the estimate and standard error may then be unreliable.
6. The package tutorial is now a Quarto book at `tutorial/` (previously `vignettes/tutorial.Rmd`); render it locally with `quarto render` from that directory.
7. `plot()` for `hbal` objects no longer calls the deprecated `ggplot2::aes_string()`; the plots are unchanged and the deprecation warning is gone.
8. `plot()`'s weight histogram now labels its x axis "Weights (log)"; it has always shown the natural log of the weights, and the plot itself is unchanged.
9. `?hbal` now lists the elements the returned object actually has, and describes what the default settings do, including that `cv`'s default (`NULL`) means no cross-validation.
10. The documented return values of `plot()`, `summary()`, `covarExpand()` and `crossValidate()` now match what those functions return.
11. `src/*.o` and `src/*.so` are no longer tracked in the repository, so installing from a GitHub clone always recompiles the C++ code.
12. `hbal()` no longer sets a random seed by default: the `seed` argument now defaults to `NULL`; pass `seed = 94035` to reproduce cross-validated results from earlier versions.
13. `plot()` for `hbal` objects reports the weight normalization with `message()` instead of `cat()`, so it can be silenced.
14. The sources and licenses of the bundled `lalonde` and `contenderJudges` datasets are now documented.
15. `att(dr = FALSE)` is no longer silently ignored by `method = "lm_lin"` and `method = "elnet"`; both now return the weighted difference in means, as `"abw"` and `"lm_robust"` already did. The undocumented value `method = "lin"`, which the old guard accepted when `dr = FALSE`, is now rejected like any other unrecognized method. (#8, thanks `@soodoku`)
16. `hbal()` no longer balances the wrong covariates when `data` already holds a column named `X1`, `X2`, ... that is not among `X`, or when `X` is not supplied in data-frame order; covariate labels in `bal.tab`, `mat` and the balance plot are unchanged. (#7, thanks `@soodoku`)
17. New exported function `balanceData()` returns an `hbal` object's balance statistics as a long data frame with a covariate-group column — the data `plot()` draws, ready for a custom `ggplot2` figure. (#4, thanks `@joshuafayallen`)

# hbal 1.2.16
1. Fix `att()` to work with estimatr >= 2.0.0, whose `tidy()` now returns a tibble, without warnings, while remaining identical under estimatr < 2.0.0.

# hbal 1.2.13
1. In `att()`, add the method "elnet" by Athey (2018)
2. Fix a bug with incorrect variable names

# hbal 1.2.12
1. In `att()`, only display the treatment effect by default
2. Fix a bug with truncated variable names
3. Show a warning when the algorithm is not converged

# hbal 1.2.11
1. Allow approximate balance on the level terms (penalty tuned via cross-validation)
2. Fix a bug with truncated variable names
   
# hbal 1.2.9
Fix a bug which occurs when there is only one covariate

# hbal 1.2.8
* First CRAN release

# hbal 1.2.2
* beta version
