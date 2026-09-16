# hbal 1.3.0
1. **`att()`'s default estimator has changed.** The default is now a cross-fitted, Neyman-orthogonal augmented balancing weights estimator (`method = "abw"`); the previous default is unchanged and available via `att(out, method = "lm_robust")`.
2. The `"abw"` outcome model is a ridge regression on the standardized columns of `mat`; its penalty is chosen once from the control sample by generalized cross-validation on a fixed grid with a cap on effective degrees of freedom (deterministic, no random numbers).
3. `att()` gains `seed` and `nfolds` arguments for `method = "abw"` (ignored by the other methods); `method = "abw"` does not accept `...` arguments such as `se_type` or `clusters`, which still work with `"lm_robust"` and `"lm_lin"`.
4. `att()` now errors on an unrecognized `method` value instead of failing with an opaque internal error.
5. `att()` with `method = "abw"` now warns when `hbalobject$converged` is 0 (the entropy-balancing weights did not converge), since the estimate and standard error may then be unreliable.
6. The package tutorial is now a Quarto book at `tutorial/` (previously `vignettes/tutorial.Rmd`); render it locally with `quarto render` from that directory.
7. `plot()` for `hbal` objects no longer calls the deprecated `ggplot2::aes_string()`; the plots are unchanged and the deprecation warning is gone.

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
