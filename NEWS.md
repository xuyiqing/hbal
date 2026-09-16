# hbal 1.3.0
1. **`att()`'s default estimator has changed.** The default is now a cross-fitted, Neyman-orthogonal augmented balancing weights estimator (`method = "abw"`); the previous default is unchanged and available via `att(out, method = "lm_robust")`.
2. The `"abw"` outcome model is a ridge regression on the standardized columns of `mat`, with a deterministic penalty chosen once from the control sample by generalized cross-validation on a fixed grid, capped in effective degrees of freedom; no data splitting and no random numbers are involved.
3. `att()` gains `seed` and `nfolds` arguments for the new default (`method = "abw"`); both are silently ignored by `"lm_robust"`, `"lm_lin"`, and `"elnet"`.
4. `att()` now errors on an unrecognized `method` value instead of failing with an opaque internal error.
5. `att()` with the default `method = "abw"` now warns if called on an `hbalobject` whose entropy-balancing weights did not converge (`hbalobject$converged` is `FALSE`), since the estimate and standard error may then be unreliable.

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
