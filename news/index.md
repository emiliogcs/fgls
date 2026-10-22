# Changelog

## fgls 0.1.0

First release. A port of the Python package ‘pyfgls’ (0.1.0), validated
against it (golden fixtures and a cross-language comparison on 4,000
simulated data sets).

### Estimators

- [`fgls()`](https://emiliogcs.github.io/fgls/reference/fgls.md):
  feasible GLS correcting for serial correlation (Perron and
  González-Coya, 2026), heteroskedasticity of unknown form
  (González-Coya and Perron, 2024) or both (Supplement S-4, “F-CH”),
  with a formula interface and `ts`/`zoo`/`xts`/time-indexed data
  frames.
- [`serial_fgls()`](https://emiliogcs.github.io/fgls/reference/serial_fgls.md):
  Durbin-regression AR filter with BIC/AIC lag selection on a common
  sample, quasi-differencing, robust (eqs. 28-29) and naive (eq. 27)
  variances, HC variants, Wallis recovery for lagged dependent
  variables, and the GMA estimator for MA(1) errors.
- [`lasso_fgls()`](https://emiliogcs.github.io/fgls/reference/lasso_fgls.md):
  adaptive-Lasso skedastic function (ridge first step, joint
  cross-validation of the penalty and the weight exponent with common
  folds, exact Lasso path), SVR and Romano-Wolf WLS-S1/S2 alternatives,
  known variances, and the bias-corrected sandwich of eqs. (6)-(8).
- Bootstrap covariances (sieve and wild) as a cross-check.

### Benchmarks and tools

- [`ols_robust()`](https://emiliogcs.github.io/fgls/reference/ols_robust.md)
  (Andrews quadratic-spectral HAC, Newey-West, HC0-HC3),
  [`cochrane_orcutt()`](https://emiliogcs.github.io/fgls/reference/ols_robust.md),
  [`durbin_estimator()`](https://emiliogcs.github.io/fgls/reference/ols_robust.md),
  [`infeasible_gls()`](https://emiliogcs.github.io/fgls/reference/ols_robust.md),
  [`hybrid_gls()`](https://emiliogcs.github.io/fgls/reference/ols_robust.md).
- [`compare_fgls()`](https://emiliogcs.github.io/fgls/reference/compare_fgls.md),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and the other
  usual methods; works with
  [`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html)
  and the broom generics.
- Diagnostics:
  [`test_predetermined()`](https://emiliogcs.github.io/fgls/reference/test_predetermined.md)
  (experimental),
  [`adf_test()`](https://emiliogcs.github.io/fgls/reference/adf_test.md),
  [`weight_diagnostics()`](https://emiliogcs.github.io/fgls/reference/weight_diagnostics.md),
  [`ar_roots()`](https://emiliogcs.github.io/fgls/reference/weight_diagnostics.md),
  [`nonexogeneity_share()`](https://emiliogcs.github.io/fgls/reference/weight_diagnostics.md).
- Simulation designs of both papers:
  [`simulate_serial()`](https://emiliogcs.github.io/fgls/reference/simulate_serial.md),
  [`simulate_hetero()`](https://emiliogcs.github.io/fgls/reference/simulate_serial.md),
  [`simulate_combined()`](https://emiliogcs.github.io/fgls/reference/simulate_serial.md),
  [`simulate_predictive()`](https://emiliogcs.github.io/fgls/reference/simulate_serial.md).
