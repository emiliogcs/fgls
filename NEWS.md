# fgls 0.1.0

First release. A port of the Python package 'pyfgls' (0.1.0), validated against it
(golden fixtures and a cross-language comparison on 4,000 simulated data sets).

## Estimators

* `fgls()`: feasible GLS correcting for serial correlation (Perron and
  González-Coya, 2026), heteroskedasticity of unknown form (González-Coya and
  Perron, 2024) or both (Supplement S-4, "F-CH"), with a formula interface and
  `ts`/`zoo`/`xts`/time-indexed data frames.
* `serial_fgls()`: Durbin-regression AR filter with BIC/AIC lag selection on a
  common sample, quasi-differencing, robust (eqs. 28-29) and naive (eq. 27)
  variances, HC variants, Wallis recovery for lagged dependent variables, and the
  GMA estimator for MA(1) errors.
* `lasso_fgls()`: adaptive-Lasso skedastic function (ridge first step, joint
  cross-validation of the penalty and the weight exponent with common folds,
  exact Lasso path), SVR and Romano-Wolf WLS-S1/S2 alternatives, known variances,
  and the bias-corrected sandwich of eqs. (6)-(8).
* Bootstrap covariances (sieve and wild) as a cross-check.

## Benchmarks and tools

* `ols_robust()` (Andrews quadratic-spectral HAC, Newey-West, HC0-HC3),
  `cochrane_orcutt()`, `durbin_estimator()`, `infeasible_gls()`,
  `hybrid_gls()`.
* `compare_fgls()`, `summary()`, `confint()`, `predict()`, `plot()` and the
  other usual methods; works with `lmtest::coeftest()` and the broom generics.
* Diagnostics: `test_predetermined()` (experimental), `adf_test()`,
  `weight_diagnostics()`, `ar_roots()`, `nonexogeneity_share()`.
* Simulation designs of both papers: `simulate_serial()`, `simulate_hetero()`,
  `simulate_combined()`, `simulate_predictive()`.

## Differences from pyfgls 0.1.0

Estimates and default variances agree with pyfgls (see `vignette("validation")`).
Elsewhere the R package deliberately differs:

* `nonexogeneity_share` is computed from the robust variance (PG26 eqs. 28-29)
  whatever `vcov` is requested (pyfgls reports 0 for non-robust types), and is
  `NA` for the benchmark estimators and the GMA estimator, which do not compute
  it.
* `durbin_estimator()` recovers the coefficients on regressors whose lags are
  collinear with the regressors (constant, trends, periodic dummies) as
  \((I - \sum_j \rho_j P_j)^{-1} c\), e.g. \(c / (1 - \sum_j \rho_j)\) for the
  intercept; pyfgls reports the Durbin coefficient \(c\).
* The `kmax` cap and the `ar_order` check count the Durbin regressors that are
  actually kept, so dropped lags of trends and seasonal dummies no longer shrink
  `kmax`. Without such columns the cap is the pyfgls one.
* `serial = "ma1"` rejects `ar_order`, `kmax`, `lagged_dep_vars` and
  `vcov = "bootstrap"` instead of ignoring them.
* Formula-interface behavior with no pyfgls counterpart: `subset`, `na.action`
  (including `na.exclude` padding), row-numbered error messages, and checks of
  time-series row order and time indexes.
