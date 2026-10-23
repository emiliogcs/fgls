# Package index

## Estimation

FGLS estimators and their options.

- [`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md) : Feasible
  GLS for serial correlation and/or heteroskedasticity
- [`serial_fgls()`](https://emiliogcs.com/fgls/reference/serial_fgls.md)
  : Feasible GLS for serially correlated errors (Perron and
  González-Coya, 2026)
- [`lasso_fgls()`](https://emiliogcs.com/fgls/reference/lasso_fgls.md) :
  Feasible GLS under heteroskedasticity of unknown form (González-Coya
  and Perron, 2024)
- [`fgls_control()`](https://emiliogcs.com/fgls/reference/fgls_control.md)
  : Advanced options for fgls

## Working with fits

- [`summary(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/summary.fgls.md)
  [`print(`*`<summary.fgls>`*`)`](https://emiliogcs.com/fgls/reference/summary.fgls.md)
  : Summarize an fgls fit
- [`vcov(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`confint(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`residuals(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`fitted(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`predict(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`nobs(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`df.residual(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`formula(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`model.matrix(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`model.frame(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`terms(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  [`print(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-methods.md)
  : Methods for fgls fits
- [`plot(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/plot.fgls.md)
  : Diagnostic plots for fgls fits
- [`compare_fgls()`](https://emiliogcs.com/fgls/reference/compare_fgls.md)
  : Compare fits side by side
- [`tidy(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-broom.md)
  [`glance(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-broom.md)
  [`augment(`*`<fgls>`*`)`](https://emiliogcs.com/fgls/reference/fgls-broom.md)
  : Tidy, glance and augment fgls fits

## Benchmarks

OLS with HAC/HC standard errors, Cochrane-Orcutt, Durbin, infeasible and
hybrid GLS.

- [`ols_robust()`](https://emiliogcs.com/fgls/reference/ols_robust.md)
  [`cochrane_orcutt()`](https://emiliogcs.com/fgls/reference/ols_robust.md)
  [`durbin_estimator()`](https://emiliogcs.com/fgls/reference/ols_robust.md)
  [`infeasible_gls()`](https://emiliogcs.com/fgls/reference/ols_robust.md)
  [`hybrid_gls()`](https://emiliogcs.com/fgls/reference/ols_robust.md) :
  Benchmark estimators

## Diagnostics

- [`test_predetermined()`](https://emiliogcs.com/fgls/reference/test_predetermined.md)
  [`print(`*`<fgls_predetermined>`*`)`](https://emiliogcs.com/fgls/reference/test_predetermined.md)
  : Test whether lags of the regressors are missing (pre-determinedness)
- [`adf_test()`](https://emiliogcs.com/fgls/reference/adf_test.md)
  [`print(`*`<fgls_adf>`*`)`](https://emiliogcs.com/fgls/reference/adf_test.md)
  : Augmented Dickey-Fuller test (screening)
- [`weight_diagnostics()`](https://emiliogcs.com/fgls/reference/weight_diagnostics.md)
  [`ar_roots()`](https://emiliogcs.com/fgls/reference/weight_diagnostics.md)
  [`nonexogeneity_share()`](https://emiliogcs.com/fgls/reference/weight_diagnostics.md)
  : Diagnostics of an fgls fit

## Candidate covariates

- [`make_dictionary()`](https://emiliogcs.com/fgls/reference/make_dictionary.md)
  : Candidate covariates for the skedastic function

## Simulation designs of the papers

- [`simulate_serial()`](https://emiliogcs.com/fgls/reference/simulate_serial.md)
  [`simulate_hetero()`](https://emiliogcs.com/fgls/reference/simulate_serial.md)
  [`simulate_combined()`](https://emiliogcs.com/fgls/reference/simulate_serial.md)
  [`simulate_predictive()`](https://emiliogcs.com/fgls/reference/simulate_serial.md)
  : Simulate the designs of the papers

## Package and conditions

- [`fgls-package`](https://emiliogcs.com/fgls/reference/fgls-package.md)
  : fgls: Feasible GLS for serial correlation and heteroskedasticity
- [`fgls-conditions`](https://emiliogcs.com/fgls/reference/fgls-conditions.md)
  : Conditions signaled by fgls
