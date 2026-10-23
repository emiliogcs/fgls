# Feasible GLS for serially correlated errors (Perron and González-Coya, 2026)

Shortcut for [`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md)
with `hetero = "none"`: the errors are approximated by an AR(k) process
whose coefficients come from the Durbin regression and whose order is
chosen by BIC (eq. 14), the data are quasi-differenced and OLS is
applied (eq. 15). The default variance is the exogeneity-robust
estimator of eqs. (28)–(29), valid when the regressors are only
pre-determined.

## Usage

``` r
serial_fgls(
  formula,
  data,
  subset,
  na.action = stats::na.fail,
  serial = c("auto", "fixed", "ma1"),
  vcov = c("robust", "naive", "HC0", "HC1", "HC2", "HC3", "bootstrap"),
  ar_order = NULL,
  kmax = NULL,
  ic = c("bic", "aic"),
  lagged_dep_vars = NULL,
  exog = c("predetermined", "exogenous"),
  ...
)
```

## Arguments

- formula:

  A model formula, as in
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html). Include an intercept
  unless the model truly has none.

- data:

  A data frame, a list, a `ts`/`mts` object, or a 'zoo'/'xts' object
  containing the variables in the formula. Rows must be in time order
  (unless `serial = "none"`).

- subset:

  Optional expression selecting rows, evaluated in `data`. For time
  series without a time index (`time`, or a ts/zoo `data`) the selected
  rows must be an increasing block of consecutive rows (e.g.
  `subset = 51:200`): repeated or reversed rows are an error, and so are
  skipped rows unless `allow_gaps = TRUE`.

- na.action:

  [`stats::na.fail()`](https://rdrr.io/r/stats/na.fail.html) (default)
  stops on missing values;
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html) or
  [`stats::na.exclude()`](https://rdrr.io/r/stats/na.fail.html) drop
  them. For time series only leading and trailing missing values can be
  dropped. With `na.exclude`,
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) are padded
  with `NA` to the rows selected by `subset`, as for
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html).

- serial:

  `"auto"` (default), `"fixed"` or `"ma1"` (GMA estimator for MA(1)
  errors, PG26 Section 5; its only variance is `vcov = "robust"`, the
  GMA covariance \\s^2 (X^{\*\prime} X^\*)^{-1}\\).

- vcov:

  `"robust"` (default, eq. 29), `"naive"` (eq. 27, exogenous regressors
  only), `"HC0"`–`"HC3"` (eq. 29 with heteroskedasticity weights) or
  `"bootstrap"` (sieve bootstrap).

- ar_order:

  Fixed AR order (implies `serial = "fixed"`).

- kmax:

  Largest AR order considered: `NULL` (default, \\\lfloor 12
  (T/200)^{1/4} \rfloor\\), an integer, or `"theory"` (\\\lfloor (T/\ln
  T)^{1/3} \rfloor\\). It is capped so that the Durbin regression uses
  at most half of the sample as regressors.

- ic:

  Lag-order criterion: `"bic"` (default) or `"aic"`.

- lagged_dep_vars:

  Names (or positions) of model-matrix columns that are lags of the
  dependent variable (not available with `serial = "ma1"`).

- exog:

  Assumption on the regressors: `"predetermined"` (default) or
  `"exogenous"`. Estimates do not depend on it; it sets the default of
  `first_stage_correction` and the notes of the summary.

- ...:

  Further arguments passed to
  [`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md), e.g.
  `wallis_regressor`, `use_t`, `time`, `control`.

## Value

An object of class `"fgls"`; see
[`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md).

## References

Perron, P. and González-Coya, E. (2026). Feasible GLS for time series
regression. *Journal of Econometrics*, 258, 106346.
[doi:10.1016/j.jeconom.2026.106346](https://doi.org/10.1016/j.jeconom.2026.106346)

González-Coya, E. and Perron, P. (2024). Estimation in the presence of
heteroskedasticity of unknown form: a Lasso-based approach. *Journal of
Econometric Methods*, 13(1), 29–48.
[doi:10.1515/jem-2023-0007](https://doi.org/10.1515/jem-2023-0007)

Romano, J. P. and Wolf, M. (2017). Resurrecting weighted least squares.
*Journal of Econometrics*, 197(1), 1–19.
[doi:10.1016/j.jeconom.2016.10.003](https://doi.org/10.1016/j.jeconom.2016.10.003)

## See also

Other fgls:
[`compare_fgls()`](https://emiliogcs.com/fgls/reference/compare_fgls.md),
[`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md),
[`fgls-methods`](https://emiliogcs.com/fgls/reference/fgls-methods.md),
[`fgls_control()`](https://emiliogcs.com/fgls/reference/fgls_control.md),
[`lasso_fgls()`](https://emiliogcs.com/fgls/reference/lasso_fgls.md),
[`plot.fgls()`](https://emiliogcs.com/fgls/reference/plot.fgls.md),
[`summary.fgls()`](https://emiliogcs.com/fgls/reference/summary.fgls.md)

## Examples

``` r
sim <- simulate_serial(500, seed = 42)
fit <- serial_fgls(y ~ x, data = sim)
fit
#> 
#> Call:
#> serial_fgls(formula = y ~ x, data = sim)
#> 
#> FGLS (serial); AR order k* = 1 (BIC), rho = 0.7669; covariance: robust
#> 
#> Coefficients:
#> (Intercept)            x  
#>      0.8491       0.9918  
#> 
fit$ar_order
#> [1] 1
fit$nonexogeneity_share
#> [1] 0.009554577

# Lagged dependent variable among the regressors (Wallis recovery)
df <- data.frame(y = sim$y[-1], ylag = sim$y[-500], x = sim$x[-1])
serial_fgls(y ~ ylag + x, data = df, lagged_dep_vars = "ylag")
#> 
#> Call:
#> serial_fgls(formula = y ~ ylag + x, data = df, lagged_dep_vars = "ylag")
#> 
#> FGLS (serial); AR order k* = 1 (BIC), rho = 0.7837; covariance: robust
#> 
#> Coefficients:
#> (Intercept)         ylag            x  
#>     0.83336      0.01755      0.99401  
#> 
```
