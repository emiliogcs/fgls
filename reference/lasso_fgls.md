# Feasible GLS under heteroskedasticity of unknown form (González-Coya and Perron, 2024)

Shortcut for [`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md)
with `serial = "none"`. The log of the squared OLS residuals is
regressed on candidate covariates `z` by an adaptive Lasso tuned by
cross-validation; the exponentiated fit gives the weights of a weighted
least-squares regression. Inference uses the sandwich of eq. (6) with
the bias correction of eq. (8), so it remains valid when the skedastic
model is misspecified.

## Usage

``` r
lasso_fgls(
  formula,
  data,
  subset,
  na.action = stats::na.fail,
  hetero = c("lasso", "svr", "wls_s1", "wls_s2", "known", "none"),
  vcov = c("robust", "nonrobust", "HC0", "HC1", "HC2", "HC3", "bootstrap"),
  z = "auto",
  seed = NULL,
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

- hetero:

  Skedastic-function estimator: `"lasso"` (default), `"svr"`,
  `"wls_s1"`, `"wls_s2"`, `"known"` or `"none"` (OLS with HC3 standard
  errors).

- vcov:

  `"robust"` (default; eq. 8 for Lasso/SVR, eq. 7 for WLS-S1/S2, HC3
  otherwise), `"nonrobust"`, `"HC0"`–`"HC3"` or `"bootstrap"` (wild
  bootstrap).

- z:

  Candidate covariates of the skedastic function: `"auto"` (default)
  builds \\\\x, x^2, \log x, (\log x)^2, \cos x, \cos 2x\\\\ from every
  non-constant regressor (see
  [`make_dictionary()`](https://emiliogcs.com/fgls/reference/make_dictionary.md));
  a character vector of regressor names; a one-sided formula evaluated
  in `data` (e.g. `~ x + I(x^2) + w`); or a numeric matrix / data frame
  with one row per row of `data`. With `hetero = "wls_s1"` all
  candidates must be strictly positive, so give them explicitly, e.g.
  `z = "x"` or `z = ~ x + I(x^2)` for a positive regressor `x`.

- seed:

  Optional seed for the cross-validation folds and the bootstrap. The
  session's random number generator state is restored afterwards.

- ...:

  Further arguments passed to
  [`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md), e.g.
  `variance`, `cv_folds`, `lambda_rule`, `weight_clip`, `use_t`,
  `control`.

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
[`plot.fgls()`](https://emiliogcs.com/fgls/reference/plot.fgls.md),
[`serial_fgls()`](https://emiliogcs.com/fgls/reference/serial_fgls.md),
[`summary.fgls()`](https://emiliogcs.com/fgls/reference/summary.fgls.md)

Other heteroskedasticity:
[`make_dictionary()`](https://emiliogcs.com/fgls/reference/make_dictionary.md)

## Examples

``` r
het <- simulate_hetero(200, skedastic = "x2", seed = 1)
fit <- lasso_fgls(y ~ x, data = het, seed = 1)
fit
#> 
#> Call:
#> lasso_fgls(formula = y ~ x, data = het, seed = 1)
#> 
#> FGLS (hetero); skedastic fit: lasso (selected: x^2, log(x), cos(2x)); covariance: robust
#> 
#> Coefficients:
#> (Intercept)            x  
#>      0.9988       1.0157  
#> 
fit$skedastic$selected
#> [1] "x^2"     "log(x)"  "cos(2x)"
ols <- lasso_fgls(y ~ x, data = het, hetero = "none")
compare_fgls(OLS = ols, FGLS = fit)
#>          term statistic       OLS      FGLS
#> 1 (Intercept)  estimate 0.9492108 0.9988174
#> 2 (Intercept) std.error 0.5227244 0.4828889
#> 3 (Intercept)  ci.width 2.0490420 1.8928895
#> 4           x  estimate 1.0142749 1.0157196
#> 5           x std.error 0.2284271 0.2148521
#> 6           x  ci.width 0.8954176 0.8422049
```
