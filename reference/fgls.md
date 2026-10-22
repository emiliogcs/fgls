# Feasible GLS for serial correlation and/or heteroskedasticity

`fgls()` fits a linear regression by feasible generalized least squares,
correcting for serially correlated errors (Perron and González-Coya,
2026), for heteroskedasticity of unknown form (González-Coya and Perron,
2024), or for both (the "F-CH" estimator of Supplement S-4 of the 2026
paper).
[`serial_fgls()`](https://emiliogcs.github.io/fgls/reference/serial_fgls.md)
and
[`lasso_fgls()`](https://emiliogcs.github.io/fgls/reference/lasso_fgls.md)
are shortcuts for the two single corrections.

## Usage

``` r
fgls(
  formula,
  data,
  subset,
  na.action = stats::na.fail,
  serial = c("auto", "fixed", "ma1", "none"),
  hetero = c("lasso", "svr", "wls_s1", "wls_s2", "known", "none"),
  vcov = c("robust", "naive", "nonrobust", "HC0", "HC1", "HC2", "HC3", "bootstrap"),
  z = "auto",
  ar_order = NULL,
  kmax = NULL,
  ic = c("bic", "aic"),
  lagged_dep_vars = NULL,
  wallis_regressor = "auto",
  exog = c("predetermined", "exogenous"),
  first_stage_correction = NULL,
  variance = NULL,
  delta = 0.1,
  psi_grid = c(0, 0.25, 0.5, 0.75, 1, 2),
  cv_folds = 10L,
  lambda_rule = c("min", "1se"),
  weight_clip = NULL,
  seed = NULL,
  use_t = FALSE,
  time = NULL,
  allow_gaps = FALSE,
  check_stationarity = FALSE,
  control = fgls_control(),
  contrasts = NULL,
  model = TRUE,
  weights,
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
  series only a contiguous block makes sense.

- na.action:

  [`stats::na.fail()`](https://rdrr.io/r/stats/na.fail.html) (default)
  stops on missing values;
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html) or
  [`stats::na.exclude()`](https://rdrr.io/r/stats/na.fail.html) drop
  them. For time series only leading and trailing missing values can be
  dropped.

- serial:

  Serial-correlation correction: `"auto"` (AR order by BIC, default),
  `"fixed"` (use `ar_order`), `"ma1"` (GMA for MA(1) errors) or
  `"none"`.

- hetero:

  Skedastic-function estimator: `"lasso"` (adaptive Lasso, default),
  `"svr"` (support vector regression, requires 'e1071'),
  `"wls_s1"`/`"wls_s2"` (the parametric fits of Romano and Wolf, 2017,
  on \\\log z\\ or \\z\\), `"known"` (variances given in `variance`) or
  `"none"`.

- vcov:

  Variance estimator: `"robust"` (default; PG26 eq. 29, GP24 eqs. 6+8,
  or Supplement A.4–A.5 depending on the corrections), `"naive"` (PG26
  eq. 27; serial only), `"nonrobust"` (heteroskedasticity only),
  `"HC0"`–`"HC3"`, or `"bootstrap"`. See
  [`vcov.fgls()`](https://emiliogcs.github.io/fgls/reference/fgls-methods.md).

- z:

  Candidate covariates of the skedastic function: `"auto"` (default)
  builds \\\\x, x^2, \log x, (\log x)^2, \cos x, \cos 2x\\\\ from every
  non-constant regressor (see
  [`make_dictionary()`](https://emiliogcs.github.io/fgls/reference/make_dictionary.md));
  a character vector of regressor names; a one-sided formula evaluated
  in `data` (e.g. `~ x + I(x^2) + w`); or a numeric matrix / data frame
  with one row per row of `data`.

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
  dependent variable.

- wallis_regressor:

  Regressor used for the Wallis recovery of the AR coefficients:
  `"auto"` (largest Durbin \\\|t\|\\), `"median"`, or a name.

- exog:

  Assumption on the regressors: `"predetermined"` (default) or
  `"exogenous"`. Estimates do not depend on it; it sets the default of
  `first_stage_correction` and the notes of the summary.

- first_stage_correction:

  Add the first-stage (Durbin) term of PG26 eq. (29) to the variance of
  the combined estimator. Default `TRUE` when `exog = "predetermined"`.

- variance:

  Known error variances for `hetero = "known"`: a numeric vector with
  one value per observation, or a function of the model matrix.

- delta:

  Floor in \\\log\max(\hat u_t^2, \delta^2)\\, in the units of the
  response (default 0.1, as in the paper).

- psi_grid:

  Exponents \\\psi\\ of the adaptive weights; 0 is the plain Lasso.

- cv_folds:

  Number of cross-validation folds (default 10).

- lambda_rule:

  Penalty at the minimum of the cross-validation error (`"min"`,
  default) or the largest penalty within one standard error of it
  (`"1se"`).

- weight_clip:

  Optional quantiles `c(lo, hi)` at which the fitted log-variance is
  winsorized, e.g. `c(0.01, 0.99)`.

- seed:

  Optional seed for the cross-validation folds and the bootstrap. The
  session's random number generator state is restored afterwards.

- use_t:

  Use Student-t (`n_eff - q` degrees of freedom) instead of normal
  inference.

- time:

  Time index for data frames: the name of a `Date`, `POSIXct` or numeric
  column of `data`, or a vector with one value per row. It is used to
  check that rows are sorted and have no gaps.

- allow_gaps:

  Accept a time index with gaps (lags are then positional).

- check_stationarity:

  Screen the response and the regressors with augmented Dickey-Fuller
  tests and warn when a unit root is not rejected.

- control:

  Advanced options, see
  [`fgls_control()`](https://emiliogcs.github.io/fgls/reference/fgls_control.md).

- contrasts:

  Optional list of contrasts for factors, as in
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html).

- model:

  Keep the model frame in the result.

- weights:

  Not supported: FGLS estimates its own weights. Use `hetero = "known"`
  with `variance` for known variances.

- ...:

  Not used; unknown arguments are an error.

## Value

An object of class `"fgls"`, a list with components including

- `coefficients`, `vcov`:

  Estimates and their covariance matrix.

- `residuals`, `fitted.values`:

  \\y_t - x_t'\hat\beta\\ and \\x_t'\hat\beta\\ on the full sample.

- `innovations`:

  Residuals of the final (filtered and/or weighted) regression.

- `ar_order`, `rho`, `bic_path`, `nonexogeneity_share`:

  Serial block (`NULL` without serial correction); more in `serial`.

- `skedastic`, `weights`:

  Heteroskedasticity block: selected covariates, coefficients, \\\psi\\,
  \\\lambda\\, degrees of freedom, fitted variances; `weights` are
  \\1/\hat v_t\\.

- `stages`:

  Intermediate estimators (`stages$serial` is F-C).

- `nobs`, `n_eff`:

  Number of observations and of observations in the final regression
  (\\T - k^\*\\).

Use the methods [`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`compare_fgls()`](https://emiliogcs.github.io/fgls/reference/compare_fgls.md).

## Serial correlation (PG26)

The errors are approximated by an AR(k) process. Its coefficients are
the coefficients on the lagged dependent variable in the Durbin
regression of \\y_t\\ on \\y\_{t-1}, \dots, y\_{t-k}\\, \\x_t, \dots,
x\_{t-k}\\ (eq. 14), and \\k\\ is chosen by BIC over \\0, \dots,
k\_{max}\\ on a common sample. The data are then quasi-differenced and
the model is estimated by OLS on the filtered data (eq. 15), dropping
the first \\k\\ observations. Unlike Cochrane-Orcutt, this is consistent
when the regressors are only *pre-determined* (uncorrelated with current
and future innovations) rather than strictly exogenous. The default
variance (`vcov = "robust"`, eqs. 28–29) accounts for the estimation of
the filter and is valid in both cases; `vcov = "naive"` (eq. 27) is
valid only for exogenous regressors. The share by which the robust
variance exceeds the naive one is reported as `nonexogeneity_share`. If
the selected order is 0, FGLS equals OLS.

## Heteroskedasticity (GP24)

The log of the squared OLS residuals (floored at `delta^2`) is regressed
on candidate covariates \\z_t\\ by an adaptive Lasso: a ridge first step
gives weights \\\|\hat\gamma_j\|^{-\psi}\\, and the penalty \\\lambda\\
and the exponent \\\psi\\ are chosen jointly by 10-fold cross-validation
with the same folds for every \\\psi\\. The exponentiated fit gives the
weights of a weighted least-squares regression. Inference uses a
sandwich variance with the bias correction of eq. (8), so it stays valid
when the skedastic model is wrong. If the Lasso selects no covariate the
weights are constant and FGLS equals OLS (no loss under
homoskedasticity).

## Both (Supplement S-4)

With `serial != "none"` and `hetero != "none"` the serial correction is
applied first ("F-C"), the skedastic function is fitted to the residuals
of the quasi-differenced regression with candidate covariates built from
the *unfiltered* regressors, and the quasi-differenced regression is
re-estimated by weighted least squares ("F-CH"). The F-C stage is kept
in `fit$stages$serial`.

## Choosing options

- Time series with possibly autocorrelated errors and no concern about
  heteroskedasticity: `hetero = "none"` (or
  [`serial_fgls()`](https://emiliogcs.github.io/fgls/reference/serial_fgls.md)).

- Cross-sectional data: `serial = "none"` (or
  [`lasso_fgls()`](https://emiliogcs.github.io/fgls/reference/lasso_fgls.md)).
  Rows can be in any order.

- Lagged dependent variables among the regressors must be declared in
  `lagged_dep_vars`; the AR coefficients are then recovered by the
  Wallis ratio (PG26 Section 4.2).

- Errors known to be MA(1), e.g. multi-step forecast errors:
  `serial = "ma1"` (GMA estimator, PG26 Section 5).

- With several stochastic regressors the BIC that penalizes only the AR
  order can select `kmax`; use
  `control = fgls_control(ic_penalty = "params")`.

- Results involving cross-validation depend on the random folds: set
  `seed` for reproducibility.

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
[`compare_fgls()`](https://emiliogcs.github.io/fgls/reference/compare_fgls.md),
[`fgls-methods`](https://emiliogcs.github.io/fgls/reference/fgls-methods.md),
[`fgls_control()`](https://emiliogcs.github.io/fgls/reference/fgls_control.md),
[`lasso_fgls()`](https://emiliogcs.github.io/fgls/reference/lasso_fgls.md),
[`plot.fgls()`](https://emiliogcs.github.io/fgls/reference/plot.fgls.md),
[`serial_fgls()`](https://emiliogcs.github.io/fgls/reference/serial_fgls.md),
[`summary.fgls()`](https://emiliogcs.github.io/fgls/reference/summary.fgls.md)

## Examples

``` r
# Serial correlation with a non-exogenous regressor (PG26 Table 1 design)
sim <- simulate_serial(300, seed = 1)
fit <- fgls(y ~ x, data = sim, hetero = "none")
summary(fit)
#> 
#> Call:
#> fgls(formula = y ~ x, data = sim, hetero = "none")
#> 
#> Method:               FGLS (serial)
#> Covariance:           robust
#> Observations:         300 (effective 299)
#> Inference:            normal
#> AR order k*:          1 [BIC, kmax = 13]
#> AR coefficients:      0.7654
#> Non-exogeneity share: 0.012
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.18739    0.23736   5.002 5.66e-07 ***
#> x            0.98640    0.03202  30.809  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Notes:
#> [1] Covariance: PG26 eq. (29), robust to non-exogenous (pre-determined)
#>     regressors.
#> [2] Collinear lagged columns dropped from the Durbin regression:
#>     L1.(Intercept).
#> 

# Heteroskedasticity of unknown form (GP24)
het <- simulate_hetero(200, skedastic = "log2", seed = 2)
fit2 <- fgls(y ~ x, data = het, serial = "none", seed = 1)
fit2$skedastic$selected
#> [1] "log(x)" "cos(x)"

# Both (Supplement S-4)
both <- simulate_combined(200, ar = 0.5, skedastic = "x2", seed = 3)
fit3 <- fgls(y ~ x, data = both, seed = 1)
compare_fgls(fit3$stages$serial, fit3)
#>          term statistic F-C (serial) F-CH (serial+hetero)
#> 1 (Intercept)  estimate    1.3723275            1.5934094
#> 2 (Intercept) std.error    0.5413481            0.4430785
#> 3 (Intercept)  ci.width    2.1220455            1.7368357
#> 4           x  estimate    0.9130718            0.8100948
#> 5           x std.error    0.2611294            0.2123379
#> 6           x  ci.width    1.0236086            0.8323494
```
