# Benchmark estimators

Estimators used as comparisons in the papers and in
[`compare_fgls()`](https://emiliogcs.github.io/fgls/reference/compare_fgls.md):

## Usage

``` r
ols_robust(
  formula,
  data,
  subset,
  na.action = stats::na.fail,
  vcov = c("HAC", "NW", "HC0", "HC1", "HC2", "HC3", "nonrobust"),
  kernel = NULL,
  bandwidth = NULL,
  small_sample = FALSE,
  use_t = FALSE,
  time = NULL,
  allow_gaps = FALSE,
  model = TRUE
)

cochrane_orcutt(
  formula,
  data,
  subset,
  na.action = stats::na.fail,
  ar_order = 1L,
  iterate = FALSE,
  max_iter = 100L,
  tol = 1e-08,
  use_t = FALSE,
  time = NULL,
  allow_gaps = FALSE,
  model = TRUE
)

durbin_estimator(
  formula,
  data,
  subset,
  na.action = stats::na.fail,
  ar_order = NULL,
  kmax = NULL,
  use_t = FALSE,
  time = NULL,
  allow_gaps = FALSE,
  model = TRUE
)

infeasible_gls(
  formula,
  data,
  subset,
  na.action = stats::na.fail,
  rho = NULL,
  omega = NULL,
  use_t = FALSE,
  time = NULL,
  allow_gaps = FALSE,
  model = TRUE
)

hybrid_gls(
  formula,
  data,
  subset,
  na.action = stats::na.fail,
  alpha = 0.05,
  lags = 1L,
  seed = NULL,
  use_t = FALSE,
  time = NULL,
  allow_gaps = FALSE,
  model = TRUE
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

- vcov:

  Covariance estimator for `ols_robust()`.

- kernel, bandwidth:

  Optional overrides for the HAC estimator: kernel `"qs"` or
  `"bartlett"`, bandwidth `"andrews"`, `"nw"` or a positive number.

- small_sample:

  Multiply HAC covariances by \\T/(T-q)\\.

- use_t:

  Use Student-t (`n_eff - q` degrees of freedom) instead of normal
  inference.

- time:

  Time index for data frames: the name of a `Date`, `POSIXct` or numeric
  column of `data`, or a vector with one value per row. It is used to
  check that rows are sorted and have no gaps.

- allow_gaps:

  Accept a time index with gaps (lags are then positional).

- model:

  Keep the model frame in the result.

- ar_order:

  AR order (`cochrane_orcutt()`, `durbin_estimator()`).

- iterate:

  Iterate Cochrane-Orcutt to convergence.

- max_iter, tol:

  Iteration limit and convergence tolerance.

- kmax:

  Largest AR order considered: `NULL` (default, \\\lfloor 12
  (T/200)^{1/4} \rfloor\\), an integer, or `"theory"` (\\\lfloor (T/\ln
  T)^{1/3} \rfloor\\). It is capped so that the Durbin regression uses
  at most half of the sample as regressors.

- rho, omega:

  Known AR coefficients or known error covariance matrix.

- alpha:

  Level of the Breusch-Godfrey test (`hybrid_gls()`).

- lags:

  Number of lags of the Breusch-Godfrey test.

- seed:

  Optional seed for the cross-validation folds and the bootstrap. The
  session's random number generator state is restored afterwards.

## Value

An object of class `"fgls"` (see
[`fgls()`](https://emiliogcs.github.io/fgls/reference/fgls.md)).

## Details

- `ols_robust()`: ordinary least squares with heteroskedasticity and
  autocorrelation consistent standard errors (`vcov = "HAC"`: Andrews'
  1991 quadratic-spectral kernel with automatic bandwidth, the benchmark
  of PG26; `"NW"`: Newey-West Bartlett kernel with \\\lfloor 4
  (T/100)^{2/9} \rfloor\\ lags) or White-type standard errors
  (`"HC0"`–`"HC3"`).

- `cochrane_orcutt()`: two-step (or iterated) Cochrane-Orcutt FGLS. The
  two-step version estimates the AR coefficients from OLS residuals and
  is **not valid with non-exogenous regressors** (PG26 Section 3); use
  [`serial_fgls()`](https://emiliogcs.github.io/fgls/reference/serial_fgls.md)
  instead.

- `durbin_estimator()`: the coefficients on the contemporaneous
  regressors of the Durbin regression (PG26 eq. 14), AR order by BIC
  unless given.

- `infeasible_gls()`: GLS with known AR coefficients (`rho`) or a known
  error covariance matrix (`omega`, up to scale); a simulation
  benchmark.

- `hybrid_gls()`: the hybrid OLS/GLS estimator of PG26 Supplement S-1
  (exogenous regressors only): OLS when the Breusch-Godfrey test does
  not reject, otherwise quasi-differencing with a random AR(1)
  coefficient of the sign of the first-order residual autocorrelation.

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

## Examples

``` r
sim <- simulate_serial(300, seed = 1)
ols_robust(y ~ x, data = sim)
#> 
#> Call:
#> ols_robust(formula = y ~ x, data = sim)
#> 
#> OLS; covariance: HAC (QS, bw = 2.03)
#> 
#> Coefficients:
#> (Intercept)            x  
#>       1.163        1.341  
#> 
cochrane_orcutt(y ~ x, data = sim)
#> 
#> Call:
#> cochrane_orcutt(formula = y ~ x, data = sim)
#> 
#> Cochrane-Orcutt (two-step); AR order k* = 1 (cochrane-orcutt), rho = 0.6133; covariance: nonrobust
#> 
#> Coefficients:
#> (Intercept)            x  
#>       1.178        1.013  
#> 
durbin_estimator(y ~ x, data = sim)
#> 
#> Call:
#> durbin_estimator(formula = y ~ x, data = sim)
#> 
#> Durbin; AR order k* = 1 (durbin), rho = 0.7654; covariance: nonrobust
#> 
#> Coefficients:
#> (Intercept)            x  
#>      0.2782       0.9938  
#> 
infeasible_gls(y ~ x, data = sim, rho = 0.8)
#> 
#> Call:
#> infeasible_gls(formula = y ~ x, data = sim, rho = 0.8)
#> 
#> GLS (infeasible); AR order k* = 1 (known), rho =   0.8; covariance: nonrobust
#> 
#> Coefficients:
#> (Intercept)            x  
#>      1.1915       0.9827  
#> 
```
