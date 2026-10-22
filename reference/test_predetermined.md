# Test whether lags of the regressors are missing (pre-determinedness)

**Experimental.** A variable-addition LM test in the spirit of Remark 9
of Perron and González-Coya (2026): FGLS needs the regressors to be
pre-determined (uncorrelated with current and future innovations), which
can fail when lags of the regressors belong in the model. The steps are
(a) the Durbin regression (AR order by BIC) giving \\\hat\beta^D\\; (b)
the residuals \\\hat u^D = y - X\hat\beta^D\\, demeaned; (c) for \\l =
1, \dots,\\ `max_lag`, the statistic \\LM_l = n (SSR_0 - SSR_1)/SSR_0\\,
asymptotically \\\chi^2\\ with as many degrees of freedom as added lags,
from the auxiliary regressions without and with lags \\1..l\\ of the
stochastic regressors.

## Usage

``` r
test_predetermined(
  formula,
  data,
  max_lag = 4L,
  ar_order = NULL,
  kmax = NULL,
  regressors = NULL,
  alpha = 0.05,
  scale = c("innovations", "errors"),
  na.action = stats::na.fail,
  time = NULL,
  allow_gaps = FALSE
)

# S3 method for class 'fgls_predetermined'
print(x, ...)
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

- max_lag:

  Largest number of lags added.

- ar_order:

  Fixed AR order (implies `serial = "fixed"`).

- kmax:

  Largest AR order considered: `NULL` (default, \\\lfloor 12
  (T/200)^{1/4} \rfloor\\), an integer, or `"theory"` (\\\lfloor (T/\ln
  T)^{1/3} \rfloor\\). It is capped so that the Durbin regression uses
  at most half of the sample as regressors.

- regressors:

  Names of the regressors whose lags are tested (default: all
  non-constant regressors).

- alpha:

  Level used for the advice.

- scale:

  `"innovations"` (default) or `"errors"`.

- na.action:

  [`stats::na.fail()`](https://rdrr.io/r/stats/na.fail.html) (default)
  stops on missing values;
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html) or
  [`stats::na.exclude()`](https://rdrr.io/r/stats/na.fail.html) drop
  them. For time series only leading and trailing missing values can be
  dropped.

- time:

  Time index for data frames: the name of a `Date`, `POSIXct` or numeric
  column of `data`, or a vector with one value per row. It is used to
  check that rows are sorted and have no gaps.

- allow_gaps:

  Accept a time index with gaps (lags are then positional).

- x:

  A `"fgls_predetermined"` object.

- ...:

  Not used.

## Value

An object of class `"fgls_predetermined"`: a list with `table` (a data
frame with columns `lags`, `LM`, `df`, `p.value`), `ar_order`,
`regressors`, `alpha`, `rejected` and `advice`.

## Details

With `scale = "innovations"` (default) the auxiliary regression is run
on the innovation scale: \\\hat u^D\\ and all regressors are filtered
with the Durbin AR filter and lagged \\\hat u^D\\ enter as nuisance
regressors, which keeps the \\\chi^2\\ approximation valid with
pre-determined but non-exogenous regressors. `scale = "errors"` uses
\\\hat u^D\\ and the raw regressors (the literal reading of Remark 9).

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

Other diagnostics:
[`adf_test()`](https://emiliogcs.github.io/fgls/reference/adf_test.md),
[`weight_diagnostics()`](https://emiliogcs.github.io/fgls/reference/weight_diagnostics.md)

## Examples

``` r
sim <- simulate_serial(300, seed = 1)
test_predetermined(y ~ x, data = sim)
#> Warning: test_predetermined() is experimental; its interface and details may change.
#> Pre-determinedness LM test (PG26 Remark 9, experimental)
#> Durbin AR order: 1; regressors tested: x; scale: innovations
#> 
#>  lags      LM df p.value
#>     1 0.05382  1  0.8165
#>     2 0.94706  2  0.6228
#>     3 1.16658  3  0.7610
#>     4 1.18574  4  0.8804
#> 
#> No rejection at 5%: no evidence of missing regressor lags.
```
