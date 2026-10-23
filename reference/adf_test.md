# Augmented Dickey-Fuller test (screening)

The regression \\\Delta x_t = a + g x\_{t-1} + \sum_j c_j \Delta
x\_{t-j} + e_t\\ with the number of lags chosen by BIC (at most
\\\lfloor 12 (T/100)^{1/4} \rfloor\\); returns the t statistic on \\g\\
with MacKinnon's asymptotic critical values. Used by
`fgls(check_stationarity = TRUE)` to screen the response and the
regressors; it is a screening device, not a replacement for a careful
unit-root analysis.

## Usage

``` r
adf_test(x, max_lags = NULL)

# S3 method for class 'fgls_adf'
print(x, ...)
```

## Arguments

- x:

  A numeric vector (a time series).

- max_lags:

  Largest number of lagged differences.

- ...:

  Not used.

## Value

A list of class `"fgls_adf"` with `statistic`, `lags`, `nobs`,
`critical_values` and `unit_root_rejected` (at 5%). The print method
returns `x` invisibly.

## See also

Other diagnostics:
[`test_predetermined()`](https://emiliogcs.com/fgls/reference/test_predetermined.md),
[`weight_diagnostics()`](https://emiliogcs.com/fgls/reference/weight_diagnostics.md)

## Examples

``` r
set.seed(1)
adf_test(cumsum(rnorm(200)))$unit_root_rejected
#> [1] FALSE
adf_test(rnorm(200))$unit_root_rejected
#> [1] TRUE
```
