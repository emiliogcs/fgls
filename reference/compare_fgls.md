# Compare fits side by side

A table of coefficients, standard errors and confidence-interval widths
for several fits, e.g. FGLS against OLS with HAC standard errors (PG26)
or OLS with HC3 standard errors (GP24), or the F-C and F-CH stages of a
combined fit.

## Usage

``` r
compare_fgls(..., level = 0.95, ols_vcov = NULL)
```

## Arguments

- ...:

  `"fgls"` (or `"lm"`) objects, optionally named; the names label the
  columns.

- level:

  Confidence level used for the interval widths.

- ols_vcov:

  Covariance type for `lm` fits: `NULL` (automatic), `"HAC"`, `"NW"` or
  `"HC0"`–`"HC3"`.

## Value

A data frame with columns `term`, `statistic` (`"estimate"`,
`"std.error"`, `"ci.width"`) and one column per fit.

## Details

Plain [`stats::lm()`](https://rdrr.io/r/stats/lm.html) fits are accepted
and converted to OLS with HAC (Andrews quadratic-spectral) standard
errors when any other fit corrects for serial correlation, and with HC3
standard errors otherwise; override with `ols_vcov`.

## See also

Other fgls:
[`fgls()`](https://emiliogcs.github.io/fgls/reference/fgls.md),
[`fgls-methods`](https://emiliogcs.github.io/fgls/reference/fgls-methods.md),
[`fgls_control()`](https://emiliogcs.github.io/fgls/reference/fgls_control.md),
[`lasso_fgls()`](https://emiliogcs.github.io/fgls/reference/lasso_fgls.md),
[`plot.fgls()`](https://emiliogcs.github.io/fgls/reference/plot.fgls.md),
[`serial_fgls()`](https://emiliogcs.github.io/fgls/reference/serial_fgls.md),
[`summary.fgls()`](https://emiliogcs.github.io/fgls/reference/summary.fgls.md)

## Examples

``` r
sim <- simulate_serial(300, seed = 1)
compare_fgls(
  OLS = stats::lm(y ~ x, data = sim),
  FGLS = serial_fgls(y ~ x, data = sim)
)
#>          term statistic        OLS       FGLS
#> 1 (Intercept)  estimate 1.16317984 1.18739202
#> 2 (Intercept) std.error 0.11500951 0.23736465
#> 3 (Intercept)  ci.width 0.45082901 0.93045231
#> 4           x  estimate 1.34140597 0.98639665
#> 5           x std.error 0.06851676 0.03201629
#> 6           x  ci.width 0.26858077 0.12550157
```
