# Summarize an fgls fit

The summary contains the coefficient table (estimates, standard errors,
z or t statistics, p-values and confidence intervals) and the main
by-products of the estimator: the selected AR order \\k^\*\\ and AR
coefficients with the non-exogeneity share for serial correlation, and
the selected candidate covariates, the tuning parameters and the degrees
of freedom of the skedastic fit for heteroskedasticity.

## Usage

``` r
# S3 method for class 'fgls'
summary(object, level = 0.95, ...)

# S3 method for class 'summary.fgls'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- object:

  An `"fgls"` object.

- level:

  Confidence level of the intervals.

- ...:

  Not used.

- x:

  A `"summary.fgls"` object.

- digits:

  Number of significant digits.

- signif.stars:

  Print significance stars.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.fgls"`: a list with the call, the coefficient table
(`coefficients`), the confidence intervals (`conf.int`), a named list of
header items (`info`) and the notes. Its
[`print()`](https://rdrr.io/r/base/print.html) method returns its
argument invisibly.

## See also

Other fgls:
[`compare_fgls()`](https://emiliogcs.github.io/fgls/reference/compare_fgls.md),
[`fgls()`](https://emiliogcs.github.io/fgls/reference/fgls.md),
[`fgls-methods`](https://emiliogcs.github.io/fgls/reference/fgls-methods.md),
[`fgls_control()`](https://emiliogcs.github.io/fgls/reference/fgls_control.md),
[`lasso_fgls()`](https://emiliogcs.github.io/fgls/reference/lasso_fgls.md),
[`plot.fgls()`](https://emiliogcs.github.io/fgls/reference/plot.fgls.md),
[`serial_fgls()`](https://emiliogcs.github.io/fgls/reference/serial_fgls.md)

## Examples

``` r
sim <- simulate_serial(200, seed = 1)
s <- summary(serial_fgls(y ~ x, data = sim))
s
#> 
#> Call:
#> serial_fgls(formula = y ~ x, data = sim)
#> 
#> Method:               FGLS (serial)
#> Covariance:           robust
#> Observations:         200 (effective 199)
#> Inference:            normal
#> AR order k*:          1 [BIC, kmax = 12]
#> AR coefficients:      0.8163
#> Non-exogeneity share: 0.010
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   1.2004     0.3582   3.352 0.000803 ***
#> x             0.9804     0.0370  26.494  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Notes:
#> [1] Covariance: PG26 eq. (29), robust to non-exogenous (pre-determined)
#>     regressors.
#> [2] Collinear lagged columns dropped from the Durbin regression:
#>     L1.(Intercept).
#> 
s$coefficients
#>              Estimate Std. Error   z value      Pr(>|z|)
#> (Intercept) 1.2004371 0.35816116  3.351667  8.032643e-04
#> x           0.9803696 0.03700277 26.494493 1.121729e-154
```
