# Methods for fgls fits

Standard accessors for objects returned by
[`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md),
[`serial_fgls()`](https://emiliogcs.com/fgls/reference/serial_fgls.md),
[`lasso_fgls()`](https://emiliogcs.com/fgls/reference/lasso_fgls.md) and
the benchmark estimators. They make the fits work with generic tools
such as
[`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html) and
[`lmtest::coefci()`](https://rdrr.io/pkg/lmtest/man/coeftest.html).

## Usage

``` r
# S3 method for class 'fgls'
vcov(object, type = NULL, ...)

# S3 method for class 'fgls'
confint(object, parm, level = 0.95, type = NULL, ...)

# S3 method for class 'fgls'
residuals(object, type = c("response", "innovations"), ...)

# S3 method for class 'fgls'
fitted(object, ...)

# S3 method for class 'fgls'
predict(
  object,
  newdata,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  type = "response",
  na.action = stats::na.pass,
  ...
)

# S3 method for class 'fgls'
nobs(object, ...)

# S3 method for class 'fgls'
df.residual(object, ...)

# S3 method for class 'fgls'
formula(x, ...)

# S3 method for class 'fgls'
model.matrix(object, ...)

# S3 method for class 'fgls'
model.frame(formula, ...)

# S3 method for class 'fgls'
terms(x, ...)

# S3 method for class 'fgls'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object, x, formula:

  An `"fgls"` object.

- type:

  For [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html): `NULL` (the
  variance chosen when fitting) or another variance type to compute from
  the stored estimates: `"robust"`, `"naive"`, `"nonrobust"`,
  `"HC0"`–`"HC3"` (where available) or `"bootstrap"` (re-runs the
  estimator); see the section *Variance types*. For
  [`residuals()`](https://rdrr.io/r/stats/residuals.html): `"response"`
  (\\y_t - x_t'\hat\beta\\, full sample) or `"innovations"` (residuals
  of the final filtered and/or weighted regression). For
  [`predict()`](https://rdrr.io/r/stats/predict.html): only `"response"`
  (the mean of the response) is supported.

- ...:

  Not used.

- parm:

  Coefficients (names or positions); all by default.

- level:

  Confidence level.

- newdata:

  Optional data frame in which to look for the regressors. Without it,
  predictions are for the observations of the fit, padded with `NA` for
  the rows dropped by
  [`stats::na.exclude()`](https://rdrr.io/r/stats/na.fail.html).

- se.fit:

  Return standard errors of the predicted means.

- interval:

  `"none"` or `"confidence"` (for the mean of the response).

- na.action:

  Missing-value handling for `newdata`.

- digits:

  Number of significant digits to print.

## Value

The usual values: a coefficient covariance matrix (`vcov`), a matrix of
confidence limits (`confint`), vectors of residuals or fitted values,
predictions (a vector, or a list/matrix with `se.fit`/`interval`), or
the corresponding model components. The
[`print()`](https://rdrr.io/r/base/print.html) method returns `x`
invisibly.

## Details

Inference is asymptotically normal (as in the papers) unless the fit was
made with `use_t = TRUE`. Accordingly,
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) returns
`Inf` for normal inference, so that
[`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html)
reproduces
[`summary.fgls()`](https://emiliogcs.com/fgls/reference/summary.fgls.md),
and \\n\_{eff} - q\\ otherwise.
[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the number of
observations in the final regression, \\n\_{eff} = T - k^\*\\.

## Variances and the sandwich package

The variance estimators of the papers are available through
`vcov(fit, type = ...)` and can be passed to other tools (for example
`lmtest::coeftest(fit, vcov. = vcov(fit, type = "naive"))`). The package
provides no
[`sandwich::estfun()`](https://zeileis.codeberg.page/sandwich/reference/estfun.html)
or
[`sandwich::bread()`](https://zeileis.codeberg.page/sandwich/reference/bread.html)
methods, because the robust variance of PG26 eq. (29) and the
bias-corrected sandwich of GP24 eq. (8) are not of the generic form
those functions assume:
[`sandwich::vcovHAC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHAC.html)
and related functions would *not* reproduce them.

## Variance types

`vcov(fit, type = t)` returns the stored covariance when `t` is the type
the fit was made with (e.g. `"robust"`, or `"HAC"` for
[`ols_robust()`](https://emiliogcs.com/fgls/reference/ols_robust.md),
whose stored `vcov_type` label also shows the kernel and bandwidth).
Other types are recomputed from the stored estimates only for fits made
by [`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md),
[`serial_fgls()`](https://emiliogcs.com/fgls/reference/serial_fgls.md)
and
[`lasso_fgls()`](https://emiliogcs.com/fgls/reference/lasso_fgls.md),
and for the F-C stage of a combined fit (`fit$stages$serial`, types
`"robust"`, `"HC0"` and `"HC3"`). The benchmark estimators store only
the covariance they were fitted with; refit them with another `vcov`
instead.

## See also

Other fgls:
[`compare_fgls()`](https://emiliogcs.com/fgls/reference/compare_fgls.md),
[`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md),
[`fgls_control()`](https://emiliogcs.com/fgls/reference/fgls_control.md),
[`lasso_fgls()`](https://emiliogcs.com/fgls/reference/lasso_fgls.md),
[`plot.fgls()`](https://emiliogcs.com/fgls/reference/plot.fgls.md),
[`serial_fgls()`](https://emiliogcs.com/fgls/reference/serial_fgls.md),
[`summary.fgls()`](https://emiliogcs.com/fgls/reference/summary.fgls.md)

## Examples

``` r
sim <- simulate_serial(200, seed = 1)
fit <- serial_fgls(y ~ x, data = sim)
vcov(fit)
#>               (Intercept)             x
#> (Intercept)  1.282794e-01 -7.488239e-05
#> x           -7.488239e-05  1.369205e-03
vcov(fit, type = "naive")
#>               (Intercept)             x
#> (Intercept)  1.282722e-01 -8.898623e-05
#> x           -8.898623e-05  1.341769e-03
confint(fit)
#>                 2.5 %   97.5 %
#> (Intercept) 0.4984541 1.902420
#> x           0.9078455 1.052894
predict(fit, newdata = data.frame(x = c(-1, 0, 1)), interval = "confidence")
#>         fit        lwr       upr
#> 1 0.2200675 -0.4860593 0.9261944
#> 2 1.2004371  0.4984541 1.9024201
#> 3 2.1808067  1.4754950 2.8861183
head(residuals(fit, type = "innovations"))
#>          2          3          4          5          6          7 
#> -0.8577166  1.5184545  0.2807810 -0.8097611  0.4004779  0.7147435 
nobs(fit)
#> [1] 199
```
