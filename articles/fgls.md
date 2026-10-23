# Getting started with fgls

``` r

library(fgls)
```

**fgls** estimates linear regressions by *feasible generalized least
squares* (FGLS). It implements two recent methods:

- **Serial correlation** (Perron and González-Coya, 2026): the errors
  are approximated by an autoregression whose order is chosen by BIC and
  whose coefficients come from a *Durbin regression*. The resulting
  estimator is consistent even when the regressors are only
  *pre-determined*, a situation where OLS and Cochrane-Orcutt are not. A
  variance estimator accounts for the estimated filter.
- **Heteroskedasticity of unknown form** (González-Coya and Perron,
  2024): the log of the squared residuals is modeled by an *adaptive
  Lasso* over a dictionary of candidate covariates, and the weighted
  regression is paired with a bias-corrected sandwich variance that
  stays valid when that model is wrong.

The two corrections can be combined (Supplement S-4 of the 2026 paper).
The R package is a port of the Python package
[pyfgls](https://pypi.org/project/pyfgls/), and the two agree to machine
precision on deterministic quantities (see
[`vignette("validation")`](https://emiliogcs.com/fgls/articles/validation.md)).

## Three examples

### Serial correlation

[`simulate_serial()`](https://emiliogcs.com/fgls/reference/simulate_serial.md)
generates the design of Table 1 of the 2026 paper: AR(1) errors with
coefficient 0.8, and a regressor that responds to past innovations,
\\x_t = v_t + e\_{t-1}\\, so it is pre-determined but not exogenous.

``` r

sim <- simulate_serial(500, seed = 1)
fit <- serial_fgls(y ~ x, data = sim)
summary(fit)
#> 
#> Call:
#> serial_fgls(formula = y ~ x, data = sim)
#> 
#> Method:               FGLS (serial)
#> Covariance:           robust
#> Observations:         500 (effective 499)
#> Inference:            normal
#> AR order k*:          1 [BIC, kmax = 15]
#> AR coefficients:      0.7800
#> Non-exogeneity share: 0.011
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.11877    0.20552   5.444 5.22e-08 ***
#> x            0.98062    0.02517  38.967  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Notes:
#> [1] Covariance: PG26 eq. (29), robust to non-exogenous (pre-determined)
#>     regressors.
#> [2] Collinear lagged columns dropped from the Durbin regression:
#>     L1.(Intercept).
```

The true slope is 1. OLS is badly biased in this design; FGLS is not:

``` r

compare_fgls(OLS = lm(y ~ x, data = sim), FGLS = fit)
#>          term statistic        OLS       FGLS
#> 1 (Intercept)  estimate 1.12860678 1.11876951
#> 2 (Intercept) std.error 0.09034505 0.20551635
#> 3 (Intercept)  ci.width 0.35414609 0.80560929
#> 4           x  estimate 1.32962748 0.98062019
#> 5           x std.error 0.05108277 0.02516513
#> 6           x  ci.width 0.20024076 0.09864549
```

`fit$ar_order` is the order chosen by BIC, `fit$rho` the AR coefficients
and `fit$nonexogeneity_share` the relative size of the correction that
makes the variance robust to non-exogenous regressors.

### Heteroskedasticity of unknown form

``` r

het <- simulate_hetero(200, skedastic = "log2", seed = 2)
fit_h <- lasso_fgls(y ~ x, data = het, seed = 1)
fit_h$skedastic$selected
#> [1] "log(x)" "cos(x)"
compare_fgls(OLS = lm(y ~ x, data = het), FGLS = fit_h)
#>          term statistic        OLS       FGLS
#> 1 (Intercept)  estimate 0.98861426 0.98253777
#> 2 (Intercept) std.error 0.18311445 0.08973843
#> 3 (Intercept)  ci.width 0.71779546 0.35176818
#> 4           x  estimate 1.02852565 1.03190100
#> 5           x std.error 0.09523622 0.06032162
#> 6           x  ci.width 0.37331913 0.23645640
```

The candidate covariates default to \\x, x^2, \log x, (\log x)^2, \cos
x, \cos 2x\\ for every non-constant regressor
([`make_dictionary()`](https://emiliogcs.com/fgls/reference/make_dictionary.md));
supply your own with `z`. The penalty and the exponent of the adaptive
weights are chosen by 10-fold cross-validation, so set `seed` for
reproducible results.

### Both

``` r

both <- simulate_combined(200, ar = 0.5, skedastic = "x2", seed = 3)
fit_b <- fgls(y ~ x, data = both, seed = 1)
fit_b
#> 
#> Call:
#> fgls(formula = y ~ x, data = both, seed = 1)
#> 
#> F-CH (serial+hetero); AR order k* = 1 (BIC), rho = 0.4976; skedastic fit: lasso (selected: x, log(x)); covariance: robust
#> 
#> Coefficients:
#> (Intercept)            x  
#>      1.5934       0.8101
compare_fgls(`F-C` = fit_b$stages$serial, `F-CH` = fit_b)
#>          term statistic       F-C      F-CH
#> 1 (Intercept)  estimate 1.3723275 1.5934094
#> 2 (Intercept) std.error 0.5413481 0.4430785
#> 3 (Intercept)  ci.width 2.1220455 1.7368357
#> 4           x  estimate 0.9130718 0.8100948
#> 5           x std.error 0.2611294 0.2123379
#> 6           x  ci.width 1.0236086 0.8323494
```

[`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md) is the general
function: `serial = "none"` gives the heteroskedasticity-only estimator
and `hetero = "none"` the serial-only one.

## Working with fits

Fits behave like `lm` objects:

``` r

coef(fit)
#> (Intercept)           x 
#>   1.1187695   0.9806202
confint(fit, level = 0.9)
#>                   5 %     95 %
#> (Intercept) 0.7807252 1.456814
#> x           0.9392272 1.022013
predict(fit, newdata = data.frame(x = c(-1, 0, 1)), interval = "confidence")
#>         fit        lwr      upr
#> 1 0.1381493 -0.2675824 0.543881
#> 2 1.1187695  0.7159649 1.521574
#> 3 2.0993897  1.6934951 2.505284
vcov(fit, type = "naive") # PG26 eq. (27), valid only for exogenous regressors
#>              (Intercept)            x
#> (Intercept) 4.223673e-02 6.788463e-06
#> x           6.788463e-06 6.197699e-04
```

They work with
[`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html),
[`lmtest::coefci()`](https://rdrr.io/pkg/lmtest/man/coeftest.html) and
the broom generics `tidy()`, `glance()` and `augment()`. `plot(fit)`
draws diagnostics: the autocorrelation function of the innovations
(residuals of the final filtered regression) should look like white
noise.

## Time series inputs

Rows must be in time order. Besides data frames,
[`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md) accepts `ts`,
`zoo` and `xts` objects; a data frame can carry a time index through
`time =`, which is checked for sorting and gaps:

``` r

dates <- seq(as.Date("2000-01-01"), by = "month", length.out = 500)
df <- data.frame(date = dates, sim)
serial_fgls(y ~ x, data = df, time = "date")$ar_order
#> [1] 1
```

## Choosing options

- With several stochastic regressors, the BIC that penalizes only the AR
  order can select the largest order; use
  `control = fgls_control(ic_penalty = "params")`.
- Declare lagged dependent variables among the regressors in
  `lagged_dep_vars`.
- For errors known to be MA(1), such as multi-step forecast errors, use
  `serial = "ma1"`.
- `exog = "exogenous"` never changes the estimates. In the serial-only
  estimator it only changes the notes, and the robust variance is valid
  either way. In the combined estimator it also drops the first-stage
  correction from the variance by default
  (`first_stage_correction = FALSE`); that correction is needed when the
  regressors are only pre-determined, so keep the default
  `exog = "predetermined"` unless they are known to be exogenous.

See
[`vignette("serial-correlation")`](https://emiliogcs.com/fgls/articles/serial-correlation.md),
[`vignette("heteroskedasticity")`](https://emiliogcs.com/fgls/articles/heteroskedasticity.md)
and
[`vignette("combined")`](https://emiliogcs.com/fgls/articles/combined.md)
for the methods in more detail.

## References

Perron, P. and González-Coya, E. (2026). Feasible GLS for time series
regression. *Journal of Econometrics*, 258, 106346.
<https://doi.org/10.1016/j.jeconom.2026.106346>

González-Coya, E. and Perron, P. (2024). Estimation in the presence of
heteroskedasticity of unknown form: a Lasso-based approach. *Journal of
Econometric Methods*, 13(1), 29–48.
<https://doi.org/10.1515/jem-2023-0007>
