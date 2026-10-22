
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fgls

<!-- badges: start -->

[![R-CMD-check](https://github.com/emiliogcs/fgls/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/emiliogcs/fgls/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**fgls** estimates linear regressions by feasible generalized least
squares when the errors are serially correlated, heteroskedastic of
unknown form, or both:

- **Serial correlation** (Perron and González-Coya, 2026). An AR filter
  is estimated from a Durbin regression, with its order chosen by BIC.
  This is consistent when the regressors are only *pre-determined* (they
  react to past shocks), where OLS and Cochrane-Orcutt are not. A
  variance estimator accounts for the estimated filter.
- **Heteroskedasticity of unknown form** (González-Coya and Perron,
  2024). An adaptive Lasso of the log squared residuals on candidate
  covariates gives the weights, and a bias-corrected sandwich variance
  keeps inference valid when the variance model is wrong. If nothing is
  selected, FGLS equals OLS.
- **Both together** (Supplement S-4 of the 2026 paper).

It is the R counterpart of the Python package
[pyfgls](https://pypi.org/project/pyfgls/), and agrees with it to
machine precision on deterministic quantities. More code by the author:
<https://emiliogcs.com/code/>.

## Installation

``` r
# development version
# install.packages("pak")
pak::pak("emiliogcs/fgls")
# from CRAN, once released
install.packages("fgls")
```

## Serial correlation with a non-exogenous regressor

``` r
library(fgls)
sim <- simulate_serial(500, seed = 1) # design of Table 1 of the 2026 paper; true slope 1
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
compare_fgls(OLS = lm(y ~ x, data = sim), FGLS = fit)
#>          term statistic        OLS       FGLS
#> 1 (Intercept)  estimate 1.12860678 1.11876951
#> 2 (Intercept) std.error 0.09034505 0.20551635
#> 3 (Intercept)  ci.width 0.35414609 0.80560929
#> 4           x  estimate 1.32962748 0.98062019
#> 5           x std.error 0.05108277 0.02516513
#> 6           x  ci.width 0.20024076 0.09864549
```

## Heteroskedasticity of unknown form

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

## Both

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
```

See `vignette("fgls-intro")` to get started, and the articles on [serial
correlation](https://emiliogcs.github.io/fgls/articles/serial-correlation.html),
[heteroskedasticity](https://emiliogcs.github.io/fgls/articles/heteroskedasticity.html),
the [combined
estimator](https://emiliogcs.github.io/fgls/articles/combined.html) and
the
[validation](https://emiliogcs.github.io/fgls/articles/validation.html)
against the papers and pyfgls.

## Citation

Perron P, González-Coya E (2026). “Feasible GLS for Time Series
Regression.” *Journal of Econometrics*, *258*, 106346.
<doi:10.1016/j.jeconom.2026.106346>
<https://doi.org/10.1016/j.jeconom.2026.106346>.

González-Coya E, Perron P (2024). “Estimation in the Presence of
Heteroskedasticity of Unknown Form: A Lasso-based Approach.” *Journal of
Econometric Methods*, *13*(1), 29-48. <doi:10.1515/jem-2023-0007>
<https://doi.org/10.1515/jem-2023-0007>.

González-Coya E (2026). *fgls: Feasible GLS for Serial Correlation and
Heteroskedasticity*. R package version 0.1.0,
<https://github.com/emiliogcs/fgls>.
