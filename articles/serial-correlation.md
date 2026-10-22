# Serial correlation (Perron and González-Coya, 2026)

``` r

library(fgls)
```

## The problem

In \\y_t = x_t'\beta + u_t\\ with autocorrelated errors, GLS (or its
feasible version) is more efficient than OLS. The classic feasible
version, Cochrane-Orcutt, estimates the error autocorrelation from OLS
residuals. That is fine when the regressors are *exogenous*,
uncorrelated with the errors at all leads and lags. Many time-series
regressors are only *pre-determined*: they are uncorrelated with current
and future innovations, but respond to past ones, as with policy
variables, expectations and lagged prices. OLS is then inconsistent, and
so is Cochrane-Orcutt, which inherits its residuals.

## The estimator

The errors are approximated by an AR(k) process, \\u_t = \rho_1
u\_{t-1} + \dots + \rho_k u\_{t-k} + e_t\\. Substituting \\u_t = y_t -
x_t'\beta\\ gives the **Durbin regression**

\\y_t = \sum\_{j=1}^k \rho_j y\_{t-j} + x_t'\beta - \sum\_{j=1}^k
x\_{t-j}'\delta_j + e_t ,\\

estimated by OLS. Every regressor of this equation is uncorrelated with
the innovation \\e_t\\ when \\x_t\\ is pre-determined, so the
coefficients on the lagged \\y\\ estimate the \\\rho_j\\ consistently.
In `fgls`:

1.  the order \\k\\ minimizes BIC over \\0, \dots, k\_{max}\\, all
    orders being fitted on the same sample (eq. 14);
2.  the data are quasi-differenced with the estimated filter, every
    column including the constant, and the first \\k\\ observations are
    dropped; OLS on the filtered data gives \\\hat\beta\\ (eq. 15).

Lagged columns that duplicate other regressors (lags of the constant, of
a linear trend, of periodic dummies) are dropped automatically from the
Durbin regression.

## Table 1 of the paper

``` r

sim <- simulate_serial(500, seed = 7) # AR(1) 0.8, x_t = v_t + e_{t-1}
fits <- list(
  OLS = ols_robust(y ~ x, data = sim),
  `Cochrane-Orcutt` = cochrane_orcutt(y ~ x, data = sim),
  Durbin = durbin_estimator(y ~ x, data = sim),
  FGLS = serial_fgls(y ~ x, data = sim),
  `GLS (known rho)` = infeasible_gls(y ~ x, data = sim, rho = 0.8)
)
sapply(fits, function(f) coef(f)[["x"]])
#>             OLS Cochrane-Orcutt          Durbin            FGLS GLS (known rho) 
#>        1.465718        1.065165        1.026262        1.024933        1.025309
```

The true slope is 1. In the paper’s Monte Carlo (10,000 replications)
OLS has a bias of about 0.40, Cochrane-Orcutt about 0.03, and FGLS is as
precise as the infeasible GLS;
[`vignette("validation")`](https://emiliogcs.github.io/fgls/articles/validation.md)
reproduces these numbers with this package.

## Variance and inference

Because \\\hat\rho\\ is estimated, and the regressors may react to past
innovations, the usual GLS variance (`vcov = "naive"`, eq. 27)
understates the uncertainty of \\\hat\beta\\. The default
`vcov = "robust"` (eqs. 28–29) adds a correction built from the
covariance between the filtered regressors and past residuals. It is
valid whether the regressors are exogenous or merely pre-determined, and
it is positive semidefinite by construction.

``` r

fit <- fits$FGLS
cbind(robust = sqrt(diag(vcov(fit))), naive = sqrt(diag(vcov(fit, type = "naive"))))
#>                 robust      naive
#> (Intercept) 0.22499434 0.22497588
#> x           0.02591955 0.02558591
fit$nonexogeneity_share
#> [1] 0.01320714
```

`nonexogeneity_share` averages the relative increase of the robust over
the naive variance across coefficients. Values near zero mean the
regressors behave as exogenous.

## Lag selection

``` r

plot(fit, which = 3)
```

![](serial-correlation_files/figure-html/ic-1.png)

``` r

fit$bic_path
#>            0            1            2            3            4            5 
#>  0.814073672 -0.027959929 -0.015479053 -0.007432489  0.005286595  0.016958147 
#>            6            7            8            9           10           11 
#>  0.026916513  0.037226391  0.049823030  0.057844341  0.061138051  0.066065650 
#>           12           13           14           15 
#>  0.076799479  0.075244282  0.084468361  0.091024190
```

`kmax` defaults to \\\lfloor 12 (T/200)^{1/4} \rfloor\\ (12 at \\T =
200\\, as in the paper’s simulations) and is capped so that the Durbin
regression uses at most half of the sample. `ar_order` fixes the order;
`ic = "aic"` and `fgls_control(min_lags = "adjusted")` are alternatives.
With several stochastic regressors,
`fgls_control(ic_penalty = "params")` penalizes every coefficient a lag
adds.

## Lagged dependent variables

When lags of \\y\\ are regressors, the Durbin coefficients on lagged
\\y\\ mix the dynamics of the model with those of the errors. Under the
common-factor structure the AR coefficients can be recovered from the
ratio of a lagged regressor’s coefficient to its contemporaneous one
(Wallis, 1967; Section 4.2 of the paper). Declare the lagged dependent
variables:

``` r

set.seed(3)
n <- 400
w <- rnorm(n + 1)
u <- as.numeric(stats::filter(rnorm(n + 1), 0.5, method = "recursive"))
yy <- numeric(n + 1)
for (t in 2:(n + 1)) yy[t] <- 0.4 * yy[t - 1] + 1 + w[t] + u[t]
dyn <- data.frame(y = yy[-1], ylag = yy[-(n + 1)], w = w[-1])
fit_dyn <- serial_fgls(y ~ ylag + w, data = dyn, lagged_dep_vars = "ylag")
coef(fit_dyn)
#> (Intercept)        ylag           w 
#>   0.8010813   0.4491766   1.0389021
fit_dyn$rho
#>      rho1 
#> 0.4586741
```

Without the declaration,
[`fgls()`](https://emiliogcs.github.io/fgls/reference/fgls.md) warns
that a regressor looks like a lagged dependent variable.

## MA(1) errors

For errors known to be MA(1), such as two-step-ahead forecast errors,
the GMA estimator (Section 5) inverts the MA filter with the
maximum-likelihood estimate of \\\theta\\:

``` r

pred <- simulate_predictive(300, theta = -0.4, seed = 4)
serial_fgls(y ~ x, data = pred, serial = "ma1")
#> 
#> Call:
#> serial_fgls(formula = y ~ x, data = pred, serial = "ma1")
#> 
#> GMA (MA(1)); MA(1) theta = -0.4643; covariance: nonrobust
#> 
#> Coefficients:
#> (Intercept)            x  
#>      0.9581       0.9370
```

## Checking pre-determinedness

[`test_predetermined()`](https://emiliogcs.github.io/fgls/reference/test_predetermined.md)
(experimental) is a variable-addition LM test for missing lags of the
regressors (Remark 9 of the paper). A rejection suggests adding lags of
the regressors to the model:

``` r

suppressWarnings(test_predetermined(y ~ x, data = sim))
#> Pre-determinedness LM test (PG26 Remark 9, experimental)
#> Durbin AR order: 1; regressors tested: x; scale: innovations
#> 
#>  lags      LM df p.value
#>     1 0.01491  1  0.9028
#>     2 0.31255  2  0.8553
#>     3 0.85698  3  0.8358
#>     4 1.58321  4  0.8118
#> 
#> No rejection at 5%: no evidence of missing regressor lags.
```
