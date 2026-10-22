# Validation: replicating the papers and agreement with pyfgls

The results below were **precomputed** (the Monte Carlo simulations take
tens of minutes) by `tools/replicate_tables.R` and
`tools/compare_python_r.R` in the source repository, and are stored in
the package:

``` r

rep <- readRDS(system.file("extdata", "replication.rds", package = "fgls"))
xl <- readRDS(system.file("extdata", "cross_language.rds", package = "fgls"))
rep$replications
#>  pg26_table1 pg26_general  gp24_table1           s9 
#>         2000         1000         1000         1000
```

The papers use 10,000 replications; with 1,000 the Monte Carlo standard
error of a coverage rate near 0.95 is about 0.007. The same designs are
run, with tolerances, by the package’s acceptance tests
(`FGLS_SLOW_TESTS=true`).

## Perron and González-Coya (2026), Table 1

AR(1) errors (0.8), a regressor that responds to past innovations, \\T =
500\\; slope estimates (true value 1):

``` r

knitr::kable(rep$pg26_table1, digits = 3)
```

| estimator        |  rmse |   bias | paper_rmse | paper_bias |
|:-----------------|------:|-------:|-----------:|-----------:|
| OLS              | 0.399 |  0.395 |      0.399 |      0.395 |
| Durbin           | 0.038 |  0.002 |      0.036 |         NA |
| GLS (infeasible) | 0.024 |  0.000 |      0.025 |         NA |
| FGLS             | 0.027 |  0.001 |      0.025 |      0.001 |
| Cochrane-Orcutt  | 0.041 |  0.031 |      0.041 |      0.031 |
| rho (FGLS)       | 0.034 | -0.008 |      0.034 |     -0.008 |
| rho (CO)         | 0.174 | -0.170 |         NA |     -0.171 |

``` r

rep$pg26_table1_median_k
#> [1] 1
```

OLS is inconsistent, Cochrane-Orcutt is biased (its estimate of \\\rho\\
is biased downwards), and FGLS matches the infeasible GLS.

## Perron and González-Coya (2026), Tables 2–5 (selected cells)

\\T = 200\\, \\x_t = 0.8 x\_{t-1} + v_t + \gamma e\_{t-1}\\; MSE
\\\times 100\\ and coverage of 95% intervals with the robust (eq. 29)
and naive (eq. 27) variances:

``` r

knitr::kable(rep$pg26_general, digits = 3)
```

| errors | gamma | mse_ols | mse_fgls | mse_gls | coverage_robust | coverage_naive | median_k |
|:---|---:|---:|---:|---:|---:|---:|---:|
| AR(1) 0.5 | 0.00 | 0.569 | 0.418 | 0.409 | 0.943 | 0.943 | 1 |
| AR(1) 0.5 | 0.25 | 0.981 | 0.406 | 0.392 | 0.947 | 0.941 | 1 |
| AR(1) 0.5 | 0.50 | 1.722 | 0.393 | 0.337 | 0.950 | 0.924 | 1 |
| AR(1) 0.8 | 0.00 | 2.221 | 0.512 | 0.504 | 0.944 | 0.944 | 1 |
| AR(1) 0.8 | 0.25 | 5.375 | 0.506 | 0.478 | 0.937 | 0.935 | 1 |
| AR(1) 0.8 | 0.50 | 11.162 | 0.479 | 0.406 | 0.929 | 0.919 | 1 |
| AR(2) (1.34, -0.42) | 0.00 | 11.401 | 0.401 | 0.400 | 0.948 | 0.947 | 2 |
| AR(2) (1.34, -0.42) | 0.25 | 27.826 | 0.394 | 0.380 | 0.946 | 0.943 | 2 |
| AR(2) (1.34, -0.42) | 0.50 | 58.056 | 0.374 | 0.324 | 0.947 | 0.934 | 2 |
| AR(2) (0.5, -0.3) | 0.00 | 0.323 | 0.283 | 0.280 | 0.948 | 0.948 | 2 |
| AR(2) (0.5, -0.3) | 0.25 | 0.351 | 0.275 | 0.270 | 0.953 | 0.950 | 2 |
| AR(2) (0.5, -0.3) | 0.50 | 0.398 | 0.261 | 0.237 | 0.955 | 0.934 | 2 |

The robust intervals have close to nominal coverage for all \\\gamma\\;
the naive ones under-cover when the regressor is not exogenous (\\\gamma
\> 0\\).

## González-Coya and Perron (2024), Table 1

Empirical MSE of the adaptive-Lasso FGLS slope relative to OLS, \\T =
400\\, with the paper’s values, and coverage of the eq. (6)+(8)
intervals:

``` r

knitr::kable(rep$gp24_table1, digits = 3)
```

| skedastic     |   T | emse_ratio | paper_ratio | coverage |
|:--------------|----:|-----------:|------------:|---------:|
| homoskedastic | 400 |      1.010 |       1.000 |    0.957 |
| x             | 400 |      0.905 |       0.920 |    0.953 |
| x2            | 400 |      0.687 |       0.700 |    0.955 |
| log2          | 400 |      0.449 |       0.460 |    0.945 |
| exp           | 400 |      0.443 |       0.450 |    0.955 |
| step          | 400 |      0.728 |       0.735 |    0.954 |

At \\T = 100\\–\\200\\ the ratios of this implementation (and of pyfgls)
are a few hundredths above the paper’s; the T = 400 column is matched.

## Supplement, Table S.9: F-C and F-CH

\\T = 200\\; MSE \\\times 100\\, coverage and average standard errors:

``` r

knitr::kable(rep$s9, digits = 3)
```

| design | mse_ols | mse_fc | mse_fch | coverage_fc | coverage_fch | se_fc | se_fch |
|:---|---:|---:|---:|---:|---:|---:|---:|
| AR(1) 0.5, x^2 | 6.802 | 3.973 | 2.825 | 0.966 | 0.962 | 0.212 | 0.174 |
| AR(1) 0.5, \[log x\]^2 | 0.809 | 0.469 | 0.212 | 0.963 | 0.968 | 0.073 | 0.048 |
| AR(1) 0.5, exp | 13.954 | 8.504 | 3.751 | 0.966 | 0.963 | 0.306 | 0.201 |
| AR(2) (1.34, -0.42), x^2 | 48.770 | 1.535 | 1.043 | 0.964 | 0.958 | 0.135 | 0.108 |
| MA(1) 0.5, x^2 | 6.339 | 3.942 | 2.849 | 0.966 | 0.958 | 0.210 | 0.175 |
| ARMA(1,1) (0.8, -0.4), x^2 | 7.042 | 4.287 | 3.120 | 0.956 | 0.956 | 0.216 | 0.180 |
| AR(1) 0.5, x^2, mixture | 5.925 | 3.813 | 2.515 | 0.965 | 0.951 | 0.211 | 0.162 |
| AR(1) 0.5, homoskedastic | 0.937 | 0.535 | 0.542 | 0.963 | 0.958 | 0.076 | 0.075 |

F-CH improves on F-C under heteroskedasticity and loses nothing under
homoskedasticity. With 1,000 replications some F-CH coverage rates are
slightly above 0.95; the intervals are conservative, not too short.

## Agreement with pyfgls

`fgls` is a port of the Python package pyfgls. Both packages fitted the
same 500 data sets per design with the same cross-validation folds
(pyfgls 0.1.0):

``` r

s <- xl$summary
s <- s[s$quantity %in% c("beta", "se", "kstar", "selected", "coverage_95"), ]
knitr::kable(s, digits = 15, row.names = FALSE)
```

| design | quantity | max_abs | median_abs | max_rel | share_different |
|:---|:---|---:|---:|---:|---:|
| serial_ar1_nonexog | beta | 2.000000e-15 | 0e+00 | 2.000000e-15 | 0.000 |
| serial_ar1_nonexog | se | 0.000000e+00 | 0e+00 | 1.000000e-15 | 0.000 |
| serial_ar1_nonexog | kstar | 0.000000e+00 | 0e+00 | 0.000000e+00 | 0.000 |
| serial_ar1_nonexog | coverage_95 | 0.000000e+00 | NA | NA | NA |
| serial_ar2_exog | beta | 2.000000e-15 | 0e+00 | 2.000000e-15 | 0.000 |
| serial_ar2_exog | se | 0.000000e+00 | 0e+00 | 2.000000e-15 | 0.000 |
| serial_ar2_exog | kstar | 0.000000e+00 | 0e+00 | 0.000000e+00 | 0.000 |
| serial_ar2_exog | coverage_95 | 0.000000e+00 | NA | NA | NA |
| serial_ma1 | beta | 2.000000e-15 | 0e+00 | 2.000000e-15 | 0.000 |
| serial_ma1 | se | 0.000000e+00 | 0e+00 | 2.000000e-15 | 0.000 |
| serial_ma1 | kstar | 0.000000e+00 | 0e+00 | 0.000000e+00 | 0.000 |
| serial_ma1 | coverage_95 | 0.000000e+00 | NA | NA | NA |
| hetero_log2 | beta | 4.110000e-13 | 1e-15 | 4.070000e-13 | 0.000 |
| hetero_log2 | se | 5.798303e-04 | 0e+00 | 1.056005e-02 | 0.012 |
| hetero_log2 | selected | NA | NA | NA | 0.012 |
| hetero_log2 | coverage_95 | 0.000000e+00 | NA | NA | NA |
| hetero_exp | beta | 6.590000e-13 | 1e-15 | 6.340000e-13 | 0.000 |
| hetero_exp | se | 3.609292e-03 | 1e-15 | 1.260608e-02 | 0.010 |
| hetero_exp | selected | NA | NA | NA | 0.010 |
| hetero_exp | coverage_95 | 0.000000e+00 | NA | NA | NA |
| combined_ar1_log2 | beta | 7.020000e-13 | 1e-15 | 7.120000e-13 | 0.000 |
| combined_ar1_log2 | se | 4.827348e-04 | 0e+00 | 1.164951e-02 | 0.004 |
| combined_ar1_log2 | kstar | 0.000000e+00 | 0e+00 | 0.000000e+00 | 0.000 |
| combined_ar1_log2 | selected | NA | NA | NA | 0.004 |
| combined_ar1_log2 | coverage_95 | 0.000000e+00 | NA | NA | NA |
| combined_ar1_exp | beta | 1.541700e-11 | 1e-15 | 1.338100e-11 | 0.000 |
| combined_ar1_exp | se | 2.956611e-03 | 0e+00 | 1.320481e-02 | 0.006 |
| combined_ar1_exp | kstar | 0.000000e+00 | 0e+00 | 0.000000e+00 | 0.000 |
| combined_ar1_exp | selected | NA | NA | NA | 0.006 |
| combined_ar1_exp | coverage_95 | 0.000000e+00 | NA | NA | NA |
| combined_ar2_x2 | beta | 1.760000e-12 | 1e-15 | 1.933000e-12 | 0.000 |
| combined_ar2_x2 | se | 1.355601e-03 | 0e+00 | 1.325978e-02 | 0.002 |
| combined_ar2_x2 | kstar | 0.000000e+00 | 0e+00 | 0.000000e+00 | 0.000 |
| combined_ar2_x2 | selected | NA | NA | NA | 0.002 |
| combined_ar2_x2 | coverage_95 | 0.000000e+00 | NA | NA | NA |

Serial-correlation quantities agree to rounding error, and the selected
AR order is identical in every data set. In about 1% of the Lasso data
sets the degrees of freedom differ by one, because pyfgls 0.1.0 counts
coefficients of order \\10^{-17}\\ (left by interpolating its Lasso path
at a knot) as selected. With the fixes proposed to pyfgls, everything
agrees and the selected sets are identical:

``` r

xl$pyfgls_pr2$version
#> [1] "0.1.0 + emiliogcs/pyfgls#2"
knitr::kable(xl$pyfgls_pr2$summary, digits = 15, row.names = FALSE)
```

| design             |      max_rel | share_support_different |
|:-------------------|-------------:|------------------------:|
| serial_ar1_nonexog | 2.300000e-14 |                      NA |
| serial_ar2_exog    | 1.673000e-12 |                      NA |
| serial_ma1         | 5.974939e-08 |                      NA |
| hetero_log2        | 1.621900e-11 |                       0 |
| hetero_exp         | 1.521700e-11 |                       0 |
| combined_ar1_log2  | 3.056000e-12 |                       0 |
| combined_ar1_exp   | 3.445600e-11 |                       0 |
| combined_ar2_x2    | 7.462000e-12 |                       0 |

The deterministic parts are also checked against golden fixtures
generated by pyfgls (`tests/testthat/test-golden-*.R`), to relative
precision \\10^{-8}\\ or better.
