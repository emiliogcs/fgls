# Diagnostics of an fgls fit

- `weight_diagnostics()`: how concentrated the FGLS weights \\1/\hat
  v_t\\ are: the ratio of the largest to the smallest fitted variance,
  the Kish effective sample size \\(\sum w)^2 / \sum w^2\\ as a share of
  the sample, and the share of the total weight carried by the top 1% of
  observations.

- `ar_roots()`: roots of the estimated lag polynomial \\1 - \hat\rho_1
  z - \dots - \hat\rho_k z^k\\; the filter is invertible when all lie
  outside the unit circle.

- `nonexogeneity_share()`: the relative size of the first-stage
  correction of PG26 eq. (29), averaged over coefficients. Values near
  zero mean the regressors behave as exogenous; large values mean the
  naive variance (eq. 27) would understate the uncertainty.

## Usage

``` r
weight_diagnostics(fit)

ar_roots(fit)

nonexogeneity_share(fit)
```

## Arguments

- fit:

  An `"fgls"` object (with a heteroskedasticity correction for
  `weight_diagnostics()`, a serial-correlation correction otherwise).

## Value

A named numeric vector, a complex vector, or a number.

## See also

Other diagnostics:
[`adf_test()`](https://emiliogcs.github.io/fgls/reference/adf_test.md),
[`test_predetermined()`](https://emiliogcs.github.io/fgls/reference/test_predetermined.md)

## Examples

``` r
fit <- fgls(y ~ x, data = simulate_combined(200, seed = 1), seed = 1)
weight_diagnostics(fit)
#>       variance_ratio    effective_n_share top1pct_weight_share 
#>           5.57109452           0.81246809           0.02399232 
Mod(ar_roots(fit))
#> [1] 1.897991
nonexogeneity_share(fit)
#> [1] 0.000919948
```
