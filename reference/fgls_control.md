# Advanced options for fgls

Tuning options of
[`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md) that most users
never change. The defaults follow the papers and the Python package
'pyfgls'.

## Usage

``` r
fgls_control(
  ic_penalty = c("order", "params"),
  min_lags = NULL,
  fold_scheme = c("random", "blocked"),
  foldid = NULL,
  cv_weights = c("fold", "full"),
  ridge_scale = c("standardized", "original"),
  n_lambda = 100L,
  lambda_min_ratio = NULL,
  lasso_engine = c("lars", "glmnet"),
  glmnet_thresh = 1e-07,
  delta_scale = FALSE,
  transforms = c("linear", "poly2", "log", "trig"),
  svr_cost = c(0.1, 1, 10, 100),
  svr_epsilon = c(0.05, 0.25, 0.5, 1),
  n_boot = 999L,
  bootstrap = c("wild", "iid"),
  wild_dist = c("rademacher", "mammen"),
  boot_seed = NULL,
  cores = 1L
)
```

## Arguments

- ic_penalty:

  Penalty of the lag-order criterion: `"order"` (default, PG26 eq. 14)
  penalizes the AR order; `"params"` penalizes every coefficient a lag
  adds to the Durbin regression, which avoids over-selection with many
  stochastic regressors.

- min_lags:

  Lower bound on the selected AR order: `NULL` (default), a non-negative
  integer, or `"adjusted"` for \\\lceil(\ln T)^{3/2}\rceil\\ (PG26
  Supplement, Lemma A.2).

- fold_scheme:

  Cross-validation folds of the skedastic fit: `"random"` (default) or
  `"blocked"` (contiguous blocks).

- foldid:

  Optional integer vector of fold labels (one per observation entering
  the skedastic fit), overriding `fold_scheme` and the seed.

- cv_weights:

  `"fold"` (default) re-estimates the ridge first step, and hence the
  adaptive weights, on every training fold; `"full"` computes the
  weights once on the full sample (the common practice with
  [`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
  and `penalty.factor`).

- ridge_scale:

  Scale on which the ridge coefficients define the adaptive weights:
  `"standardized"` (default) or `"original"` units of `z`.

- n_lambda:

  Length of the Lasso penalty grid.

- lambda_min_ratio:

  Smallest penalty relative to the largest; default `1e-4` when there
  are more observations than candidate covariates and `1e-2` otherwise
  (as in 'glmnet').

- lasso_engine:

  Solver of the Lasso problems: `"lars"` (default) computes the exact
  solution path by the LARS homotopy with the Lasso modification, as
  'pyfgls' does; `"glmnet"` uses
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
  (coordinate descent) followed by an exact active-set refinement. Both
  give the same solutions; `"lars"` is faster for the moderate numbers
  of candidate covariates typical here.

- glmnet_thresh:

  Convergence threshold passed to
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
  when `lasso_engine = "glmnet"`.

- delta_scale:

  Use the scale-aware floor `delta * sigma_OLS` in \\\log\max(\hat
  u_t^2, \delta^2)\\.

- transforms:

  Transform families used by `z = "auto"`; see
  [`make_dictionary()`](https://emiliogcs.com/fgls/reference/make_dictionary.md).

- svr_cost, svr_epsilon:

  Grids for the SVR skedastic fit, multiplied by the standard deviation
  of the response.

- n_boot:

  Number of bootstrap replications for `vcov = "bootstrap"`.

- bootstrap:

  Innovation scheme of the sieve bootstrap: `"wild"` (default) or
  `"iid"` resampling.

- wild_dist:

  Wild-bootstrap multipliers: `"rademacher"` (default) or `"mammen"`.

- boot_seed:

  Seed of the bootstrap; default the `seed` of the fit.

- cores:

  Number of cores for the bootstrap (forked processes via 'parallel';
  ignored on Windows). Defaults to 1.

## Value

A list of class `"fgls_control"`.

## See also

Other fgls:
[`compare_fgls()`](https://emiliogcs.com/fgls/reference/compare_fgls.md),
[`fgls()`](https://emiliogcs.com/fgls/reference/fgls.md),
[`fgls-methods`](https://emiliogcs.com/fgls/reference/fgls-methods.md),
[`lasso_fgls()`](https://emiliogcs.com/fgls/reference/lasso_fgls.md),
[`plot.fgls()`](https://emiliogcs.com/fgls/reference/plot.fgls.md),
[`serial_fgls()`](https://emiliogcs.com/fgls/reference/serial_fgls.md),
[`summary.fgls()`](https://emiliogcs.com/fgls/reference/summary.fgls.md)

## Examples

``` r
fgls_control(ic_penalty = "params", n_boot = 199)
#> $ic_penalty
#> [1] "params"
#> 
#> $min_lags
#> NULL
#> 
#> $fold_scheme
#> [1] "random"
#> 
#> $foldid
#> NULL
#> 
#> $cv_weights
#> [1] "fold"
#> 
#> $ridge_scale
#> [1] "standardized"
#> 
#> $n_lambda
#> [1] 100
#> 
#> $lambda_min_ratio
#> NULL
#> 
#> $lasso_engine
#> [1] "lars"
#> 
#> $glmnet_thresh
#> [1] 1e-07
#> 
#> $delta_scale
#> [1] FALSE
#> 
#> $transforms
#> [1] "linear" "poly2"  "log"    "trig"  
#> 
#> $svr_cost
#> [1]   0.1   1.0  10.0 100.0
#> 
#> $svr_epsilon
#> [1] 0.05 0.25 0.50 1.00
#> 
#> $n_boot
#> [1] 199
#> 
#> $bootstrap
#> [1] "wild"
#> 
#> $wild_dist
#> [1] "rademacher"
#> 
#> $boot_seed
#> NULL
#> 
#> $cores
#> [1] 1
#> 
#> attr(,"class")
#> [1] "fgls_control"
```
