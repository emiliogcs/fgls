# Simulate the designs of the papers

Data-generating processes used in the Monte Carlo studies of Perron and
González-Coya (2026) (`simulate_serial()`, `simulate_combined()`,
`simulate_predictive()`) and González-Coya and Perron (2024)
(`simulate_hetero()`). They are handy for examples, teaching and
replication.

## Usage

``` r
simulate_serial(
  n,
  ar = 0.8,
  ma = numeric(0),
  rho_x = 0,
  gamma = 1,
  beta = c(1, 1),
  burn = 0,
  seed = NULL
)

simulate_hetero(
  n,
  skedastic = "x2",
  beta = c(1, 1),
  n_extra = 0L,
  extra_corr = 0,
  seed = NULL
)

simulate_combined(
  n,
  ar = 0.5,
  ma = numeric(0),
  skedastic = "x2",
  beta = c(1, 1),
  driver = c("x", "mixture"),
  phi_prob = 0.75,
  gamma = 0,
  burn = 100L,
  seed = NULL
)

simulate_predictive(
  n,
  theta = 0.5,
  rho_x = 0,
  gamma = 0,
  beta = c(1, 1),
  burn = 100L,
  seed = NULL
)
```

## Arguments

- n:

  Number of observations returned.

- ar, ma:

  Autoregressive and moving-average coefficients of the errors.

- rho_x:

  Autoregressive coefficient of the regressor.

- gamma:

  Feedback from past innovations to the regressor; the regressor is not
  exogenous when `gamma != 0`.

- beta:

  Intercept and slope.

- burn:

  Number of initial observations discarded.

- seed:

  Optional seed. The random number generator state of the session is
  restored afterwards.

- skedastic:

  Skedastic function: one of `"homoskedastic"`, `"x"`, `"x2"`, `"log2"`
  (\\\[\log x\]^2\\), `"exp"` (\\\exp(0.2x + 0.2x^2)\\), `"step"`, or a
  function of `x`.

- n_extra:

  Number of extra candidate covariates `w1, w2, ...` (columns of the
  result, not regressors).

- extra_corr:

  Weight \\c\\ of `x` in the extra covariates.

- driver:

  Variance driver of `simulate_combined()`: `"x"` or `"mixture"`.

- phi_prob:

  Probability that the variance driver equals `x` when
  `driver = "mixture"`.

- theta:

  MA(1) coefficient of `simulate_predictive()`.

## Value

A data frame with columns `y` and `x` (and `w` or `w1, w2, ...` where
relevant), with attributes `beta` (true coefficients), `u` (errors), `e`
(innovations) and `info` (design parameters and, for the heteroskedastic
designs, the true variances `v`).

## Details

- `simulate_serial()`: \\y_t = \beta_0 + \beta_1 x_t + u_t\\ with ARMA
  errors driven by standard normal innovations \\e_t\\ and a regressor
  \\x_t = \rho_x x\_{t-1} + v_t + \gamma e\_{t-1}\\, which is
  pre-determined and exogenous only when \\\gamma = 0\\. The defaults
  (\\\rho_x = 0\\, \\\gamma = 1\\, AR(1) coefficient 0.8, no burn-in,
  \\u_0 = 0\\) give the design of Table 1 of the 2026 paper.

- `simulate_hetero()`: \\y_t = \alpha + \beta x_t + \sqrt{v(x_t)}
  \varepsilon_t\\, \\x_t \sim U\[1, 4\]\\ (eq. 9 of the 2024 paper),
  optionally with extra candidate covariates \\w_j = c x + (1 - c)
  q_j\\, \\q_j \sim U\[1, 4\]\\.

- `simulate_combined()`: ARMA errors whose innovations are
  heteroskedastic, \\e_t = \sqrt{v(w_t)} \varepsilon_t\\, with \\w_t =
  x_t\\ (`driver = "x"`) or \\w_t = \phi_t x_t + (1 - \phi_t) q_t\\,
  \\\phi_t \sim\\ Bernoulli(`phi_prob`) (`driver = "mixture"`)
  (Supplement S-4).

- `simulate_predictive()`: the predictive regression \\y\_{t+2} =
  \alpha + \beta x_t + u\_{t+2}\\ with MA(1) errors (Supplement S-3),
  returned aligned so that row \\t\\ holds \\y\_{t+2}\\ and \\x_t\\.

## References

Perron, P. and González-Coya, E. (2026). Feasible GLS for time series
regression. *Journal of Econometrics*, 258, 106346.
[doi:10.1016/j.jeconom.2026.106346](https://doi.org/10.1016/j.jeconom.2026.106346)

González-Coya, E. and Perron, P. (2024). Estimation in the presence of
heteroskedasticity of unknown form: a Lasso-based approach. *Journal of
Econometric Methods*, 13(1), 29–48.
[doi:10.1515/jem-2023-0007](https://doi.org/10.1515/jem-2023-0007)

## Examples

``` r
sim <- simulate_serial(200, seed = 1)
head(sim)
#>           y          x
#> 1 2.2460628  1.0624195
#> 2 2.0815178  1.7702318
#> 3 0.8777732 -1.1665364
#> 4 1.4750008 -0.6899547
#> 5 3.9386654  2.8271694
#> 6 1.4232236 -0.1534022
attr(sim, "beta")
#> [1] 1 1

het <- simulate_hetero(100, skedastic = "log2", seed = 2)
summary(attr(het, "info")$v)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.0004453 0.2643021 0.7135063 0.8469443 1.4413330 1.8987080 
```
