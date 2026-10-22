# Candidate covariates for the skedastic function

Builds the default dictionary of candidate covariates used by
[`fgls()`](https://emiliogcs.github.io/fgls/reference/fgls.md) and
[`lasso_fgls()`](https://emiliogcs.github.io/fgls/reference/lasso_fgls.md)
with `z = "auto"`: for every non-constant regressor \\x\\, the columns
\\x\\, \\x^2\\, \\\log x\\ and \\(\log x)^2\\ (only when \\x \> 0\\),
\\\cos x\\ and \\\cos 2x\\. This contains the specification \\z = (1, x,
\[\log x\]^2, x^2, \cos x, \cos 2x)\\ of the simulations of
González-Coya and Perron (2024); the intercept of the skedastic
regression is always included and never penalized. Binary (dummy)
columns are passed through untransformed, undefined or constant
transforms are skipped and duplicated columns are removed.

## Usage

``` r
make_dictionary(
  x,
  transforms = c("linear", "poly2", "log", "trig"),
  skip = NULL
)
```

## Arguments

- x:

  A numeric matrix or data frame of regressors (constant columns are
  ignored), or a numeric vector.

- transforms:

  Transform families: `"linear"` (\\x\\), `"poly2"` (\\x^2\\), `"log"`
  (\\\log x\\, \\(\log x)^2\\ for positive \\x\\) and `"trig"` (\\\cos
  x\\, \\\cos 2x\\).

- skip:

  Columns (names or positions) to leave out, e.g. lagged dependent
  variables.

## Value

A numeric matrix with one column per candidate covariate and informative
column names.

## References

Perron, P. and González-Coya, E. (2026). Feasible GLS for time series
regression. *Journal of Econometrics*, 258, 106346.
[doi:10.1016/j.jeconom.2026.106346](https://doi.org/10.1016/j.jeconom.2026.106346)

González-Coya, E. and Perron, P. (2024). Estimation in the presence of
heteroskedasticity of unknown form: a Lasso-based approach. *Journal of
Econometric Methods*, 13(1), 29–48.
[doi:10.1515/jem-2023-0007](https://doi.org/10.1515/jem-2023-0007)

Romano, J. P. and Wolf, M. (2017). Resurrecting weighted least squares.
*Journal of Econometrics*, 197(1), 1–19.
[doi:10.1016/j.jeconom.2016.10.003](https://doi.org/10.1016/j.jeconom.2016.10.003)

## See also

Other heteroskedasticity:
[`lasso_fgls()`](https://emiliogcs.github.io/fgls/reference/lasso_fgls.md)

## Examples

``` r
X <- cbind(const = 1, x = seq(1, 4, length.out = 10), d = rep(0:1, 5))
colnames(make_dictionary(X))
#> [1] "x"        "x^2"      "log(x)"   "log(x)^2" "cos(x)"   "cos(2x)"  "d"       
colnames(make_dictionary(X, transforms = c("linear", "poly2")))
#> [1] "x"   "x^2" "d"  
```
