# Tidy, glance and augment fgls fits

Methods for the 'broom' generics (from the 'generics' package),
registered when 'generics' is loaded.

## Usage

``` r
# S3 method for class 'fgls'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'fgls'
glance(x, ...)

# S3 method for class 'fgls'
augment(
  x,
  data = NULL,
  newdata = NULL,
  se_fit = FALSE,
  interval = c("none", "confidence"),
  conf.level = 0.95,
  ...
)
```

## Arguments

- x:

  An `"fgls"` object.

- conf.int:

  Include confidence intervals.

- conf.level:

  Confidence level.

- ...:

  Not used.

- data:

  Data to augment: by default the model frame. It must have one row per
  observation used in the fit or, for a fit with
  `na.action = na.exclude`, one row per row of the original data (the
  excluded rows then get `NA`).

- newdata:

  Optional new data; fitted values (and, if requested, standard errors
  and confidence limits) are added, but no residuals.

- se_fit:

  Add the standard errors of the fitted means (`.se.fit`).

- interval:

  `"none"` or `"confidence"`: add confidence limits for the mean of the
  response (`.lower`, `.upper`) at level `conf.level`.

## Value

A data frame (a tibble if 'tibble' is installed): one row per
coefficient for `tidy()`, one row for `glance()`, one row per
observation for `augment()` (columns `.fitted`, optionally `.lower`,
`.upper` and `.se.fit`, then `.resid`, `.innovation` and, with a
heteroskedasticity correction, `.weight`). `nonexogeneity.share` in
`glance()` is `NA` when the share was not computed.

## Examples

``` r
if (requireNamespace("broom", quietly = TRUE)) {
  fit <- serial_fgls(y ~ x, data = simulate_serial(200, seed = 1))
  broom::tidy(fit, conf.int = TRUE)
  broom::glance(fit)
  head(broom::augment(fit))
}
#> # A tibble: 6 × 5
#>       y      x .fitted   .resid .innovation
#>   <dbl>  <dbl>   <dbl>    <dbl>       <dbl>
#> 1 2.25   1.06   2.24    0.00406      NA    
#> 2 2.08   1.77   2.94   -0.854        -0.858
#> 3 0.878 -1.17   0.0568  0.821         1.52 
#> 4 1.48  -0.690  0.524   0.951         0.281
#> 5 3.94   2.83   3.97   -0.0334       -0.810
#> 6 1.42  -0.153  1.05    0.373         0.400
```
