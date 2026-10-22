# Conditions signaled by fgls

All warnings raised by the package inherit from class `"fgls_warning"`
and all errors from `"fgls_error"`, so they can be caught or muffled as
a group with [`tryCatch()`](https://rdrr.io/r/base/conditions.html),
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) or
[`suppressWarnings()`](https://rdrr.io/r/base/warning.html) (R \>= 4.3
accepts a `classes` argument).

## Warning classes

- `fgls_warning_nonstationary`: the estimated autoregressive filter has
  a root on, inside or close to the unit circle.

- `fgls_warning_collinearity`: columns were dropped or a matrix was not
  positive definite.

- `fgls_warning_convergence`: an iterative routine did not fully
  converge.

- `fgls_warning_weight`: the estimated weights are extreme.

- `fgls_warning_sample_size`: the sample is small relative to the number
  of parameters, lags or folds.

- `fgls_warning_specification`: the specification looks unusual (e.g. no
  intercept, the selected order hits `kmax`).

- `fgls_warning_lagged_dependent`: a regressor looks like a lag of the
  dependent variable (also a specification warning).

- `fgls_warning_experimental`: the feature is experimental.

## Error classes

`fgls_error_input` (invalid arguments) and its subclasses
`fgls_error_missing`, `fgls_error_time_index`, `fgls_error_collinearity`
and `fgls_error_insufficient_data`.

## Examples

``` r
# an interior missing value in a time series is an error of class fgls_error_missing
df <- data.frame(y = c(1, 2, NA, 4, 5, 3, 2), x = c(2, 1, 3, 5, 4, 2, 1))
tryCatch(serial_fgls(y ~ x, data = df), fgls_error = function(e) class(e)[1])
#> [1] "fgls_error_missing"

# warnings can be caught by class too
withCallingHandlers(
  serial_fgls(y ~ x - 1, data = simulate_serial(100, seed = 1)),
  fgls_warning_specification = function(w) invokeRestart("muffleWarning")
)$ar_order
#> [1] 1
```
