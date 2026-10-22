# Conditions (warnings and errors) with classes users can catch.
#
# Every warning inherits from "fgls_warning" and every error from "fgls_error",
# mirroring the PyFGLSWarning / PyFGLSError hierarchy of pyfgls:
#
#   warnings: fgls_warning_nonstationary, fgls_warning_collinearity,
#             fgls_warning_convergence, fgls_warning_weight,
#             fgls_warning_sample_size, fgls_warning_specification,
#             fgls_warning_lagged_dependent (also a specification warning),
#             fgls_warning_experimental
#   errors:   fgls_error_input, fgls_error_missing, fgls_error_time_index,
#             fgls_error_collinearity, fgls_error_insufficient_data
#             (all of them also inherit from fgls_error_input)

fgls_condition <- function(message, class, call = NULL) {
  structure(class = class, list(message = message, call = call))
}

# Signal a classed warning. `class` is the specific subclass, without prefix.
fgls_warn <- function(message, class = NULL) {
  cls <- c(if (!is.null(class)) paste0("fgls_warning_", class), "fgls_warning")
  if (identical(class, "lagged_dependent")) {
    cls <- append(cls, "fgls_warning_specification", after = 1L)
  }
  warning(fgls_condition(message, c(cls, "warning", "condition")))
}

# Signal a classed error. All fgls errors are input errors in the pyfgls sense.
fgls_abort <- function(message, class = "input") {
  cls <- unique(c(paste0("fgls_error_", class), "fgls_error_input", "fgls_error"))
  stop(fgls_condition(message, c(cls, "error", "condition")))
}

#' Conditions signaled by fgls
#'
#' All warnings raised by the package inherit from class `"fgls_warning"` and
#' all errors from `"fgls_error"`, so they can be caught or muffled as a group
#' with [tryCatch()], [withCallingHandlers()] or [suppressWarnings()] (R >= 4.3
#' accepts a `classes` argument).
#'
#' @section Warning classes:
#' * `fgls_warning_nonstationary`: the estimated autoregressive filter has a
#'   root on, inside or close to the unit circle.
#' * `fgls_warning_collinearity`: columns were dropped or a matrix was not
#'   positive definite.
#' * `fgls_warning_convergence`: an iterative routine did not fully converge.
#' * `fgls_warning_weight`: the estimated weights are extreme.
#' * `fgls_warning_sample_size`: the sample is small relative to the number of
#'   parameters, lags or folds.
#' * `fgls_warning_specification`: the specification looks unusual (e.g. no
#'   intercept, the selected order hits `kmax`).
#' * `fgls_warning_lagged_dependent`: a regressor looks like a lag of the
#'   dependent variable (also a specification warning).
#' * `fgls_warning_experimental`: the feature is experimental.
#'
#' @section Error classes:
#' `fgls_error_input` (invalid arguments) and its subclasses
#' `fgls_error_missing`, `fgls_error_time_index`, `fgls_error_collinearity` and
#' `fgls_error_insufficient_data`.
#'
#' @name fgls-conditions
#' @examples
#' # an interior missing value in a time series is an error of class fgls_error_missing
#' df <- data.frame(y = c(1, 2, NA, 4, 5, 3, 2), x = c(2, 1, 3, 5, 4, 2, 1))
#' tryCatch(serial_fgls(y ~ x, data = df), fgls_error = function(e) class(e)[1])
#'
#' # warnings can be caught by class too
#' withCallingHandlers(
#'   serial_fgls(y ~ x - 1, data = simulate_serial(100, seed = 1)),
#'   fgls_warning_specification = function(w) invokeRestart("muffleWarning")
#' )$ar_order
NULL
