#' Feasible GLS for serially correlated errors (Perron and González-Coya, 2026)
#'
#' Shortcut for [fgls()] with `hetero = "none"`: the errors are approximated by
#' an AR(k) process whose coefficients come from the Durbin regression and
#' whose order is chosen by BIC (eq. 14), the data are quasi-differenced and
#' OLS is applied (eq. 15). The default variance is the exogeneity-robust
#' estimator of eqs. (28)--(29), valid when the regressors are only
#' pre-determined.
#'
#' @inheritParams fgls
#' @param serial `"auto"` (default), `"fixed"` or `"ma1"` (GMA estimator for
#'   MA(1) errors, PG26 Section 5; its only variance is `vcov = "robust"`, the
#'   GMA covariance \eqn{s^2 (X^{*\prime} X^*)^{-1}}).
#' @param vcov `"robust"` (default, eq. 29), `"naive"` (eq. 27, exogenous
#'   regressors only), `"HC0"`--`"HC3"` (eq. 29 with heteroskedasticity
#'   weights) or `"bootstrap"` (sieve bootstrap).
#' @param ... Further arguments passed to [fgls()], e.g. `wallis_regressor`,
#'   `use_t`, `time`, `control`.
#' @return An object of class `"fgls"`; see [fgls()].
#' @inherit fgls references
#' @family fgls
#' @export
#' @examples
#' sim <- simulate_serial(500, seed = 42)
#' fit <- serial_fgls(y ~ x, data = sim)
#' fit
#' fit$ar_order
#' fit$nonexogeneity_share
#'
#' # Lagged dependent variable among the regressors (Wallis recovery)
#' df <- data.frame(y = sim$y[-1], ylag = sim$y[-500], x = sim$x[-1])
#' serial_fgls(y ~ ylag + x, data = df, lagged_dep_vars = "ylag")
serial_fgls <- function(formula, data, subset, na.action = stats::na.fail,
                        serial = c("auto", "fixed", "ma1"),
                        vcov = c("robust", "naive", "HC0", "HC1", "HC2", "HC3", "bootstrap"),
                        ar_order = NULL, kmax = NULL, ic = c("bic", "aic"),
                        lagged_dep_vars = NULL, exog = c("predetermined", "exogenous"), ...) {
  cl <- match.call()
  cl[[1L]] <- quote(fgls::fgls)
  cl$serial <- arg_default(serial, c("auto", "fixed", "ma1"), "serial")
  cl$hetero <- "none"
  res <- eval(cl, parent.frame())
  res$call <- match.call()
  res
}
