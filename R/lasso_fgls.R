#' Feasible GLS under heteroskedasticity of unknown form (González-Coya and Perron, 2024)
#'
#' Shortcut for [fgls()] with `serial = "none"`. The log of the squared OLS
#' residuals is regressed on candidate covariates `z` by an adaptive Lasso tuned
#' by cross-validation; the exponentiated fit gives the weights of a weighted
#' least-squares regression. Inference uses the sandwich of eq. (6) with the
#' bias correction of eq. (8), so it remains valid when the skedastic model is
#' misspecified.
#'
#' @inheritParams fgls
#' @param hetero Skedastic-function estimator: `"lasso"` (default), `"svr"`,
#'   `"wls_s1"`, `"wls_s2"`, `"known"` or `"none"` (OLS with HC3 standard
#'   errors).
#' @param vcov `"robust"` (default; eq. 8 for Lasso/SVR, eq. 7 for WLS-S1/S2,
#'   HC3 otherwise), `"nonrobust"`, `"HC0"`--`"HC3"` or `"bootstrap"` (wild
#'   bootstrap).
#' @param ... Further arguments passed to [fgls()], e.g. `variance`,
#'   `cv_folds`, `lambda_rule`, `weight_clip`, `use_t`, `control`.
#' @return An object of class `"fgls"`; see [fgls()].
#' @inherit fgls references
#' @family fgls
#' @family heteroskedasticity
#' @export
#' @examples
#' het <- simulate_hetero(200, skedastic = "x2", seed = 1)
#' fit <- lasso_fgls(y ~ x, data = het, seed = 1)
#' fit
#' fit$skedastic$selected
#' ols <- lasso_fgls(y ~ x, data = het, hetero = "none")
#' compare_fgls(OLS = ols, FGLS = fit)
lasso_fgls <- function(formula, data, subset, na.action = stats::na.fail,
                       hetero = c("lasso", "svr", "wls_s1", "wls_s2", "known", "none"),
                       vcov = c("robust", "nonrobust", "HC0", "HC1", "HC2", "HC3", "bootstrap"),
                       z = "auto", seed = NULL, ...) {
  cl <- match.call()
  cl[[1L]] <- quote(fgls::fgls)
  cl$hetero <- arg_default(hetero, c("lasso", "svr", "wls_s1", "wls_s2", "known", "none"), "hetero")
  cl$serial <- "none"
  res <- eval(cl, parent.frame())
  res$call <- match.call()
  res
}
