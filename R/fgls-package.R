#' fgls: Feasible GLS for serial correlation and heteroskedasticity
#'
#' Feasible generalized least squares (FGLS) for linear regressions with
#' serially correlated errors (Perron and González-Coya, 2026), heteroskedastic
#' errors of unknown form (González-Coya and Perron, 2024), or both
#' (Supplement S-4 of the 2026 paper).
#'
#' @references
#' Perron, P. and González-Coya, E. (2026). Feasible GLS for time series
#' regression. *Journal of Econometrics*, 258, 106346.
#' \doi{10.1016/j.jeconom.2026.106346}
#'
#' González-Coya, E. and Perron, P. (2024). Estimation in the presence of
#' heteroskedasticity of unknown form: a Lasso-based approach. *Journal of
#' Econometric Methods*, 13(1), 29--48. \doi{10.1515/jem-2023-0007}
#'
#' @keywords internal
#' @importFrom glmnet glmnet
#' @importFrom graphics par
#' @importFrom grDevices dev.interactive
#' @importFrom stats coef confint df.residual fitted formula model.frame model.matrix
#' @importFrom stats nobs predict residuals terms vcov
"_PACKAGE"
