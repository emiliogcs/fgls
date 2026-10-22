#' Feasible GLS for serial correlation and/or heteroskedasticity
#'
#' `fgls()` fits a linear regression by feasible generalized least squares,
#' correcting for serially correlated errors (Perron and González-Coya, 2026),
#' for heteroskedasticity of unknown form (González-Coya and Perron, 2024), or
#' for both (the "F-CH" estimator of Supplement S-4 of the 2026 paper).
#' [serial_fgls()] and [lasso_fgls()] are shortcuts for the two single
#' corrections.
#'
#' @section Serial correlation (PG26):
#' The errors are approximated by an AR(k) process. Its coefficients are the
#' coefficients on the lagged dependent variable in the Durbin regression of
#' \eqn{y_t} on \eqn{y_{t-1}, \dots, y_{t-k}}, \eqn{x_t, \dots, x_{t-k}}
#' (eq. 14), and \eqn{k} is chosen by BIC over \eqn{0, \dots, k_{max}} on a
#' common sample. The data are then quasi-differenced and the model is
#' estimated by OLS on the filtered data (eq. 15), dropping the first \eqn{k}
#' observations. Unlike Cochrane-Orcutt, this is consistent when the
#' regressors are only *pre-determined* (uncorrelated with current and future
#' innovations) rather than strictly exogenous. The default variance
#' (`vcov = "robust"`, eqs. 28--29) accounts for the estimation of the filter
#' and is valid in both cases; `vcov = "naive"` (eq. 27) is valid only for
#' exogenous regressors. The share by which the robust variance exceeds the
#' naive one is reported as `nonexogeneity_share`. If the selected order is 0,
#' FGLS equals OLS.
#'
#' @section Heteroskedasticity (GP24):
#' The log of the squared OLS residuals (floored at `delta^2`) is regressed on
#' candidate covariates \eqn{z_t} by an adaptive Lasso: a ridge first step
#' gives weights \eqn{|\hat\gamma_j|^{-\psi}}, and the penalty \eqn{\lambda} and
#' the exponent \eqn{\psi} are chosen jointly by 10-fold cross-validation with
#' the same folds for every \eqn{\psi}. The exponentiated fit gives the
#' weights of a weighted least-squares regression. Inference uses a sandwich
#' variance with the bias correction of eq. (8), so it stays valid when the
#' skedastic model is wrong. If the Lasso selects no covariate the weights are
#' constant and FGLS equals OLS (no loss under homoskedasticity).
#'
#' @section Both (Supplement S-4):
#' With `serial != "none"` and `hetero != "none"` the serial correction is
#' applied first ("F-C"), the skedastic function is fitted to the residuals of
#' the quasi-differenced regression with candidate covariates built from the
#' *unfiltered* regressors, and the quasi-differenced regression is
#' re-estimated by weighted least squares ("F-CH"). The F-C stage is kept in
#' `fit$stages$serial`.
#'
#' @section Choosing options:
#' * Time series with possibly autocorrelated errors and no concern about
#'   heteroskedasticity: `hetero = "none"` (or [serial_fgls()]).
#' * Cross-sectional data: `serial = "none"` (or [lasso_fgls()]). Rows can be
#'   in any order.
#' * Lagged dependent variables among the regressors must be declared in
#'   `lagged_dep_vars`; the AR coefficients are then recovered by the Wallis
#'   ratio (PG26 Section 4.2).
#' * Errors known to be MA(1), e.g. multi-step forecast errors: `serial =
#'   "ma1"` (GMA estimator, PG26 Section 5).
#' * With several stochastic regressors the BIC that penalizes only the AR
#'   order can select `kmax`; use `control = fgls_control(ic_penalty =
#'   "params")`.
#' * Results involving cross-validation depend on the random folds: set `seed`
#'   for reproducibility.
#'
#' @param formula A model formula, as in [stats::lm()]. Include an intercept
#'   unless the model truly has none.
#' @param data A data frame, a list, a `ts`/`mts` object, or a 'zoo'/'xts'
#'   object containing the variables in the formula. Rows must be in time
#'   order (unless `serial = "none"`).
#' @param subset Optional expression selecting rows, evaluated in `data`.
#'   For time series only a contiguous block makes sense.
#' @param na.action [stats::na.fail()] (default) stops on missing values;
#'   [stats::na.omit()] or [stats::na.exclude()] drop them. For time series only
#'   leading and trailing missing values can be dropped.
#' @param serial Serial-correlation correction: `"auto"` (AR order by BIC,
#'   default), `"fixed"` (use `ar_order`), `"ma1"` (GMA for MA(1) errors) or
#'   `"none"`.
#' @param hetero Skedastic-function estimator: `"lasso"` (adaptive Lasso,
#'   default), `"svr"` (support vector regression, requires 'e1071'),
#'   `"wls_s1"`/`"wls_s2"` (the parametric fits of Romano and Wolf, 2017, on
#'   \eqn{\log z} or \eqn{z}), `"known"` (variances given in `variance`) or
#'   `"none"`.
#' @param vcov Variance estimator: `"robust"` (default; PG26 eq. 29, GP24
#'   eqs. 6+8, or Supplement A.4--A.5 depending on the corrections), `"naive"`
#'   (PG26 eq. 27; serial only), `"nonrobust"` (heteroskedasticity only),
#'   `"HC0"`--`"HC3"`, or `"bootstrap"`. See [vcov.fgls()].
#' @param z Candidate covariates of the skedastic function: `"auto"` (default)
#'   builds \eqn{\{x, x^2, \log x, (\log x)^2, \cos x, \cos 2x\}} from every
#'   non-constant regressor (see [make_dictionary()]); a character vector of
#'   regressor names; a one-sided formula evaluated in `data` (e.g. `~ x +
#'   I(x^2) + w`); or a numeric matrix / data frame with one row per row of
#'   `data`.
#' @param ar_order Fixed AR order (implies `serial = "fixed"`).
#' @param kmax Largest AR order considered: `NULL` (default,
#'   \eqn{\lfloor 12 (T/200)^{1/4} \rfloor}), an integer, or `"theory"`
#'   (\eqn{\lfloor (T/\ln T)^{1/3} \rfloor}). It is capped so that the Durbin
#'   regression uses at most half of the sample as regressors.
#' @param ic Lag-order criterion: `"bic"` (default) or `"aic"`.
#' @param lagged_dep_vars Names (or positions) of model-matrix columns that are
#'   lags of the dependent variable.
#' @param wallis_regressor Regressor used for the Wallis recovery of the AR
#'   coefficients: `"auto"` (largest Durbin \eqn{|t|}), `"median"`, or a name.
#' @param exog Assumption on the regressors: `"predetermined"` (default) or
#'   `"exogenous"`. Estimates do not depend on it; it sets the default of
#'   `first_stage_correction` and the notes of the summary.
#' @param first_stage_correction Add the first-stage (Durbin) term of PG26
#'   eq. (29) to the variance of the combined estimator. Default `TRUE` when
#'   `exog = "predetermined"`.
#' @param variance Known error variances for `hetero = "known"`: a numeric
#'   vector with one value per observation, or a function of the model matrix.
#' @param delta Floor in \eqn{\log\max(\hat u_t^2, \delta^2)}, in the units of
#'   the response (default 0.1, as in the paper).
#' @param psi_grid Exponents \eqn{\psi} of the adaptive weights; 0 is the
#'   plain Lasso.
#' @param cv_folds Number of cross-validation folds (default 10).
#' @param lambda_rule Penalty at the minimum of the cross-validation error
#'   (`"min"`, default) or the largest penalty within one standard error of it
#'   (`"1se"`).
#' @param weight_clip Optional quantiles `c(lo, hi)` at which the fitted
#'   log-variance is winsorized, e.g. `c(0.01, 0.99)`.
#' @param seed Optional seed for the cross-validation folds and the bootstrap.
#'   The session's random number generator state is restored afterwards.
#' @param use_t Use Student-t (`n_eff - q` degrees of freedom) instead of
#'   normal inference.
#' @param time Time index for data frames: the name of a `Date`, `POSIXct` or
#'   numeric column of `data`, or a vector with one value per row. It is used
#'   to check that rows are sorted and have no gaps.
#' @param allow_gaps Accept a time index with gaps (lags are then positional).
#' @param check_stationarity Screen the response and the regressors with
#'   augmented Dickey-Fuller tests and warn when a unit root is not rejected.
#' @param control Advanced options, see [fgls_control()].
#' @param contrasts Optional list of contrasts for factors, as in
#'   [stats::lm()].
#' @param model Keep the model frame in the result.
#' @param weights Not supported: FGLS estimates its own weights. Use `hetero
#'   = "known"` with `variance` for known variances.
#' @param ... Not used; unknown arguments are an error.
#'
#' @return An object of class `"fgls"`, a list with components including
#' \describe{
#'   \item{`coefficients`, `vcov`}{Estimates and their covariance matrix.}
#'   \item{`residuals`, `fitted.values`}{\eqn{y_t - x_t'\hat\beta} and
#'     \eqn{x_t'\hat\beta} on the full sample.}
#'   \item{`innovations`}{Residuals of the final (filtered and/or weighted)
#'     regression.}
#'   \item{`ar_order`, `rho`, `bic_path`, `nonexogeneity_share`}{Serial block
#'     (`NULL` without serial correction); more in `serial`.}
#'   \item{`skedastic`, `weights`}{Heteroskedasticity block: selected
#'     covariates, coefficients, \eqn{\psi}, \eqn{\lambda}, degrees of freedom,
#'     fitted variances; `weights` are \eqn{1/\hat v_t}. For the Lasso,
#'     `skedastic$lambda` is the selected penalty on the scale of 'glmnet',
#'     whose objective is that of GP24 eq. (4) divided by the number of
#'     observations \eqn{n} of the skedastic regression, and
#'     `skedastic$lambda_paper` \eqn{= n \lambda} is the same penalty on the
#'     scale of eq. (4). Both refer to the standardized candidate covariates.}
#'   \item{`stages`}{Intermediate estimators (`stages$serial` is F-C). They
#'     are `"fgls"` objects with the same methods as the final fit;
#'     `vcov(fit$stages$serial, type =)` accepts `"robust"`, `"HC0"` and
#'     `"HC3"`.}
#'   \item{`nobs`, `n_eff`}{Number of observations and of observations in the
#'     final regression (\eqn{T - k^*}).}
#' }
#' Use the methods [summary()], [coef()], [vcov()], [confint()],
#' [predict()], [residuals()], [plot()] and [compare_fgls()].
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
#' Romano, J. P. and Wolf, M. (2017). Resurrecting weighted least squares.
#' *Journal of Econometrics*, 197(1), 1--19.
#' \doi{10.1016/j.jeconom.2016.10.003}
#' @family fgls
#' @export
#' @examples
#' # Serial correlation with a non-exogenous regressor (PG26 Table 1 design)
#' sim <- simulate_serial(300, seed = 1)
#' fit <- fgls(y ~ x, data = sim, hetero = "none")
#' summary(fit)
#'
#' # Heteroskedasticity of unknown form (GP24)
#' het <- simulate_hetero(200, skedastic = "log2", seed = 2)
#' fit2 <- fgls(y ~ x, data = het, serial = "none", seed = 1)
#' fit2$skedastic$selected
#'
#' # Both (Supplement S-4)
#' both <- simulate_combined(200, ar = 0.5, skedastic = "x2", seed = 3)
#' fit3 <- fgls(y ~ x, data = both, seed = 1)
#' compare_fgls(fit3$stages$serial, fit3)
fgls <- function(formula, data, subset, na.action = stats::na.fail,
                 serial = c("auto", "fixed", "ma1", "none"),
                 hetero = c("lasso", "svr", "wls_s1", "wls_s2", "known", "none"),
                 vcov = c("robust", "naive", "nonrobust", "HC0", "HC1", "HC2", "HC3", "bootstrap"),
                 z = "auto", ar_order = NULL, kmax = NULL, ic = c("bic", "aic"),
                 lagged_dep_vars = NULL, wallis_regressor = "auto",
                 exog = c("predetermined", "exogenous"), first_stage_correction = NULL,
                 variance = NULL, delta = 0.1, psi_grid = c(0, 0.25, 0.5, 0.75, 1, 2),
                 cv_folds = 10L, lambda_rule = c("min", "1se"), weight_clip = NULL,
                 seed = NULL, use_t = FALSE, time = NULL, allow_gaps = FALSE,
                 check_stationarity = FALSE, control = fgls_control(), contrasts = NULL,
                 model = TRUE, weights, ...) {
  cl <- match.call()
  if (...length()) {
    nms <- names(list(...))
    nms <- if (is.null(nms)) "(unnamed)" else ifelse(nms == "", "(unnamed)", nms)
    fgls_abort(sprintf("Unknown argument(s): %s.", paste(nms, collapse = ", ")))
  }
  if (!missing(weights)) {
    fgls_abort(paste0(
      "`weights` is not supported: FGLS estimates its own weights. For known error ",
      "variances use hetero = \"known\" with `variance`."
    ))
  }
  serial <- arg_default(serial, c("auto", "fixed", "ma1", "none"), "serial")
  hetero <- arg_default(hetero, c("lasso", "svr", "wls_s1", "wls_s2", "known", "none"), "hetero")
  vcov_type <- arg_default(
    vcov, c("robust", "naive", "nonrobust", "HC0", "HC1", "HC2", "HC3", "bootstrap"), "vcov"
  )
  ic <- arg_default(ic, c("bic", "aic"), "ic")
  exog <- arg_default(exog, c("predetermined", "exogenous"), "exog")
  lambda_rule <- arg_default(lambda_rule, c("min", "1se"), "lambda_rule")
  if (!inherits(control, "fgls_control")) {
    fgls_abort("`control` must be created by fgls_control().")
  }
  check_seed(seed)
  if (!is.null(ar_order) && serial == "auto") serial <- "fixed"
  if (serial == "fixed" && is.null(ar_order)) {
    fgls_abort("serial = \"fixed\" needs ar_order = <integer>.")
  }
  if (serial == "ma1" && hetero != "none") {
    fgls_abort(paste0(
      "serial = \"ma1\" cannot be combined with a heteroskedasticity correction; ",
      "use hetero = \"none\"."
    ))
  }
  if (hetero == "known" && is.null(variance)) {
    fgls_abort("hetero = \"known\" needs the error variances in `variance`.")
  }
  if (!(is.numeric(delta) && length(delta) == 1L && delta > 0)) {
    fgls_abort("`delta` must be a positive number.")
  }
  if (!is_count(cv_folds) || cv_folds < 2L) fgls_abort("`cv_folds` must be an integer >= 2.")
  if (!is.numeric(psi_grid) || !length(psi_grid) || any(psi_grid < 0)) {
    fgls_abort("`psi_grid` must contain non-negative values.")
  }
  if (!is.null(weight_clip) &&
    !(is.numeric(weight_clip) && length(weight_clip) == 2L &&
      weight_clip[1L] >= 0 && weight_clip[1L] < weight_clip[2L] && weight_clip[2L] <= 1)) {
    fgls_abort("`weight_clip` must be c(lo, hi) with 0 <= lo < hi <= 1.")
  }
  mode <- if (serial == "none") "hetero" else if (hetero == "none") "serial" else "combined"
  allowed <- allowed_vcov(mode, serial)
  if (!(vcov_type %in% allowed)) {
    fgls_abort(sprintf(
      "vcov = \"%s\" is not available for this model; choose from %s.", vcov_type,
      paste0("\"", allowed, "\"", collapse = ", ")
    ))
  }
  time_series <- serial != "none"
  use_z <- hetero %in% c("lasso", "svr", "wls_s1", "wls_s2")

  frame <- fgls_frame(
    formula,
    data = if (missing(data)) NULL else data,
    subset_expr = if (missing(subset)) NULL else substitute(subset),
    env = parent.frame(), time = time, policy = na_policy(na.action),
    time_series = time_series, z = if (use_z) z else "auto", contrasts = contrasts
  )
  d <- prepare_data(frame$y, frame$X,
    y_name = frame$y_name, index = frame$index,
    time_series = time_series, allow_gaps = allow_gaps
  )
  lagged_dep_idx <- resolve_columns(lagged_dep_vars, d$x_names, "lagged_dep_vars")
  if (use_z) {
    cand <- resolve_candidates(z, d, frame, transforms = control$transforms, skip = lagged_dep_idx)
  } else {
    cand <- list(Z = matrix(0, d$nobs, 0L), names = character(0))
  }
  if (time_series && length(lagged_dep_idx) == 0L) {
    warn_lagged_dependent(detect_lagged_dependent(d$y, d$X, d$x_names))
  }
  if (isTRUE(check_stationarity) && time_series) screen_stationarity(d)
  if (hetero == "known") variance <- resolve_known_variance(variance, d)

  spec <- list(
    mode = mode, serial = serial, hetero = hetero, X = d$X, x_names = d$x_names,
    y_name = d$y_name, has_const = d$has_const, Z = cand$Z, z_names = cand$names,
    lagged_dep_idx = lagged_dep_idx,
    options = list(
      ar_order = ar_order, kmax = kmax, ic = ic, wallis_regressor = wallis_regressor,
      exog = exog,
      first_stage_correction = if (is.null(first_stage_correction)) {
        exog == "predetermined"
      } else {
        isTRUE(first_stage_correction)
      },
      variance = variance, delta = delta, psi_grid = sort(unique(as.numeric(psi_grid))),
      cv_folds = as.integer(cv_folds), lambda_rule = lambda_rule, weight_clip = weight_clip,
      seed = seed, use_t = isTRUE(use_t), control = control, z = z
    )
  )
  base_cov <- if (vcov_type == "bootstrap") "robust" else vcov_type
  est <- fgls_estimate(spec, d$y, base_cov)
  res <- fgls_build(spec, est, d$y, base_cov, frame$rownames)
  if (vcov_type == "bootstrap") res <- fgls_bootstrap(spec, est, res, d$y)

  res$spec <- spec
  res$internal <- est
  attach_model_info(res, frame, d, cl, model)
}

# Model bookkeeping shared by fgls() and the benchmark estimators. The stages of
# a combined fit (stages$serial is F-C) receive the same model information, and
# the F-C stage the estimates its variance needs, so that model.frame(),
# augment(), predict() and vcov(type =) work on them as on the final fit.
attach_model_info <- function(res, frame, d, cl, model = TRUE) {
  res$call <- cl
  res$formula <- stats::formula(frame$terms)
  res$terms <- frame$terms
  if (isTRUE(model)) res$model <- frame$mf
  res$x <- d$X
  res$y <- stats::setNames(d$y, frame$rownames)
  res$contrasts <- frame$contrasts
  res$xlevels <- frame$xlevels
  res$na.action <- frame$na.action
  res$index <- d$index
  res$vcov_code <- vcov_code_of(res)
  shared <- c(
    "call", "formula", "terms", "model", "x", "y", "contrasts", "xlevels", "na.action", "index"
  )
  shared <- intersect(shared, names(res))
  for (nm in names(res$stages)) {
    res$stages[[nm]][shared] <- res[shared]
    res$stages[[nm]]$stage_name <- nm
    res$stages[[nm]]$vcov_code <- vcov_code_of(res$stages[[nm]])
  }
  if (!is.null(res$stages[["serial"]]) && identical(res$spec$mode, "combined")) {
    res$stages$serial$stage_internal <- fc_stage_internal(res$internal)
  }
  res
}

allowed_vcov <- function(mode, serial = "auto") {
  if (mode == "serial" && serial == "ma1") {
    return(c("robust", "bootstrap"))
  }
  switch(mode,
    serial = c(serial_cov_types, "bootstrap"),
    hetero = c(hetero_cov_types, "bootstrap"),
    combined = c(combined_cov_types, "bootstrap")
  )
}

# Run the estimator selected by `spec` on the response y.
fgls_estimate <- function(spec, y, cov_type) {
  o <- spec$options
  ctl <- o$control
  serial_args <- list(
    x_names = spec$x_names, y_name = spec$y_name, has_const = spec$has_const,
    ar_order = o$ar_order, kmax = o$kmax, ic = o$ic, min_lags = ctl$min_lags,
    lagged_dep_idx = spec$lagged_dep_idx, wallis_regressor = o$wallis_regressor,
    ic_penalty = ctl$ic_penalty
  )
  hetero_args <- list(
    delta = o$delta, delta_scale = ctl$delta_scale, weight_clip = o$weight_clip,
    psi_grid = o$psi_grid, cv_folds = o$cv_folds, lambda_rule = o$lambda_rule,
    control = ctl, seed = o$seed
  )
  if (spec$mode == "serial") {
    if (spec$serial == "ma1") {
      return(estimate_gma(y, spec$X, spec$x_names))
    }
    return(do.call(estimate_serial, c(list(y = y, X = spec$X, cov_type = cov_type), serial_args)))
  }
  if (spec$mode == "hetero") {
    return(do.call(estimate_hetero, c(
      list(
        y = y, X = spec$X, Z = spec$Z, z_names = spec$z_names, method = spec$hetero,
        known = o$variance, cov_type = cov_type
      ),
      hetero_args
    )))
  }
  do.call(estimate_combined, c(
    list(
      y = y, X = spec$X, Z = spec$Z, z_names = spec$z_names, method = spec$hetero,
      known = o$variance, first_stage_correction = o$first_stage_correction,
      cov_type = cov_type
    ),
    serial_args, hetero_args
  ))
}

# ------------------------------------------------------------------- model frame
as_fgls_data <- function(data) {
  index <- NULL
  if (inherits(data, "zoo")) {
    if (!requireNamespace("zoo", quietly = TRUE)) {
      fgls_abort("Package 'zoo' is needed for zoo/xts data.")
    }
    index <- zoo::index(data)
    cd <- zoo::coredata(data)
    if (is.null(dim(cd)) || is.null(colnames(cd))) {
      fgls_abort("A zoo/xts `data` object must have named columns.")
    }
    df <- as.data.frame(cd)
  } else if (stats::is.ts(data)) {
    if (is.null(dim(data)) || is.null(colnames(data))) {
      fgls_abort("A ts `data` object must be a multivariate series with named columns.")
    }
    index <- as.numeric(stats::time(data))
    df <- as.data.frame(matrix(as.numeric(data), nrow(data), dimnames = list(NULL, colnames(data))))
  } else if (is.data.frame(data)) {
    df <- data
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.list(data) || is.environment(data)) {
    df <- data
  } else {
    fgls_abort("`data` must be a data frame, a list, a ts object or a zoo/xts object.")
  }
  list(df = df, index = index)
}

index_from_env <- function(formula) {
  lhs <- tryCatch(eval(formula[[2L]], environment(formula)), error = function(e) NULL)
  if (inherits(lhs, "zoo") && requireNamespace("zoo", quietly = TRUE)) {
    return(zoo::index(lhs))
  }
  if (stats::is.ts(lhs)) {
    return(as.numeric(stats::time(lhs)))
  }
  NULL
}

fgls_frame <- function(formula, data, subset_expr, env, time, policy, time_series, z,
                       contrasts = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    fgls_abort("`formula` must be a two-sided formula such as y ~ x.")
  }
  index <- NULL
  df <- NULL
  if (!is.null(data)) {
    conv <- as_fgls_data(data)
    df <- conv$df
    index <- conv$index
  }
  mf <- if (is.null(df)) {
    stats::model.frame(formula, na.action = stats::na.pass, drop.unused.levels = TRUE)
  } else {
    stats::model.frame(formula, data = df, na.action = stats::na.pass, drop.unused.levels = TRUE)
  }
  mt <- attr(mf, "terms")
  n0 <- nrow(mf)
  if (!is.null(time)) {
    if (is.character(time) && length(time) == 1L) {
      if (is.null(df) || !(time %in% names(df))) {
        fgls_abort(sprintf("`time` column '%s' not found in `data`.", time))
      }
      index <- df[[time]]
    } else {
      index <- time
    }
  } else if (is.null(index) && is.null(df)) {
    index <- index_from_env(formula)
  }
  if (!is.null(index) && length(index) != n0) {
    fgls_abort(sprintf(
      "The time index has %d values but the data have %d rows.", length(index), n0
    ))
  }
  # candidate covariates given as a formula or a matrix share the row handling
  zmat <- NULL
  if (inherits(z, "formula")) {
    zf <- if (is.null(df)) {
      stats::model.frame(z, na.action = stats::na.pass)
    } else {
      stats::model.frame(z, data = df, na.action = stats::na.pass)
    }
    zmat <- stats::model.matrix(attr(zf, "terms"), zf)
    zmat <- zmat[, colnames(zmat) != "(Intercept)", drop = FALSE]
  } else if (is.matrix(z) || is.data.frame(z)) {
    zmat <- as.matrix(z)
    if (!is.numeric(zmat)) fgls_abort("A matrix or data frame `z` must be numeric.")
    if (is.null(colnames(zmat))) colnames(zmat) <- sprintf("z%d", seq_len(ncol(zmat)))
  }
  if (!is.null(zmat) && nrow(zmat) != n0) {
    fgls_abort(sprintf("`z` has %d rows but the data have %d rows.", nrow(zmat), n0))
  }
  rows <- seq_len(n0)
  if (!is.null(subset_expr)) {
    s <- if (is.null(df)) eval(subset_expr, envir = env) else eval(subset_expr, df, env)
    if (is.logical(s)) {
      if (length(s) != n0) fgls_abort("`subset` must have one value per row of the data.")
      rows <- which(s & !is.na(s))
    } else if (is.numeric(s)) {
      rows <- seq_len(n0)[s]
    } else {
      fgls_abort("`subset` must be a logical or numeric vector.")
    }
  }
  rn <- rownames(mf)
  mf <- mf[rows, , drop = FALSE]
  if (!is.null(index)) index <- index[rows]
  if (!is.null(zmat)) zmat <- zmat[rows, , drop = FALSE]

  na_rows <- !stats::complete.cases(mf)
  if (!is.null(zmat)) na_rows <- na_rows | !stats::complete.cases(zmat)
  keep <- missing_rows_keep(na_rows, policy, time_series, index)
  na_act <- NULL
  if (!all(keep)) {
    omitted <- rows[!keep]
    names(omitted) <- rn[omitted]
    na_act <- structure(omitted, class = if (policy == "exclude") "exclude" else "omit")
    mf <- mf[keep, , drop = FALSE]
    if (!is.null(index)) index <- index[keep]
    if (!is.null(zmat)) zmat <- zmat[keep, , drop = FALSE]
  }
  attr(mf, "terms") <- mt
  if (!is.null(stats::model.offset(mf))) {
    fgls_abort("Offsets are not supported; subtract the offset from the response instead.")
  }
  yraw <- stats::model.response(mf)
  if (is.null(yraw)) fgls_abort("The formula has no response.")
  if (is.factor(yraw) || is.character(yraw) || is.logical(yraw)) {
    fgls_abort("The response must be numeric.")
  }
  if (is.matrix(yraw) && ncol(yraw) > 1L) fgls_abort("The response must be a single column.")
  y <- as.numeric(yraw)
  X <- stats::model.matrix(mt, mf, contrasts.arg = contrasts)
  if (nrow(X) != length(y)) fgls_abort("The model matrix and the response do not match.")
  list(
    y = y, X = X, index = index, mf = mf, terms = mt,
    xlevels = stats::.getXlevels(mt, mf), contrasts = attr(X, "contrasts"),
    zmat = zmat, na.action = na_act, rownames = rownames(mf),
    y_name = paste(deparse(formula[[2L]], width.cutoff = 500L), collapse = " ")
  )
}

resolve_known_variance <- function(variance, d) {
  vals <- if (is.function(variance)) variance(d$X) else variance
  vals <- as.numeric(vals)
  if (length(vals) != d$nobs) {
    fgls_abort(sprintf("Known variances have length %d; expected %d.", length(vals), d$nobs))
  }
  if (any(!is.finite(vals)) || any(vals <= 0)) {
    fgls_abort("Known variances must be positive and finite.")
  }
  vals
}
