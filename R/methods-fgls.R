# S3 methods for "fgls" objects.

#' Methods for fgls fits
#'
#' Standard accessors for objects returned by [fgls()], [serial_fgls()],
#' [lasso_fgls()] and the benchmark estimators. They make the fits work with
#' generic tools such as `lmtest::coeftest()` and `lmtest::coefci()`.
#'
#' Inference is asymptotically normal (as in the papers) unless the fit was made
#' with `use_t = TRUE`. Accordingly, `df.residual()` returns `Inf` for normal
#' inference, so that `lmtest::coeftest()` reproduces [summary.fgls()], and
#' \eqn{n_{eff} - q} otherwise. `nobs()` returns the number of observations in
#' the final regression, \eqn{n_{eff} = T - k^*}.
#'
#' @section Variances and the sandwich package:
#' The variance estimators of the papers are available through
#' `vcov(fit, type = ...)` and can be passed to other tools (for example
#' `lmtest::coeftest(fit, vcov. = vcov(fit, type = "naive"))`). The package
#' provides no `sandwich::estfun()` or `sandwich::bread()` methods, because the
#' robust variance of PG26 eq. (29) and the bias-corrected sandwich of GP24 eq. (8)
#' are not of the generic form those functions assume: `sandwich::vcovHAC()` and
#' related functions would *not* reproduce them.
#'
#' @section Variance types:
#' `vcov(fit, type = t)` returns the stored covariance when `t` is the type the
#' fit was made with (e.g. `"robust"`, or `"HAC"` for [ols_robust()], whose
#' stored `vcov_type` label also shows the kernel and bandwidth). Other types
#' are recomputed from the stored estimates only for fits made by [fgls()],
#' [serial_fgls()] and [lasso_fgls()], and for the F-C stage of a combined fit
#' (`fit$stages$serial`, types `"robust"`, `"HC0"` and `"HC3"`). The benchmark
#' estimators store only the covariance they were fitted with; refit them with
#' another `vcov` instead.
#'
#' @param object,x,formula An `"fgls"` object.
#' @param type For `vcov()` and `confint()`: `NULL` (the variance chosen when
#'   fitting) or another variance type to compute from the stored estimates:
#'   `"robust"`, `"naive"`, `"nonrobust"`, `"HC0"`--`"HC3"` (where
#'   available) or `"bootstrap"` (re-runs the estimator); see the section
#'   *Variance types*. For `residuals()`: `"response"`
#'   (\eqn{y_t - x_t'\hat\beta}, full sample) or `"innovations"` (residuals
#'   of the final filtered and/or weighted regression). For `predict()`: only
#'   `"response"` (the mean of the response) is supported.
#' @param parm Coefficients (names or positions); all by default.
#' @param level Confidence level.
#' @param newdata Optional data frame in which to look for the regressors.
#'   Without it, predictions are for the observations of the fit, padded with
#'   `NA` for the rows dropped by [stats::na.exclude()].
#' @param se.fit Return standard errors of the predicted means.
#' @param interval `"none"` or `"confidence"` (for the mean of the response).
#' @param na.action Missing-value handling for `newdata`.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @return The usual values: a coefficient covariance matrix (`vcov`), a
#'   matrix of confidence limits (`confint`), vectors of residuals or fitted
#'   values, predictions (a vector, or a list/matrix with `se.fit`/`interval`),
#'   or the corresponding model components. The `print()` method returns `x`
#'   invisibly.
#' @name fgls-methods
#' @family fgls
#' @examples
#' sim <- simulate_serial(200, seed = 1)
#' fit <- serial_fgls(y ~ x, data = sim)
#' vcov(fit)
#' vcov(fit, type = "naive")
#' confint(fit)
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)), interval = "confidence")
#' head(residuals(fit, type = "innovations"))
#' nobs(fit)
NULL

#' @rdname fgls-methods
#' @export
vcov.fgls <- function(object, type = NULL, ...) {
  # the stored matrix for its own type: compare with the canonical code, not
  # only with the display label (e.g. "HAC (QS, bw = 2.34)" has code "HAC")
  if (is.null(type) || identical(type, vcov_code_of(object)) ||
    identical(type, object$vcov_type)) {
    return(object$vcov)
  }
  fgls_recompute_vcov(object, type)
}

fgls_recompute_vcov <- function(object, type) {
  # `[[` matches exactly ($ would partially match "stages" on the final fit)
  if (!is.null(object[["stage_name"]])) {
    return(stage_recompute_vcov(object, type))
  }
  spec <- object$spec
  est <- object$internal
  if (is.null(spec) || is.null(est)) {
    fgls_abort(sprintf(
      paste0(
        "vcov(type = \"%s\") is not available: this fit stores only the covariance it was ",
        "fitted with (type = \"%s\"). Other variance types can be recomputed only for fits ",
        "made by fgls(), serial_fgls() or lasso_fgls(); refit with the variance you need."
      ),
      paste(format(type), collapse = ", "), vcov_code_of(object)
    ))
  }
  arg_choice(type, allowed_vcov(spec$mode, spec$serial), "type")
  if (type == "bootstrap") {
    return(fgls_bootstrap(spec, est, object, object$y)$vcov)
  }
  nm <- names(object$coefficients)
  V <- switch(spec$mode,
    serial = {
      if (spec$serial == "ma1") {
        est$cov
      } else {
        fr <- est$filtered
        serial_variance(fr$Xf, fr$ehat, fr$uhat, est$design, J = est$J, cov_type = type)$cov
      }
    },
    hetero = hetero_covariance_from_estimate(est, spec$X, type),
    combined = combined_covariance_from_estimate(est, spec, object$y, type)$cov_fch
  )
  dimnames(V) <- list(nm, nm)
  V
}

# Variance of an intermediate stage of a combined fit. The F-C stage
# (stages$serial) uses the Supplement (A.4)-(A.5) sandwich with unit weights,
# computed by the same code as when fitting from the estimates in
# `stage_internal` (see fc_stage_internal()).
stage_recompute_vcov <- function(object, type) {
  stage <- object[["stage_name"]]
  est <- object[["stage_internal"]]
  if (!identical(stage, "serial") || is.null(est)) {
    fgls_abort(sprintf("vcov(type =) cannot recompute the variance of the \"%s\" stage.", stage))
  }
  if (identical(type, "bootstrap")) {
    fgls_abort(paste0(
      "The bootstrap variance is not available for the F-C stage (fit$stages$serial); ",
      "use serial_fgls(..., vcov = \"bootstrap\") for a serial-correction-only bootstrap."
    ))
  }
  arg_choice(type, combined_cov_types, "type")
  V <- combined_covariance_from_estimate(est, list(X = object$x), object$y, type)$cov_fc
  nm <- names(object$coefficients)
  dimnames(V) <- list(nm, nm)
  V
}

# The estimates of a combined fit that combined_covariances() reads, kept with
# the F-C stage instead of a second copy of all the internals of the fit.
fc_stage_internal <- function(est) {
  fr <- est$serial$filtered
  list(
    serial = list(filtered = fr[c("yf", "Xf", "ehat", "uhat", "n_eff")], order = est$serial$order),
    hetero = est$hetero[c("beta", "v", "df")],
    C = est$C, h_fc = est$h_fc
  )
}

fgls_df <- function(object) {
  if (isTRUE(object$use_t)) object$df_resid else Inf
}

#' @rdname fgls-methods
#' @export
confint.fgls <- function(object, parm, level = 0.95, type = NULL, ...) {
  cf <- object$coefficients
  pnames <- names(cf)
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  a <- (1 - level) / 2
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
  crit <- critical_value(1 - level, fgls_df(object))
  se <- sqrt(pmax(diag(stats::vcov(object, type = type)), 0))[parm]
  ci <- cbind(cf[parm] - crit * se, cf[parm] + crit * se)
  dimnames(ci) <- list(parm, pct)
  ci
}

#' @rdname fgls-methods
#' @export
residuals.fgls <- function(object, type = c("response", "innovations"), ...) {
  type <- arg_default(type, c("response", "innovations"), "type")
  if (type == "innovations") {
    return(object$innovations)
  }
  stats::naresid(object$na.action, object$residuals)
}

#' @rdname fgls-methods
#' @export
fitted.fgls <- function(object, ...) {
  stats::napredict(object$na.action, object$fitted.values)
}

#' @rdname fgls-methods
#' @export
predict.fgls <- function(object, newdata, se.fit = FALSE, interval = c("none", "confidence"),
                         level = 0.95, type = "response", na.action = stats::na.pass, ...) {
  if (!identical(type, "response")) {
    fgls_abort(sprintf(
      "predict() supports only type = \"response\" for fgls fits; got %s.",
      paste0("\"", format(type), "\"", collapse = ", ")
    ))
  }
  interval <- arg_default(interval, c("none", "confidence"), "interval")
  if (missing(newdata) || is.null(newdata)) {
    X <- object$x
    offset_na <- object$na.action
  } else {
    tt <- stats::delete.response(object$terms)
    m <- stats::model.frame(tt, newdata, na.action = na.action, xlev = object$xlevels)
    X <- stats::model.matrix(tt, m, contrasts.arg = object$contrasts)
    offset_na <- NULL
  }
  if (!identical(colnames(X), names(object$coefficients))) {
    fgls_abort("The regressors in `newdata` do not match the fitted model.")
  }
  # without newdata, pad every component for the rows dropped by na.exclude,
  # as predict.lm() does
  pad <- function(v) if (is.null(offset_na)) v else stats::napredict(offset_na, v)
  if (!se.fit && interval == "none") {
    return(pad(drop(X %*% object$coefficients)))
  }
  p <- predict_mean(object, X, level)
  df <- fgls_df(object)
  if (interval == "confidence") {
    fit <- pad(cbind(fit = p$fit, lwr = p$lwr, upr = p$upr))
    if (!se.fit) {
      return(fit)
    }
    return(list(fit = fit, se.fit = pad(p$se), df = df))
  }
  list(fit = pad(p$fit), se.fit = pad(p$se), df = df)
}

# Predicted means x'b, their standard errors (x' V x)^{1/2} and confidence
# limits (normal, or Student-t with use_t = TRUE) for the rows of X.
predict_mean <- function(object, X, level = 0.95) {
  fit <- drop(X %*% object$coefficients)
  se <- sqrt(pmax(rowSums((X %*% object$vcov) * X), 0))
  crit <- critical_value(1 - level, fgls_df(object))
  list(fit = fit, se = se, lwr = fit - crit * se, upr = fit + crit * se)
}

#' @rdname fgls-methods
#' @export
nobs.fgls <- function(object, ...) {
  object$n_eff
}

#' @rdname fgls-methods
#' @export
df.residual.fgls <- function(object, ...) {
  fgls_df(object)
}

#' @rdname fgls-methods
#' @export
formula.fgls <- function(x, ...) {
  x$formula
}

#' @rdname fgls-methods
#' @export
model.matrix.fgls <- function(object, ...) {
  object$x
}

#' @rdname fgls-methods
#' @export
model.frame.fgls <- function(formula, ...) {
  if (!is.null(formula$model)) {
    return(formula$model)
  }
  fgls_abort("The model frame was not kept; refit with model = TRUE.")
}

#' @rdname fgls-methods
#' @export
terms.fgls <- function(x, ...) {
  x$terms
}

#' @rdname fgls-methods
#' @export
print.fgls <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  }
  cat("\n", fgls_headline(x), "\n", sep = "")
  cat("\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

fgls_headline <- function(x) {
  parts <- x$method
  s <- x$serial
  if (!is.null(s)) {
    if (identical(s$method, "ma1")) {
      parts <- c(parts, sprintf("MA(1) theta = %.4g", s$extra$theta))
    } else {
      how <- if (s$method %in% c("bic", "aic")) toupper(s$method) else s$method
      rho <- if (length(s$rho)) {
        paste0(", rho = ", paste(formatC(s$rho, digits = 4, format = "fg"), collapse = ", "))
      } else {
        ""
      }
      parts <- c(parts, sprintf("AR order k* = %d (%s)%s", s$ar_order, how, rho))
    }
  }
  h <- x$skedastic
  if (!is.null(h)) {
    sel <- if (length(h$selected)) paste(h$selected, collapse = ", ") else "none"
    parts <- c(parts, sprintf("skedastic fit: %s (selected: %s)", h$method, sel))
  }
  paste0(paste(parts, collapse = "; "), "; covariance: ", x$vcov_type)
}
