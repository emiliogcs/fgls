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
#' @param object,x,formula An `"fgls"` object.
#' @param type For `vcov()` and `confint()`: `NULL` (the variance chosen when
#'   fitting) or another variance type to compute from the stored estimates:
#'   `"robust"`, `"naive"`, `"nonrobust"`, `"HC0"`--`"HC3"` (where
#'   available) or `"bootstrap"` (re-runs the estimator). For `residuals()`:
#'   `"response"` (\eqn{y_t - x_t'\hat\beta}, full sample) or `"innovations"`
#'   (residuals of the final filtered and/or weighted regression).
#' @param parm Coefficients (names or positions); all by default.
#' @param level Confidence level.
#' @param newdata Optional data frame in which to look for the regressors.
#' @param se.fit Return standard errors of the predicted means.
#' @param interval `"none"` or `"confidence"` (for the mean of the response).
#' @param na.action Missing-value handling for `newdata`.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @return The usual values: a coefficient covariance matrix (`vcov`), a
#'   matrix of confidence limits (`confint`), vectors of residuals or fitted
#'   values, predictions (a vector, or a list/matrix with `se.fit`/`interval`),
#'   or the corresponding model components.
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
  if (is.null(type) || identical(type, object$vcov_type)) {
    return(object$vcov)
  }
  fgls_recompute_vcov(object, type)
}

fgls_recompute_vcov <- function(object, type) {
  spec <- object$spec
  est <- object$internal
  if (is.null(spec) || is.null(est)) {
    fgls_abort("This fit does not store what is needed to recompute the variance.")
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
                         level = 0.95, na.action = stats::na.pass, ...) {
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
  beta <- object$coefficients
  if (!identical(colnames(X), names(beta))) {
    fgls_abort("The regressors in `newdata` do not match the fitted model.")
  }
  pred <- drop(X %*% beta)
  if (!se.fit && interval == "none") {
    return(if (is.null(offset_na)) pred else stats::napredict(offset_na, pred))
  }
  se <- sqrt(pmax(rowSums((X %*% object$vcov) * X), 0))
  if (interval == "confidence") {
    crit <- critical_value(1 - level, fgls_df(object))
    fit <- cbind(fit = pred, lwr = pred - crit * se, upr = pred + crit * se)
    if (!se.fit) {
      return(fit)
    }
    return(list(fit = fit, se.fit = se, df = fgls_df(object)))
  }
  list(fit = pred, se.fit = se, df = fgls_df(object))
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
