#' Tidy, glance and augment fgls fits
#'
#' Methods for the 'broom' generics (from the 'generics' package), registered
#' when 'generics' is loaded.
#'
#' @param x An `"fgls"` object.
#' @param conf.int Include confidence intervals.
#' @param conf.level Confidence level.
#' @param data Data to augment: by default the model frame.
#' @param newdata Optional new data; only fitted values are added.
#' @param ... Not used.
#' @return A data frame (a tibble if 'tibble' is installed): one row per
#'   coefficient for `tidy()`, one row for `glance()`, one row per observation
#'   for `augment()` (columns `.fitted`, `.resid`, `.innovation` and, with a
#'   heteroskedasticity correction, `.weight`).
#' @name fgls-broom
#' @examples
#' if (requireNamespace("broom", quietly = TRUE)) {
#'   fit <- serial_fgls(y ~ x, data = simulate_serial(200, seed = 1))
#'   broom::tidy(fit, conf.int = TRUE)
#'   broom::glance(fit)
#'   head(broom::augment(fit))
#' }
NULL

as_tidy_frame <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' @rdname fgls-broom
#' @exportS3Method generics::tidy
tidy.fgls <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  s <- summary(x, level = conf.level)$coefficients
  out <- data.frame(
    term = rownames(s), estimate = s[, 1L], std.error = s[, 2L], statistic = s[, 3L],
    p.value = s[, 4L], row.names = NULL, stringsAsFactors = FALSE
  )
  if (isTRUE(conf.int)) {
    ci <- stats::confint(x, level = conf.level)
    out$conf.low <- unname(ci[, 1L])
    out$conf.high <- unname(ci[, 2L])
  }
  as_tidy_frame(out)
}

#' @rdname fgls-broom
#' @exportS3Method generics::glance
glance.fgls <- function(x, ...) {
  s <- x$serial
  h <- x$skedastic
  out <- data.frame(
    method = x$method, vcov.type = x$vcov_type, nobs = x$nobs, n.eff = x$n_eff,
    ar.order = if (is.null(s)) NA_integer_ else s$ar_order,
    nonexogeneity.share = if (is.null(s)) NA_real_ else s$nonexogeneity_share,
    skedastic = if (is.null(h)) NA_character_ else h$method,
    n.selected = if (is.null(h)) NA_integer_ else length(h$selected),
    df.hat = if (is.null(h)) NA_real_ else h$df_hat,
    sigma = sqrt(mean(x$innovations^2)), df.residual = fgls_df(x),
    stringsAsFactors = FALSE
  )
  as_tidy_frame(out)
}

#' @rdname fgls-broom
#' @exportS3Method generics::augment
augment.fgls <- function(x, data = NULL, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    out <- as.data.frame(newdata)
    out$.fitted <- unname(stats::predict(x, newdata = newdata))
    return(as_tidy_frame(out))
  }
  out <- if (is.null(data)) as.data.frame(x$model) else as.data.frame(data)
  rn <- names(x$residuals)
  if (nrow(out) != length(rn)) {
    fgls_abort("`data` must have one row per observation used in the fit.")
  }
  out$.fitted <- unname(x$fitted.values)
  out$.resid <- unname(x$residuals)
  out$.innovation <- unname(x$innovations[rn])
  if (!is.null(x$weights)) {
    w <- rep(NA_real_, length(rn))
    w[match(names(x$innovations), rn)] <- x$weights
    out$.weight <- w
  }
  as_tidy_frame(out)
}
