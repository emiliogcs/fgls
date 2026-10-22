#' Tidy, glance and augment fgls fits
#'
#' Methods for the 'broom' generics (from the 'generics' package), registered
#' when 'generics' is loaded.
#'
#' @param x An `"fgls"` object.
#' @param conf.int Include confidence intervals.
#' @param conf.level Confidence level.
#' @param data Data to augment: by default the model frame. It must have one
#'   row per observation used in the fit or, for a fit with
#'   `na.action = na.exclude`, one row per row of the original data (the
#'   excluded rows then get `NA`).
#' @param newdata Optional new data; fitted values (and, if requested,
#'   standard errors and confidence limits) are added, but no residuals.
#' @param se_fit Add the standard errors of the fitted means (`.se.fit`).
#' @param interval `"none"` or `"confidence"`: add confidence limits for the
#'   mean of the response (`.lower`, `.upper`) at level `conf.level`.
#' @param ... Not used.
#' @return A data frame (a tibble if 'tibble' is installed): one row per
#'   coefficient for `tidy()`, one row for `glance()`, one row per observation
#'   for `augment()` (columns `.fitted`, optionally `.lower`, `.upper` and
#'   `.se.fit`, then `.resid`, `.innovation` and, with a heteroskedasticity
#'   correction, `.weight`). `nonexogeneity.share` in `glance()` is `NA` when
#'   the share was not computed.
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
    nonexogeneity.share = if (is.null(s$nonexogeneity_share)) {
      NA_real_
    } else {
      as.numeric(s$nonexogeneity_share)
    },
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
augment.fgls <- function(x, data = NULL, newdata = NULL, se_fit = FALSE,
                         interval = c("none", "confidence"), conf.level = 0.95, ...) {
  interval <- arg_default(interval, c("none", "confidence"), "interval")
  if (!is.null(newdata)) {
    out <- as.data.frame(newdata)
    if (!isTRUE(se_fit) && interval == "none") {
      out$.fitted <- unname(stats::predict(x, newdata = newdata))
      return(as_tidy_frame(out))
    }
    p <- stats::predict(x,
      newdata = newdata, se.fit = TRUE, interval = "confidence", level = conf.level
    )
    p <- list(fit = p$fit[, "fit"], lwr = p$fit[, "lwr"], upr = p$fit[, "upr"], se = p$se.fit)
    return(as_tidy_frame(augment_prediction_columns(out, p, se_fit, interval, identity)))
  }
  if (is.null(data)) {
    if (is.null(x$model)) {
      fgls_abort(paste0(
        "The model frame was not kept (model = FALSE): pass the data of the fit as ",
        "`data`, or refit with model = TRUE."
      ))
    }
    out <- as.data.frame(x$model)
  } else {
    out <- as.data.frame(data)
  }
  rn <- names(x$residuals)
  n_fit <- length(rn)
  na <- x$na.action
  pad <- identity
  if (nrow(out) != n_fit) {
    if (inherits(na, "exclude") && nrow(out) == n_fit + length(na)) {
      # one row per row of the original data: pad as residuals() and fitted() do
      pad <- function(v) stats::naresid(na, v)
    } else {
      fgls_abort(sprintf(
        "`data` must have one row per observation used in the fit (%d)%s; it has %d.",
        n_fit,
        if (inherits(na, "exclude")) {
          sprintf(" or one row per row of the original data (%d)", n_fit + length(na))
        } else {
          ""
        },
        nrow(out)
      ))
    }
  }
  p <- list(fit = unname(x$fitted.values))
  if (isTRUE(se_fit) || interval != "none") {
    p <- predict_mean(x, x$x, level = conf.level)
    p$fit <- unname(x$fitted.values)
  }
  out <- augment_prediction_columns(out, p, se_fit, interval, pad)
  out$.resid <- pad(unname(x$residuals))
  out$.innovation <- pad(unname(x$innovations[rn]))
  if (!is.null(x$weights)) {
    w <- rep(NA_real_, n_fit)
    w[match(names(x$innovations), rn)] <- x$weights
    out$.weight <- pad(w)
  }
  as_tidy_frame(out)
}

# Add .fitted and, as requested, .lower/.upper and .se.fit (the column order of
# broom's augment.lm()); `pad` aligns the values with the rows of `out`.
augment_prediction_columns <- function(out, p, se_fit, interval, pad) {
  out$.fitted <- pad(unname(p$fit))
  if (interval == "confidence") {
    out$.lower <- pad(unname(p$lwr))
    out$.upper <- pad(unname(p$upr))
  }
  if (isTRUE(se_fit)) out$.se.fit <- pad(unname(p$se))
  out
}
