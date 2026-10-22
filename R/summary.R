#' Summarize an fgls fit
#'
#' The summary contains the coefficient table (estimates, standard errors,
#' z or t statistics, p-values and confidence intervals) and the main
#' by-products of the estimator: the selected AR order \eqn{k^*} and AR
#' coefficients with the non-exogeneity share for serial correlation, and the
#' selected candidate covariates, the tuning parameters and the degrees of
#' freedom of the skedastic fit for heteroskedasticity.
#'
#' @param object An `"fgls"` object.
#' @param level Confidence level of the intervals.
#' @param x A `"summary.fgls"` object.
#' @param digits Number of significant digits.
#' @param signif.stars Print significance stars.
#' @param ... Not used.
#'
#' @return `summary()` returns an object of class `"summary.fgls"`: a list
#'   with the call, the coefficient table (`coefficients`), the confidence
#'   intervals (`conf.int`), a named list of header items (`info`) and the
#'   notes. Its `print()` method returns its argument invisibly.
#' @family fgls
#' @export
#' @examples
#' sim <- simulate_serial(200, seed = 1)
#' s <- summary(serial_fgls(y ~ x, data = sim))
#' s
#' s$coefficients
summary.fgls <- function(object, level = 0.95, ...) {
  cf <- object$coefficients
  se <- sqrt(pmax(diag(object$vcov), 0))
  stat <- ifelse(se > 0, cf / se, NA_real_)
  df <- fgls_df(object)
  pv <- two_sided_pvalues(stat, df)
  lab <- if (is.finite(df)) "t" else "z"
  tab <- cbind(cf, se, stat, pv)
  dimnames(tab) <- list(names(cf), c(
    "Estimate", "Std. Error", paste(lab, "value"),
    paste0("Pr(>|", lab, "|)")
  ))
  out <- list(
    call = object$call, method = object$method, vcov_type = object$vcov_type,
    coefficients = tab, conf.int = stats::confint(object, level = level),
    info = summary_info(object), notes = object$notes, df = df, level = level
  )
  class(out) <- "summary.fgls"
  out
}

summary_info <- function(x) {
  info <- list(
    "Method" = x$method,
    "Covariance" = x$vcov_type,
    "Observations" = if (x$n_eff != x$nobs) {
      sprintf("%d (effective %d)", x$nobs, x$n_eff)
    } else {
      as.character(x$nobs)
    },
    "Inference" = if (x$use_t) sprintf("Student-t (%d df)", x$df_resid) else "normal"
  )
  s <- x$serial
  if (!is.null(s)) {
    if (identical(s$method, "ma1")) {
      info[["MA(1) theta"]] <- sprintf("%.4f", s$extra$theta)
      info[["Durbin order k_T"]] <- as.character(s$extra$durbin_order)
    } else {
      how <- if (s$method %in% c("bic", "aic")) toupper(s$method) else s$method
      if (!is.na(s$kmax)) how <- sprintf("%s, kmax = %d", how, s$kmax)
      info[["AR order k*"]] <- sprintf("%d [%s]", s$ar_order, how)
      info[["AR coefficients"]] <- if (length(s$rho)) {
        paste(sprintf("%.4f", s$rho), collapse = ", ")
      } else {
        "-"
      }
      share <- s$nonexogeneity_share
      info[["Non-exogeneity share"]] <- if (is.na(share)) "-" else sprintf("%.3f", share)
      if (length(s$lagged_dep_vars)) {
        info[["Lagged dep. vars"]] <- paste(s$lagged_dep_vars, collapse = ", ")
      }
    }
  }
  h <- x$skedastic
  if (!is.null(h)) {
    info[["Skedastic fit"]] <- h$method
    if (h$method %in% c("lasso", "wls_s1", "wls_s2", "svr")) {
      info[["Selected z"]] <- if (length(h$selected)) {
        paste(h$selected, collapse = ", ")
      } else {
        "(none)"
      }
    }
    if (identical(h$method, "lasso") && !is.null(h$psi)) {
      info[["Lasso psi, lambda"]] <- sprintf("%g, %.4g", h$psi, h$lambda)
    }
    info[["df_hat"]] <- format(h$df_hat)
    v <- h$variance
    if (length(v)) info[["Variance ratio max/min"]] <- format(signif(max(v) / min(v), 3))
  }
  info
}

#' @rdname summary.fgls
#' @export
print.summary.fgls <- function(x, digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"), ...) {
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  }
  cat("\n")
  keys <- names(x$info)
  w <- max(nchar(keys)) + 2L
  for (k in keys) cat(formatC(paste0(k, ":"), width = -w), x$info[[k]], "\n", sep = "")
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients,
    digits = digits, signif.stars = signif.stars,
    na.print = "NA", ...
  )
  if (length(x$notes)) {
    cat("\nNotes:\n")
    for (i in seq_along(x$notes)) {
      cat(strwrap(sprintf("[%d] %s", i, x$notes[i]), exdent = 4L), sep = "\n")
    }
  }
  cat("\n")
  invisible(x)
}
