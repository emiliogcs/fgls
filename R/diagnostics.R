# Diagnostics (port of pyfgls/diagnostics.py): the variable-addition LM test of
# PG26 Remark 9 (experimental), augmented Dickey-Fuller screening, weight
# diagnostics and accessors for the AR roots and the non-exogeneity share.

#' Test whether lags of the regressors are missing (pre-determinedness)
#'
#' **Experimental.** A variable-addition LM test in the spirit of Remark 9 of
#' Perron and González-Coya (2026): FGLS needs the regressors to be
#' pre-determined (uncorrelated with current and future innovations), which can
#' fail when lags of the regressors belong in the model. The steps are (a) the
#' Durbin regression (AR order by BIC) giving \eqn{\hat\beta^D}; (b) the
#' residuals \eqn{\hat u^D = y - X\hat\beta^D}, demeaned; (c) for
#' \eqn{l = 1, \dots,} `max_lag`, the statistic
#' \eqn{LM_l = n (SSR_0 - SSR_1)/SSR_0}, asymptotically \eqn{\chi^2} with as
#' many degrees of freedom as added lags, from the auxiliary regressions without
#' and with lags \eqn{1..l} of the stochastic regressors.
#'
#' With `scale = "innovations"` (default) the auxiliary regression is run on
#' the innovation scale: \eqn{\hat u^D} and all regressors are filtered with
#' the Durbin AR filter and lagged \eqn{\hat u^D} enter as nuisance regressors,
#' which keeps the \eqn{\chi^2} approximation valid with pre-determined but
#' non-exogenous regressors. `scale = "errors"` uses \eqn{\hat u^D} and the raw
#' regressors (the literal reading of Remark 9).
#'
#' @inheritParams fgls
#' @param data A data frame, a list, a `ts`/`mts` object, or a 'zoo'/'xts'
#'   object containing the variables in the formula. Rows must be in time
#'   order.
#' @param ar_order Fixed AR order of the Durbin regression; `NULL` (default)
#'   selects it by BIC up to `kmax`.
#' @param max_lag Largest number of lags added.
#' @param regressors Names of the regressors whose lags are tested (default:
#'   all non-constant regressors).
#' @param alpha Level used for the advice.
#' @param scale `"innovations"` (default) or `"errors"`.
#' @param x A `"fgls_predetermined"` object.
#' @param ... Not used.
#' @return An object of class `"fgls_predetermined"`: a list with
#'   * `table`: a data frame with columns `lags`, `LM`, `df` and `p.value`,
#'     one row per number of added lags;
#'   * `ar_order`: the AR order of the Durbin regression;
#'   * `regressors`: the names of the regressors whose lags were tested;
#'   * `alpha`: the level used for the advice;
#'   * `rejected`: `TRUE` if some test rejects at level `alpha`;
#'   * `advice`: a sentence summarizing the result;
#'   * `scale`: the `scale` used, `"innovations"` or `"errors"`.
#'
#'   The print method returns `x` invisibly.
#' @inherit fgls references
#' @family diagnostics
#' @export
#' @examples
#' sim <- simulate_serial(300, seed = 1)
#' test_predetermined(y ~ x, data = sim)
test_predetermined <- function(formula, data, max_lag = 4L, ar_order = NULL, kmax = NULL,
                               regressors = NULL, alpha = 0.05,
                               scale = c("innovations", "errors"), na.action = stats::na.fail,
                               time = NULL, allow_gaps = FALSE) {
  fgls_warn(
    "test_predetermined() is experimental; its interface and details may change.",
    "experimental"
  )
  scale <- arg_default(scale, c("innovations", "errors"), "scale")
  if (!is_count(max_lag) || max_lag < 1L) fgls_abort("`max_lag` must be a positive integer.")
  bd <- bench_data(formula, if (missing(data)) NULL else data, NULL, parent.frame(),
    na.action, time, allow_gaps,
    time_series = TRUE
  )
  d <- bd$d
  n_all <- d$nobs
  X <- d$X
  if (is.null(ar_order)) {
    km <- resolve_kmax(kmax, n_all, d$k_vars, has_const = d$has_const)
    k <- select_ar_order(d$y, X, d$x_names, km)$order
  } else {
    if (!is_count(ar_order)) fgls_abort("`ar_order` must be a non-negative integer.")
    k <- as.integer(ar_order)
  }
  fit <- fit_durbin(durbin_design(d$y, X, d$x_names, k, y_name = d$y_name))
  rho <- if (k > 0L) unname(fit$coef[design_ylag_columns(fit$design, k)]) else numeric(0)
  uD <- drop(d$y - X %*% fit$beta)
  uD <- uD - mean(uD)
  if (is.null(regressors)) {
    idx <- which(apply(X, 2L, function(v) diff(range(v)) > 0))
  } else {
    idx <- match(regressors, d$x_names)
    if (anyNA(idx)) {
      fgls_abort(sprintf("Unknown regressors %s.", paste(regressors[is.na(idx)], collapse = ", ")))
    }
  }
  m <- length(idx)
  lagged <- lapply(seq_len(max_lag), function(j) {
    rbind(matrix(NA_real_, j, m), X[seq_len(n_all - j), idx, drop = FALSE])
  })
  if (scale == "innovations") {
    e <- ar_filter(uD, rho)
    base <- ar_filter(X, rho)
    if (k > 0L) {
      ulags <- vapply(seq_len(k), function(j) uD[(k - j + 1L):(n_all - j)], numeric(n_all - k))
      base <- cbind(base, ulags)
    }
    cands <- lapply(lagged, function(L) {
      L[is.na(L)] <- 0
      ar_filter(L, rho)
    })
    first <- k
  } else {
    e <- uD
    base <- X
    cands <- lagged
    first <- 0L
  }
  keep <- (max(max_lag + k, first) - first + 1L):length(e)
  u <- e[keep]
  B <- base[keep, , drop = FALSE]
  r0 <- u - drop(B %*% lsq(B, u))
  ssr0 <- sum(r0^2)
  n <- length(u)
  rows <- lapply(seq_len(max_lag), function(l) {
    Cm <- do.call(cbind, lapply(cands[seq_len(l)], function(L) L[keep, , drop = FALSE]))
    kept <- independent_columns(Cm, B, tol = 1e-8)$keep
    dof <- length(kept)
    if (dof == 0L) {
      return(c(l, 0, 0, 1))
    }
    A <- cbind(B, Cm[, kept, drop = FALSE])
    r1 <- u - drop(A %*% lsq(A, u))
    lm <- max(0, n * (ssr0 - sum(r1^2)) / ssr0)
    c(l, lm, dof, stats::pchisq(lm, dof, lower.tail = FALSE))
  })
  tab <- as.data.frame(do.call(rbind, rows))
  names(tab) <- c("lags", "LM", "df", "p.value")
  tab$lags <- as.integer(tab$lags)
  tab$df <- as.integer(tab$df)
  nm <- d$x_names[idx]
  rej <- tab$lags[tab$p.value < alpha]
  advice <- if (length(rej)) {
    sprintf(
      paste0(
        "Rejected at %g%% with %d lag(s): the regressors may not be pre-determined. Add ",
        "lag(s) 1..%d of %s as regressors, then refit."
      ),
      100 * alpha, rej[1L], rej[1L], paste(nm, collapse = ", ")
    )
  } else {
    sprintf("No rejection at %g%%: no evidence of missing regressor lags.", 100 * alpha)
  }
  structure(
    list(
      table = tab, ar_order = k, regressors = nm, alpha = alpha,
      rejected = length(rej) > 0L, advice = advice, scale = scale
    ),
    class = "fgls_predetermined"
  )
}

#' @rdname test_predetermined
#' @export
print.fgls_predetermined <- function(x, ...) {
  cat("Pre-determinedness LM test (PG26 Remark 9, experimental)\n")
  cat(sprintf(
    "Durbin AR order: %d; regressors tested: %s; scale: %s\n\n", x$ar_order,
    paste(x$regressors, collapse = ", "), x$scale
  ))
  print(format(x$table, digits = 4), row.names = FALSE)
  cat("\n", x$advice, "\n", sep = "")
  invisible(x)
}

# MacKinnon (2010) asymptotic critical values of the ADF test with a constant.
.adf_critical <- c("0.01" = -3.43, "0.05" = -2.86, "0.10" = -2.57)

#' Augmented Dickey-Fuller test (screening)
#'
#' The regression \eqn{\Delta x_t = a + g x_{t-1} + \sum_j c_j \Delta x_{t-j}
#' + e_t} with the number of lags chosen by BIC (at most \eqn{\lfloor 12
#' (T/100)^{1/4} \rfloor}); returns the t statistic on \eqn{g} with MacKinnon's
#' asymptotic critical values. Used by `fgls(check_stationarity = TRUE)` to
#' screen the response and the regressors; it is a screening device, not a
#' replacement for a careful unit-root analysis.
#'
#' @param x A numeric vector (a time series).
#' @param max_lags Largest number of lagged differences.
#' @return A list of class `"fgls_adf"` with `statistic`, `lags`, `nobs`,
#'   `critical_values` and `unit_root_rejected` (at 5%). The print method
#'   returns `x` invisibly.
#' @family diagnostics
#' @export
#' @examples
#' set.seed(1)
#' adf_test(cumsum(rnorm(200)))$unit_root_rejected
#' adf_test(rnorm(200))$unit_root_rejected
adf_test <- function(x, max_lags = NULL) {
  v <- as.numeric(x)
  n <- length(v)
  pmax_ <- if (is.null(max_lags)) floor(12 * (n / 100)^0.25) else as.integer(max_lags)
  pmax_ <- max(0L, min(pmax_, n %/% 4L))
  dx <- diff(v)
  design <- function(p, start) {
    cols <- list(rep(1, n - start), v[start:(n - 1L)])
    for (j in seq_len(p)) cols[[length(cols) + 1L]] <- dx[(start - j):(n - 1L - j)]
    list(A = do.call(cbind, cols), y = dx[start:(n - 1L)])
  }
  start <- pmax_ + 1L
  bics <- vapply(0:pmax_, function(p) {
    ds <- design(p, start)
    e <- ds$y - drop(ds$A %*% lsq(ds$A, ds$y))
    m <- length(ds$y)
    log(sum(e^2) / m) + ncol(ds$A) * log(m) / m
  }, numeric(1))
  p <- which.min(bics) - 1L
  ds <- design(p, p + 1L)
  b <- lsq(ds$A, ds$y)
  e <- ds$y - drop(ds$A %*% b)
  m <- nrow(ds$A)
  s2 <- sum(e^2) / (m - ncol(ds$A))
  se <- sqrt(s2 * solve(crossprod(ds$A))[2L, 2L])
  stat <- b[2L] / se
  structure(
    list(
      statistic = stat, lags = p, nobs = m, critical_values = .adf_critical,
      unit_root_rejected = stat < .adf_critical[["0.05"]]
    ),
    class = "fgls_adf"
  )
}

#' @rdname adf_test
#' @param ... Not used.
#' @export
print.fgls_adf <- function(x, ...) {
  cat(sprintf(
    paste(
      "ADF test (constant, %d lag(s), %d obs.): statistic %.3f; 5%% critical value",
      "%.2f; unit root %s\n"
    ),
    x$lags, x$nobs, x$statistic, x$critical_values[["0.05"]],
    if (x$unit_root_rejected) "rejected" else "not rejected"
  ))
  invisible(x)
}

# ADF screening of y and every non-constant, non-binary regressor (warnings).
screen_stationarity <- function(d) {
  series <- c(list(d$y), lapply(seq_len(d$k_vars), function(j) d$X[, j]))
  names(series) <- c(d$y_name, d$x_names)
  flagged <- character(0)
  out <- list()
  for (nm in names(series)) {
    s <- series[[nm]]
    if (length(unique(s)) <= 2L) next
    r <- adf_test(s)
    out[[nm]] <- r
    if (!r$unit_root_rejected) flagged <- c(flagged, sprintf("%s (ADF %.2f)", nm, r$statistic))
  }
  if (length(flagged)) {
    fgls_warn(
      paste0(
        "ADF tests do not reject a unit root at 5% for: ", paste(flagged, collapse = ", "),
        ". FGLS assumes stationary data; consider differencing (or cointegration methods)."
      ),
      "nonstationary"
    )
  }
  invisible(out)
}

#' Diagnostics of an fgls fit
#'
#' * `weight_diagnostics()`: how concentrated the FGLS weights \eqn{1/\hat v_t}
#'   are: the ratio of the largest to the smallest fitted variance, the Kish
#'   effective sample size \eqn{(\sum w)^2 / \sum w^2} as a share of the
#'   sample, and the share of the total weight carried by the top 1% of
#'   observations.
#' * `ar_roots()`: roots of the estimated lag polynomial
#'   \eqn{1 - \hat\rho_1 z - \dots - \hat\rho_k z^k}; the filter is invertible
#'   when all lie outside the unit circle.
#' * `nonexogeneity_share()`: the relative size of the first-stage correction
#'   of PG26 eq. (29), averaged over coefficients. Values near zero mean the
#'   regressors behave as exogenous; large values mean the naive variance (eq.
#'   27) would understate the uncertainty.
#'
#' @param fit An `"fgls"` object (with a heteroskedasticity correction for
#'   `weight_diagnostics()`, a serial-correlation correction otherwise).
#' @return
#' * `weight_diagnostics()`: a named numeric vector with elements
#'   `variance_ratio`, `effective_n_share` and `top1pct_weight_share`.
#' * `ar_roots()`: a complex vector with the \eqn{k^*} roots (empty when
#'   \eqn{k^* = 0}).
#' * `nonexogeneity_share()`: a single number.
#' @family diagnostics
#' @export
#' @examples
#' fit <- fgls(y ~ x, data = simulate_combined(200, seed = 1), seed = 1)
#' weight_diagnostics(fit)
#' Mod(ar_roots(fit))
#' nonexogeneity_share(fit)
weight_diagnostics <- function(fit) {
  if (!inherits(fit, "fgls") || is.null(fit$skedastic)) {
    fgls_abort("This fit has no heteroskedasticity weights.")
  }
  v <- fit$skedastic$variance
  w <- 1 / v
  ws <- sort(w, decreasing = TRUE)
  top <- max(1L, ceiling(0.01 * length(w)))
  c(
    variance_ratio = max(v) / min(v),
    effective_n_share = sum(w)^2 / sum(w^2) / length(w),
    top1pct_weight_share = sum(ws[seq_len(top)]) / sum(w)
  )
}

#' @rdname weight_diagnostics
#' @export
ar_roots <- function(fit) {
  if (!inherits(fit, "fgls") || is.null(fit$serial)) {
    fgls_abort("This fit has no serial-correlation correction.")
  }
  ar_poly_roots(fit$serial$rho)
}

#' @rdname weight_diagnostics
#' @export
nonexogeneity_share <- function(fit) {
  if (!inherits(fit, "fgls") || is.null(fit$serial)) {
    fgls_abort("This fit has no serial-correlation correction.")
  }
  fit$serial$nonexogeneity_share
}
