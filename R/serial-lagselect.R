# Selection of the autoregressive order by an information criterion, PG26
# eq. (14) (port of pyfgls/serial/lagselect.py).
#
# Every candidate order k = 0, ..., kmax is estimated on the *common sample*
# t = kmax + 1, ..., T so the criteria are comparable (Ng and Perron, 2005):
#
#   BIC(k) = ln(sigma2_k) + k ln(T - kmax) / (T - kmax),
#   AIC(k) = ln(sigma2_k) + 2 k / (T - kmax),
#
# with sigma2_k the mean squared residual of the order-k Durbin regression. The
# penalty counts the AR order k, as in the paper (ic_penalty = "order");
# ic_penalty = "params" penalizes every coefficient a lag adds instead.
#
# The order-k design is a prefix of the order-kmax design, so one QR
# factorization Z = QR gives every residual sum of squares: with r = Q'y,
# RSS_k = RSS_kmax + sum_{i > p_k} r_i^2.

# Upper bound on the number of Durbin regressors relative to the sample.
.max_regressor_share <- 0.5

# Conservative count of Durbin regressors for order k: every lag of y and of
# every regressor is kept except the lagged constant (the formula of pyfgls).
n_durbin_regressors <- function(q, k, has_const) {
  q + k * (1 + q - as.integer(has_const))
}

# Number of regressors of the order-k Durbin regression on its own sample
# t = k + 1, ..., T once the collinear lagged columns are dropped (lags of the
# constant, of trends and of periodic dummies; lags of y that are regressors):
# the number of columns durbin_design() keeps. It never exceeds the
# conservative count, and equals it when only lagged constants are dropped.
# Returns the conservative count when the sample is too short to build the
# design (the order is then infeasible anyway).
durbin_regressor_count <- function(y, X, k, has_const = TRUE) {
  q <- ncol(X)
  if (k == 0L) {
    return(q)
  }
  conservative <- n_durbin_regressors(q, k, has_const)
  if (nrow(X) - k <= q) {
    return(conservative)
  }
  des <- tryCatch(
    durbin_design(y, X, sprintf("x%d", seq_len(q)), k),
    error = function(e) NULL
  )
  if (is.null(des)) conservative else ncol(des$Z)
}

is_count <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 0 && x == round(x)
}

# Turn the kmax option into an order the sample can support. The default is
# int(12 (T/200)^(1/4)); "theory" gives floor((T / ln T)^(1/3)). The order is
# capped at the largest k whose Durbin regression has at most half as many
# regressors as observations (T - k), with a warning when the cap binds.
#
# The conservative count n_durbin_regressors() gives a first cap, exactly as
# pyfgls. When y and X are given, the orders above it are then checked with
# the number of regressors the Durbin design actually keeps (lags of trends and
# of periodic dummies are collinear with the contemporaneous block and dropped),
# so that, e.g., a monthly model with seasonal dummies is not capped by columns
# the regression never contains. Without such columns the two counts agree and
# the cap is the pyfgls one.
resolve_kmax <- function(kmax, n, q, has_const, extra_lags = 0L, y = NULL, X = NULL) {
  if (is.null(kmax) || identical(kmax, "default")) {
    k <- default_kmax(n)
  } else if (identical(kmax, "theory")) {
    k <- theory_kmax(n)
  } else if (is_count(kmax)) {
    k <- as.integer(kmax)
  } else {
    fgls_abort(sprintf(
      "`kmax` must be a non-negative integer, \"default\" or \"theory\"; got %s.",
      paste(format(kmax), collapse = ", ")
    ))
  }
  limit <- function(order) .max_regressor_share * (n - order)
  if (n_durbin_regressors(q, extra_lags, has_const) > limit(extra_lags)) {
    fgls_abort(
      sprintf(
        paste0(
          "T = %d observations are too few for %d regressors: the Durbin regression ",
          "may have at most half as many regressors as observations."
        ),
        n, q
      ),
      "insufficient_data"
    )
  }
  cap <- k
  while (cap > 0L &&
    n_durbin_regressors(q, cap + extra_lags, has_const) > limit(cap + extra_lags)) {
    cap <- cap - 1L
  }
  rejected <- n_durbin_regressors(q, cap + 1L + extra_lags, has_const)
  if (!is.null(X)) {
    # raise the cap while the regressors actually kept fit
    while (cap < k) {
      order <- cap + 1L + extra_lags
      rejected <- durbin_regressor_count(y, X, order, has_const)
      if (rejected > limit(order)) break
      cap <- cap + 1L
    }
  }
  if (cap < k) {
    order <- cap + 1L + extra_lags
    fgls_warn(
      sprintf(
        paste0(
          "kmax was reduced from %d to %d: with %d lags the Durbin regression would have ",
          "%d regressors (collinear lagged columns excluded) for %d observations, more ",
          "than half of them."
        ),
        k, cap, order, as.integer(rejected), n - order
      ),
      "sample_size"
    )
  }
  as.integer(cap)
}

# Abort when a fixed AR order leaves the Durbin regression of PG26 eq. (14)
# without residual degrees of freedom (regressors counted after dropping the
# collinear lagged columns).
check_ar_order <- function(k, y, X, has_const) {
  n <- nrow(X)
  m <- n_durbin_regressors(ncol(X), k, has_const)
  if (m >= n - k) m <- durbin_regressor_count(y, X, k, has_const)
  if (m >= n - k) {
    fgls_abort(
      sprintf(
        paste0(
          "ar_order = %d is too large for T = %d observations: the Durbin regression ",
          "would have %d regressors (collinear lagged columns excluded) for %d ",
          "observations."
        ),
        k, n, as.integer(m), max(n - k, 0L)
      ),
      "insufficient_data"
    )
  }
  invisible(m)
}

resolve_min_lags <- function(min_lags, n) {
  if (is.null(min_lags)) {
    return(NULL)
  }
  if (identical(min_lags, "adjusted")) {
    return(adjusted_min_lags(n))
  }
  if (is_count(min_lags)) {
    return(as.integer(min_lags))
  }
  fgls_abort("`min_lags` must be NULL, \"adjusted\" or a non-negative integer.")
}

# Residual sums of squares of the nested Durbin regressions of orders 0..kmax.
nested_rss <- function(design, kmax) {
  Z <- design$Z
  yv <- design$y
  p <- ncol(Z)
  Q <- design$Q
  if (!is.null(Q) && ncol(Q) == p) {
    # orthonormal basis of the design, columns in design order (Gram-Schmidt)
    r <- drop(crossprod(Q, yv))
    res <- yv - drop(Q %*% r)
    res <- res - drop(Q %*% crossprod(Q, res))
    rss_full <- sum(res^2)
    tail_ss <- rev(cumsum(rev(r^2)))
    return(vapply(0:kmax, function(k) {
      pk <- design$block_end[k + 1L]
      rss_full + if (pk < p) tail_ss[pk + 1L] else 0
    }, numeric(1)))
  }
  qz <- qr(Z, tol = 1e-11)
  if (qz$rank == p && identical(qz$pivot, seq_len(p))) {
    r <- qr.qty(qz, yv)[seq_len(p)]
    rss_full <- sum(qr.resid(qz, yv)^2)
    tail_ss <- rev(cumsum(rev(r^2)))
    vapply(0:kmax, function(k) {
      pk <- design$block_end[k + 1L]
      rss_full + if (pk < p) tail_ss[pk + 1L] else 0
    }, numeric(1))
  } else { # pragma: fallback when QR pivoted (near-collinear columns)
    vapply(0:kmax, function(k) {
      sum(qr.resid(qr(design_prefix(design, k), tol = 1e-11), yv)^2)
    }, numeric(1))
  }
}

# Choose the AR order minimizing BIC/AIC over k = 0..kmax on the common sample.
select_ar_order <- function(y, X, x_names, kmax, ic = "bic", min_lags = NULL,
                            y_name = "y", ic_penalty = "order") {
  arg_choice(ic, c("bic", "aic"), "ic")
  arg_choice(ic_penalty, c("order", "params"), "ic_penalty")
  design <- durbin_design(y, X, x_names, kmax, start = kmax, y_name = y_name)
  n <- length(design$y)
  rss <- nested_rss(design, kmax)
  sigma2 <- rss / n
  npar <- if (ic_penalty == "order") {
    as.numeric(0:kmax)
  } else {
    as.numeric(design$block_end - design$block_end[1L])
  }
  penalty <- if (ic == "bic") log(n) / n else 2 / n
  values <- suppressWarnings(log(sigma2)) + npar * penalty
  k_ic <- which.min(values) - 1L
  k <- k_ic
  if (!is.null(min_lags)) k <- max(k, min_lags)
  list(
    order = as.integer(k), selected_by_ic = as.integer(k_ic), kmax = as.integer(kmax),
    criterion = ic, values = stats::setNames(values, 0:kmax),
    sigma2 = sigma2, n_common = n, dropped = design$dropped
  )
}
