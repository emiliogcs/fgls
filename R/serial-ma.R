# GMA: feasible GLS for MA(1) errors, PG26 Section 5 (port of
# pyfgls/serial/ma.py).
#
# For u_t = e_t + theta e_{t-1} the GLS transformation inverts the MA filter,
# y*_t = y_t - theta y*_{t-1} (and likewise every column of X), and
# y*_t = x*_t' beta + e_t is estimated by OLS with Var = s^2 (X*'X*)^{-1}.
# theta is the exact Gaussian maximum-likelihood estimate (innovations
# algorithm, sigma^2 profiled out) from u_t = y_t - x_t' beta_tilde, where
# beta_tilde is the Durbin-regression estimate with k_T = int(4 (T/100)^(2/9))
# lags, consistent whether or not the regressors are exogenous.

.theta_bound <- 0.995

gma_durbin_order <- function(n) {
  as.integer(floor(4 * (n / 100)^(2 / 9)))
}

# Profile negative Gaussian log-likelihood of an MA(1) (innovations algorithm).
ma1_negloglik <- function(theta, u) {
  n <- length(u)
  g0 <- 1 + theta^2
  g1 <- theta
  r <- numeric(n)
  innov <- numeric(n)
  r[1L] <- g0
  innov[1L] <- u[1L]
  for (t in seq_len(n)[-1L]) {
    th <- g1 / r[t - 1L]
    r[t] <- g0 - th^2 * r[t - 1L]
    innov[t] <- u[t] - th * innov[t - 1L]
  }
  s2 <- mean(innov^2 / r)
  0.5 * n * log(s2) + 0.5 * sum(log(r))
}

# Exact Gaussian MLE of theta in u_t = e_t + theta e_{t-1} (Brent's method on
# [-0.995, 0.995], as scipy's bounded scalar minimiser).
ma1_mle <- function(u) {
  uc <- u - mean(u)
  theta <- stats::optimize(ma1_negloglik, c(-.theta_bound, .theta_bound), u = uc, tol = 1e-8)$minimum
  if (abs(theta) > 0.98) {
    fgls_warn(
      sprintf(
        "Estimated MA coefficient %.3f is close to the non-invertibility boundary.",
        theta
      ),
      "nonstationary"
    )
  }
  theta
}

# a*_t = a_t - theta a*_{t-1} with a*_0 = 0, column by column.
inverse_ma_filter <- function(a, theta) {
  f <- function(v) as.numeric(stats::filter(v, -theta, method = "recursive"))
  if (is.matrix(a)) {
    out <- apply(a, 2L, f)
    dim(out) <- dim(a)
    dimnames(out) <- dimnames(a)
    out
  } else {
    f(a)
  }
}

estimate_gma <- function(y, X, x_names = colnames(X), durbin_order = NULL) {
  n <- nrow(X)
  q <- ncol(X)
  kT <- if (is.null(durbin_order)) gma_durbin_order(n) else as.integer(durbin_order)
  beta_tilde <- fit_durbin(durbin_design(y, X, x_names, kT))$beta
  theta <- ma1_mle(as.numeric(y - X %*% beta_tilde))
  ys <- inverse_ma_filter(y, theta)
  Xs <- inverse_ma_filter(X, theta)
  beta <- lsq(Xs, ys)
  names(beta) <- x_names
  e <- as.numeric(ys - Xs %*% beta)
  s2 <- sum(e^2) / (n - q)
  cov <- s2 * inv_psd(crossprod(Xs), name = "GMA transformed X'X")
  dimnames(cov) <- list(x_names, x_names)
  list(
    beta = beta, cov = cov, theta = theta, durbin_order = kT, ystar = ys, Xstar = Xs,
    ehat = e, uhat = as.numeric(y - X %*% beta), n_eff = n
  )
}
