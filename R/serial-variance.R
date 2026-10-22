# Variance estimators of the serial-correlation FGLS estimator, PG26 Theorem 2
# (port of pyfgls/serial/variance.py).
#
# Notation (effective sample t = k+1..T, n = T - k): xf_t the filtered
# regressors, u_t the *unfiltered* FGLS residuals, e_t the residuals of the
# filtered regression, sigma2 = n^{-1} sum e_t^2 and Q = n^{-1} sum xf_t xf_t'.
#
# Naive variance, PG26 eq. (27), valid only for exogenous regressors:
#
#   AVar_N = sigma2 Q^{-1}
#
# Exogeneity-robust variance, PG26 eqs. (28)-(29), the default:
#
#   g_j    = n^{-1} sum_t xf_t u_{t-j},   G = [g_1, ..., g_k]            (q x k)
#   psi_t  = residual of (y_{t-1}, ..., y_{t-k})' on the other Durbin
#            regressors W_t = (x_t', x_{t-1}', ..., x_{t-k}')'           (k x 1)
#   Sigma  = n^{-1} sum_t psi_t psi_t'
#   AVar_R = sigma2 (Q^{-1} + Q^{-1} G Sigma^{-1} G' Q^{-1})
#
# and Var(beta) = AVar / n. The correction is positive semidefinite.
#
# The same correction written through the influence function of the estimator,
# a_t = xf_t - G C_t with C = J (Z'Z/n)^{-1} Z' (Z the Durbin design, J the
# Jacobian of rho with respect to the Durbin coefficients), gives
# n^{-1} sum a_t a_t' = Q + G Sigma^{-1} G' when rho are the coefficients on
# the lagged y. This general form serves the Wallis recovery (lagged dependent
# variables, PG26 Section 4.2) and the HC0-HC3 variants, which replace sigma2
# by e_t^2 inside the sum.

serial_cov_types <- c("robust", "naive", "HC0", "HC1", "HC2", "HC3")

# G = [g_1..g_k], g_j = n^{-1} sum_t xf_t u_{t-j} (PG26 eq. (28)). Xf holds
# observations start+1..T; uhat the unfiltered residuals on the full sample.
# The raw lagged residual u_{t-j} is used.
g_matrix <- function(Xf, uhat, start, k) {
  n <- nrow(Xf)
  n_all <- length(uhat)
  G <- matrix(0, ncol(Xf), k)
  for (j in seq_len(k)) {
    G[, j] <- crossprod(Xf, uhat[(start - j + 1L):(n_all - j)]) / n
  }
  G
}

# psi_t of PG26 eq. (29): the part of (y_{t-1}, ..., y_{t-k}) not explained by
# W_t = (x_t', x_{t-1}', ..., x_{t-k}')' (collinear lagged columns removed).
psi_residuals <- function(design, k) {
  Z <- design_prefix(design, k)
  ycols <- design_ylag_columns(design, k)
  other <- setdiff(seq_len(ncol(Z)), ycols)
  Yl <- Z[, ycols, drop = FALSE]
  W <- Z[, other, drop = FALSE]
  Yl - W %*% lsq(W, Yl)
}

# C = J (Z'Z/n)^{-1} Z' (k x n): first-stage influence of the Durbin regression.
influence_weights <- function(Z, J) {
  n <- nrow(Z)
  qz <- qr(Z, tol = 1e-11)
  R <- qr.R(qz)
  Q <- qr.Q(qz)
  Jp <- J[, qz$pivot, drop = FALSE]
  JRinv <- t(forwardsolve(t(R), t(Jp))) # J R^{-1}
  n * JRinv %*% t(Q)
}

# First-stage influence weights C (k x n) of the AR coefficients:
# C_t = Sigma^{-1} psi_t in the standard case (J = NULL) and
# C = J (Z'Z/n)^{-1} Z' in general.
first_stage_influence <- function(design, J = NULL) {
  k <- design$order
  n <- length(design$y)
  if (k == 0L) {
    return(matrix(0, 0L, n))
  }
  if (is.null(J)) {
    Psi <- psi_residuals(design, k)
    Sigma <- crossprod(Psi) / n
    return(inv_psd(Sigma, name = "variance of psi_t (eq. 29)") %*% t(Psi))
  }
  influence_weights(design_prefix(design, k), J)
}

hc_weights <- function(cov_type, ehat, Xf) {
  n <- nrow(Xf)
  q <- ncol(Xf)
  e2 <- ehat^2
  if (cov_type == "HC0") {
    return(e2)
  }
  if (cov_type == "HC1") {
    return(e2 * n / max(n - q, 1L))
  }
  h <- pmin(pmax(hat_diag(Xf), 0), 1 - 1e-12)
  if (cov_type == "HC2") e2 / (1 - h) else e2 / (1 - h)^2
}

# Variance of the FGLS estimator with an estimated AR(k) filter.
#
# Xf: filtered regressors (n x q); ehat: filtered-regression residuals; uhat:
# unfiltered FGLS residuals (T); design: Durbin design of order k with
# design$start == k; J: Jacobian of rho with respect to the Durbin
# coefficients, NULL when rho are the coefficients on the lagged y (eq. 29).
serial_variance <- function(Xf, ehat, uhat, design, J = NULL, cov_type = "robust") {
  n <- nrow(Xf)
  q <- ncol(Xf)
  k <- design$order
  start <- design$start
  sigma2 <- sum(ehat^2) / n
  Qhat <- crossprod(Xf) / n
  Qinv <- inv_psd(Qhat, name = "filtered regressor moment matrix")
  naive <- sigma2 * Qinv
  G <- g_matrix(Xf, uhat, start, k)

  if (k == 0L) {
    C <- matrix(0, 0L, n)
    Sigma <- matrix(0, 0L, 0L)
  } else if (is.null(J)) {
    Psi <- psi_residuals(design, k)
    Sigma <- crossprod(Psi) / n
    C <- inv_psd(Sigma, name = "variance of psi_t (eq. 29)") %*% t(Psi)
  } else {
    Sigma <- NULL
    C <- influence_weights(design_prefix(design, k), J)
  }

  if (cov_type == "naive") {
    avar <- naive
  } else if (cov_type == "robust" && is.null(J)) {
    if (k == 0L) {
      avar <- naive
    } else {
      # PG26 eq. (29), literally
      corr <- Qinv %*% G %*% inv_psd(Sigma, name = "variance of psi_t (eq. 29)") %*%
        t(G) %*% Qinv
      avar <- sigma2 * (Qinv + make_psd(corr, name = "eq. (29) correction"))
    }
  } else {
    A <- Xf - crossprod(C, t(G)) # influence a_t, rows t
    if (cov_type == "robust") {
      meat <- sigma2 * crossprod(A) / n
    } else if (cov_type %in% c("HC0", "HC1", "HC2", "HC3")) {
      w <- hc_weights(cov_type, ehat, Xf)
      meat <- crossprod(A * w, A) / n
    } else {
      stop("unknown cov_type ", cov_type)
    }
    avar <- Qinv %*% meat %*% Qinv
  }
  avar <- make_psd(avar, name = "FGLS covariance matrix")
  nm <- colnames(Xf)
  dimnames(avar) <- dimnames(naive) <- list(nm, nm)
  correction <- if (cov_type == "robust") avar - naive else matrix(0, q, q, dimnames = list(nm, nm))
  d <- diag(naive)
  by_coef <- ifelse(d > 0, diag(correction) / d, 0)
  names(by_coef) <- nm
  list(
    cov = avar / n, avar = avar, naive_avar = naive, correction = correction,
    G = G, Sigma = Sigma, Q = Qhat, sigma2 = sigma2, n = n, cov_type = cov_type,
    nonexogeneity_by_coef = by_coef,
    nonexogeneity_share = if (length(by_coef)) mean(by_coef) else 0
  )
}
