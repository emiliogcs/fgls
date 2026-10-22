# Weights and sandwich variances for FGLS with an estimated skedastic function,
# GP24 Section 4 (port of pyfgls/hetero/weights.py and variance.py).
#
# Because the skedastic model is only an approximation, inference uses the
# sandwich (GP24 eq. (6))
#
#   Var(b) = (T^{-1} X'W^{-1}X)^{-1} Omega (T^{-1} X'W^{-1}X)^{-1} / T,
#   Omega  = T^{-1} sum_t Sigma_tt x_t x_t',
#
# with W = diag(v_t) and, for the Lasso and SVR fits, the Miller-Startz type
# bias correction of eq. (8):
#
#   Sigma_tt = (u_t^2 / v_t^2) (1/(1 - h_t^FGLS)^2 + 4 (h_t^OLS / k) df),
#
# u_t = y_t - x_t' b_FGLS, h^FGLS the leverages of the weighted regression,
# h^OLS those of OLS, k the number of regressors and df the degrees of freedom
# of the skedastic fit (non-zero Lasso coefficients including the intercept,
# or the SVR count). For the Romano-Wolf fits (eq. (7)),
# Sigma_tt = e_t^2 / v_t^2 with the HC3-type OLS residual
# e_t^2 = u_t^{OLS 2} / (1 - h_t^OLS)^2 (open item B of pyfgls).
#
# The same code serves the combined estimator (PG26 Supplement (A.4)-(A.5)),
# where X holds the filtered regressors and "OLS" is the unweighted filtered
# (F-C) regression.

hetero_cov_types <- c("robust", "nonrobust", "HC0", "HC1", "HC2", "HC3")

# Warn when the largest fitted variance exceeds the smallest by this factor.
.max_variance_ratio <- 1e6
.log_clip <- 700

# Response of the skedastic regression, log(max(u_t^2, delta^2)).
log_squared_residuals <- function(u, delta) {
  if (!(delta > 0)) fgls_abort(sprintf("delta must be positive, got %s.", format(delta)))
  log(pmax(u^2, delta^2))
}

# Fitted log-variances to variances normalised to mean one (optional quantile
# winsorisation). Normalising by a constant changes neither the estimates nor
# the sandwich variance.
variances_from_log <- function(h, clip = NULL) {
  if (!is.null(clip)) {
    q <- stats::quantile(h, clip, names = FALSE, type = 7)
    h <- pmin(pmax(h, q[1L]), q[2L])
  }
  h <- pmin(pmax(h - max(h), -.log_clip), 0)
  v <- exp(h)
  v <- v / mean(v)
  ratio <- max(v) / min(v)
  if (ratio > .max_variance_ratio) {
    fgls_warn(
      sprintf(
        paste0(
          "The fitted variances span a factor of %.3g; a few observations get almost ",
          "all the weight. Consider weight_clip = c(0.01, 0.99), a smaller candidate set ",
          "z, or fgls_control(delta_scale = TRUE)."
        ),
        ratio
      ),
      "weight"
    )
  }
  v
}

known_variances <- function(v) {
  v <- as.numeric(v)
  if (any(!is.finite(v)) || any(v <= 0)) fgls_abort("Known variances must be positive and finite.")
  v / mean(v)
}

# Leverages of the weighted regression: diag of X (X'W^{-1}X)^{-1} X'W^{-1}.
weighted_leverage <- function(X, v) {
  hat_diag(X / sqrt(v))
}

# (T^{-1}X'W^{-1}X)^{-1} Omega (T^{-1}X'W^{-1}X)^{-1} / T, Omega = T^{-1} X'diag(s)X.
sandwich_cov <- function(X, v, sigma_tt) {
  n <- nrow(X)
  A <- inv_psd(crossprod(X / v, X) / n, name = "weighted moment matrix X'W^-1 X")
  Omega <- crossprod(X * sigma_tt, X) / n
  make_psd(A %*% Omega %*% A / n, name = "FGLS sandwich covariance")
}

# Diagonal Sigma_tt of the middle matrix for each covariance type.
hetero_sigma <- function(cov_type, X, v, u_fgls, u_ols, h_ols, df, method) {
  n <- nrow(X)
  k <- ncol(X)
  h_w <- pmin(pmax(weighted_leverage(X, v), 0), 1 - 1e-12)
  base <- u_fgls^2 / v^2
  if (cov_type == "robust") {
    if (method %in% c("lasso", "svr")) {
      # GP24 eq. (8): Miller-Startz / Rothenberg bias correction
      return(base * (1 / (1 - h_w)^2 + 4 * (h_ols / k) * df))
    }
    if (method %in% c("wls_s1", "wls_s2")) {
      # GP24 eq. (7): HC3-type OLS residuals
      h_o <- pmin(pmax(h_ols, 0), 1 - 1e-12)
      return(u_ols^2 / (1 - h_o)^2 / v^2)
    }
    cov_type <- "HC3"
  }
  switch(cov_type,
    HC0 = base,
    HC1 = base * n / max(n - k, 1L),
    HC2 = base / (1 - h_w),
    HC3 = base / (1 - h_w)^2,
    stop("unknown cov_type ", cov_type)
  )
}

# Covariance of the weighted (FGLS) coefficients. "nonrobust" assumes the
# variance model is correct: s^2 (X'W^{-1}X)^{-1}.
hetero_covariance <- function(cov_type, X, v, u_fgls, u_ols, h_ols, df, method) {
  n <- nrow(X)
  k <- ncol(X)
  if (cov_type == "nonrobust") {
    s2 <- sum(u_fgls^2 / v) / max(n - k, 1L)
    return(s2 * inv_psd(crossprod(X / v, X), name = "weighted moment matrix"))
  }
  sandwich_cov(X, v, hetero_sigma(cov_type, X, v, u_fgls, u_ols, h_ols, df, method))
}
