# Quasi-differencing and the second-step FGLS regression, PG26 eq. (15)
# (port of pyfgls/serial/filter.py).

# Apply A(L) = 1 - sum_j rho_j L^j to y and every column of X, dropping the
# first k = length(rho) observations (no Prais-Winsten transformation).
quasi_difference <- function(y, X, rho) {
  list(yf = ar_filter(y, rho), Xf = ar_filter(X, rho))
}

# Second FGLS step: OLS of the quasi-differenced y on the quasi-differenced X.
# With known rho this is exactly GLS for AR(k) errors conditional on the first
# k observations. Because every column (including a constant) is filtered, the
# intercept is reported on the original scale.
filtered_regression <- function(y, X, rho) {
  qd <- quasi_difference(y, X, rho)
  beta <- lsq(qd$Xf, qd$yf)
  names(beta) <- colnames(X)
  ehat <- as.numeric(qd$yf - qd$Xf %*% beta)
  list(
    beta = beta, yf = qd$yf, Xf = qd$Xf, ehat = ehat,
    uhat = as.numeric(y - X %*% beta), start = length(rho),
    n_eff = length(qd$yf), sigma2 = sum(ehat^2) / length(qd$yf)
  )
}
