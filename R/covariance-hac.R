# Heteroskedasticity and autocorrelation consistent (HAC) covariance for OLS
# (port of pyfgls/covariance/hac.py):
#
# * Andrews (1991) quadratic-spectral kernel k(x) = (3/z^2)(sin z / z - cos z),
#   z = 6 pi x / 5, with the automatic bandwidth 1.3221 (alpha(2) T)^(1/5) from
#   AR(1) approximations of the score components (the benchmark of PG26);
# * Newey and West (1987) Bartlett kernel with floor(4 (T/100)^(2/9)) lags, or
#   Andrews' automatic bandwidth 1.1447 (alpha(1) T)^(1/3).

qs_kernel <- function(x) {
  out <- rep(1, length(x))
  nz <- x != 0
  z <- 6 * pi * x[nz] / 5
  out[nz] <- 3 / z^2 * (sin(z) / z - cos(z))
  out
}

bartlett_kernel <- function(x) {
  pmax(1 - abs(x), 0)
}

# Andrews (1991) automatic bandwidth from AR(1) fits to each score column;
# `skip` excludes one column (conventionally the intercept score).
andrews_bandwidth <- function(scores, kernel = "qs", skip = NULL) {
  n <- nrow(scores)
  q <- ncol(scores)
  num <- 0
  den <- 0
  for (a in seq_len(q)) {
    if (!is.null(skip) && !is.na(skip) && a == skip && q > 1L) next
    s <- scores[, a]
    s0 <- s[-n]
    s1 <- s[-1L]
    denom <- sum(s0 * s0)
    rho <- if (denom > 0) sum(s0 * s1) / denom else 0
    rho <- min(max(rho, -0.97), 0.97)
    sig2 <- mean((s1 - rho * s0)^2)
    num <- num + if (kernel == "qs") {
      4 * rho^2 * sig2^2 / (1 - rho)^8
    } else {
      4 * rho^2 * sig2^2 / ((1 - rho)^6 * (1 + rho)^2)
    }
    den <- den + sig2^2 / (1 - rho)^4
  }
  if (den <= 0) {
    return(1)
  }
  alpha <- num / den
  if (kernel == "qs") 1.3221 * (alpha * n)^0.2 else 1.1447 * (alpha * n)^(1 / 3)
}

newey_west_lags <- function(n) {
  as.integer(floor(4 * (n / 100)^(2 / 9)))
}

# Kernel estimate of lim Var(T^{-1/2} sum_t s_t) and the bandwidth used.
long_run_covariance <- function(scores, kernel = "qs", bandwidth = "andrews", skip = NULL) {
  n <- nrow(scores)
  if (identical(bandwidth, "andrews")) {
    bw <- andrews_bandwidth(scores, kernel, skip = skip)
  } else if (identical(bandwidth, "nw")) {
    bw <- if (kernel == "bartlett") newey_west_lags(n) + 1 else andrews_bandwidth(scores, kernel, skip = skip)
  } else {
    bw <- as.numeric(bandwidth)
    if (!(bw > 0)) fgls_abort("`bandwidth` must be positive.")
  }
  S <- crossprod(scores) / n
  lags <- seq_len(n - 1L)
  w <- if (kernel == "qs") qs_kernel(lags / bw) else bartlett_kernel(lags / bw)
  for (j in lags) {
    wj <- w[j]
    if (kernel == "bartlett" && wj == 0) break
    if (kernel == "qs" && abs(wj) < 1e-12 && j > 50 * bw) break
    Gj <- crossprod(scores[(j + 1L):n, , drop = FALSE], scores[seq_len(n - j), , drop = FALSE]) / n
    S <- S + wj * (Gj + t(Gj))
  }
  list(S = make_psd(S, name = "long-run covariance"), bandwidth = bw)
}

# HAC covariance (X'X)^{-1} T S (X'X)^{-1} of OLS coefficients.
hac_covariance <- function(X, resid, kernel = "qs", bandwidth = "andrews",
                           small_sample = FALSE, const_col = NULL) {
  n <- nrow(X)
  q <- ncol(X)
  lr <- long_run_covariance(X * resid, kernel = kernel, bandwidth = bandwidth, skip = const_col)
  A <- inv_psd(crossprod(X) / n, name = "X'X")
  V <- A %*% lr$S %*% A / n
  if (small_sample) V <- V * n / (n - q)
  list(cov = make_psd(V, name = "HAC covariance"), bandwidth = lr$bandwidth)
}

# White-type covariances (X'X)^{-1} X' diag(w_t u_t^2) X (X'X)^{-1}.
hc_covariance <- function(X, resid, cov_type = "HC3") {
  n <- nrow(X)
  q <- ncol(X)
  A <- inv_psd(crossprod(X), name = "X'X")
  u2 <- resid^2
  if (cov_type == "nonrobust") {
    return(sum(u2) / (n - q) * A)
  }
  w <- switch(cov_type,
    HC0 = u2,
    HC1 = u2 * n / (n - q),
    HC2 = ,
    HC3 = {
      h <- pmin(pmax(hat_diag(X), 0), 1 - 1e-12)
      if (cov_type == "HC2") u2 / (1 - h) else u2 / (1 - h)^2
    },
    fgls_abort(sprintf("Unknown covariance type '%s' (HC0-HC3 or nonrobust).", cov_type))
  )
  make_psd(A %*% crossprod(X * w, X) %*% A, name = "HC covariance")
}
