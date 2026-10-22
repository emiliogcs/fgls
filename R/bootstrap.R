# Bootstrap covariance matrices, an independent cross-check of the analytic
# variances (port of pyfgls/covariance/bootstrap.py).
#
# * Heteroskedasticity only (cross section): wild bootstrap
#   y*_t = x_t' b + u_t eta_t with Rademacher or Mammen multipliers.
# * Serial correlation (with or without heteroskedasticity): residual-based
#   sieve bootstrap with the estimated AR(k*) filter,
#   u*_t = sum_j rho_j u*_{t-j} + e*_t, with e*_t = e_t eta_t (wild, default) or
#   resampled with replacement ("iid"); the first k* errors are the observed
#   residuals.
#
# The whole estimator (lag selection, skedastic fit with cross-validation) is
# re-run on every bootstrap sample. The regressors are held fixed, so with
# non-exogenous regressors the bootstrap does not reproduce the feedback from
# past errors to the regressors.

wild_multipliers <- function(n, dist) {
  if (dist == "rademacher") {
    return(sample(c(-1, 1), n, replace = TRUE))
  }
  s5 <- sqrt(5)
  p <- (s5 + 1) / (2 * s5)
  ifelse(stats::runif(n) < p, (1 - s5) / 2, (1 + s5) / 2)
}

wild_sample <- function(fitted, resid, dist) {
  fitted + resid * wild_multipliers(length(resid), dist)
}

sieve_sample <- function(fitted, uhat, innovations, rho, scheme, dist) {
  k <- length(rho)
  e <- innovations - mean(innovations)
  e_star <- if (scheme == "wild") {
    e * wild_multipliers(length(e), dist)
  } else {
    sample(e, length(e), replace = TRUE)
  }
  if (k == 0L) {
    return(fitted + e_star)
  }
  # u*_t = sum_j rho_j u*_{t-j} + e*_t, started from the observed u_1..u_k
  tail_u <- stats::filter(e_star, rho, method = "recursive", init = rev(uhat[seq_len(k)]))
  fitted + c(uhat[seq_len(k)], as.numeric(tail_u))
}

# Bootstrap covariance of refit(y*) over n_boot samples y* = draw(). Every
# replication has its own seed, so results do not depend on `cores`.
bootstrap_covariance <- function(refit, draw, n_boot, seed = NULL, cores = 1L) {
  seeds <- with_seed(seed, sample.int(.Machine$integer.max, n_boot))
  one <- function(s) with_seed(s, refit(draw()))
  if (cores > 1L && .Platform$OS.type != "windows" && requireNamespace("parallel", quietly = TRUE)) {
    draws <- parallel::mclapply(seeds, one, mc.cores = cores)
  } else {
    draws <- lapply(seeds, one)
  }
  D <- do.call(rbind, draws)
  list(cov = stats::cov(D), draws = D, n_boot = n_boot)
}

fgls_bootstrap <- function(spec, est, res, y) {
  ctl <- spec$options$control
  X <- spec$X
  beta <- res$coefficients
  fitted <- drop(X %*% beta)
  if (spec$mode == "hetero") {
    resid <- y - fitted
    draw <- function() wild_sample(fitted, resid, ctl$wild_dist)
    label <- sprintf("bootstrap (wild, %s)", ctl$wild_dist)
  } else {
    rho <- if (spec$serial == "ma1") numeric(0) else unname(res$serial$rho)
    k <- length(rho)
    uhat <- y - fitted
    innov <- ar_filter(uhat, rho)
    draw <- function() sieve_sample(fitted, uhat, innov, rho, ctl$bootstrap, ctl$wild_dist)
    label <- sprintf("bootstrap (sieve AR(%d), %s)", k, ctl$bootstrap)
  }
  refit <- function(ystar) {
    b <- suppressWarnings(fgls_estimate(spec, ystar, "robust"))$beta
    as.numeric(b)
  }
  bseed <- if (is.null(ctl$boot_seed)) spec$options$seed else ctl$boot_seed
  boot <- bootstrap_covariance(refit, draw, ctl$n_boot, seed = bseed, cores = ctl$cores)
  nm <- names(beta)
  V <- boot$cov
  dimnames(V) <- list(nm, nm)
  colnames(boot$draws) <- nm
  res$vcov <- V
  res$vcov_type <- label
  res$notes <- sprintf("Covariance: %s, %d replications.", label, ctl$n_boot)
  res$bootstrap <- list(draws = boot$draws, n_boot = ctl$n_boot, scheme = label)
  res
}
