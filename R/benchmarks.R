# Benchmark estimators used for comparisons, tests and the papers' simulations
# (port of pyfgls/benchmarks.py).

bench_data <- function(formula, data, subset_expr, env, na.action, time, allow_gaps,
                       time_series) {
  frame <- fgls_frame(formula,
    data = data, subset_expr = subset_expr, env = env, time = time,
    policy = na_policy(na.action), time_series = time_series, z = "auto"
  )
  d <- prepare_data(frame$y, frame$X,
    y_name = frame$y_name, index = frame$index,
    time_series = time_series, allow_gaps = allow_gaps
  )
  list(frame = frame, d = d)
}

simple_serial_info <- function(rho, method, names, ...) {
  k <- length(rho)
  list(
    ar_order = k, rho = stats::setNames(rho, if (k) paste0("rho", seq_len(k)) else character(0)),
    method = method, kmax = NA_integer_, ic_path = NULL, durbin_coef = NULL,
    dropped_lags = character(0),
    nonexogeneity_by_coef = stats::setNames(numeric(length(names)), names),
    nonexogeneity_share = 0, lagged_dep_vars = character(0), wallis_regressors = character(0),
    ar_roots = ar_poly_roots(rho), extra = list(...)
  )
}

bench_result <- function(bd, beta, cov, method, vcov_type, use_t, cl, model, n_eff = NULL,
                         innovations = NULL, start = 0L, serial = NULL, notes = character(0),
                         options = list()) {
  d <- bd$d
  rn <- bd$frame$rownames
  names(beta) <- d$x_names
  u <- drop(d$y - d$X %*% beta)
  inn <- if (is.null(innovations)) u else innovations
  inn_rows <- if (start > 0L) rn[-seq_len(start)] else rn
  res <- new_fgls(
    coefficients = beta, vcov = cov, method = method, vcov_type = vcov_type,
    nobs = d$nobs, n_eff = if (is.null(n_eff)) d$nobs else n_eff, use_t = use_t,
    residuals = stats::setNames(u, rn), innovations = stats::setNames(inn, inn_rows),
    y = d$y, serial = serial, notes = notes, options = options
  )
  attach_model_info(res, bd$frame, d, cl, model)
}

ols_core <- function(X, y, vcov_type, kernel = NULL, bandwidth = NULL, small_sample = FALSE,
                     const_col = NA_integer_) {
  b <- lsq(X, y)
  u <- drop(y - X %*% b)
  if (vcov_type %in% c("HAC", "NW")) {
    kern <- if (is.null(kernel)) (if (vcov_type == "HAC") "qs" else "bartlett") else kernel
    bw <- if (is.null(bandwidth)) (if (vcov_type == "HAC") "andrews" else "nw") else bandwidth
    h <- hac_covariance(X, u,
      kernel = kern, bandwidth = bw, small_sample = small_sample,
      const_col = if (is.na(const_col)) NULL else const_col
    )
    label <- sprintf("HAC (%s, bw = %.2f)", if (kern == "qs") "QS" else "Bartlett", h$bandwidth)
    list(beta = b, cov = h$cov, label = label, notes = sprintf("Covariance: %s.", label))
  } else {
    list(beta = b, cov = hc_covariance(X, u, vcov_type), label = vcov_type, notes = character(0))
  }
}

#' Benchmark estimators
#'
#' Estimators used as comparisons in the papers and in [compare_fgls()]:
#'
#' * `ols_robust()`: ordinary least squares with heteroskedasticity and
#'   autocorrelation consistent standard errors (`vcov = "HAC"`: Andrews' 1991
#'   quadratic-spectral kernel with automatic bandwidth, the benchmark of PG26;
#'   `"NW"`: Newey-West Bartlett kernel with \eqn{\lfloor 4 (T/100)^{2/9}
#'   \rfloor} lags) or White-type standard errors (`"HC0"`--`"HC3"`).
#' * `cochrane_orcutt()`: two-step (or iterated) Cochrane-Orcutt FGLS. The
#'   two-step version estimates the AR coefficients from OLS residuals and is
#'   **not valid with non-exogenous regressors** (PG26 Section 3); use
#'   [serial_fgls()] instead.
#' * `durbin_estimator()`: the coefficients on the contemporaneous regressors
#'   of the Durbin regression (PG26 eq. 14), AR order by BIC unless given.
#' * `infeasible_gls()`: GLS with known AR coefficients (`rho`) or a known
#'   error covariance matrix (`omega`, up to scale); a simulation benchmark.
#' * `hybrid_gls()`: the hybrid OLS/GLS estimator of PG26 Supplement S-1
#'   (exogenous regressors only): OLS when the Breusch-Godfrey test does not
#'   reject, otherwise quasi-differencing with a random AR(1) coefficient of the
#'   sign of the first-order residual autocorrelation.
#'
#' @inheritParams fgls
#' @param vcov Covariance estimator for `ols_robust()`.
#' @param kernel,bandwidth Optional overrides for the HAC estimator: kernel
#'   `"qs"` or `"bartlett"`, bandwidth `"andrews"`, `"nw"` or a positive number.
#' @param small_sample Multiply HAC covariances by \eqn{T/(T-q)}.
#' @param ar_order AR order (`cochrane_orcutt()`, `durbin_estimator()`).
#' @param iterate Iterate Cochrane-Orcutt to convergence.
#' @param max_iter,tol Iteration limit and convergence tolerance.
#' @param rho,omega Known AR coefficients or known error covariance matrix.
#' @param alpha Level of the Breusch-Godfrey test (`hybrid_gls()`).
#' @param lags Number of lags of the Breusch-Godfrey test.
#' @return An object of class `"fgls"` (see [fgls()]).
#' @inherit fgls references
#' @family benchmarks
#' @export
#' @examples
#' sim <- simulate_serial(300, seed = 1)
#' ols_robust(y ~ x, data = sim)
#' cochrane_orcutt(y ~ x, data = sim)
#' durbin_estimator(y ~ x, data = sim)
#' infeasible_gls(y ~ x, data = sim, rho = 0.8)
ols_robust <- function(formula, data, subset, na.action = stats::na.fail,
                       vcov = c("HAC", "NW", "HC0", "HC1", "HC2", "HC3", "nonrobust"),
                       kernel = NULL, bandwidth = NULL, small_sample = FALSE, use_t = FALSE,
                       time = NULL, allow_gaps = FALSE, model = TRUE) {
  cl <- match.call()
  vtype <- arg_default(vcov, c("HAC", "NW", "HC0", "HC1", "HC2", "HC3", "nonrobust"), "vcov")
  if (!is.null(kernel)) arg_choice(kernel, c("qs", "bartlett"), "kernel")
  bd <- bench_data(formula, if (missing(data)) NULL else data,
    if (missing(subset)) NULL else substitute(subset), parent.frame(), na.action, time,
    allow_gaps, time_series = vtype %in% c("HAC", "NW")
  )
  o <- ols_core(bd$d$X, bd$d$y, vtype, kernel, bandwidth, small_sample, bd$d$const_col)
  bench_result(bd, o$beta, o$cov, "OLS", o$label, use_t, cl, model, notes = o$notes)
}

#' @rdname ols_robust
#' @export
cochrane_orcutt <- function(formula, data, subset, na.action = stats::na.fail, ar_order = 1L,
                            iterate = FALSE, max_iter = 100L, tol = 1e-8, use_t = FALSE,
                            time = NULL, allow_gaps = FALSE, model = TRUE) {
  cl <- match.call()
  if (!is_count(ar_order) || ar_order < 1L) fgls_abort("`ar_order` must be an integer >= 1.")
  bd <- bench_data(formula, if (missing(data)) NULL else data,
    if (missing(subset)) NULL else substitute(subset), parent.frame(), na.action, time,
    allow_gaps,
    time_series = TRUE
  )
  X <- bd$d$X
  y <- bd$d$y
  p <- as.integer(ar_order)
  n <- nrow(X)
  b <- lsq(X, y)
  rho <- numeric(p)
  n_iter <- 0L
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    n_iter <- it
    u <- drop(y - X %*% b)
    L <- lag_matrix(u, p)
    rho_new <- lsq(L, u[(p + 1L):n])
    b <- lsq(ar_filter(X, rho_new), ar_filter(y, rho_new))
    done <- max(abs(rho_new - rho)) < tol
    rho <- rho_new
    if (!iterate || done) {
      converged <- TRUE
      break
    }
  }
  if (!converged) fgls_warn("Cochrane-Orcutt did not converge.", "convergence")
  Xf <- ar_filter(X, rho)
  yf <- ar_filter(y, rho)
  e <- drop(yf - Xf %*% b)
  nf <- nrow(Xf)
  V <- sum(e^2) / (nf - ncol(X)) * inv_psd(crossprod(Xf), name = "filtered X'X")
  info <- simple_serial_info(rho, "cochrane-orcutt", bd$d$x_names, iterations = n_iter)
  note <- if (iterate) {
    "Iterated Cochrane-Orcutt: NLS of the AR(p) common-factor model."
  } else {
    "Two-step Cochrane-Orcutt is not valid with non-exogenous regressors (PG26 Section 3)."
  }
  bench_result(bd, b, V, paste0("Cochrane-Orcutt", if (iterate) " (iterated)" else " (two-step)"),
    "nonrobust", use_t, cl, model,
    n_eff = nf, innovations = e, start = p, serial = info,
    notes = note
  )
}

#' @rdname ols_robust
#' @export
durbin_estimator <- function(formula, data, subset, na.action = stats::na.fail, ar_order = NULL,
                             kmax = NULL, use_t = FALSE, time = NULL, allow_gaps = FALSE,
                             model = TRUE) {
  cl <- match.call()
  bd <- bench_data(formula, if (missing(data)) NULL else data,
    if (missing(subset)) NULL else substitute(subset), parent.frame(), na.action, time,
    allow_gaps,
    time_series = TRUE
  )
  d <- bd$d
  if (is.null(ar_order)) {
    km <- resolve_kmax(kmax, d$nobs, d$k_vars, has_const = d$has_const)
    k <- select_ar_order(d$y, d$X, d$x_names, km)$order
  } else {
    if (!is_count(ar_order)) fgls_abort("`ar_order` must be a non-negative integer.")
    k <- as.integer(ar_order)
  }
  fit <- fit_durbin(durbin_design(d$y, d$X, d$x_names, k, y_name = d$y_name))
  q <- d$k_vars
  Z <- fit$design$Z
  s2 <- sum(fit$resid^2) / (nrow(Z) - ncol(Z))
  V <- s2 * inv_psd(crossprod(Z), name = "Durbin design")[seq_len(q), seq_len(q), drop = FALSE]
  rho <- unname(fit$coef[design_ylag_columns(fit$design, k)])
  bench_result(bd, unname(fit$beta), V, "Durbin", "nonrobust", use_t, cl, model,
    n_eff = nrow(Z), serial = simple_serial_info(rho, "durbin", d$x_names)
  )
}

#' @rdname ols_robust
#' @export
infeasible_gls <- function(formula, data, subset, na.action = stats::na.fail, rho = NULL,
                           omega = NULL, use_t = FALSE, time = NULL, allow_gaps = FALSE,
                           model = TRUE) {
  cl <- match.call()
  if (is.null(rho) == is.null(omega)) fgls_abort("Give exactly one of `rho` or `omega`.")
  bd <- bench_data(formula, if (missing(data)) NULL else data,
    if (missing(subset)) NULL else substitute(subset), parent.frame(), na.action, time,
    allow_gaps,
    time_series = !is.null(rho)
  )
  d <- bd$d
  if (!is.null(rho)) {
    Xw <- ar_filter(d$X, rho)
    yw <- ar_filter(d$y, rho)
    start <- length(rho)
    serial <- simple_serial_info(as.numeric(rho), "known", d$x_names)
  } else {
    Om <- as.matrix(omega)
    if (!identical(dim(Om), c(d$nobs, d$nobs))) {
      fgls_abort(sprintf("`omega` must be %d x %d.", d$nobs, d$nobs))
    }
    Lc <- t(chol(Om))
    Xw <- forwardsolve(Lc, d$X)
    yw <- drop(forwardsolve(Lc, d$y))
    start <- 0L
    serial <- NULL
  }
  b <- lsq(Xw, yw)
  e <- drop(yw - Xw %*% b)
  V <- sum(e^2) / (nrow(Xw) - ncol(Xw)) * inv_psd(crossprod(Xw), name = "whitened X'X")
  bench_result(bd, b, V, "GLS (infeasible)", "nonrobust", use_t, cl, model,
    n_eff = nrow(Xw), innovations = e, start = start, serial = serial
  )
}

#' @rdname ols_robust
#' @export
hybrid_gls <- function(formula, data, subset, na.action = stats::na.fail, alpha = 0.05,
                       lags = 1L, seed = NULL, use_t = FALSE, time = NULL, allow_gaps = FALSE,
                       model = TRUE) {
  cl <- match.call()
  check_seed(seed)
  bd <- bench_data(formula, if (missing(data)) NULL else data,
    if (missing(subset)) NULL else substitute(subset), parent.frame(), na.action, time,
    allow_gaps,
    time_series = TRUE
  )
  d <- bd$d
  X <- d$X
  y <- d$y
  n <- d$nobs
  b <- lsq(X, y)
  u <- drop(y - X %*% b)
  L <- vapply(seq_len(lags), function(j) c(rep(0, j), u[seq_len(n - j)]), numeric(n))
  aux <- cbind(X, L)
  e <- drop(u - aux %*% lsq(aux, u))
  r2 <- 1 - sum(e^2) / sum((u - mean(u))^2)
  pval <- stats::pchisq(n * r2, lags, lower.tail = FALSE)
  if (pval >= alpha) {
    o <- ols_core(X, y, "HAC", const_col = d$const_col)
    return(bench_result(bd, o$beta, o$cov, "Hybrid (OLS)", o$label, use_t, cl, model,
      notes = sprintf("Breusch-Godfrey LM p-value %.3g.", pval)
    ))
  }
  r1 <- sum(u[-1L] * u[-n]) / sum(u^2)
  rho <- with_seed(seed, stats::runif(1, 0.1, 0.9)) * if (r1 >= 0) 1 else -1
  Xf <- ar_filter(X, rho)
  yf <- ar_filter(y, rho)
  bh <- lsq(Xf, yf)
  eh <- drop(yf - Xf %*% bh)
  h <- hac_covariance(Xf, eh, const_col = if (is.na(d$const_col)) NULL else d$const_col)
  bench_result(bd, bh, h$cov, "Hybrid (GLS)", "HAC", use_t, cl, model,
    n_eff = n - 1L, innovations = eh, start = 1L,
    serial = simple_serial_info(rho, "hybrid", d$x_names, lm_pvalue = pval),
    notes = sprintf(
      "Breusch-Godfrey LM p-value %.3g; quasi-differenced with rho = %.3f.", pval, rho
    )
  )
}
