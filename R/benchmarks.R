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

# Serial block of a benchmark fit. The non-exogeneity share of PG26 eq. (29)
# is not computed for these estimators, so it is NA (not 0, which would read as
# evidence of exogenous regressors).
simple_serial_info <- function(rho, method, names, ...) {
  k <- length(rho)
  list(
    ar_order = k, rho = stats::setNames(rho, if (k) paste0("rho", seq_len(k)) else character(0)),
    method = method, kmax = NA_integer_, ic_path = NULL, durbin_coef = NULL,
    dropped_lags = character(0),
    nonexogeneity_by_coef = stats::setNames(rep(NA_real_, length(names)), names),
    nonexogeneity_share = NA_real_, lagged_dep_vars = character(0),
    wallis_regressors = character(0),
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
#' * `durbin_estimator()`: the Durbin regression (PG26 eq. 14) of \eqn{y_t}
#'   on \eqn{y_{t-j}}, \eqn{x_t} and \eqn{x_{t-j}}, \eqn{j = 1, \dots, k},
#'   with the AR order \eqn{k} chosen by BIC unless given. The coefficients
#'   on the stochastic regressors are those on \eqn{x_t}. Lags of the
#'   constant, of trends and of periodic dummies are collinear with the
#'   regressors and dropped, so the Durbin coefficients \eqn{c} on these
#'   regressors estimate \eqn{(I - \sum_j \rho_j P_j)\beta}, where the
#'   dropped lags are \eqn{x_{t-j} = x_t P_j} (plus kept lags). They are
#'   converted back to \eqn{\hat\beta = (I - \sum_j \hat\rho_j
#'   P_j)^{-1} c}, e.g. \eqn{c / (1 - \sum_j \hat\rho_j)} for the
#'   intercept, with a delta-method covariance, and are `NA` (with a note)
#'   when that is not possible. Residuals, fitted values and predictions use
#'   these coefficients; `innovations` are the residuals of the Durbin
#'   regression. (The Python package 'pyfgls' reports the untransformed \eqn{c}.)
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
#' @param alpha Level of the Breusch-Godfrey test (`hybrid_gls()`), strictly
#'   between 0 and 1.
#' @param lags Number of lags of the Breusch-Godfrey test, a positive integer.
#' @param seed Optional seed for the random AR(1) coefficient that
#'   `hybrid_gls()` draws when the Breusch-Godfrey test rejects. The session's
#'   random number generator state is restored afterwards.
#' @return An object of class `"fgls"` (see [fgls()]). The non-exogeneity
#'   share of PG26 eq. (29) is not computed for these estimators and is `NA`.
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
  # the AR(p) regression of the OLS residuals needs n - p > p rows and the
  # quasi-differenced regression n - p > q rows (positive degrees of freedom)
  if (n - p <= max(p, ncol(X))) {
    fgls_abort(
      sprintf(
        paste0(
          "Cochrane-Orcutt with ar_order = %d needs more than %d observations ",
          "for %d regressors; got %d."
        ),
        p, p + max(p, ncol(X)), ncol(X), n
      ),
      "insufficient_data"
    )
  }
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
    km <- resolve_kmax(kmax, d$nobs, d$k_vars, has_const = d$has_const, y = d$y, X = d$X)
    k <- select_ar_order(d$y, d$X, d$x_names, km)$order
  } else {
    if (!is_count(ar_order)) fgls_abort("`ar_order` must be a non-negative integer.")
    k <- as.integer(ar_order)
    check_ar_order(k, d$y, d$X, d$has_const)
  }
  fit <- fit_durbin(durbin_design(d$y, d$X, d$x_names, k, y_name = d$y_name))
  Z <- fit$design$Z
  s2 <- sum(fit$resid^2) / (nrow(Z) - ncol(Z))
  V <- s2 * inv_psd(crossprod(Z), name = "Durbin design")
  rho <- unname(fit$coef[design_ylag_columns(fit$design, k)])
  st <- durbin_structural(fit, d$X, V)
  bench_result(bd, st$beta, st$cov, "Durbin", "nonrobust", use_t, cl, model,
    n_eff = nrow(Z), innovations = fit$resid, start = k,
    serial = simple_serial_info(rho, "durbin", d$x_names), notes = st$notes
  )
}

# Structural coefficients from the Durbin regression of PG26 eq. (14).
#
# Under A(L) u_t = e_t the Durbin regression has coefficients beta on x_t,
# rho_j on y_{t-j} and -rho_j beta on x_{t-j}. When a lag x_{t-j,i} is dropped
# because it is a linear combination of the kept columns (lags of the
# constant, of trends and of periodic dummies), its term -rho_j beta_i x_{t-j,i}
# is absorbed by those columns. Writing the dropped lag as
# x_{t-j,i} = x_{t,E}' r_{ji} + (kept lags), with E the regressors that have
# dropped lags, the coefficients c on x_{t,E} estimate
#
#   c = M beta_E,   M = I - sum_j rho_j P_j,   P_j = sum_i r_{ji} e_i',
#
# (for D_{t-j} = D_t P_j this is the P_j of the lag relation; for a constant,
# c = (1 - sum_j rho_j) beta_0). Hence beta_E = M^{-1} c, with the
# delta-method covariance from the covariance V of all the Durbin coefficients:
# d beta_E / d c = M^{-1} and d beta_E / d rho_j = M^{-1} P_j beta_E. The
# coefficients on the other (stochastic) regressors are the Durbin ones,
# unchanged. When a dropped lag is not such a combination, a lag of y is
# collinear with the regressors or M is singular, the coefficients on E are NA.
durbin_structural <- function(fit, X, V) {
  des <- fit$design
  k <- des$order
  q <- des$n_x
  p <- ncol(des$Z)
  beta <- unname(fit$beta)
  Vq <- V[seq_len(q), seq_len(q), drop = FALSE]
  drops <- which(is.na(des$xlag_pos[seq_len(k), , drop = FALSE]), arr.ind = TRUE)
  if (k == 0L || !nrow(drops)) {
    return(list(beta = beta, cov = Vq, notes = character(0)))
  }
  E <- sort(unique(drops[, 2L]))
  e_names <- paste(des$names[E], collapse = ", ")
  not_recovered <- function(why) {
    beta[E] <- NA_real_
    Vq[E, ] <- NA_real_
    Vq[, E] <- NA_real_
    list(beta = beta, cov = Vq, notes = sprintf(
      paste0(
        "Durbin estimator: the coefficients on %s are NA because %s (their lags are ",
        "dropped from the Durbin regression; see ?durbin_estimator)."
      ),
      e_names, why
    ))
  }
  ypos <- des$ylag_pos[seq_len(k)]
  if (anyNA(ypos)) {
    return(not_recovered("lags of the response are collinear with the regressors"))
  }
  rows <- (k + 1L):nrow(X)
  # dropped lags on x_{t,E} and the kept lag columns, on the Durbin sample
  B <- cbind(des$Z[, E, drop = FALSE], des$Z[, seq_len(p - q) + q, drop = FALSE])
  A <- vapply(seq_len(nrow(drops)), function(a) X[rows - drops[a, 1L], drops[a, 2L]],
    numeric(length(rows)))
  A <- matrix(A, nrow = length(rows))
  R <- as.matrix(lsq(B, A))
  res <- A - B %*% R
  if (any(sqrt(colSums(res^2)) > 1e-8 * column_norms(A))) {
    return(not_recovered("their dropped lags are not linear combinations of the kept columns"))
  }
  m <- length(E)
  P <- replicate(k, matrix(0, m, m), simplify = FALSE)
  for (a in seq_len(nrow(drops))) {
    j <- drops[a, 1L]
    i <- match(drops[a, 2L], E)
    P[[j]][, i] <- P[[j]][, i] + R[seq_len(m), a]
  }
  rho <- fit$coef[ypos]
  M <- diag(m)
  for (j in seq_len(k)) M <- M - rho[j] * P[[j]]
  sv <- svd(M, 0L, 0L)$d
  if (min(sv) <= 1e-8 * max(1, sv[1L])) {
    return(not_recovered(
      "the estimated AR polynomial has a unit root that leaves them unidentified"
    ))
  }
  Minv <- solve(M)
  gamma <- drop(Minv %*% beta[E])
  Jac <- diag(1, q, p) # derivatives of the reported coefficients
  Jac[E, E] <- Minv
  for (j in seq_len(k)) Jac[E, ypos[j]] <- Minv %*% P[[j]] %*% gamma
  beta[E] <- gamma
  list(
    beta = beta, cov = symmetrize(Jac %*% V %*% t(Jac)),
    notes = sprintf(
      paste0(
        "Durbin estimator: the coefficients on %s, whose lags are collinear with the ",
        "regressors, are recovered as (I - sum_j rho_j P_j)^{-1} c from the Durbin ",
        "coefficients c (delta-method covariance); see ?durbin_estimator."
      ),
      e_names
    )
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
  if (nrow(Xw) <= ncol(Xw)) {
    fgls_abort(
      sprintf(
        paste0(
          "The quasi-differenced sample has %d observations for %d regressors; ",
          "need more observations than regressors."
        ),
        nrow(Xw), ncol(Xw)
      ),
      "insufficient_data"
    )
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
  if (!is_count(lags) || lags < 1L) {
    fgls_abort("`lags` must be a positive integer (the order of the Breusch-Godfrey test).")
  }
  if (!(is.numeric(alpha) && length(alpha) == 1L && !is.na(alpha) && alpha > 0 && alpha < 1)) {
    fgls_abort("`alpha` must be a number strictly between 0 and 1.")
  }
  lags <- as.integer(lags)
  bd <- bench_data(formula, if (missing(data)) NULL else data,
    if (missing(subset)) NULL else substitute(subset), parent.frame(), na.action, time,
    allow_gaps,
    time_series = TRUE
  )
  d <- bd$d
  X <- d$X
  y <- d$y
  n <- d$nobs
  if (n <= ncol(X) + lags) {
    fgls_abort(
      sprintf(
        "The Breusch-Godfrey test with %d lag(s) needs more than %d observations; got %d.",
        lags, ncol(X) + lags, n
      ),
      "insufficient_data"
    )
  }
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
