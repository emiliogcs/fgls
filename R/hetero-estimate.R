# Heteroskedasticity FGLS of GP24 Section 2: orchestration (port of
# pyfgls/hetero/fgls.py and the heteroskedasticity part of _builders.py).
#
# 1. unweighted regression (OLS) and its residuals u_t;
# 2. skedastic fit of log(max(u_t^2, delta^2)) on the candidate covariates z_t;
# 3. weighted least squares with weights 1/exp(fitted) (rows scaled by
#    exp(-fitted/2), then OLS);
# 4. sandwich variance with the bias correction of eq. (8) (eq. (7) for the WLS
#    fits).
#
# The function works on whatever (y, X) it is given; the combined estimator
# passes quasi-differenced data, so that step 1 is the F-C regression and step 3
# gives F-CH.

fit_skedastic <- function(r, Z, names, method, psi_grid = .psi_grid, cv_folds = 10L,
                          lambda_rule = "min", control = fgls_control(), seed = NULL) {
  switch(method,
    lasso = adaptive_lasso(r, Z, names,
      psi_grid = psi_grid, cv_folds = cv_folds,
      lambda_rule = lambda_rule, fold_scheme = control$fold_scheme,
      n_lambda = control$n_lambda, lambda_min_ratio = control$lambda_min_ratio,
      seed = seed, foldid = control$foldid, cv_weights = control$cv_weights,
      ridge_scale = control$ridge_scale, engine = control$lasso_engine,
      thresh = control$glmnet_thresh
    ),
    svr = svr_skedastic(r, Z, names,
      cost_grid = control$svr_cost, eps_grid = control$svr_epsilon,
      cv_folds = cv_folds, fold_scheme = control$fold_scheme, seed = seed,
      foldid = control$foldid
    ),
    wls_s1 = ,
    wls_s2 = wls_skedastic(r, Z, names, kind = method),
    fgls_abort(sprintf("Unknown skedastic method '%s'.", method))
  )
}

estimate_hetero <- function(y, X, Z, z_names = colnames(Z), method = "lasso", delta = 0.1,
                            delta_scale = FALSE, weight_clip = NULL, known = NULL,
                            cov_type = "robust", psi_grid = .psi_grid, cv_folds = 10L,
                            lambda_rule = "min", control = fgls_control(), seed = NULL) {
  n <- nrow(X)
  k <- ncol(X)
  b0 <- lsq(X, y)
  u0 <- drop(y - X %*% b0)
  h0 <- hat_diag(X)
  floor_ <- delta
  if (isTRUE(delta_scale)) {
    floor_ <- delta * sqrt(sum(u0^2) / max(n - k, 1L))
    if (!(floor_ > 0)) floor_ <- delta
  }
  fit <- NULL
  if (method == "none") {
    v <- rep(1, n)
  } else if (method == "known") {
    v <- known_variances(known)
    if (length(v) != n) fgls_abort(sprintf("Known variances have length %d; expected %d.", length(v), n))
  } else {
    r <- log_squared_residuals(u0, floor_)
    fit <- fit_skedastic(r, Z, z_names, method,
      psi_grid = psi_grid, cv_folds = cv_folds,
      lambda_rule = lambda_rule, control = control, seed = seed
    )
    fit$response <- r
    v <- variances_from_log(fit$fitted, clip = weight_clip)
  }
  s <- 1 / sqrt(v)
  beta <- lsq(X * s, y * s)
  names(beta) <- colnames(X)
  u <- drop(y - X %*% beta)
  label <- if (!is.null(fit)) method else if (method == "none") "ols" else "known"
  df <- if (is.null(fit)) 1 else fit$df
  cov <- hetero_covariance(cov_type, X, v, u, u0, h0, df, label)
  list(
    beta = beta, cov = cov, v = v, fit = fit, beta_unweighted = b0, u_unweighted = u0,
    h_unweighted = h0, u_weighted = u, delta = floor_, method = method, var_label = label,
    cov_type = cov_type, df = df
  )
}

hetero_covariance_from_estimate <- function(est, X, type) {
  hetero_covariance(
    type, X, est$v, est$u_weighted, est$u_unweighted, est$h_unweighted, est$df,
    est$var_label
  )
}

# Candidate covariates from the `z` option.
resolve_candidates <- function(z, d, frame, transforms = .transforms, skip = integer(0)) {
  if (is.character(z) && length(z) == 1L && identical(z, "auto")) {
    out <- dictionary_impl(d$X, d$x_names, transforms, skip)
  } else if (!is.null(frame$zmat)) {
    out <- list(Z = frame$zmat, names = colnames(frame$zmat))
  } else if (is.character(z)) {
    idx <- resolve_columns(z, d$x_names, "z")
    out <- list(Z = d$X[, idx, drop = FALSE], names = d$x_names[idx])
  } else {
    fgls_abort(paste0(
      "`z` must be \"auto\", a character vector of regressor names, a one-sided ",
      "formula, or a numeric matrix / data frame."
    ))
  }
  if (any(!is.finite(out$Z))) fgls_abort("The candidate covariates `z` have non-finite values.")
  out
}

hetero_info <- function(est, z_names) {
  fit <- est$fit
  info <- list(
    method = est$method, coefficients = NULL, selected = character(0),
    variance = est$v, df_hat = est$df, delta = est$delta
  )
  if (!is.null(fit)) {
    if (fit$method != "svr") {
      info$coefficients <- stats::setNames(c(fit$intercept, fit$coef), c("(Intercept)", z_names))
    }
    info$selected <- fit$selected
    info$fitted_logvar <- fit$fitted
    tun <- fit$tuning
    info$psi <- tun$psi
    info$lambda <- tun$lambda
    info$lambda_paper <- tun$lambda_paper
    info$cv_mse <- tun$cv_mse
    info$tuning <- tun
  }
  info
}

hetero_notes <- function(est, info, cov_type) {
  notes <- character(0)
  if (cov_type == "robust") {
    notes <- if (est$method %in% c("lasso", "svr")) {
      "Covariance: GP24 eqs. (6)+(8), sandwich with bias correction."
    } else if (est$method %in% c("wls_s1", "wls_s2")) {
      "Covariance: GP24 eq. (7), HC3-type sandwich (Romano-Wolf)."
    } else {
      "Covariance: HC3 sandwich."
    }
  }
  if (est$method == "lasso" && !is.null(info$coefficients) && !length(info$selected)) {
    notes <- c(notes, "The Lasso selected no covariate: constant weights, FGLS equals OLS.")
  }
  notes
}

build_hetero <- function(spec, est, y, cov_type, rn) {
  info <- hetero_info(est, spec$z_names)
  method <- switch(spec$hetero,
    none = "OLS",
    known = "GLS (known variances)",
    "FGLS (hetero)"
  )
  new_fgls(
    coefficients = est$beta, vcov = est$cov, method = method, vcov_type = cov_type,
    nobs = length(y), n_eff = length(y), use_t = spec$options$use_t,
    residuals = stats::setNames(est$u_weighted, rn),
    innovations = stats::setNames(est$u_weighted / sqrt(est$v), rn),
    y = y, skedastic = info, notes = hetero_notes(est, info, cov_type),
    options = spec$options
  )
}
