# FGLS for serial correlation and heteroskedasticity, PG26 Supplement S-4
# (port of pyfgls/combined.py).
#
# 1. F-C: the serial-correlation FGLS of PG26 Section 4 (AR order by BIC,
#    quasi-differencing, OLS); residuals e_t, t = k*+1..T.
# 2. Skedastic fit: adaptive Lasso of log(max(e_t^2, delta^2)) on candidate
#    covariates built from the *unfiltered* regressors, aligned to t = k*+1..T.
# 3. F-CH: weighted least squares on the quasi-differenced regression with
#    weights 1 / v_t (Supplement eqs. (A.6)-(A.7)).
# 4. Variance (Supplement eqs. (A.4)-(A.5)): the GP24 eq. (6)/(8) sandwich on
#    the filtered regressors xf_t,
#
#      Var = A^{-1} [n^{-1} sum_t e_t^2 adj_t a_t a_t'] A^{-1} / n,
#      A = n^{-1} sum_t xf_t xf_t' / v_t,
#      adj_t = 1/(1 - h_t^{F-CH})^2 + 4 (h_t^{F-C} / q) df,
#
#    with a_t = xf_t / v_t, h^{F-CH} the leverages of the weighted filtered
#    regression, h^{F-C} those of the unweighted one and df the Lasso degrees of
#    freedom (intercept included). F-C intervals are built the same way with
#    unit weights.
#
# First-stage correction (open item A of pyfgls): the Supplement does not write
# out the PG26 eq. (29) term for the combined estimator. With
# first_stage_correction = TRUE the influence of the estimated AR coefficients
# is added to a_t,
#
#   a_t = xf_t / v_t - G_w C_t,   G_w = [n^{-1} sum_t xf_t u_{t-j} / v_t]_j,
#
# with u the unfiltered residuals and C_t the Durbin-regression influence. With
# v_t = 1 and homoskedastic innovations this reduces to eq. (29); it is zero for
# exogenous regressors. It is on by default for pre-determined regressors.

combined_cov_types <- c("robust", "HC0", "HC3")

combined_covariance <- function(Xf, e, uhat, v, start, C, h_fc, df, cov_type = "robust") {
  n <- nrow(Xf)
  q <- ncol(Xf)
  a <- Xf / v
  if (!is.null(C) && nrow(C) > 0L) {
    Gw <- g_matrix(Xf / v, uhat, start, nrow(C))
    a <- a - crossprod(C, t(Gw))
  }
  h_w <- pmin(pmax(weighted_leverage(Xf, v), 0), 1 - 1e-12)
  adj <- switch(cov_type,
    robust = 1 / (1 - h_w)^2 + 4 * (h_fc / q) * df,
    HC3 = 1 / (1 - h_w)^2,
    HC0 = rep(1, n),
    stop("unknown cov_type ", cov_type)
  )
  meat <- crossprod(a * (e^2 * adj), a) / n
  Ainv <- inv_psd(crossprod(Xf / v, Xf) / n, name = "weighted filtered moment matrix")
  make_psd(Ainv %*% meat %*% Ainv / n, name = "combined FGLS covariance")
}

estimate_combined <- function(y, X, Z, z_names = colnames(Z), x_names = colnames(X),
                              y_name = "y", has_const = TRUE, ar_order = NULL, kmax = NULL,
                              ic = "bic", min_lags = NULL, lagged_dep_idx = integer(0),
                              wallis_regressor = "auto", method = "lasso", known = NULL,
                              delta = 0.1, delta_scale = FALSE, weight_clip = NULL,
                              psi_grid = .psi_grid, cv_folds = 10L, lambda_rule = "min",
                              control = fgls_control(), seed = NULL,
                              first_stage_correction = TRUE, cov_type = "robust",
                              ic_penalty = "order") {
  ser <- estimate_serial(y, X, x_names,
    y_name = y_name, has_const = has_const, ar_order = ar_order,
    kmax = kmax, ic = ic, min_lags = min_lags, lagged_dep_idx = lagged_dep_idx,
    wallis_regressor = wallis_regressor, cov_type = "robust", ic_penalty = ic_penalty
  )
  k <- ser$order
  fr <- ser$filtered
  rows <- seq_len(nrow(X))[-seq_len(k)]
  known_k <- if (!is.null(known) && !is.function(known)) as.numeric(known)[rows] else known
  het <- estimate_hetero(fr$yf, fr$Xf, Z[rows, , drop = FALSE], z_names,
    method = method,
    delta = delta, delta_scale = delta_scale, weight_clip = weight_clip, known = known_k,
    cov_type = if (cov_type == "robust") "robust" else cov_type, psi_grid = psi_grid,
    cv_folds = cv_folds, lambda_rule = lambda_rule, control = control, seed = seed
  )
  C <- if (first_stage_correction && k > 0L) first_stage_influence(ser$design, ser$J) else NULL
  h_fc <- hat_diag(fr$Xf)
  covs <- combined_covariances(ser, het, y, X, C, h_fc, cov_type)
  list(
    serial = ser, hetero = het, cov_fc = covs$cov_fc, cov_fch = covs$cov_fch,
    first_stage_correction = first_stage_correction, cov_type = cov_type,
    beta = het$beta, beta_fc = ser$beta, order = k, C = C, h_fc = h_fc
  )
}

combined_covariances <- function(ser, het, y, X, C, h_fc, cov_type) {
  fr <- ser$filtered
  k <- ser$order
  beta <- het$beta
  e_fch <- drop(fr$yf - fr$Xf %*% beta)
  u_fch <- drop(y - X %*% beta)
  nm <- colnames(X)
  cov_fch <- combined_covariance(fr$Xf, e_fch, u_fch, het$v, k, C, h_fc, het$df, cov_type)
  cov_fc <- combined_covariance(
    fr$Xf, fr$ehat, fr$uhat, rep(1, fr$n_eff), k, C, h_fc, het$df,
    cov_type
  )
  dimnames(cov_fch) <- dimnames(cov_fc) <- list(nm, nm)
  list(cov_fc = cov_fc, cov_fch = cov_fch)
}

combined_covariance_from_estimate <- function(est, spec, y, type) {
  combined_covariances(est$serial, est$hetero, y, spec$X, est$C, est$h_fc, type)
}

build_combined <- function(spec, est, y, cov_type, rn) {
  ser <- est$serial
  het <- est$hetero
  k <- ser$order
  s_info <- serial_info(ser, spec)
  h_info <- hetero_info(het, spec$z_names)
  fr <- ser$filtered
  corr <- if (est$first_stage_correction) "with" else "without"
  note1 <- sprintf(
    paste(
      "Covariance: PG26 Supplement (A.4)-(A.5) sandwich (%s), %s the first-stage",
      "(Durbin) correction."
    ),
    cov_type, corr
  )
  notes <- note1
  if (k == 0L) notes <- c(notes, "BIC selected k* = 0: no serial correlation detected.")
  if (het$method == "lasso" && !is.null(h_info$coefficients) && !length(h_info$selected)) {
    notes <- c(notes, "The Lasso selected no covariate: F-CH equals F-C.")
  }
  if (spec$options$exog == "exogenous" && !est$first_stage_correction) {
    notes <- c(notes, "exog = \"exogenous\": first-stage correction omitted.")
  }
  inn_rows <- rn[seq_len(length(y))[-seq_len(k)]]
  fc <- new_fgls(
    coefficients = est$beta_fc, vcov = est$cov_fc, method = "F-C (serial)",
    vcov_type = cov_type, nobs = length(y), n_eff = fr$n_eff, use_t = spec$options$use_t,
    residuals = stats::setNames(fr$uhat, rn),
    innovations = stats::setNames(fr$ehat, inn_rows), y = y, serial = s_info,
    notes = sub("sandwich", "sandwich with unit weights", note1, fixed = TRUE),
    options = spec$options
  )
  beta <- est$beta
  e <- drop(fr$yf - fr$Xf %*% beta)
  new_fgls(
    coefficients = beta, vcov = est$cov_fch, method = "F-CH (serial+hetero)",
    vcov_type = cov_type, nobs = length(y), n_eff = fr$n_eff, use_t = spec$options$use_t,
    residuals = stats::setNames(drop(y - spec$X %*% beta), rn),
    innovations = stats::setNames(e / sqrt(het$v), inn_rows), y = y,
    serial = s_info, skedastic = h_info, stages = list(serial = fc), notes = notes,
    options = spec$options
  )
}
