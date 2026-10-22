# Assembly of "fgls" result objects from estimator internals (port of
# pyfgls/_builders.py).

new_fgls <- function(coefficients, vcov, method, vcov_type, nobs, n_eff, use_t, residuals,
                     innovations, y, serial = NULL, skedastic = NULL, stages = list(),
                     notes = character(0), options = list()) {
  nm <- names(coefficients)
  dimnames(vcov) <- list(nm, nm)
  res <- list(
    coefficients = coefficients,
    vcov = vcov,
    residuals = residuals,
    fitted.values = stats::setNames(y - residuals, names(residuals)),
    innovations = innovations,
    method = method,
    vcov_type = vcov_type,
    nobs = as.integer(nobs),
    n_eff = as.integer(n_eff),
    df_resid = as.integer(n_eff - length(coefficients)),
    use_t = isTRUE(use_t),
    serial = serial,
    skedastic = skedastic,
    weights = if (!is.null(skedastic)) 1 / skedastic$variance,
    ar_order = serial$ar_order,
    rho = serial$rho,
    bic_path = serial$ic_path,
    nonexogeneity_share = serial$nonexogeneity_share,
    stages = stages,
    notes = notes,
    options = options
  )
  class(res) <- "fgls"
  res
}

fgls_build <- function(spec, est, y, cov_type, rn) {
  switch(spec$mode,
    serial = if (spec$serial == "ma1") {
      build_gma(spec, est, y, rn)
    } else {
      build_serial(spec, est, y, cov_type, rn)
    },
    hetero = build_hetero(spec, est, y, cov_type, rn),
    combined = build_combined(spec, est, y, cov_type, rn)
  )
}

serial_info <- function(est, spec) {
  k <- est$order
  sel <- est$selection
  v <- est$variance
  list(
    ar_order = k,
    rho = stats::setNames(est$rho, if (k > 0L) paste0("rho", seq_len(k)) else character(0)),
    method = if (is.null(sel)) "fixed" else sel$criterion,
    kmax = if (is.null(sel)) NA_integer_ else sel$kmax,
    ic_path = if (is.null(sel)) NULL else sel$values,
    durbin_coef = est$durbin$coef,
    dropped_lags = est$design$dropped,
    nonexogeneity_by_coef = v$nonexogeneity_by_coef,
    nonexogeneity_share = v$nonexogeneity_share,
    lagged_dep_vars = spec$x_names[spec$lagged_dep_idx],
    wallis_regressors = if (is.null(est$wallis)) character(0) else est$wallis$used,
    ar_roots = ar_poly_roots(est$rho),
    extra = list()
  )
}

serial_notes <- function(info, cov_type, exog) {
  notes <- switch(cov_type,
    robust = "Covariance: PG26 eq. (29), robust to non-exogenous (pre-determined) regressors.",
    naive = "Covariance: PG26 eq. (27), valid only for exogenous regressors.",
    sprintf("Covariance: PG26 eq. (29) with %s heteroskedasticity weights.", cov_type)
  )
  if (info$ar_order == 0L) {
    notes <- c(notes, "BIC selected k* = 0: no serial correlation detected; FGLS equals OLS.")
  }
  if (length(info$dropped_lags)) {
    notes <- c(notes, paste0(
      "Collinear lagged columns dropped from the Durbin regression: ",
      paste(info$dropped_lags, collapse = ", "), "."
    ))
  }
  if (exog == "exogenous" && cov_type == "robust" && info$nonexogeneity_share > 0.1) {
    notes <- c(notes, sprintf(
      paste0(
        "You declared exog = \"exogenous\", but the first-stage correction adds %.0f%% ",
        "to the variance on average: the regressors may not be exogenous."
      ),
      100 * info$nonexogeneity_share
    ))
  }
  if (length(info$lagged_dep_vars)) {
    notes <- c(notes, paste0(
      "AR coefficients recovered by the Wallis ratio from regressor(s) ",
      paste(info$wallis_regressors, collapse = ", "), " (PG26 Section 4.2)."
    ))
  }
  notes
}

build_serial <- function(spec, est, y, cov_type, rn, method = "FGLS (serial)") {
  info <- serial_info(est, spec)
  fr <- est$filtered
  k <- est$order
  new_fgls(
    coefficients = fr$beta, vcov = est$cov, method = method, vcov_type = cov_type,
    nobs = length(y), n_eff = fr$n_eff, use_t = spec$options$use_t,
    residuals = stats::setNames(fr$uhat, rn),
    innovations = stats::setNames(fr$ehat, rn[seq_len(length(y))[-seq_len(k)]]),
    y = y, serial = info, notes = serial_notes(info, cov_type, spec$options$exog),
    options = spec$options
  )
}

build_gma <- function(spec, est, y, rn) {
  q <- length(est$beta)
  info <- list(
    ar_order = 0L, rho = stats::setNames(numeric(0), character(0)), method = "ma1",
    kmax = NA_integer_, ic_path = NULL, durbin_coef = NULL, dropped_lags = character(0),
    nonexogeneity_by_coef = stats::setNames(rep(NA_real_, q), spec$x_names),
    nonexogeneity_share = NA_real_, lagged_dep_vars = character(0), # not computed for GMA
    wallis_regressors = character(0), ar_roots = complex(0),
    extra = list(theta = est$theta, durbin_order = est$durbin_order)
  )
  new_fgls(
    coefficients = est$beta, vcov = est$cov, method = "GMA (MA(1))", vcov_type = "nonrobust",
    nobs = length(y), n_eff = est$n_eff, use_t = spec$options$use_t,
    residuals = stats::setNames(est$uhat, rn), innovations = stats::setNames(est$ehat, rn),
    y = y, serial = info,
    notes = "GMA for MA(1) errors (PG26 Section 5); covariance s^2 (X*'X*)^{-1}.",
    options = spec$options
  )
}
