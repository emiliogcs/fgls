# Serial-correlation FGLS of PG26 Section 4: orchestration of the three steps
# (port of pyfgls/serial/fgls.py).
#
# 1. Durbin regression with the AR order chosen by BIC on a common sample
#    (eq. (14));
# 2. OLS on quasi-differenced data (eq. (15));
# 3. naive (eq. (27)) or exogeneity-robust (eqs. (28)-(29)) variance.

# Warn when the smallest root of the estimated AR polynomial has modulus
# below this.
.near_unit_root <- 1.05

# Warn if the estimated AR filter is (nearly) non-invertible (PG26 Assumption
# 2). Never modifies rho.
check_ar_roots <- function(rho, threshold = .near_unit_root) {
  if (length(rho) == 0L) {
    return(invisible(NULL))
  }
  mod <- Mod(ar_poly_roots(rho))
  if (length(mod) == 0L) {
    return(invisible(NULL))
  }
  m <- min(mod)
  if (m <= 1) {
    fgls_warn(
      sprintf(
        paste0(
          "The estimated AR filter has a root with modulus %.3f <= 1: the errors look ",
          "non-stationary (unit root or explosive). FGLS results are unreliable; ",
          "consider differencing the data or testing for a unit root."
        ),
        m
      ),
      "nonstationary"
    )
  } else if (m < threshold) {
    fgls_warn(
      sprintf(
        paste0(
          "The estimated AR filter has a root with modulus %.3f close to the unit ",
          "circle (near unit root). Inference may be unreliable; consider ",
          "differencing or a unit-root test."
        ),
        m
      ),
      "nonstationary"
    )
  }
  invisible(NULL)
}

# Run the serial-correlation FGLS estimator of PG26 on validated data.
estimate_serial <- function(y, X, x_names = colnames(X), y_name = "y", has_const = TRUE,
                            ar_order = NULL, kmax = NULL, ic = "bic", min_lags = NULL,
                            lagged_dep_idx = integer(0), wallis_regressor = "auto",
                            cov_type = "robust", ic_penalty = "order") {
  n <- nrow(X)
  q <- ncol(X)
  selection <- NULL
  if (is.null(ar_order)) {
    km <- resolve_kmax(kmax, n, q, has_const = has_const)
    selection <- select_ar_order(y, X, x_names, km,
      ic = ic, min_lags = resolve_min_lags(min_lags, n),
      y_name = y_name, ic_penalty = ic_penalty
    )
    k <- selection$order
    if (km > 0L && selection$selected_by_ic == km) {
      hint <- if (ic_penalty == "order" && q > 2L) {
        paste0(
          " With many stochastic regressors the paper's criterion, which penalizes ",
          "only the AR order, tends to over-select: consider ic_penalty = \"params\"."
        )
      } else {
        ""
      }
      fgls_warn(
        sprintf(
          paste0(
            "The %s selected the largest order considered (kmax = %d); the errors may ",
            "need a longer autoregression (increase kmax).%s"
          ),
          toupper(ic), km, hint
        ),
        "specification"
      )
    }
    if (k > selection$kmax) { # min_lags above kmax: cap by what the sample supports
      k <- resolve_kmax(k, n, q, has_const = has_const)
    }
  } else {
    if (!is_count(ar_order)) {
      fgls_abort(sprintf("`ar_order` must be a non-negative integer, got %s.", format(ar_order)))
    }
    k <- as.integer(ar_order)
    if (n_durbin_regressors(q, k, has_const) >= n - k) {
      fgls_abort(
        sprintf("ar_order = %d is too large for T = %d observations and %d regressors.", k, n, q),
        "insufficient_data"
      )
    }
  }

  design <- durbin_design(y, X, x_names, k, start = k, y_name = y_name)
  durbin <- fit_durbin(design)
  wallis <- NULL
  J <- NULL
  if (length(lagged_dep_idx)) {
    wallis <- wallis_rho(durbin, k, x_names, lagged_dep_idx, wallis_regressor)
    rho <- wallis$rho
    J <- wallis$J
  } else {
    missing_lags <- which(is.na(design$ylag_pos[seq_len(k)]))
    if (length(missing_lags)) {
      fgls_abort(sprintf(
        paste0(
          "The lag(s) %s of %s are collinear with the regressors, so the model ",
          "seems to contain lags of the dependent variable. Declare them with ",
          "lagged_dep_vars (PG26 Section 4.2)."
        ),
        paste(missing_lags, collapse = ", "), y_name
      ))
    }
    rho <- unname(durbin$coef[design_ylag_columns(design, k)])
  }
  check_ar_roots(rho)
  filtered <- filtered_regression(y, X, rho)
  variance <- serial_variance(filtered$Xf, filtered$ehat, filtered$uhat, design,
    J = J, cov_type = cov_type
  )
  list(
    rho = rho, order = as.integer(k), selection = selection, durbin = durbin,
    filtered = filtered, variance = variance, wallis = wallis, J = J,
    beta = filtered$beta, cov = variance$cov, design = design
  )
}
