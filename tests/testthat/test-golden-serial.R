# Golden tests against pyfgls for the serial-correlation estimator (scenarios 1-6).
# Tier (i): one R function with pyfgls intermediates injected; tier (ii): the
# estimator with pyfgls' k*; tier (iii): data -> fgls().

serial_ids <- c("s1", "s2", "s3", "s4", "s5", "s6")

lagged_idx <- function(g, X) {
  ld <- g$spec$lagged_dep_vars
  if (is.null(ld)) integer(0) else match(ld, colnames(X))
}

check_variance_block <- function(v, r, tol) {
  q <- nrow(v$avar)
  k <- r$kstar
  expect_rel(v$G, as_mat(r$Ghat, q, k), tol, "Ghat")
  if (!is.null(r$Sigma_psi) && length(v$Sigma)) {
    expect_rel(v$Sigma, as_mat(r$Sigma_psi, k, k), tol, "Sigma_psi")
  }
  expect_rel(v$Q, as_mat(r$Qhat, q, q), tol, "Qhat")
  expect_rel(v$sigma2, r$sigma2_e, tol, "sigma2_e")
  expect_rel(v$naive_avar, as_mat(r$AVar_naive, q, q), tol, "AVar_naive")
  expect_rel(v$avar, as_mat(r$AVar_R, q, q), tol, "AVar_R")
  expect_rel(v$cov, as_mat(r$vcov, q, q), tol, "vcov")
  expect_abs(v$nonexogeneity_share, r$nonexogeneity_share, 1e-8, "nonexogeneity_share")
}

for (id in serial_ids) {
  test_that(sprintf("golden %s tier (i): components given pyfgls intermediates", id), {
    skip_if_no_golden()
    g <- read_golden(id)
    r <- g$res
    X <- golden_X(g)
    y <- g$data$y
    k <- r$kstar
    tol <- golden_tol$serial
    # BIC path on the common sample of the pyfgls kmax
    sel <- fgls:::select_ar_order(y, X, colnames(X), r$kmax)
    expect_abs(sel$values, r$bic_path, golden_tol$exact, "bic_path")
    expect_identical(sel$selected_by_ic, as.integer(r$selected_by_ic))
    expect_identical(sel$dropped, as.character(r$dropped_kmax))
    # Durbin regression at k*
    des <- fgls:::durbin_design(y, X, colnames(X), k)
    expect_identical(des$names, as.character(r$durbin_names))
    expect_identical(des$dropped, as.character(unlist(r$dropped_kstar)))
    dfit <- fgls:::fit_durbin(des)
    expect_rel(dfit$coef, r$durbin_coef, tol, "durbin_coef")
    expect_rel(dfit$bse, r$durbin_bse, tol, "durbin_bse")
    # filter + variance given pyfgls' rho (and Wallis Jacobian)
    rho <- as.numeric(r$rho_hat)
    fr <- fgls:::filtered_regression(y, X, rho)
    expect_rel(fr$beta, r$beta_fgls, tol, "beta_fgls")
    expect_rel(fr$ehat, r$e_hat, tol, "e_hat")
    expect_rel(fr$uhat, r$u_hat, tol, "u_hat")
    J <- if (!is.null(r$wallis)) as_mat(r$wallis$J, k, length(r$durbin_coef))
    v <- fgls:::serial_variance(fr$Xf, fr$ehat, fr$uhat, des, J = J)
    check_variance_block(v, r, tol)
    for (ct in c("naive", "HC0", "HC1", "HC2", "HC3")) {
      vo <- fgls:::serial_variance(fr$Xf, fr$ehat, fr$uhat, des, J = J, cov_type = ct)
      expect_rel(vo$cov, as_mat(r$vcov_other[[ct]], ncol(X), ncol(X)), tol, ct)
    }
    expect_abs(sort(Mod(fgls:::ar_poly_roots(rho))), r$ar_roots_mod, golden_tol$exact, "roots")
  })

  test_that(sprintf("golden %s tier (ii): estimator given pyfgls' k*", id), {
    skip_if_no_golden()
    g <- read_golden(id)
    r <- g$res
    X <- golden_X(g)
    est <- suppressWarnings(fgls:::estimate_serial(g$data$y, X, colnames(X),
      ar_order = r$kstar, lagged_dep_idx = lagged_idx(g, X)
    ))
    tol <- golden_tol$serial
    expect_rel(est$rho, r$rho_hat, tol, "rho")
    expect_rel(est$durbin$beta, r$beta_durbin, tol, "beta_durbin")
    expect_rel(est$beta, r$beta_fgls, tol, "beta_fgls")
    check_variance_block(est$variance, r, tol)
    if (!is.null(r$wallis)) {
      expect_identical(est$wallis$used, as.character(r$wallis$used))
      expect_rel(est$J, as_mat(r$wallis$J, r$kstar, length(r$durbin_coef)), tol, "Wallis J")
    }
  })

  test_that(sprintf("golden %s tier (iii): end to end through fgls()", id), {
    skip_if_no_golden()
    g <- read_golden(id)
    r <- g$res
    fit <- suppressWarnings(fgls(stats::as.formula(g$spec$formula),
      data = g$data, hetero = "none",
      lagged_dep_vars = g$spec$lagged_dep_vars
    ))
    tol <- golden_tol$serial
    expect_identical(fit$ar_order, as.integer(r$kstar))
    expect_identical(fit$serial$kmax, as.integer(r$kmax))
    expect_abs(fit$bic_path, r$bic_path, golden_tol$exact, "bic_path")
    expect_rel(fit$rho, r$rho_hat, tol, "rho")
    expect_rel(coef(fit), r$beta_fgls, tol, "beta_fgls")
    expect_rel(vcov(fit), as_mat(r$vcov, ncol(fit$x), ncol(fit$x)), tol, "vcov")
    expect_abs(fit$nonexogeneity_share, r$nonexogeneity_share, 1e-8, "share")
    expect_identical(fit$n_eff, as.integer(r$n_eff))
    if (!is.null(r$public)) {
      expect_rel(coef(fit), r$public$params, tol, "public params")
      expect_rel(sqrt(diag(vcov(fit))), r$public$bse, tol, "public bse")
    }
  })
}

test_that("golden s1: fixed AR order 3 end to end", {
  skip_if_no_golden()
  g <- read_golden("s1")
  r <- g$res$fixed_order_3
  fit <- fgls(y ~ x, data = g$data, hetero = "none", ar_order = 3)
  expect_identical(fit$serial$method, "fixed")
  expect_rel(fit$rho, r$rho_hat, golden_tol$serial, "rho")
  expect_rel(coef(fit), r$beta_fgls, golden_tol$serial, "beta")
  expect_rel(vcov(fit), as_mat(r$vcov, 2, 2), golden_tol$serial, "vcov")
  for (ct in c("naive", "HC0", "HC3")) {
    expect_rel(vcov(fit, type = ct), as_mat(r$vcov_other[[ct]], 2, 2), golden_tol$serial, ct)
  }
})

test_that("golden s6: Wallis median rule with fixed order 2", {
  skip_if_no_golden()
  g <- read_golden("s6")
  r <- g$res$median
  fit <- suppressWarnings(fgls(y ~ ylag1 + x,
    data = g$data, hetero = "none", ar_order = 2,
    lagged_dep_vars = "ylag1", wallis_regressor = "median"
  ))
  expect_rel(fit$rho, r$rho_hat, golden_tol$serial, "rho")
  expect_rel(coef(fit), r$beta_fgls, golden_tol$serial, "beta")
  expect_rel(vcov(fit), as_mat(r$vcov, 3, 3), golden_tol$serial, "vcov")
})

test_that("golden s3: GMA for MA(1) errors", {
  skip_if_no_golden()
  g <- read_golden("s3")
  r <- g$res$gma
  X <- golden_X(g)
  y <- g$data$y
  # tier (i): given pyfgls' theta
  ys <- fgls:::inverse_ma_filter(y, r$theta)
  Xs <- fgls:::inverse_ma_filter(X, r$theta)
  b <- fgls:::lsq(Xs, ys)
  expect_rel(b, r$beta, golden_tol$serial, "GMA beta given theta")
  e <- ys - Xs %*% b
  V <- sum(e^2) / (nrow(X) - ncol(X)) * solve(crossprod(Xs))
  expect_rel(V, as_mat(r$vcov, 2, 2), golden_tol$serial, "GMA vcov given theta")
  expect_identical(fgls:::gma_durbin_order(nrow(X)), as.integer(r$durbin_order))
  # tier (iii): R's own MLE
  fit <- fgls(y ~ x, data = g$data, serial = "ma1", hetero = "none")
  expect_rel(fit$serial$extra$theta, r$theta, golden_tol$gma_theta, "theta")
  expect_rel(coef(fit), r$beta, golden_tol$gma_theta, "GMA beta")
  expect_rel(vcov(fit), as_mat(r$vcov, 2, 2), golden_tol$gma_theta, "GMA vcov")
})

test_that("golden s1: benchmark estimators agree with pyfgls.benchmarks", {
  skip_if_no_golden()
  g <- read_golden("s1")
  b <- g$res$benchmarks
  tol <- golden_tol$serial
  for (ct in c("HAC", "NW", "HC3")) {
    f <- ols_robust(y ~ x, data = g$data, vcov = ct)
    r <- b[[paste0("ols_", ct)]]
    expect_rel(coef(f), r$params, tol, paste("OLS", ct, "beta"))
    expect_rel(vcov(f), as_mat(r$cov, 2, 2), tol, paste("OLS", ct, "vcov"))
  }
  co <- cochrane_orcutt(y ~ x, data = g$data)
  expect_rel(coef(co), b$co$params, tol, "CO beta")
  expect_rel(vcov(co), as_mat(b$co$cov, 2, 2), tol, "CO vcov")
  expect_rel(co$rho, b$co$rho, tol, "CO rho")
  ci <- cochrane_orcutt(y ~ x, data = g$data, iterate = TRUE)
  expect_rel(coef(ci), b$co_iter$params, 1e-6, "iterated CO beta (convergence tol 1e-8)")
  du <- durbin_estimator(y ~ x, data = g$data)
  # pyfgls reports the untransformed Durbin coefficient c on the constant
  # (whose lags are dropped as collinear); fgls recovers the intercept as
  # c / (1 - sum(rho)), so the intercept is compared with the transformed pyfgls
  # value, and the slope and its variance directly
  expect_rel(coef(du)[["x"]], b$durbin$params[2], tol, "Durbin beta (x)")
  expect_rel(coef(du)[["(Intercept)"]], b$durbin$params[1] / (1 - sum(b$durbin$rho)), tol,
    "Durbin intercept, c / (1 - sum(rho))"
  )
  expect_rel(vcov(du)["x", "x"], as_mat(b$durbin$cov, 2, 2)[2, 2], tol, "Durbin var(x)")
  ig <- infeasible_gls(y ~ x, data = g$data, rho = 0.8)
  expect_rel(coef(ig), b$gls_rho$params, tol, "GLS beta")
  expect_rel(vcov(ig), as_mat(b$gls_rho$cov, 2, 2), tol, "GLS vcov")
})

test_that("golden s1: diagnostics agree with pyfgls.diagnostics", {
  skip_if_no_golden()
  g <- read_golden("s1")
  dg <- g$res$diagnostics
  for (sc in c("innovations", "errors")) {
    t <- suppressWarnings(test_predetermined(y ~ x, data = g$data, scale = sc))
    r <- dg[[paste0("predetermined_", sc)]]
    expect_identical(t$ar_order, as.integer(r$ar_order))
    expect_rel(t$table$LM, r$LM, golden_tol$serial, paste(sc, "LM"))
    expect_identical(t$table$df, as.integer(r$df))
    expect_rel(t$table$p.value, r$pvalue, golden_tol$serial, paste(sc, "p-value"))
  }
  a <- adf_test(g$data$y)
  expect_rel(a$statistic, dg$adf_y$statistic, golden_tol$serial, "ADF y")
  expect_identical(a$lags, as.integer(dg$adf_y$lags))
  expect_identical(adf_test(g$data$x)$lags, as.integer(dg$adf_x$lags))
  expect_rel(adf_test(g$data$x)$statistic, dg$adf_x$statistic, golden_tol$serial, "ADF x")
})
