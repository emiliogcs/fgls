# Golden tests against pyfgls for the combined F-C / F-CH estimator (scenario 9
# and the combined part of scenario 10).

test_that("golden s9 tier (i): F-CH given pyfgls' k* and weights", {
  skip_if_no_golden()
  g <- read_golden("s9")
  r <- g$res
  X <- golden_X(g)
  y <- g$data$y
  k <- r$kstar
  ser <- fgls:::estimate_serial(y, X, colnames(X), ar_order = k)
  expect_rel(ser$rho, r$rho_hat, golden_tol$serial, "rho")
  expect_rel(ser$beta, r$beta_FC, golden_tol$serial, "beta F-C")
  fr <- ser$filtered
  v <- as.numeric(r$weights)
  b <- fgls:::lsq(fr$Xf / sqrt(v), fr$yf / sqrt(v))
  expect_rel(b, r$beta_FCH, golden_tol$given_weights, "beta F-CH given weights")
  C <- fgls:::first_stage_influence(ser$design)
  h_fc <- fgls:::hat_diag(fr$Xf)
  expect_rel(h_fc, r$h_fc, golden_tol$given_weights, "h F-C")
  expect_rel(fgls:::weighted_leverage(fr$Xf, v), r$h_fch, golden_tol$given_weights, "h F-CH")
  e <- drop(fr$yf - fr$Xf %*% b)
  u <- drop(y - X %*% b)
  V <- fgls:::combined_covariance(fr$Xf, e, u, v, k, C, h_fc, r$df_hat)
  expect_rel(V, as_mat(r$vcov_FCH, 2, 2), golden_tol$given_weights, "vcov F-CH given weights")
  Vc <- fgls:::combined_covariance(fr$Xf, fr$ehat, fr$uhat, rep(1, fr$n_eff), k, C, h_fc, r$df_hat)
  expect_rel(Vc, as_mat(r$vcov_FC, 2, 2), golden_tol$given_weights, "vcov F-C")
  V0 <- fgls:::combined_covariance(fr$Xf, e, u, v, k, NULL, h_fc, r$df_hat)
  expect_rel(V0, as_mat(r$no_first_stage$vcov_FCH, 2, 2), golden_tol$given_weights,
    "vcov F-CH without first-stage correction"
  )
})

test_that("golden s9 tiers (ii)-(iii): combined estimator with pyfgls folds", {
  skip_if_no_golden()
  g <- read_golden("s9")
  r <- g$res
  X <- golden_X(g)
  ctl <- fgls_control(foldid = as.integer(r$foldid))
  Z <- make_dictionary(X)
  expect_identical(colnames(Z), as.character(r$z_names))
  est <- fgls:::estimate_combined(g$data$y, X, Z, colnames(Z), ar_order = r$kstar, control = ctl)
  expect_rel(est$beta, r$beta_FCH, golden_tol$given_weights, "beta F-CH (stage)")
  expect_rel(est$cov_fch, as_mat(r$vcov_FCH, 2, 2), golden_tol$given_weights, "vcov F-CH (stage)")
  expect_rel(est$hetero$v, r$weights, golden_tol$weights, "weights (stage)")
  fit <- fgls(y ~ x, data = g$data, control = ctl)
  expect_identical(fit$ar_order, as.integer(r$kstar))
  expect_identical(fit$skedastic$selected, as.character(unlist(r$nonzero_set)))
  expect_equal(fit$skedastic$psi, r$psi_star)
  expect_rel(coef(fit), r$beta_FCH, golden_tol$fch_beta_e2e, "beta F-CH (e2e)")
  expect_rel(sqrt(diag(vcov(fit))), sqrt(diag(as_mat(r$vcov_FCH, 2, 2))), golden_tol$fch_se_e2e, "se F-CH (e2e)")
  expect_rel(coef(fit$stages$serial), r$beta_FC, golden_tol$serial, "beta F-C (e2e)")
  expect_rel(vcov(fit$stages$serial), as_mat(r$vcov_FC, 2, 2), golden_tol$fch_se_e2e, "vcov F-C (e2e)")
  f0 <- fgls(y ~ x, data = g$data, control = ctl, first_stage_correction = FALSE)
  expect_rel(vcov(f0), as_mat(r$no_first_stage$vcov_FCH, 2, 2), golden_tol$fch_se_e2e, "no first stage")
})

test_that("golden s10: homoskedastic white noise, F-CH = F-C = OLS", {
  skip_if_no_golden()
  g <- read_golden("s10")
  r <- g$res$combined
  fit <- fgls(y ~ x, data = g$data, control = fgls_control(foldid = as.integer(r$foldid)))
  expect_identical(fit$ar_order, as.integer(r$kstar))
  expect_length(fit$skedastic$selected, 0L)
  expect_rel(coef(fit), r$beta_FCH, golden_tol$given_weights, "beta")
  expect_rel(vcov(fit), as_mat(r$vcov_FCH, 2, 2), golden_tol$given_weights, "vcov")
  expect_equal(coef(fit), coef(stats::lm(y ~ x, data = g$data)), tolerance = 1e-10)
})
