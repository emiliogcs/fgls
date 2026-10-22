# Golden tests against pyfgls for the heteroskedasticity estimator (scenarios 7,
# 8 and 10). R reads pyfgls' fold assignment, so cross-validation is comparable.

# Smallest penalty at which scikit-learn's LARS (used by pyfgls) is exact.
.lars_valid_lambda <- 1.2e-7

golden_Z <- function(g) {
  d <- g$data
  r <- g$res
  x <- d$x
  switch(g$spec$z,
    gp24 = {
      Z <- cbind(x, log(x)^2, x^2, cos(x), cos(2 * x))
      colnames(Z) <- r$z_names
      Z
    },
    many = {
      cols <- list(x, x^2, log(x)^2)
      for (j in 1:19) {
        w <- d[[paste0("w", j)]]
        cols <- c(cols, list(w, cos(w), cos(2 * w)))
      }
      Z <- do.call(cbind, cols)
      colnames(Z) <- r$z_names
      Z
    },
    auto = make_dictionary(golden_X(g))
  )
}

# GP24 eq. (4) objective, 1/2 RSS + lambda_paper sum_j w_j |g_j| (standardised
# coefficients, w_j = |ridge_j|^-psi), evaluated at coefficients `coef`.
lasso_objective <- function(r, Z, intercept, coef, ridge, psi, lambda_paper) {
  st <- fgls:::standardize_cols(Z)
  g_std <- coef[st$keep] * st$sd[st$keep]
  w <- ifelse(abs(ridge) > 0, abs(ridge)^(-psi), Inf)
  pen <- sum(ifelse(g_std != 0, w * abs(g_std), 0))
  0.5 * sum((r - intercept - drop(Z %*% coef))^2) + lambda_paper * pen
}

check_given_weights <- function(blk, X, y, label) {
  v <- as.numeric(blk$weights)
  b <- fgls:::lsq(X / sqrt(v), y / sqrt(v))
  expect_rel(b, blk$beta_fgls, golden_tol$given_weights, paste(label, "WLS beta given weights"))
  u <- drop(y - X %*% b)
  u0 <- drop(y - X %*% fgls:::lsq(X, y))
  h0 <- fgls:::hat_diag(X)
  lab <- if (blk$method %in% c("none", "known")) {
    if (blk$method == "none") "ols" else "known"
  } else {
    blk$method
  }
  for (ct in names(blk$vcov_types)) {
    V <- fgls:::hetero_covariance(ct, X, v, u, u0, h0, blk$df_hat, lab)
    expect_rel(V, as_mat(blk$vcov_types[[ct]], ncol(X), ncol(X)), golden_tol$given_weights,
      paste(label, ct, "given weights")
    )
  }
}

for (id in c("s7", "s8", "s10")) {
  test_that(sprintf("golden %s tier (i): ridge, CV curves, (lambda, psi), objective", id), {
    skip_if_no_golden()
    g <- read_golden(id)
    r <- g$res
    Z <- golden_Z(g)
    expect_identical(colnames(Z), as.character(r$z_names))
    X <- golden_X(g)
    y <- g$data$y
    expect_rel(fgls:::hat_diag(X), r$h_ols, golden_tol$serial, "h_ols")
    resp <- as.numeric(r$response)
    st <- fgls:::standardize_cols(Z)
    expect_identical(unname(st$keep), as.integer(r$keep))
    rid <- fgls:::ridge_cv(st$Zs, resp, as.integer(r$foldid))
    expect_rel(rid$lambda, r$lambda_ridge, golden_tol$exact, "lambda_ridge")
    expect_rel(rid$coef, r$ridge_coef, golden_tol$ridge, "ridge coefficients")
    fit <- fgls:::adaptive_lasso(resp, Z, colnames(Z), foldid = as.integer(r$foldid))
    cvm <- as_mat(r$cv_error_matrix, length(r$psi_grid), r$n_lambda)
    for (i in seq_along(r$psi_grid)) {
      p <- fit$tuning$cv_path[[i]]
      expect_equal(p$psi, r$psi_grid[i])
      expect_rel(range(p$lambda), sort(r$cv_lambda_range[i, ]), golden_tol$exact, "lambda grid")
      # pyfgls solves the Lasso by scikit-learn's LARS, which stops once the
      # penalty falls below float32 epsilon (~1.19e-7, an absolute threshold)
      # and then reuses the last knot's coefficients: its CV errors are not the
      # exact Lasso errors there (pyfgls issue, see tools/README-validation.md).
      # Compare where the reference is valid; R's solutions satisfy the KKT
      # conditions everywhere (test-hetero.R).
      ok <- p$lambda >= .lars_valid_lambda
      expect_abs(p$cv_mse[ok], cvm[i, ok], golden_tol$cv_curve_abs, sprintf("CV curve psi=%g", p$psi))
    }
    expect_identical(fit$tuning$psi, r$psi_star)
    expect_rel(fit$tuning$lambda, r$lambda_star, golden_tol$exact, "lambda*")
    obj_r <- lasso_objective(resp, Z, fit$intercept, fit$coef, rid$coef, r$psi_star, r$lambda_paper)
    obj_py <- lasso_objective(resp, Z, r$intercept, r$gamma_hat, rid$coef, r$psi_star, r$lambda_paper)
    expect_rel(obj_py, r$objective, golden_tol$lasso_objective, "objective (pyfgls, recomputed)")
    expect_rel(obj_r, r$objective, golden_tol$lasso_objective, "objective at R solution")
    expect_abs(fit$coef, r$gamma_hat, golden_tol$lasso_coef_abs, "Lasso coefficients")
    expect_identical(fit$selected, as.character(unlist(r$nonzero_set)))
    expect_equal(fit$df, r$df_hat)
    expect_rel(fit$fitted, r$fitted_logvar, golden_tol$weights, "fitted log-variance")
    expect_rel(fgls:::variances_from_log(fit$fitted), r$weights, golden_tol$weights, "weights")
    check_given_weights(r, X, y, id)
    expect_rel(fgls:::weighted_leverage(X, as.numeric(r$weights)), r$h_fgls,
      golden_tol$given_weights, "h_fgls"
    )
  })

  test_that(sprintf("golden %s tiers (ii)-(iii): estimator and lasso_fgls() with pyfgls folds", id), {
    skip_if_no_golden()
    g <- read_golden(id)
    r <- g$res
    Z <- golden_Z(g)
    X <- golden_X(g)
    ctl <- fgls_control(foldid = as.integer(r$foldid))
    est <- fgls:::estimate_hetero(g$data$y, X, Z, colnames(Z), control = ctl)
    expect_rel(est$beta, r$beta_fgls, golden_tol$weights, "beta (stage)")
    expect_rel(est$cov, as_mat(r$vcov_eq6_8, 2, 2), golden_tol$weights, "vcov (stage)")
    zarg <- if (g$spec$z == "auto") "auto" else Z
    fit <- lasso_fgls(y ~ x, data = g$data, z = zarg, control = ctl)
    expect_rel(coef(fit), r$beta_fgls, golden_tol$weights, "beta (e2e)")
    expect_rel(vcov(fit), as_mat(r$vcov_eq6_8, 2, 2), golden_tol$weights, "vcov (e2e)")
    expect_identical(fit$skedastic$selected, as.character(unlist(r$nonzero_set)))
    expect_equal(fit$skedastic$df_hat, r$df_hat)
    expect_equal(fit$skedastic$psi, r$psi_star)
    for (ct in c("nonrobust", "HC0", "HC3")) {
      expect_rel(vcov(fit, type = ct), as_mat(r$vcov_types[[ct]], 2, 2), golden_tol$weights, ct)
    }
    if (!is.null(r$public)) {
      expect_rel(coef(fit), r$public$params, golden_tol$weights, "public params")
      expect_rel(sqrt(diag(vcov(fit))), r$public$bse, golden_tol$weights, "public bse")
    }
  })
}

test_that("golden s10: the automatic dictionary matches pyfgls and nothing is selected", {
  skip_if_no_golden()
  g <- read_golden("s10")
  expect_identical(colnames(make_dictionary(golden_X(g))), as.character(g$res$dictionary_names))
  fit <- lasso_fgls(y ~ x, data = g$data, control = fgls_control(foldid = as.integer(g$res$foldid)))
  expect_length(fit$skedastic$selected, 0L)
  expect_equal(coef(fit), coef(stats::lm(y ~ x, data = g$data)), tolerance = 1e-10)
  expect_match(fit$notes, "selected no covariate", all = FALSE)
})

test_that("golden s7: 1se rule, WLS-S1/S2, OLS and SVR", {
  skip_if_no_golden()
  g <- read_golden("s7")
  r <- g$res
  Z <- golden_Z(g)
  X <- golden_X(g)
  y <- g$data$y
  ctl <- fgls_control(foldid = as.integer(r$foldid))
  # 1se rule
  o <- r$one_se
  f1 <- lasso_fgls(y ~ x, data = g$data, z = Z, lambda_rule = "1se", control = ctl)
  expect_rel(f1$skedastic$lambda, o$lambda_star, golden_tol$exact, "1se lambda")
  expect_identical(f1$skedastic$selected, as.character(unlist(o$nonzero_set)))
  expect_rel(coef(f1), o$beta_fgls, golden_tol$weights, "1se beta")
  check_given_weights(o, X, y, "1se")
  # Romano-Wolf WLS (deterministic)
  for (kind in c("wls_s1", "wls_s2")) {
    b <- r[[kind]]
    zk <- if (kind == "wls_s1") Z[, c(1, 3)] else Z
    fk <- lasso_fgls(y ~ x, data = g$data, z = zk, hetero = kind)
    expect_rel(coef(fk), b$beta_fgls, golden_tol$serial, paste(kind, "beta"))
    expect_rel(fk$skedastic$variance, b$weights, golden_tol$serial, paste(kind, "weights"))
    expect_rel(vcov(fk), as_mat(b$vcov_eq6_8, 2, 2), golden_tol$serial, paste(kind, "eq. (7) vcov"))
    expect_equal(fk$skedastic$df_hat, b$df_hat)
    check_given_weights(b, X, y, kind)
  }
  # OLS with HC3 (hetero = "none")
  fo <- lasso_fgls(y ~ x, data = g$data, hetero = "none")
  expect_rel(coef(fo), r$ols$beta_fgls, golden_tol$serial, "OLS beta")
  expect_rel(vcov(fo), as_mat(r$ols$vcov_eq6_8, 2, 2), golden_tol$serial, "OLS HC3")
  # SVR: both use libsvm; tuning is compared exactly, fits loosely
  skip_if_not_installed("e1071")
  s <- r$svr
  fs <- lasso_fgls(y ~ x, data = g$data, z = Z, hetero = "svr", control = ctl)
  expect_rel(fs$skedastic$tuning$C, s$svr$C, golden_tol$exact, "SVR C")
  expect_rel(fs$skedastic$tuning$epsilon, s$svr$epsilon, golden_tol$exact, "SVR epsilon")
  expect_rel(fs$skedastic$fitted_logvar, s$fitted_logvar, 1e-2, "SVR fitted log-variance")
  expect_rel(coef(fs), s$beta_fgls, 1e-3, "SVR beta")
  check_given_weights(s, X, y, "svr")
})
