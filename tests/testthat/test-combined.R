comb <- simulate_combined(250, ar = 0.5, skedastic = "x2", seed = 21)
paper_w <- function(x) {
  W <- cbind(x, x^2, cos(x), cos(2 * x))
  colnames(W) <- c("w", "w2", "cw", "c2w")
  W
}

test_that("unit weights make F-CH equal F-C", {
  X <- cbind("(Intercept)" = 1, x = comb$x)
  est <- fgls:::estimate_combined(comb$y, X, matrix(0, 250, 0), character(0),
    method = "known", known = rep(1, 250)
  )
  expect_equal(est$beta, est$beta_fc, tolerance = 1e-10)
  expect_equal(est$cov_fch, est$cov_fc, tolerance = 1e-10)
})

test_that("without first-stage correction the variance is GP24 eq. (8) on filtered data", {
  X <- cbind("(Intercept)" = 1, x = comb$x)
  est <- fgls:::estimate_combined(comb$y, X, paper_w(comb$x),
    seed = 1,
    first_stage_correction = FALSE
  )
  fr <- est$serial$filtered
  het <- est$hetero
  V <- fgls:::hetero_covariance(
    "robust", fr$Xf, het$v, drop(fr$yf - fr$Xf %*% het$beta),
    fr$ehat, fgls:::hat_diag(fr$Xf), het$df, "lasso"
  )
  expect_equal(unname(est$cov_fch), unname(V), tolerance = 1e-8)
})

test_that("the first-stage correction with unit weights matches the serial HC0 variance", {
  X <- cbind("(Intercept)" = 1, x = comb$x)
  ser <- fgls:::estimate_serial(comb$y, X, ar_order = 2)
  fr <- ser$filtered
  C <- fgls:::first_stage_influence(ser$design)
  V <- fgls:::combined_covariance(fr$Xf, fr$ehat, fr$uhat, rep(1, fr$n_eff), 2, C,
    fgls:::hat_diag(fr$Xf), 1,
    cov_type = "HC0"
  )
  ref <- fgls:::serial_variance(fr$Xf, fr$ehat, fr$uhat, ser$design, cov_type = "HC0")
  expect_equal(unname(V), unname(ref$cov), tolerance = 1e-8)
})

test_that("the first-stage correction increases the variance with a non-exogenous regressor", {
  s <- simulate_combined(400, ar = 0.8, skedastic = "x2", gamma = 1, seed = 1)
  a <- fgls(y ~ x, data = s, seed = 1)
  b <- fgls(y ~ x, data = s, seed = 1, first_stage_correction = FALSE)
  expect_equal(coef(a), coef(b))
  expect_gt(sqrt(vcov(a)[2, 2]), sqrt(vcov(b)[2, 2]))
  c <- fgls(y ~ x, data = s, seed = 1, exog = "exogenous")
  expect_equal(vcov(c), vcov(b))
  expect_match(c$notes, "first-stage correction omitted", all = FALSE)
})

test_that("fgls() result structure for the combined estimator", {
  fit <- fgls(y ~ x, data = comb, seed = 0)
  expect_identical(fit$method, "F-CH (serial+hetero)")
  expect_gte(fit$ar_order, 1L)
  expect_identical(fit$skedastic$method, "lasso")
  fc <- fit$stages$serial
  expect_identical(fc$method, "F-C (serial)")
  expect_equal(fc$rho, fit$rho)
  expect_equal(fit$n_eff, 250L - fit$ar_order)
  expect_length(residuals(fit, type = "innovations"), fit$n_eff)
  expect_output(print(summary(fit)), "Skedastic fit")
  expect_output(print(summary(fit)), "AR order")
  expect_output(print(fit), "F-CH")
  for (ct in c("HC0", "HC3")) {
    expect_true(all(is.finite(vcov(fgls(y ~ x, data = comb, seed = 0, vcov = ct)))))
    expect_equal(vcov(fit, type = ct), vcov(fgls(y ~ x, data = comb, seed = 0, vcov = ct)))
  }
  expect_error(fgls(y ~ x, data = comb, vcov = "naive"), class = "fgls_error_input")
  expect_error(vcov(fc, type = "HC0"), class = "fgls_error_input")
})

test_that("fgls() dispatches to the right estimator", {
  expect_identical(fgls(y ~ x, data = comb, hetero = "none")$method, "FGLS (serial)")
  expect_identical(fgls(y ~ x, data = comb, serial = "none", seed = 0)$method, "FGLS (hetero)")
  expect_identical(fgls(y ~ x, data = comb, serial = "none", hetero = "none")$method, "OLS")
  expect_equal(fgls(y ~ x, data = comb, serial = "fixed", ar_order = 2, hetero = "none")$ar_order, 2L)
  k <- fgls(y ~ x, data = comb, hetero = "known", variance = attr(comb, "info")$v)
  expect_identical(k$skedastic$method, "known")
  expect_error(fgls(y ~ x, data = comb, serial = "arma"), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = comb, hetero = "kernel"), class = "fgls_error_input")
})

test_that("known variances in the combined estimator use the effective sample", {
  v <- attr(comb, "info")$v
  fit <- fgls(y ~ x, data = comb, hetero = "known", variance = v)
  k <- fit$ar_order
  X <- model.matrix(fit)
  Xf <- fgls:::ar_filter(X, fit$rho)
  yf <- fgls:::ar_filter(comb$y, fit$rho)
  vk <- v[-seq_len(k)]
  b <- fgls:::lsq(Xf / sqrt(vk), yf / sqrt(vk))
  expect_equal(unname(coef(fit)), unname(b), tolerance = 1e-10)
})

test_that("wild multipliers and the sieve recursion", {
  set.seed(2)
  r <- fgls:::wild_multipliers(2e5, "rademacher")
  expect_setequal(unique(r), c(-1, 1))
  m <- fgls:::wild_multipliers(2e5, "mammen")
  expect_lt(abs(mean(m)), 0.01)
  expect_lt(abs(mean(m^2) - 1), 0.01)
  expect_lt(abs(mean(m^3) - 1), 0.03)
  out <- fgls:::sieve_sample(numeric(6), c(1, 2, 0, 0, 0, 0), numeric(4), c(0.5, 0.25), "iid", "rademacher")
  expect_equal(out, c(1, 2, 1.25, 1.125, 0.875, 0.71875))
  w <- fgls:::wild_sample(rep(1, 4), c(1, -1, 2, 0), "rademacher")
  expect_equal(abs(w - 1), c(1, 1, 2, 0))
})

test_that("the bootstrap is reproducible and does not depend on cores", {
  refit <- function(y) c(mean(y), stats::sd(y))
  draw <- function() stats::rnorm(50)
  a <- fgls:::bootstrap_covariance(refit, draw, 40, seed = 5)
  b <- fgls:::bootstrap_covariance(refit, draw, 40, seed = 5, cores = 2L)
  expect_equal(a$cov, b$cov)
  expect_equal(a$draws, b$draws)
})

test_that("bootstrap for the combined estimator", {
  fb <- fgls(y ~ x, data = comb, seed = 0, vcov = "bootstrap", cv_folds = 5,
    control = fgls_control(n_boot = 20))
  expect_true(all(is.finite(sqrt(diag(vcov(fb))))))
  expect_false(is.null(fb$stages$serial))
  expect_match(fb$vcov_type, sprintf("sieve AR\\(%d\\)", fb$ar_order))
})
