het <- simulate_hetero(300, skedastic = "x2", seed = 11)
paper_z <- function(x) {
  Z <- cbind(x, log(x)^2, x^2, cos(x), cos(2 * x))
  colnames(Z) <- c("x", "log(x)^2", "x^2", "cos(x)", "cos(2x)")
  Z
}

test_that("the dictionary transforms regressors and passes dummies through", {
  X <- cbind(const = 1, x = seq(1, 4, length.out = 10), d = rep(0:1, 5))
  expect_equal(
    colnames(make_dictionary(X)),
    c("x", "x^2", "log(x)", "log(x)^2", "cos(x)", "cos(2x)", "d")
  )
  expect_equal(colnames(make_dictionary(X, transforms = "linear")), c("x", "d"))
  neg <- cbind(z = seq(-2, 2, length.out = 10))
  expect_false(any(grepl("log", colnames(make_dictionary(neg)))))
  expect_equal(ncol(make_dictionary(cbind(1, 1))), 0L)
  expect_equal(colnames(make_dictionary(X, skip = "d")), colnames(make_dictionary(X))[1:6])
  dup <- cbind(a = 1:10 + 0.5, b = 2 * (1:10 + 0.5))
  expect_equal(colnames(make_dictionary(dup, transforms = "linear")), "a")
  v <- seq(1, 2, length.out = 8)
  expect_equal(ncol(make_dictionary(v)), 6L)
  expect_error(make_dictionary(letters), class = "fgls_error_input")
  expect_error(fgls:::dictionary_impl(X, colnames(X), "sqrt"), class = "fgls_error_input")
  big <- matrix(stats::runif(40, 1, 4), 4, 10)
  expect_warning(make_dictionary(big), class = "fgls_warning_sample_size")
})

test_that("the LARS homotopy and glmnet give the exact Lasso solution", {
  set.seed(1)
  n <- 120
  Z <- matrix(rnorm(n * 6), n)
  Z[, 2] <- Z[, 1] + 0.1 * rnorm(n)
  r <- drop(Z[, 1:3] %*% c(1, -1, 0.5)) + rnorm(n)
  Zc <- sweep(Z, 2, colMeans(Z))
  G <- crossprod(Zc)
  cv <- drop(crossprod(Zc, r - mean(r)))
  lam <- max(abs(cv)) / n * 10^seq(0, -4, length.out = 50)
  B1 <- fgls:::lars_lasso_coefs(G, cv, n, lam)
  B2 <- fgls:::lasso_path(Z, r, lam)$coef
  expect_equal(B1, B2, tolerance = 1e-10)
  expect_true(all(fgls:::kkt_gaps(G, cv, n, lam, B1) < 1e-12))
  # against glmnet without refinement (coordinate-descent tolerance)
  raw <- fgls:::lasso_path(Z, r, lam, polish = FALSE)$coef
  expect_lt(max(abs(raw - B1)), 1e-3)
  # one column and a zero column
  one <- fgls:::lasso_path(Z[, 1, drop = FALSE], r, lam[1:5])$coef
  expect_equal(dim(one), c(1L, 5L))
  z0 <- fgls:::lasso_path(cbind(Z[, 1], 0), r, lam[1:5])$coef
  expect_true(all(z0[2, ] == 0))
  expect_equal(fgls:::lars_lasso_coefs(G, cv * 0, n, lam), matrix(0, 6, 50))
})

test_that("feature-sign search reaches the exact solution from any start", {
  set.seed(2)
  n <- 80
  Z <- matrix(rnorm(n * 5), n)
  r <- Z[, 1] - Z[, 2] + rnorm(n)
  Zc <- sweep(Z, 2, colMeans(Z))
  G <- crossprod(Zc)
  cv <- drop(crossprod(Zc, r - mean(r)))
  lam <- 0.05
  ref <- fgls:::lars_lasso_coefs(G, cv, n, lam)
  for (b0 in list(numeric(5), c(5, -5, 5, -5, 5), drop(ref) * 1.1)) {
    b <- fgls:::feature_sign(G, cv, n, lam, b0)
    expect_equal(b, drop(ref), tolerance = 1e-10)
  }
})

test_that("psi = 0 is the plain Lasso and the adaptive Lasso matches a reweighted brute force", {
  x <- het$x
  Z <- paper_z(x)
  r <- log(pmax((het$y - stats::fitted(stats::lm(y ~ x, het)))^2, 0.01))
  fit0 <- fgls:::adaptive_lasso(r, Z, psi_grid = 0, seed = 1)
  st <- fgls:::standardize_cols(Z)
  lam <- fit0$tuning$lambda
  g <- glmnet::glmnet(st$Zs, r, lambda = lam, standardize = FALSE, thresh = 1e-14, maxit = 1e7)
  expect_equal(unname(fit0$coef[st$keep] * st$sd[st$keep]), as.numeric(g$beta), tolerance = 1e-5)
  fit1 <- fgls:::adaptive_lasso(r, Z, psi_grid = 1, seed = 1)
  w <- abs(fit1$tuning$ridge_coef)
  Zt <- sweep(st$Zs, 2, w, "*")
  b <- fgls:::lars_lasso_coefs(
    crossprod(sweep(Zt, 2, colMeans(Zt))),
    drop(crossprod(sweep(Zt, 2, colMeans(Zt)), r - mean(r))), length(r), fit1$tuning$lambda
  )
  expect_equal(unname(fit1$coef * st$sd), unname(drop(b) * w), tolerance = 1e-10)
})

test_that("df counts the intercept, and a constant response gives a constant fit", {
  x <- het$x
  r <- log(pmax(stats::rnorm(300)^2 * x^2, 0.01))
  fit <- fgls:::adaptive_lasso(r, paper_z(x), seed = 2)
  expect_equal(fit$df, sum(fit$coef != 0) + 1)
  cst <- fgls:::adaptive_lasso(rep(1, 300), paper_z(x), seed = 2)
  expect_equal(cst$df, 1)
  expect_equal(cst$fitted, rep(1, 300))
  none <- fgls:::adaptive_lasso(r, matrix(0, 300, 0))
  expect_equal(none$df, 1)
  allc <- fgls:::adaptive_lasso(r, matrix(1, 300, 2, dimnames = list(NULL, c("a", "b"))))
  expect_equal(allc$tuning$reason, "all candidate covariates are constant")
})

test_that("ridge handles more covariates than observations", {
  set.seed(3)
  Z <- matrix(rnorm(40 * 80), 40)
  r <- Z[, 1] + rnorm(40)
  st <- fgls:::standardize_cols(Z)
  rid <- fgls:::ridge_cv(st$Zs, r, fgls:::make_folds(40, 5, seed = 1))
  expect_length(rid$coef, 80L)
  lam <- rid$lambda
  Zc <- sweep(st$Zs, 2, colMeans(st$Zs))
  direct <- solve(crossprod(Zc) + lam * diag(80), crossprod(Zc, r - mean(r)))
  expect_equal(rid$coef, drop(direct), tolerance = 1e-8)
  expect_warning(fit <- fgls:::adaptive_lasso(r, Z, cv_folds = 5, seed = 1), NA)
  expect_equal(length(fit$coef), 80L)
})

test_that("folds are balanced, reproducible and validated", {
  f <- fgls:::make_folds(103, 10, seed = 5)
  expect_true(all(table(f) %in% c(10L, 11L)))
  expect_identical(f, fgls:::make_folds(103, 10, seed = 5))
  b <- fgls:::make_folds(10, 5, "blocked")
  expect_equal(b, rep(1:5, each = 2))
  expect_error(fgls:::make_folds(10, 1), class = "fgls_error_input")
  expect_error(fgls:::make_folds(3, 5), class = "fgls_error_input")
  expect_error(fgls:::check_foldid(1:3, 4), class = "fgls_error_input")
  expect_error(fgls:::check_foldid(rep(1, 4), 4), class = "fgls_error_input")
  expect_equal(fgls:::check_foldid(c(7, 7, 9, 9), 4), c(1L, 1L, 2L, 2L))
  expect_warning(fgls:::adaptive_lasso(rnorm(15), matrix(rnorm(30), 15), cv_folds = 10),
    class = "fgls_warning_sample_size"
  )
})

test_that("the 1se rule is at least as sparse as the minimum rule", {
  f_min <- lasso_fgls(y ~ x, data = het, seed = 3)
  f_1se <- lasso_fgls(y ~ x, data = het, seed = 3, lambda_rule = "1se")
  expect_gte(f_1se$skedastic$lambda, f_min$skedastic$lambda)
  expect_lte(length(f_1se$skedastic$selected), length(f_min$skedastic$selected) + 1L)
})

test_that("the two Lasso engines agree end to end", {
  a <- lasso_fgls(y ~ x, data = het, seed = 4)
  b <- lasso_fgls(y ~ x, data = het, seed = 4, control = fgls_control(lasso_engine = "glmnet"))
  expect_equal(coef(a), coef(b), tolerance = 1e-8)
  expect_equal(vcov(a), vcov(b), tolerance = 1e-8)
  expect_identical(a$skedastic$selected, b$skedastic$selected)
  c2 <- lasso_fgls(y ~ x, data = het, seed = 4, control = fgls_control(cv_weights = "full"))
  expect_identical(c2$skedastic$tuning$cv_weights, "full")
  c3 <- lasso_fgls(y ~ x, data = het, seed = 4, control = fgls_control(ridge_scale = "original"))
  expect_true(is.finite(c3$skedastic$lambda))
})

test_that("WLS-S1/S2 and SVR skedastic fits", {
  w1 <- lasso_fgls(y ~ x, data = het, hetero = "wls_s1", z = "x")
  expect_equal(w1$skedastic$df_hat, 2)
  expect_match(w1$notes, "eq. \\(7\\)", all = FALSE)
  w2 <- lasso_fgls(y ~ x, data = het, hetero = "wls_s2")
  expect_equal(w2$skedastic$df_hat, 1 + length(w2$skedastic$selected))
  df <- het
  df$xn <- het$x - 2.5
  expect_error(lasso_fgls(y ~ x, data = df, hetero = "wls_s1", z = ~xn), class = "fgls_error_input")
  skip_if_not_installed("e1071")
  sv <- lasso_fgls(y ~ x, data = het, hetero = "svr", seed = 1)
  expect_gte(sv$skedastic$df_hat, 1)
  expect_null(sv$skedastic$coefficients)
  expect_true(all(is.finite(vcov(sv))))
})

test_that("variances are normalised, clipped and checked", {
  v <- fgls:::variances_from_log(c(0, 1, 2, 3))
  expect_equal(mean(v), 1)
  vc <- fgls:::variances_from_log(1:100, clip = c(0.1, 0.9))
  expect_lt(max(vc) / min(vc), exp(99))
  expect_warning(fgls:::variances_from_log(c(0, 20)), class = "fgls_warning_weight")
  expect_error(fgls:::log_squared_residuals(1, -1), class = "fgls_error_input")
  expect_error(fgls:::known_variances(c(1, -1)), class = "fgls_error_input")
})

test_that("hetero = 'none' is OLS with HC3 and known variances give closed-form GLS", {
  ols <- lasso_fgls(y ~ x, data = het, hetero = "none")
  lmf <- stats::lm(y ~ x, data = het)
  expect_equal(coef(ols), coef(lmf))
  X <- stats::model.matrix(lmf)
  h <- stats::hatvalues(lmf)
  u <- residuals(lmf)
  A <- solve(crossprod(X))
  hc3 <- A %*% crossprod(X * (u / (1 - h))) %*% A
  expect_equal(unname(vcov(ols)), unname(hc3), tolerance = 1e-10)
  v <- attr(het, "info")$v
  gls <- lasso_fgls(y ~ x, data = het, hetero = "known", variance = v)
  W <- diag(1 / v)
  b <- solve(t(X) %*% W %*% X, t(X) %*% W %*% het$y)
  expect_equal(unname(coef(gls)), unname(drop(b)), tolerance = 1e-10)
  gf <- lasso_fgls(y ~ x, data = het, hetero = "known", variance = function(X) X[, 2]^2)
  expect_equal(coef(gf), coef(gls), tolerance = 1e-10)
  expect_error(lasso_fgls(y ~ x, data = het, hetero = "known", variance = 1:3), class = "fgls_error_input")
  expect_error(lasso_fgls(y ~ x, data = het, hetero = "known", variance = rep(-1, 300)),
    class = "fgls_error_input"
  )
  expect_identical(gls$method, "GLS (known variances)")
})

test_that("eq. (8) and eq. (7) variances by hand", {
  fit <- lasso_fgls(y ~ x, data = het, seed = 5)
  X <- model.matrix(fit)
  v <- fit$skedastic$variance
  u <- residuals(fit)
  k <- ncol(X)
  n <- nrow(X)
  hw <- fgls:::hat_diag(X / sqrt(v))
  ho <- fgls:::hat_diag(X)
  sig <- u^2 / v^2 * (1 / (1 - hw)^2 + 4 * ho / k * fit$skedastic$df_hat)
  A <- solve(crossprod(X / v, X) / n)
  V <- A %*% (crossprod(X * sig, X) / n) %*% A / n
  expect_equal(unname(vcov(fit)), unname(V), tolerance = 1e-10)
  w1 <- lasso_fgls(y ~ x, data = het, hetero = "wls_s1", z = "x")
  v1 <- w1$skedastic$variance
  u0 <- residuals(stats::lm(y ~ x, het))
  sig7 <- u0^2 / (1 - ho)^2 / v1^2
  A1 <- solve(crossprod(X / v1, X) / n)
  expect_equal(unname(vcov(w1)), unname(A1 %*% (crossprod(X * sig7, X) / n) %*% A1 / n),
    tolerance = 1e-10
  )
})

test_that("the variance is invariant to the normalisation of the weights", {
  v <- attr(het, "info")$v
  a <- lasso_fgls(y ~ x, data = het, hetero = "known", variance = v)
  b <- lasso_fgls(y ~ x, data = het, hetero = "known", variance = 7 * v)
  expect_equal(vcov(a), vcov(b), tolerance = 1e-12)
})

test_that("other covariance types are available", {
  fit <- lasso_fgls(y ~ x, data = het, seed = 6)
  for (ct in c("nonrobust", "HC0", "HC1", "HC2", "HC3")) {
    V <- vcov(fit, type = ct)
    expect_true(all(eigen(V)$values >= 0), info = ct)
  }
  f2 <- lasso_fgls(y ~ x, data = het, seed = 6, vcov = "HC0")
  expect_equal(vcov(f2), vcov(fit, type = "HC0"))
  expect_error(vcov(fit, type = "naive"), class = "fgls_error_input")
})

test_that("lasso_fgls end to end beats OLS under strong heteroskedasticity", {
  fit <- lasso_fgls(y ~ x, data = het, seed = 7)
  ols <- lasso_fgls(y ~ x, data = het, hetero = "none")
  expect_lt(sqrt(vcov(fit)[2, 2]), sqrt(vcov(ols)[2, 2]))
  expect_true(any(c("x", "x^2", "log(x)", "log(x)^2") %in% fit$skedastic$selected))
  expect_equal(mean(fit$skedastic$variance), 1)
  expect_equal(fit$weights, 1 / fit$skedastic$variance)
  expect_output(print(summary(fit)), "Skedastic fit")
  expect_equal(fit$skedastic$lambda_paper, 300 * fit$skedastic$lambda)
})

test_that("user-supplied z: names, formula and matrix", {
  df <- het
  df$w <- stats::runif(300, 1, 4)
  f1 <- lasso_fgls(y ~ x, data = df, z = "x", seed = 1)
  expect_equal(names(f1$skedastic$coefficients), c("(Intercept)", "x"))
  f2 <- lasso_fgls(y ~ x, data = df, z = ~ x + I(x^2) + w, seed = 1)
  expect_equal(names(f2$skedastic$coefficients), c("(Intercept)", "x", "I(x^2)", "w"))
  f3 <- lasso_fgls(y ~ x, data = df, z = cbind(a = df$x, b = df$w), seed = 1)
  expect_equal(names(f3$skedastic$coefficients), c("(Intercept)", "a", "b"))
  f4 <- lasso_fgls(y ~ x, data = df, z = data.frame(df$x, df$w), seed = 1)
  expect_length(f4$skedastic$coefficients, 3L)
  expect_error(lasso_fgls(y ~ x, data = df, z = cbind(1:3)), class = "fgls_error_input")
  expect_error(lasso_fgls(y ~ x, data = df, z = list(1)), class = "fgls_error_input")
  expect_error(lasso_fgls(y ~ x, data = df, z = "nope"), class = "fgls_error_input")
  bad <- cbind(a = df$x)
  bad[3] <- Inf
  expect_error(lasso_fgls(y ~ x, data = df, z = bad), class = "fgls_error_input")
})

test_that("missing values are dropped jointly with z", {
  df <- het
  df$w <- stats::runif(300, 1, 4)
  df$w[c(5, 50)] <- NA
  expect_error(lasso_fgls(y ~ x, data = df, z = ~w), class = "fgls_error_missing")
  fit <- lasso_fgls(y ~ x, data = df, z = ~w, na.action = stats::na.omit, seed = 1)
  expect_equal(fit$nobs, 298L)
})

test_that("delta_scale and weight_clip", {
  a <- lasso_fgls(y ~ x, data = het, seed = 8, control = fgls_control(delta_scale = TRUE))
  expect_gt(a$skedastic$delta, 0.1)
  b <- lasso_fgls(y ~ x, data = het, seed = 8, weight_clip = c(0.05, 0.95))
  v <- b$skedastic$variance
  expect_true(all(is.finite(v)))
})

test_that("scale invariance with a scaled delta and invariance to row order", {
  a <- lasso_fgls(y ~ x, data = het, seed = 9)
  het5 <- transform(het, y = 5 * y)
  b <- lasso_fgls(y ~ x, data = het5, seed = 9, delta = 0.5)
  expect_equal(unname(coef(b)), unname(5 * coef(a)), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(b)))), unname(5 * sqrt(diag(vcov(a)))), tolerance = 1e-8)
  o <- sample(300)
  f_perm <- lasso_fgls(y ~ x, data = het[o, ], hetero = "wls_s2")
  expect_equal(coef(f_perm), coef(lasso_fgls(y ~ x, data = het, hetero = "wls_s2")), tolerance = 1e-10)
})

test_that("wild bootstrap for heteroskedasticity", {
  fit <- lasso_fgls(y ~ x, data = het, hetero = "wls_s2", vcov = "bootstrap", seed = 1,
    control = fgls_control(n_boot = 99))
  ref <- lasso_fgls(y ~ x, data = het, hetero = "wls_s2")
  ratio <- sqrt(diag(vcov(fit))) / sqrt(diag(vcov(ref)))
  expect_true(all(ratio > 0.6 & ratio < 1.6))
  expect_match(fit$vcov_type, "wild")
})
