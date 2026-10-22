ar_data <- function(n, ar, seed, gamma = 0, rho_x = 0.5, burn = 100L) {
  simulate_serial(n, ar = ar, rho_x = rho_x, gamma = gamma, burn = burn, seed = seed)
}

test_that("AR(1) quasi-differencing equals explicit D-matrix GLS", {
  set.seed(1)
  n <- 30
  X <- cbind(1, rnorm(n))
  y <- drop(X %*% c(1, 2)) + rnorm(n)
  rho <- 0.6
  D <- matrix(0, n - 1, n)
  for (i in seq_len(n - 1)) D[i, c(i, i + 1)] <- c(-rho, 1)
  b_gls <- solve(crossprod(D %*% X), crossprod(D %*% X, D %*% y))
  fr <- fgls:::filtered_regression(y, X, rho)
  expect_equal(unname(fr$beta), drop(b_gls), tolerance = 1e-10)
})

test_that("known AR(2) filter equals closed-form GLS conditional on initial values", {
  set.seed(2)
  n <- 40
  X <- cbind(1, rnorm(n), rnorm(n))
  y <- drop(X %*% c(1, -1, 0.5)) + rnorm(n)
  rho <- c(0.5, -0.3)
  D <- matrix(0, n - 2, n)
  for (i in seq_len(n - 2)) D[i, i:(i + 2)] <- c(-rho[2], -rho[1], 1)
  Oinv <- crossprod(D)
  b <- solve(t(X) %*% Oinv %*% X, t(X) %*% Oinv %*% y)
  expect_equal(unname(fgls:::filtered_regression(y, X, rho)$beta), drop(b), tolerance = 1e-10)
  expect_equal(fgls:::filtered_regression(y, X, rho)$Xf[, 1], rep(1 - sum(rho), n - 2))
})

test_that("Durbin designs are nested and the lag-selection RSS matches direct fits", {
  sim <- ar_data(200, c(0.6, -0.2), seed = 3)
  X <- cbind("(Intercept)" = 1, x = sim$x)
  d <- fgls:::durbin_design(sim$y, X, colnames(X), 4, start = 4)
  expect_true("L1.(Intercept)" %in% d$dropped)
  rss <- fgls:::nested_rss(d, 4)
  direct <- vapply(0:4, function(k) {
    sum(stats::lm.fit(fgls:::design_prefix(d, k), d$y)$residuals^2)
  }, numeric(1))
  expect_equal(rss, direct, tolerance = 1e-10)
})

test_that("deterministic regressors neither crash the Durbin step nor change rho", {
  set.seed(6)
  n <- 240
  tt <- seq_len(n)
  quarter <- (seq_len(n) - 1) %% 4
  Dq <- sapply(1:3, function(j) as.numeric(quarter == j))
  x <- rnorm(n)
  X <- cbind(1, tt, Dq, x)
  colnames(X) <- c("(Intercept)", "trend", "q2", "q3", "q4", "x")
  u <- fgls:::arma_errors(rnorm(n), 0.6)
  y <- drop(X %*% c(1, 0.01, 0.5, -0.5, 0.2, 1)) + u
  d <- fgls:::durbin_design(y, X, colnames(X), 2)
  expect_true(all(c("L1.(Intercept)", "L1.trend", "L1.q2", "L1.q3", "L1.q4") %in% d$dropped))
  rho <- fgls:::fit_durbin(d)$coef[fgls:::design_ylag_columns(d)]
  full <- cbind(X[3:n, ], y[2:(n - 1)], X[2:(n - 1), ], y[1:(n - 2)], X[1:(n - 2), ])
  cf <- qr.coef(qr(full), y[3:n])
  expect_equal(unname(rho), unname(cf[c(7, 14)]), tolerance = 1e-6)
  df <- data.frame(y = y, X[, -1])
  fit <- fgls(y ~ trend + q2 + q3 + q4 + x, data = df, hetero = "none", ar_order = 2)
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
})

eq29_by_hand <- function(y, X, k) {
  n_all <- nrow(X)
  rows <- (k + 1):n_all
  Zy <- sapply(seq_len(k), function(j) y[rows - j])
  W <- cbind(X[rows, ], do.call(cbind, lapply(seq_len(k), function(j) X[rows - j, -1, drop = FALSE])))
  cf <- qr.coef(qr(cbind(Zy, W)), y[rows])
  rho <- cf[seq_len(k)]
  Xf <- X[rows, ]
  yf <- y[rows]
  for (j in seq_len(k)) {
    Xf <- Xf - rho[j] * X[rows - j, ]
    yf <- yf - rho[j] * y[rows - j]
  }
  beta <- solve(crossprod(Xf), crossprod(Xf, yf))
  e <- yf - Xf %*% beta
  u <- y - X %*% beta
  n <- n_all - k
  s2 <- sum(e^2) / n
  Q <- crossprod(Xf) / n
  G <- sapply(seq_len(k), function(j) crossprod(Xf, u[rows - j]) / n)
  G <- matrix(G, ncol(X), k)
  psi <- Zy - W %*% qr.coef(qr(W), Zy)
  Sig <- crossprod(psi) / n
  Qi <- solve(Q)
  avar <- s2 * (Qi + Qi %*% G %*% solve(Sig) %*% t(G) %*% Qi)
  list(beta = drop(beta), cov = avar / n, naive = s2 * Qi / n)
}

test_that("robust variance matches an independent implementation of eqs. (28)-(29)", {
  sim <- ar_data(300, c(0.5, 0.2), seed = 7, gamma = 0.5)
  X <- cbind("(Intercept)" = 1, x = sim$x)
  for (k in 1:3) {
    est <- fgls:::estimate_serial(sim$y, X, ar_order = k)
    ref <- eq29_by_hand(sim$y, X, k)
    expect_equal(unname(est$beta), unname(ref$beta), tolerance = 1e-8)
    expect_equal(unname(est$cov), unname(ref$cov), tolerance = 1e-8)
    naive <- fgls:::serial_variance(est$filtered$Xf, est$filtered$ehat, est$filtered$uhat,
      est$design,
      cov_type = "naive"
    )
    expect_equal(unname(naive$cov), unname(ref$naive), tolerance = 1e-8)
  }
})

test_that("the influence form equals eq. (29) in the standard case", {
  sim <- ar_data(250, 0.7, seed = 8, gamma = 0.5)
  X <- cbind("(Intercept)" = 1, x = sim$x)
  est <- fgls:::estimate_serial(sim$y, X, ar_order = 2)
  fr <- est$filtered
  k <- 2
  p <- ncol(est$design$Z)
  J <- matrix(0, k, p)
  J[cbind(1:k, fgls:::design_ylag_columns(est$design))] <- 1
  v_inf <- fgls:::serial_variance(fr$Xf, fr$ehat, fr$uhat, est$design, J = J)
  expect_equal(v_inf$cov, est$cov, tolerance = 1e-8)
  C1 <- fgls:::first_stage_influence(est$design)
  C2 <- fgls:::first_stage_influence(est$design, J)
  expect_equal(C1, C2, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("robust variance is PSD and the correction is non-negative", {
  sim <- ar_data(200, 0.8, seed = 9, gamma = 1)
  fit <- serial_fgls(y ~ x, data = sim)
  ev <- eigen(vcov(fit), symmetric = TRUE)$values
  expect_true(all(ev >= 0))
  expect_true(all(diag(fit$serial$nonexogeneity_by_coef) >= 0))
  expect_true(all(diag(vcov(fit)) >= diag(vcov(fit, type = "naive")) - 1e-12))
  expect_gte(fit$nonexogeneity_share, 0)
})

test_that("k* = 0 reduces FGLS to OLS", {
  set.seed(10)
  df <- data.frame(x = rnorm(200))
  df$y <- 1 + df$x + rnorm(200)
  fit <- fgls(y ~ x, data = df, hetero = "none", ar_order = 0)
  ols <- stats::lm(y ~ x, data = df)
  expect_equal(coef(fit), coef(ols), tolerance = 1e-10)
  expect_equal(fit$n_eff, 200L)
  expect_equal(unname(vcov(fit)), unname(vcov(ols) * 198 / 200), tolerance = 1e-10)
  expect_match(fit$notes, "k\\* = 0", all = FALSE)
})

test_that("HC variants are close to robust under homoskedasticity", {
  sim <- ar_data(400, 0.5, seed = 11)
  fit <- serial_fgls(y ~ x, data = sim)
  for (ct in c("HC0", "HC1", "HC2", "HC3")) {
    r <- sqrt(diag(vcov(fit, type = ct))) / sqrt(diag(vcov(fit)))
    expect_true(all(abs(r - 1) < 0.2), info = ct)
  }
})

test_that("BIC picks k = 0 for white noise and the true order in large samples", {
  ks <- vapply(1:40, function(s) {
    sim <- ar_data(200, numeric(0), seed = 100 + s)
    X <- cbind(1, sim$x)
    fgls:::estimate_serial(sim$y, X, c("c", "x"))$order
  }, integer(1))
  expect_gte(mean(ks == 0L), 0.8)
  sim <- ar_data(2000, c(0.5, -0.3), seed = 12)
  expect_equal(serial_fgls(y ~ x, data = sim)$ar_order, 2L)
})

test_that("AIC, min_lags and a fixed order", {
  sim <- ar_data(300, 0.5, seed = 13)
  fa <- serial_fgls(y ~ x, data = sim, ic = "aic")
  expect_identical(fa$serial$method, "aic")
  fm <- serial_fgls(y ~ x, data = sim, control = fgls_control(min_lags = 4))
  expect_gte(fm$ar_order, 4L)
  fadj <- serial_fgls(y ~ x, data = sim, control = fgls_control(min_lags = "adjusted"))
  expect_gte(fadj$ar_order, fgls:::adjusted_min_lags(300))
  ff <- serial_fgls(y ~ x, data = sim, ar_order = 2)
  expect_identical(ff$serial$method, "fixed")
  expect_null(ff$bic_path)
  fp <- serial_fgls(y ~ x, data = sim, control = fgls_control(ic_penalty = "params"))
  expect_identical(fp$ar_order, serial_fgls(y ~ x, data = sim)$ar_order)
  expect_error(serial_fgls(y ~ x, data = sim, ar_order = 200), class = "fgls_error_insufficient_data")
  expect_error(serial_fgls(y ~ x, data = sim, ar_order = -1), class = "fgls_error_input")
  expect_error(serial_fgls(y ~ x, data = sim, kmax = "huge"), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, serial = "fixed"), class = "fgls_error_input")
  expect_equal(serial_fgls(y ~ x, data = sim, kmax = "theory")$serial$kmax, fgls:::theory_kmax(300))
})

test_that("the kmax cap warns for short samples", {
  sim <- ar_data(40, 0.5, seed = 14)
  expect_warning(serial_fgls(y ~ x, data = sim, kmax = 30), class = "fgls_warning_sample_size")
  sim2 <- ar_data(5, 0.5, seed = 14)
  expect_error(
    suppressWarnings(fgls:::resolve_kmax(NULL, 5, 4, TRUE)),
    class = "fgls_error_insufficient_data"
  )
})

test_that("selecting kmax warns with a hint", {
  set.seed(15)
  n <- 300
  u <- fgls:::arma_errors(rnorm(n), ma = c(0.9, 0.8, 0.7, 0.6, 0.5))
  df <- data.frame(x = rnorm(n))
  df$y <- 1 + df$x + u
  expect_warning(serial_fgls(y ~ x, data = df, kmax = 1), class = "fgls_warning_specification")
})

test_that("scaling invariance", {
  sim <- ar_data(300, 0.6, seed = 15, gamma = 0.5)
  base <- serial_fgls(y ~ x, data = sim)
  s1 <- serial_fgls(y ~ I(10 * x), data = sim)
  expect_equal(unname(coef(s1)[2]), unname(coef(base)[2] / 10), tolerance = 1e-8)
  tv <- function(f) coef(f) / sqrt(diag(vcov(f)))
  expect_equal(unname(tv(s1)), unname(tv(base)), tolerance = 1e-6)
  s2 <- serial_fgls(I(5 * y) ~ x, data = sim)
  expect_identical(s2$ar_order, base$ar_order)
  expect_equal(unname(coef(s2)), unname(5 * coef(base)), tolerance = 1e-8)
  expect_equal(unname(tv(s2)), unname(tv(base)), tolerance = 1e-6)
  s3 <- serial_fgls(I(2 * y + 3 * x) ~ x, data = sim)
  expect_equal(unname(coef(s3)), unname(2 * coef(base) + c(0, 3)), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(s3)))), unname(2 * sqrt(diag(vcov(base)))), tolerance = 1e-6)
})

test_that("a near unit root warns", {
  set.seed(17)
  n <- 400
  df <- data.frame(x = rnorm(n))
  df$y <- 1 + df$x + fgls:::arma_errors(rnorm(n), 0.995)
  expect_warning(serial_fgls(y ~ x, data = df), class = "fgls_warning_nonstationary")
  expect_warning(fgls:::check_ar_roots(1.2), class = "fgls_warning_nonstationary")
  expect_silent(fgls:::check_ar_roots(numeric(0)))
})

test_that("an undeclared lagged dependent variable is detected", {
  sim <- ar_data(200, 0.5, seed = 18)
  df <- data.frame(y = sim$y[-1], x = sim$x[-1], ylag = sim$y[-200])
  expect_error(
    expect_warning(serial_fgls(y ~ x + ylag, data = df), class = "fgls_warning_lagged_dependent"),
    "lagged_dep_vars"
  )
})

test_that("the lagged-dependent-variable path recovers rho", {
  set.seed(19)
  n <- 2000
  w <- rnorm(n + 1)
  u <- fgls:::arma_errors(rnorm(n + 1), 0.5)
  yy <- numeric(n + 1)
  for (t in 2:(n + 1)) yy[t] <- 0.3 * yy[t - 1] + 1 + w[t] + u[t]
  df <- data.frame(y = yy[-1], ylag = yy[-(n + 1)], w = w[-1])
  fit <- serial_fgls(y ~ ylag + w, data = df, lagged_dep_vars = "ylag")
  expect_equal(fit$ar_order, 1L)
  expect_lt(abs(unname(fit$rho) - 0.5), 0.1)
  expect_lt(abs(unname(coef(fit)["ylag"]) - 0.3), 0.1)
  expect_identical(fit$serial$wallis_regressors, "w")
  fm <- serial_fgls(y ~ ylag + w, data = df, lagged_dep_vars = "ylag", wallis_regressor = "w")
  expect_equal(coef(fm), coef(fit))
  expect_error(
    serial_fgls(y ~ ylag + w, data = df, lagged_dep_vars = "ylag", wallis_regressor = "ylag"),
    class = "fgls_error_input"
  )
  expect_error(
    serial_fgls(y ~ ylag + w, data = df, lagged_dep_vars = "ylag", wallis_regressor = "nope"),
    class = "fgls_error_input"
  )
  df2 <- data.frame(y = yy[-1], ylag = yy[-(n + 1)])
  expect_error(
    serial_fgls(y ~ ylag, data = df2, lagged_dep_vars = "ylag", ar_order = 1),
    "stochastic regressor"
  )
})

test_that("an exogeneity declaration adds a note when violated", {
  sim <- ar_data(500, 0.8, seed = 20, gamma = 1, rho_x = 0.8)
  fit <- serial_fgls(y ~ x, data = sim, exog = "exogenous")
  expect_gt(fit$nonexogeneity_share, 0.1)
  expect_match(fit$notes, "exog = \"exogenous\"", all = FALSE, fixed = TRUE)
  fit2 <- serial_fgls(y ~ x, data = sim)
  expect_equal(coef(fit), coef(fit2))
})

test_that("filtered innovations match the quasi-difference", {
  sim <- ar_data(200, 0.6, seed = 21)
  fit <- serial_fgls(y ~ x, data = sim)
  k <- fit$ar_order
  e <- fgls:::ar_filter(residuals(fit), fit$rho)
  expect_equal(unname(residuals(fit, type = "innovations")), unname(e), tolerance = 1e-10)
  expect_length(residuals(fit, type = "innovations"), 200L - k)
})

test_that("GMA recovers an MA(1) coefficient", {
  sim <- simulate_serial(400, ar = numeric(0), ma = 0.6, gamma = 0, burn = 50, seed = 1)
  fit <- serial_fgls(y ~ x, data = sim, serial = "ma1")
  expect_equal(fit$serial$extra$theta, 0.6, tolerance = 0.15)
  expect_identical(fit$n_eff, 400L)
  expect_error(serial_fgls(y ~ x, data = sim, serial = "ma1", vcov = "naive"), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, serial = "ma1"), class = "fgls_error_input")
  expect_warning(fgls:::ma1_mle(c(1, -1, 1, -1, 1, -1, 1, -1)), class = "fgls_warning_nonstationary")
})

test_that("sieve bootstrap gives standard errors close to the analytic ones", {
  sim <- ar_data(300, 0.5, seed = 22)
  fit <- serial_fgls(y ~ x, data = sim)
  fb <- serial_fgls(y ~ x,
    data = sim, vcov = "bootstrap", seed = 1,
    control = fgls_control(n_boot = 99)
  )
  ratio <- sqrt(diag(vcov(fb))) / sqrt(diag(vcov(fit)))
  expect_true(all(ratio > 0.6 & ratio < 1.6))
  expect_match(fb$vcov_type, "sieve AR")
  expect_equal(dim(fb$bootstrap$draws), c(99L, 2L))
  fb2 <- serial_fgls(y ~ x,
    data = sim, vcov = "bootstrap", seed = 1,
    control = fgls_control(n_boot = 99)
  )
  expect_identical(vcov(fb), vcov(fb2))
  fi <- serial_fgls(y ~ x,
    data = sim, vcov = "bootstrap", seed = 1,
    control = fgls_control(n_boot = 20, bootstrap = "iid", wild_dist = "mammen")
  )
  expect_match(fi$vcov_type, "iid")
  vb <- vcov(fit, type = "bootstrap")
  expect_equal(dim(vb), c(2L, 2L))
})
