# Regression tests for the pre-submission audit fixes in the serial-correlation
# internals, the benchmark estimators and the bootstrap (audit findings 2, 5,
# 6, 7, 12, 18, 26, 27, 31, 32).

nonexog_data <- function(n = 500, seed = 2) {
  simulate_serial(n, ar = 0.8, rho_x = 0.8, gamma = 1, burn = 100, seed = seed)
}

monthly_data <- function(n = 120, seed = 5) {
  sim <- simulate_serial(n, ar = c(0.5, 0.3), seed = seed)
  sim$tr <- seq_len(n)
  sim$m <- factor(rep(1:12, length.out = n))
  sim
}

# ------------------------------------------------------------ finding 5: share
test_that("the non-exogeneity share does not depend on the variance type", {
  sim <- nonexog_data()
  rob <- serial_fgls(y ~ x, data = sim)
  expect_gt(rob$nonexogeneity_share, 0.05)
  for (ct in c("naive", "HC0", "HC1", "HC2", "HC3")) {
    f <- serial_fgls(y ~ x, data = sim, vcov = ct)
    expect_equal(f$nonexogeneity_share, rob$nonexogeneity_share, tolerance = 1e-12)
    expect_equal(f$serial$nonexogeneity_by_coef, rob$serial$nonexogeneity_by_coef,
      tolerance = 1e-12
    )
    expect_equal(nonexogeneity_share(f), rob$nonexogeneity_share, tolerance = 1e-12)
  }
  hc3 <- serial_fgls(y ~ x, data = sim, vcov = "HC3")
  expect_identical(
    summary(hc3)$info[["Non-exogeneity share"]],
    sprintf("%.3f", rob$nonexogeneity_share)
  )
  # the variance itself still follows the requested type
  expect_false(isTRUE(all.equal(vcov(hc3), vcov(rob))))
})

test_that("the share of a non-robust variance equals the robust share with a Wallis Jacobian", {
  sim <- simulate_serial(400, seed = 3)
  df <- data.frame(y = sim$y[-1], ylag = sim$y[-400], x = sim$x[-1])
  fits <- lapply(c("robust", "naive", "HC3"), function(ct) {
    suppressWarnings(serial_fgls(y ~ ylag + x,
      data = df, lagged_dep_vars = "ylag", ar_order = 1, vcov = ct
    ))
  })
  expect_false(is.null(fits[[1]]$internal$J))
  expect_equal(fits[[2]]$nonexogeneity_share, fits[[1]]$nonexogeneity_share, tolerance = 1e-12)
  expect_equal(fits[[3]]$nonexogeneity_share, fits[[1]]$nonexogeneity_share, tolerance = 1e-12)
})

test_that("golden: the share of every variance type is the pyfgls robust share", {
  skip_if_no_golden()
  for (id in c("s1", "s2", "s5", "s6")) {
    g <- read_golden(id)
    r <- g$res
    X <- golden_X(g)
    des <- fgls:::durbin_design(g$data$y, X, colnames(X), r$kstar)
    fr <- fgls:::filtered_regression(g$data$y, X, as.numeric(r$rho_hat))
    J <- if (!is.null(r$wallis)) as_mat(r$wallis$J, r$kstar, length(r$durbin_coef))
    for (ct in c("naive", "HC0", "HC3")) {
      v <- fgls:::serial_variance(fr$Xf, fr$ehat, fr$uhat, des, J = J, cov_type = ct)
      # pyfgls reports 0 for the non-robust types; fgls reports the robust share
      expect_abs(v$nonexogeneity_share, r$nonexogeneity_share, 1e-8, paste(id, ct, "share"))
    }
  }
})

# ------------------------------------------------------- finding 6: ma1 bootstrap
test_that("the bootstrap is not available for the GMA estimator", {
  pred <- simulate_predictive(300, theta = 0.8, seed = 4)
  expect_identical(fgls:::allowed_vcov("serial", "ma1"), "robust")
  expect_error(
    serial_fgls(y ~ x, data = pred, serial = "ma1", vcov = "bootstrap"),
    class = "fgls_error_input"
  )
  fit <- serial_fgls(y ~ x, data = pred, serial = "ma1")
  expect_error(vcov(fit, type = "bootstrap"), class = "fgls_error_input")
  expect_error(
    fgls:::fgls_bootstrap(fit$spec, fit$internal, fit, fit$y),
    "MA\\(1\\) structure", class = "fgls_error_input"
  )
  expect_true("bootstrap" %in% fgls:::allowed_vcov("serial", "auto"))
})

# -------------------------------------------------------- finding 7: kmax cap
test_that("the kmax cap counts only the Durbin regressors that are kept", {
  sim <- monthly_data()
  fit <- suppressWarnings(serial_fgls(y ~ x + tr + m, data = sim))
  expect_identical(fit$serial$kmax, fgls:::default_kmax(120))
  expect_no_warning(
    suppressWarnings(serial_fgls(y ~ x + tr + m, data = sim),
      classes = c("fgls_warning_specification", "fgls_warning_nonstationary")
    ),
    class = "fgls_warning_sample_size"
  )
  # an explicit kmax and a fixed order the real design supports
  f9 <- suppressWarnings(serial_fgls(y ~ x + tr + m, data = sim, kmax = 9))
  expect_identical(f9$serial$kmax, 9L)
  f8 <- suppressWarnings(serial_fgls(y ~ x + tr + m, data = sim, ar_order = 8))
  expect_identical(f8$ar_order, 8L)
  expect_identical(f8$n_eff, 112L)
  # the cap still binds, with the actual count in the message
  X <- fit$x
  expect_warning(
    k <- fgls:::resolve_kmax(40, 120, ncol(X), TRUE, y = sim$y, X = X),
    "regressors \\(collinear lagged columns excluded\\)", class = "fgls_warning_sample_size"
  )
  cnt <- fgls:::durbin_regressor_count(sim$y, X, k, TRUE)
  expect_lte(cnt, 0.5 * (120 - k))
  expect_gt(fgls:::durbin_regressor_count(sim$y, X, k + 1L, TRUE), 0.5 * (120 - k - 1L))
  # without the data the conservative pyfgls count is used
  expect_warning(
    k0 <- fgls:::resolve_kmax(NULL, 120, ncol(X), TRUE),
    class = "fgls_warning_sample_size"
  )
  expect_identical(k0, 3L)
  expect_error(
    serial_fgls(y ~ x + tr + m, data = sim, ar_order = 60),
    "collinear lagged columns excluded", class = "fgls_error_insufficient_data"
  )
})

test_that("without deterministic lags the cap is the pyfgls one", {
  set.seed(21)
  n <- 60
  X <- cbind(1, matrix(rnorm(n * 4), n, 4))
  colnames(X) <- c("(Intercept)", paste0("x", 1:4))
  y <- drop(X %*% rep(1, 5)) + rnorm(n)
  for (k in 0:6) {
    expect_identical(
      fgls:::durbin_regressor_count(y, X, k, TRUE),
      as.integer(fgls:::n_durbin_regressors(5, k, TRUE))
    )
  }
  expect_warning(
    k1 <- fgls:::resolve_kmax(20, n, 5, TRUE, y = y, X = X),
    class = "fgls_warning_sample_size"
  )
  expect_warning(k2 <- fgls:::resolve_kmax(20, n, 5, TRUE), class = "fgls_warning_sample_size")
  expect_identical(k1, k2)
  # no constant: nothing is dropped
  X0 <- X[, -1]
  expect_identical(
    fgls:::durbin_regressor_count(y, X0, 3, FALSE),
    as.integer(fgls:::n_durbin_regressors(4, 3, FALSE))
  )
})

test_that("durbin_estimator() checks a fixed order against the kept regressors", {
  sim <- monthly_data()
  du <- durbin_estimator(y ~ x + tr + m, data = sim, ar_order = 8)
  expect_identical(du$serial$ar_order, 8L)
  expect_error(
    durbin_estimator(y ~ x, data = simulate_serial(40, seed = 1), ar_order = 30),
    class = "fgls_error_insufficient_data"
  )
})

# ---------------------------------------- findings 18 and 27: Durbin estimator
test_that("durbin_estimator() recovers the intercept and the trend coefficient", {
  s2 <- simulate_serial(2000, seed = 1)
  du <- durbin_estimator(y ~ x, data = s2)
  fg <- serial_fgls(y ~ x, data = s2)
  X <- cbind("(Intercept)" = 1, x = s2$x)
  raw <- fgls:::fit_durbin(fgls:::durbin_design(s2$y, X, colnames(X), du$ar_order))
  rho <- unname(du$rho)
  expect_identical(coef(du)[["x"]], unname(raw$beta[2]))
  expect_equal(coef(du)[["(Intercept)"]], unname(raw$beta[1]) / (1 - sum(rho)), tolerance = 1e-12)
  expect_equal(coef(du)[["(Intercept)"]], coef(fg)[["(Intercept)"]], tolerance = 0.01)
  expect_lt(abs(mean(residuals(du))), 0.01)
  expect_equal(unname(predict(du, data.frame(x = 0))), coef(du)[["(Intercept)"]])
  expect_match(du$notes, "recovered", all = FALSE)

  s2$t <- seq_len(2000)
  s2$y2 <- s2$y + 0.5 * s2$t
  dt <- durbin_estimator(y2 ~ x + t, data = s2)
  ft <- serial_fgls(y2 ~ x + t, data = s2)
  expect_equal(coef(dt)[["t"]], coef(ft)[["t"]], tolerance = 1e-3)
  expect_equal(coef(dt)[["(Intercept)"]], coef(ft)[["(Intercept)"]], tolerance = 0.05)
})

test_that("the recovered coefficients have the delta-method covariance", {
  sim <- simulate_serial(400, ar = c(0.5, 0.2), seed = 9)
  sim$t <- seq_len(400) / 400
  sim$q <- as.numeric((seq_len(400) - 1) %% 4 == 0)
  du <- durbin_estimator(y ~ x + t + q, data = sim, ar_order = 5)
  X <- du$x
  des <- fgls:::durbin_design(du$y, X, colnames(X), 5)
  fit <- fgls:::fit_durbin(des)
  Z <- des$Z
  V <- sum(fit$resid^2) / (nrow(Z) - ncol(Z)) * solve(crossprod(Z))
  st <- function(theta) {
    f <- fit
    f$coef <- theta
    f$beta <- theta[seq_len(ncol(X))]
    fgls:::durbin_structural(f, X, V)$beta
  }
  theta <- fit$coef
  expect_equal(unname(coef(du)), st(theta), tolerance = 1e-12)
  h <- 1e-6
  Jn <- vapply(seq_along(theta), function(i) {
    e <- replace(numeric(length(theta)), i, h)
    (st(theta + e) - st(theta - e)) / (2 * h)
  }, numeric(ncol(X)))
  expect_equal(unname(vcov(du)), Jn %*% V %*% t(Jn), tolerance = 1e-6)
  # the stochastic regressor is untouched
  expect_identical(coef(du)[["x"]], unname(fit$beta[["x"]]))
  expect_equal(unname(vcov(du)["x", "x"]), V[2, 2], tolerance = 1e-12)
})

test_that("seasonal dummies and partial dummies are recovered", {
  set.seed(11)
  n <- 1200
  x <- rnorm(n)
  m <- factor(rep(1:12, length.out = n))
  seas <- c(0, 1, -1, 0.5, 2, -0.5, 1.5, -2, 0.3, -0.3, 1, -1)
  u <- fgls:::arma_errors(rnorm(n), c(0.5, 0.2))
  df <- data.frame(y = 1 + x + seas[as.integer(m)] + u, x = x, m = m)
  du <- durbin_estimator(y ~ x + m, data = df)
  fg <- serial_fgls(y ~ x + m, data = df)
  expect_equal(unname(coef(du)), unname(coef(fg)), tolerance = 0.05)
  # quarterly dummy for the first quarter only: L1 and L2 are kept, L3 = 1 -
  # q - L1 - L2 is dropped and absorbed by the constant and by q
  q <- as.numeric((seq_len(n) - 1) %% 4 == 0)
  d2 <- data.frame(y = 1 + x + 0.8 * q + u, x = x, q = q)
  dq <- durbin_estimator(y ~ x + q, data = d2, ar_order = 3)
  fq <- serial_fgls(y ~ x + q, data = d2, ar_order = 3)
  expect_true("L3.q" %in% fgls:::durbin_design(d2$y, dq$x, colnames(dq$x), 3)$dropped)
  expect_equal(unname(coef(dq)), unname(coef(fq)), tolerance = 0.1)
})

test_that("coefficients that cannot be recovered are NA with a note", {
  sim <- simulate_serial(300, seed = 4)
  df <- data.frame(y = sim$y[-1], ylag = sim$y[-300], x = sim$x[-1])
  du <- suppressWarnings(durbin_estimator(y ~ ylag + x, data = df, ar_order = 1))
  expect_true(is.na(coef(du)[["(Intercept)"]]))
  expect_false(anyNA(coef(du)[c("ylag", "x")]))
  expect_match(du$notes, "are NA", all = FALSE)
  expect_output(print(summary(du)), "Durbin")
  # a unit root in the estimated filter leaves the intercept unidentified
  X <- cbind("(Intercept)" = 1, x = sim$x)
  des <- fgls:::durbin_design(sim$y, X, colnames(X), 1)
  fit <- fgls:::fit_durbin(des)
  fit$coef[des$ylag_pos[1]] <- 1
  V <- diag(ncol(des$Z))
  st <- fgls:::durbin_structural(fit, X, V)
  expect_true(is.na(st$beta[1]))
  expect_match(st$notes, "unit root")
  # nothing to recover without a constant
  d0 <- suppressWarnings(durbin_estimator(y ~ x - 1, data = sim, ar_order = 1))
  r0 <- fgls:::fit_durbin(fgls:::durbin_design(sim$y, cbind(x = sim$x), "x", 1))
  expect_identical(unname(coef(d0)), unname(r0$beta))
  expect_length(d0$notes, 0L)
})

test_that("durbin_estimator() innovations are the Durbin-regression residuals", {
  sim <- simulate_serial(200, seed = 1)
  du <- durbin_estimator(y ~ x, data = sim)
  k <- du$ar_order
  X <- cbind("(Intercept)" = 1, x = sim$x)
  raw <- fgls:::fit_durbin(fgls:::durbin_design(sim$y, X, colnames(X), k))
  expect_length(du$innovations, 200L - k)
  expect_equal(unname(du$innovations), raw$resid, tolerance = 1e-12)
  expect_identical(names(du$innovations), as.character((k + 1L):200L))
  expect_lt(abs(stats::acf(du$innovations, plot = FALSE)$acf[2]), 0.2)
})

# ---------------------------------------------------- finding 26: share is NA
test_that("benchmark fits and GMA report no non-exogeneity share", {
  sim <- simulate_serial(200, seed = 1)
  fits <- list(
    co = cochrane_orcutt(y ~ x, data = sim),
    du = durbin_estimator(y ~ x, data = sim),
    gls = infeasible_gls(y ~ x, data = sim, rho = 0.8),
    hy = hybrid_gls(y ~ x, data = sim, seed = 1),
    ma1 = serial_fgls(y ~ x, data = sim, serial = "ma1")
  )
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_true(is.na(nonexogeneity_share(f)), label = nm)
    expect_true(is.na(f$nonexogeneity_share), label = nm)
  }
  expect_identical(summary(fits$co)$info[["Non-exogeneity share"]], "-")
  expect_output(print(summary(fits$co)), "Non-exogeneity share: -")
  skip_if_not_installed("broom")
  expect_true(is.na(broom::glance(fits$co)$nonexogeneity.share))
  expect_true(is.na(broom::glance(fits$ma1)$nonexogeneity.share))
  expect_false(is.na(broom::glance(serial_fgls(y ~ x, data = sim))$nonexogeneity.share))
})

# ------------------------------------------------- finding 31: hybrid_gls args
test_that("hybrid_gls() validates lags and alpha", {
  sim <- simulate_serial(200, seed = 1)
  for (bad in list(0, -1, 1.5, NA, c(1, 2), "1")) {
    expect_error(hybrid_gls(y ~ x, data = sim, lags = bad), "lags", class = "fgls_error_input")
  }
  for (bad in list(0, 1, 2, -0.1, NA, c(0.05, 0.1), "0.05")) {
    expect_error(hybrid_gls(y ~ x, data = sim, alpha = bad), "alpha", class = "fgls_error_input")
  }
  expect_s3_class(hybrid_gls(y ~ x, data = sim, lags = 2, alpha = 0.1, seed = 1), "fgls")
  expect_error(
    hybrid_gls(y ~ x, data = sim[1:4, ], lags = 2),
    class = "fgls_error_insufficient_data"
  )
})

# ------------------------------------- finding 32: zero residual degrees of freedom
test_that("Cochrane-Orcutt and infeasible GLS need positive degrees of freedom", {
  d <- data.frame(y = c(1, 3, 2), x = c(0.5, 1, 3))
  expect_error(cochrane_orcutt(y ~ x, data = d), class = "fgls_error_insufficient_data")
  d4 <- data.frame(y = c(1, 3, 2, 5), x = c(0.5, 1, 3, 2))
  expect_error(cochrane_orcutt(y ~ x, data = d4, ar_order = 2),
    class = "fgls_error_insufficient_data"
  )
  expect_error(infeasible_gls(y ~ x, data = d, rho = 0.5), class = "fgls_error_insufficient_data")
  expect_error(infeasible_gls(y ~ x, data = d4, rho = c(0.5, 0.1)),
    class = "fgls_error_insufficient_data"
  )
  ok <- suppressWarnings(cochrane_orcutt(y ~ x, data = d4))
  expect_true(all(is.finite(vcov(ok))))
  expect_true(all(is.finite(vcov(infeasible_gls(y ~ x, data = d4, rho = 0.5)))))
})

# -------------------------------------------- findings 2 and 12: hybrid seed
test_that("the hybrid_gls() seed drives the random AR coefficient only", {
  sim <- simulate_serial(200, seed = 1)
  h1 <- hybrid_gls(y ~ x, data = sim, seed = 1)
  h2 <- hybrid_gls(y ~ x, data = sim, seed = 2)
  expect_identical(h1$method, "Hybrid (GLS)")
  expect_false(identical(h1$rho, h2$rho))
  expect_identical(h1$rho, hybrid_gls(y ~ x, data = sim, seed = 1)$rho)
  set.seed(99)
  before <- .Random.seed
  hybrid_gls(y ~ x, data = sim, seed = 5)
  expect_identical(.Random.seed, before)
})
