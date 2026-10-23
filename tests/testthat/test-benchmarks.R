sim <- simulate_serial(300, ar = 0.6, gamma = 0, rho_x = 0.5, burn = 50, seed = 5)

test_that("HAC, Newey-West and HC covariances agree with the sandwich package", {
  skip_if_not_installed("sandwich")
  lmf <- stats::lm(y ~ x, data = sim)
  qs <- ols_robust(y ~ x, data = sim)
  X <- cbind(1, sim$x)
  bw <- fgls:::andrews_bandwidth(X * residuals(lmf), "qs", skip = 1L)
  ref <- sandwich::kernHAC(lmf,
    kernel = "Quadratic Spectral", bw = bw, prewhite = FALSE, adjust = FALSE
  )
  expect_equal(unname(vcov(qs)), unname(ref), tolerance = 1e-10)
  # sandwich's bwAndrews() demeans the scores in its AR(1) fits; pyfgls does not
  bw_sw <- sandwich::bwAndrews(lmf, kernel = "Quadratic Spectral", prewhite = FALSE)
  expect_equal(bw, bw_sw, tolerance = 0.02)
  nw <- ols_robust(y ~ x, data = sim, vcov = "NW")
  lag <- floor(4 * (300 / 100)^(2 / 9))
  ref_nw <- sandwich::NeweyWest(lmf, lag = lag, prewhite = FALSE, adjust = FALSE)
  expect_equal(unname(vcov(nw)), unname(ref_nw), tolerance = 1e-10)
  for (ct in c("HC0", "HC1", "HC2", "HC3")) {
    h <- ols_robust(y ~ x, data = sim, vcov = ct)
    expect_equal(unname(vcov(h)), unname(sandwich::vcovHC(lmf, type = ct)), tolerance = 1e-10)
  }
  nr <- ols_robust(y ~ x, data = sim, vcov = "nonrobust")
  expect_equal(unname(vcov(nr)), unname(vcov(lmf)), tolerance = 1e-10)
  fixed <- ols_robust(y ~ x, data = sim, kernel = "bartlett", bandwidth = 5, small_sample = TRUE)
  ref_b <- sandwich::kernHAC(lmf,
    kernel = "Bartlett", bw = 5, prewhite = FALSE,
    adjust = TRUE
  )
  expect_equal(unname(vcov(fixed)), unname(ref_b), tolerance = 1e-10)
  expect_error(ols_robust(y ~ x, data = sim, kernel = "parzen"), class = "fgls_error_input")
  expect_error(ols_robust(y ~ x, data = sim, bandwidth = -1), class = "fgls_error_input")
})

test_that("the QS kernel and the Andrews bandwidth", {
  expect_equal(fgls:::qs_kernel(0), 1)
  x <- c(0.5, 1, 2)
  z <- 6 * pi * x / 5
  expect_equal(fgls:::qs_kernel(x), 3 / z^2 * (sin(z) / z - cos(z)))
  expect_equal(fgls:::bartlett_kernel(c(0, 0.5, 2)), c(1, 0.5, 0))
  expect_equal(fgls:::andrews_bandwidth(matrix(0, 10, 2)), 1)
})

test_that("Cochrane-Orcutt two-step and iterated", {
  co <- cochrane_orcutt(y ~ x, data = sim)
  u <- residuals(stats::lm(y ~ x, data = sim))
  rho <- sum(u[-1] * u[-300]) / sum(u[-300]^2)
  expect_equal(unname(co$rho), rho, tolerance = 1e-10)
  X <- cbind(1, sim$x)
  b <- stats::lm.fit(fgls:::ar_filter(X, rho), fgls:::ar_filter(sim$y, rho))$coefficients
  expect_equal(unname(coef(co)), unname(b), tolerance = 1e-10)
  expect_match(co$notes, "not valid", all = FALSE)
  it <- cochrane_orcutt(y ~ x, data = sim, iterate = TRUE)
  expect_identical(it$method, "Cochrane-Orcutt (iterated)")
  expect_gte(it$serial$extra$iterations, 2)
  co2 <- cochrane_orcutt(y ~ x, data = sim, ar_order = 2)
  expect_length(co2$rho, 2L)
  expect_error(cochrane_orcutt(y ~ x, data = sim, ar_order = 0), class = "fgls_error_input")
  expect_warning(cochrane_orcutt(y ~ x, data = sim, iterate = TRUE, max_iter = 1, tol = 0),
    class = "fgls_warning_convergence"
  )
})

test_that("the Durbin estimator equals OLS on the Durbin design", {
  du <- durbin_estimator(y ~ x, data = sim, ar_order = 2)
  n <- 300
  Z <- cbind(1, sim$x[3:n], sim$y[2:(n - 1)], sim$x[2:(n - 1)], sim$y[1:(n - 2)], sim$x[1:(n - 2)])
  lmd <- stats::lm.fit(Z, sim$y[3:n])
  cf <- unname(lmd$coefficients)
  # the slope is the Durbin coefficient on x_t; the lags of the constant are
  # collinear with it and dropped, so the intercept is recovered from the
  # Durbin coefficient c on the constant as c / (1 - rho_1 - rho_2)
  m <- 1 - cf[3] - cf[5]
  expect_equal(unname(coef(du)), c(cf[1] / m, cf[2]), tolerance = 1e-10)
  expect_equal(unname(du$rho), cf[c(3, 5)], tolerance = 1e-10)
  ref <- unname(vcov(stats::lm(sim$y[3:n] ~ Z - 1)))
  jac <- rbind(c(1 / m, 0, cf[1] / m^2, 0, cf[1] / m^2, 0), c(0, 1, 0, 0, 0, 0))
  expect_equal(unname(vcov(du)), jac %*% ref %*% t(jac), tolerance = 1e-10)
  expect_equal(durbin_estimator(y ~ x, data = sim)$serial$ar_order, serial_fgls(y ~ x, data = sim)$ar_order)
  expect_error(durbin_estimator(y ~ x, data = sim, ar_order = -1), class = "fgls_error_input")
})

test_that("infeasible GLS with known rho or a known covariance matrix", {
  g <- infeasible_gls(y ~ x, data = sim, rho = 0.6)
  expect_equal(coef(g), coef(serial_fgls(y ~ x, data = sim, ar_order = 1)), tolerance = 0.2)
  X <- cbind(1, sim$x)
  b <- stats::lm.fit(fgls:::ar_filter(X, 0.6), fgls:::ar_filter(sim$y, 0.6))$coefficients
  expect_equal(unname(coef(g)), unname(b), tolerance = 1e-10)
  # A small T x T omega: with an OpenMP OpenBLAS and OMP_THREAD_LIMIT below the
  # number of cores (as on some check machines) a multithreaded Cholesky
  # factorization of a large matrix can deadlock.
  s50 <- sim[1:50, ]
  Om <- 0.6^abs(outer(1:50, 1:50, "-")) / (1 - 0.36)
  go <- infeasible_gls(y ~ x, data = s50, omega = Om)
  X50 <- cbind(1, s50$x)
  Oi <- solve(Om)
  bo <- solve(t(X50) %*% Oi %*% X50, t(X50) %*% Oi %*% s50$y)
  expect_equal(unname(coef(go)), drop(bo), tolerance = 1e-8)
  expect_error(infeasible_gls(y ~ x, data = sim), class = "fgls_error_input")
  expect_error(infeasible_gls(y ~ x, data = sim, rho = 0.5, omega = Om), class = "fgls_error_input")
  expect_error(infeasible_gls(y ~ x, data = sim, omega = diag(3)), class = "fgls_error_input")
})

test_that("hybrid GLS: OLS without serial correlation, quasi-differencing with it", {
  wn <- simulate_serial(300, ar = numeric(0), gamma = 0, seed = 8)
  h0 <- hybrid_gls(y ~ x, data = wn, seed = 1)
  if (h0$method == "Hybrid (OLS)") expect_equal(coef(h0), coef(stats::lm(y ~ x, data = wn)))
  h1 <- hybrid_gls(y ~ x, data = sim, seed = 1)
  expect_identical(h1$method, "Hybrid (GLS)")
  expect_true(h1$rho > 0.1 && h1$rho < 0.9)
  expect_identical(h1$rho, hybrid_gls(y ~ x, data = sim, seed = 1)$rho)
})

sandwich_free_hc3 <- function(d) {
  X <- cbind(1, d$x)
  u <- residuals(stats::lm(y ~ x, data = d))
  h <- fgls:::hat_diag(X)
  A <- solve(crossprod(X))
  A %*% crossprod(X * (u / (1 - h))) %*% A
}

test_that("compare_fgls() builds a side-by-side table", {
  fit <- serial_fgls(y ~ x, data = sim)
  lmf <- stats::lm(y ~ x, data = sim)
  tab <- compare_fgls(OLS = lmf, FGLS = fit)
  expect_named(tab, c("term", "statistic", "OLS", "FGLS"))
  expect_equal(nrow(tab), 6L)
  se_hac <- sqrt(diag(vcov(ols_robust(y ~ x, data = sim))))
  expect_equal(tab$OLS[tab$statistic == "std.error"], unname(se_hac), tolerance = 1e-10)
  expect_equal(tab$FGLS[tab$statistic == "ci.width"], unname(2 * qnorm(0.975) * sqrt(diag(vcov(fit)))))
  t2 <- compare_fgls(fit, ols_robust(y ~ x, data = sim), level = 0.9)
  expect_named(t2, c("term", "statistic", "FGLS (serial)", "OLS"))
  t3 <- compare_fgls(fit, serial_fgls(y ~ x, data = sim, vcov = "naive"))
  expect_true(any(grepl("\\[naive\\]", names(t3))))
  het <- simulate_hetero(100, seed = 1)
  t4 <- compare_fgls(stats::lm(y ~ x, data = het), lasso_fgls(y ~ x, data = het, seed = 1))
  se_hc3 <- sqrt(diag(sandwich_free_hc3(het)))
  expect_equal(t4$OLS[t4$statistic == "std.error"], unname(se_hc3), tolerance = 1e-10)
  expect_error(compare_fgls(), class = "fgls_error_input")
  expect_error(compare_fgls(1), class = "fgls_error_input")
  t5 <- compare_fgls(fit, lmf, ols_vcov = "HC3")
  expect_equal(ncol(t5), 4L)
})

test_that("broom tidiers", {
  skip_if_not_installed("broom")
  fit <- fgls(y ~ x, data = simulate_combined(200, seed = 2), seed = 1)
  td <- broom::tidy(fit, conf.int = TRUE)
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$std.error, unname(sqrt(diag(vcov(fit)))))
  expect_equal(td$conf.low, unname(confint(fit)[, 1]))
  gl <- broom::glance(fit)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$ar.order, fit$ar_order)
  au <- broom::augment(fit)
  expect_equal(nrow(au), 200L)
  expect_equal(sum(is.na(au$.innovation)), fit$ar_order)
  expect_equal(sum(is.na(au$.weight)), fit$ar_order)
  an <- broom::augment(fit, newdata = data.frame(x = c(1, 2)))
  expect_equal(an$.fitted, unname(predict(fit, data.frame(x = c(1, 2)))))
  expect_error(broom::augment(fit, data = data.frame(a = 1)), class = "fgls_error_input")
  s <- broom::tidy(serial_fgls(y ~ x, data = sim))
  expect_false("conf.low" %in% names(s))
})

test_that("plot() draws the diagnostics", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fit <- fgls(y ~ x, data = simulate_combined(200, seed = 3), seed = 1)
  expect_invisible(plot(fit))
  expect_invisible(plot(serial_fgls(y ~ x, data = sim), which = 2:4))
  expect_invisible(plot(serial_fgls(y ~ x, data = sim, ar_order = 1), which = 3))
  expect_error(plot(fit, which = 7), class = "fgls_error_input")
})
