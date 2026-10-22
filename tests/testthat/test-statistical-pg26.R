# Acceptance tests replicating the PG26 simulation results (serial correlation).
# Fewer replications than the paper (10,000); tolerances cover the Monte Carlo
# noise (relative 15% on RMSE/MSE, absolute 0.01 on biases).

ols_slope <- function(y, X) fgls:::lsq(X, y)[2]

cochrane_orcutt_two_step <- function(y, X) {
  u <- drop(y - X %*% fgls:::lsq(X, y))
  n <- length(u)
  rho <- sum(u[-1] * u[-n]) / sum(u[-n]^2)
  list(beta = fgls:::lsq(fgls:::ar_filter(X, rho), fgls:::ar_filter(y, rho)), rho = rho)
}

test_that("PG26 Table 1: FGLS is consistent with a non-exogenous regressor", {
  skip_if_not_slow()
  r <- run_mc(2000, 20260101, function() {
    s <- simulate_serial(500) # AR(1) 0.8, x_t = v_t + e_{t-1}, u_0 = 0
    X <- cbind(1, s$x)
    est <- fgls:::estimate_serial(s$y, X, c("c", "x"), kmax = 12)
    co <- cochrane_orcutt_two_step(s$y, X)
    list(
      ols = ols_slope(s$y, X) - 1, fgls = est$beta[2] - 1, durbin = est$durbin$beta[2] - 1,
      gls = fgls:::lsq(fgls:::ar_filter(X, 0.8), fgls:::ar_filter(s$y, 0.8))[2] - 1,
      rho = sum(est$rho) - 0.8, k = est$order, co = co$beta[2] - 1, co_rho = co$rho - 0.8
    )
  })
  expect_equal(rmse(r$ols), 0.399, tolerance = 0.15)
  expect_lt(abs(mean(r$ols) - 0.395), 0.01)
  expect_equal(rmse(r$fgls), 0.025, tolerance = 0.15)
  expect_lt(abs(mean(r$fgls) - 0.001), 0.01)
  expect_equal(rmse(r$gls), 0.025, tolerance = 0.15)
  expect_equal(rmse(r$durbin), 0.036, tolerance = 0.15)
  expect_lt(rmse(r$fgls), rmse(r$durbin))
  expect_equal(rmse(r$rho), 0.034, tolerance = 0.15)
  expect_lt(abs(mean(r$rho) + 0.008), 0.01)
  expect_equal(stats::median(r$k), 1)
  expect_equal(rmse(r$co), 0.041, tolerance = 0.15)
  expect_lt(abs(mean(r$co) - 0.031), 0.01)
  expect_lt(abs(mean(r$co_rho) + 0.171), 0.01)
})

pg26_cells <- list(
  "AR1(0.5)" = 0.5, "AR1(0.8)" = 0.8,
  "AR2(1.34,-0.42)" = c(1.34, -0.42), "AR2(0.5,-0.3)" = c(0.5, -0.3)
)

pg26_general <- function(ar, gamma, R = 1000, n = 200) {
  run_mc(R, 4242, function() {
    s <- simulate_serial(n, ar = ar, rho_x = 0.8, gamma = gamma, burn = 100)
    X <- cbind(1, s$x)
    est <- fgls:::estimate_serial(s$y, X, c("c", "x"), kmax = 12)
    b <- est$beta[2]
    v <- est$variance
    list(
      ols = ols_slope(s$y, X) - 1, fgls = b - 1,
      gls = fgls:::lsq(fgls:::ar_filter(X, ar), fgls:::ar_filter(s$y, ar))[2] - 1,
      cov_r = covered(b, sqrt(est$cov[2, 2]), 1),
      cov_n = covered(b, sqrt(v$naive_avar[2, 2] / v$n), 1), k = est$order
    )
  })
}

test_that("PG26 Tables 2-5: robust coverage, naive under-coverage, MSE, median k*", {
  skip_if_not_slow()
  for (cell in names(pg26_cells)) {
    for (gamma in c(0, 0.25, 0.5)) {
      ar <- pg26_cells[[cell]]
      r <- pg26_general(ar, gamma)
      lab <- sprintf("%s, gamma = %g", cell, gamma)
      cov_r <- mean(r$cov_r)
      cov_n <- mean(r$cov_n)
      expect_true(cov_r >= 0.92 && cov_r <= 0.96, label = paste(lab, "robust coverage", cov_r))
      if (gamma == 0.5) {
        expect_gte(cov_r - cov_n, 0.01, label = paste(lab, "naive under-coverage"))
      } else if (gamma == 0) {
        expect_lte(abs(cov_r - cov_n), 0.01, label = paste(lab, "robust vs naive"))
      } else {
        expect_gte(cov_r, cov_n, label = lab)
      }
      expect_lte(mse(r$fgls), mse(r$ols), label = paste(lab, "FGLS vs OLS MSE"))
      expect_lte(mse(r$fgls), 1.35 * mse(r$gls), label = paste(lab, "FGLS vs GLS MSE"))
      expect_equal(stats::median(r$k), length(ar), label = paste(lab, "median k*"))
    }
  }
})

test_that("PG26: AR(2) (1.34, -0.42) errors, ~96% MSE reduction over OLS", {
  skip_if_not_slow()
  r <- pg26_general(c(1.34, -0.42), 0)
  ols <- 100 * mse(r$ols)
  fg <- 100 * mse(r$fgls)
  expect_equal(ols, 11.5, tolerance = 0.15)
  expect_equal(fg, 0.41, tolerance = 0.15)
  expect_lt(fg / ols, 0.05)
})

test_that("PG26: FGLS beats OLS with MA(1) errors", {
  skip_if_not_slow()
  for (theta in c(-0.4, 0.5)) {
    r <- run_mc(500, 77, function() {
      s <- simulate_serial(200, ar = numeric(0), ma = theta, rho_x = 0.8, gamma = 0.25, burn = 100)
      X <- cbind(1, s$x)
      est <- fgls:::estimate_serial(s$y, X, c("c", "x"), kmax = 12)
      list(ols = ols_slope(s$y, X) - 1, fgls = est$beta[2] - 1)
    })
    expect_lt(mse(r$fgls), mse(r$ols))
  }
})

test_that("PG26 Section 4.2: lagged dependent variable (Wallis) path", {
  skip_if_not_slow()
  n_obs <- 500
  r <- run_mc(600, 99, function() {
    n <- n_obs + 201
    e <- stats::rnorm(n)
    w <- as.numeric(stats::filter(stats::rnorm(n), 0.5, method = "recursive"))
    u <- as.numeric(stats::filter(e, 0.5, method = "recursive"))
    y <- numeric(n)
    for (t in 2:n) y[t] <- 0.5 * y[t - 1] + 1 + w[t] + u[t]
    Y <- y[202:n]
    X <- cbind(1, y[201:(n - 1)], w[202:n])
    est <- fgls:::estimate_serial(Y, X, c("c", "ylag", "w"), kmax = 12, lagged_dep_idx = 2L)
    se <- sqrt(diag(est$cov))
    list(
      a = est$beta[2] - 0.5, b = est$beta[3] - 1, cov_a = covered(est$beta[2], se[2], 0.5),
      cov_b = covered(est$beta[3], se[3], 1), ols_a = fgls:::lsq(X, Y)[2] - 0.5, k = est$order
    )
  })
  expect_lt(abs(mean(r$a)), 0.01)
  expect_lt(abs(mean(r$b)), 0.01)
  expect_gt(mean(r$ols_a), 0.1)
  expect_true(mean(r$cov_a) >= 0.92 && mean(r$cov_a) <= 0.97)
  expect_true(mean(r$cov_b) >= 0.92 && mean(r$cov_b) <= 0.97)
  expect_equal(stats::median(r$k), 1)
})
