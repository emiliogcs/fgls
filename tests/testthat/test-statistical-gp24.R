# Acceptance tests replicating the GP24 simulation results (heteroskedasticity).
#
# DGP (eq. 9): y = 1 + x + sqrt(v(x)) eps, x ~ U[1, 4]; candidate covariates
# z = (x, [log x]^2, x^2, cos x, cos 2x) plus the unpenalised intercept. The
# efficiency measure is the empirical MSE of the slope relative to OLS (GP24
# Table 1). As in pyfgls, the T = 400 column is matched (within 0.05); at
# T = 100-200 the ratios of the adaptive-Lasso FGLS are a few hundredths above
# the paper's (an open item inherited from pyfgls).

paper_z <- function(x) cbind(x, log(x)^2, x^2, cos(x), cos(2 * x))

ols_hc3 <- function(X, y) {
  b <- fgls:::lsq(X, y)
  u <- drop(y - X %*% b)
  h <- fgls:::hat_diag(X)
  A <- solve(crossprod(X))
  V <- A %*% crossprod(X * (u / (1 - h))) %*% A
  list(b = b, se = sqrt(diag(V)))
}

gp24_table1 <- function(skedastic, n, R, seed) {
  run_mc(R, seed, function() {
    s <- simulate_hetero(n, skedastic = skedastic)
    X <- cbind(1, s$x)
    est <- fgls:::estimate_hetero(s$y, X, paper_z(s$x), letters[1:5],
      seed = sample.int(.Machine$integer.max, 1)
    )
    se <- sqrt(diag(est$cov))
    o <- ols_hc3(X, s$y)
    list(
      fgls = est$beta[2] - 1, ols = o$b[2] - 1,
      cover = abs(est$beta[2] - 1) <= qnorm(0.975) * se[2],
      len = 2 * qnorm(0.975) * se[2], len_ols = 2 * qnorm(0.975) * o$se[2]
    )
  })
}

emse_ratio <- function(r) mean(r$fgls^2) / mean(r$ols^2)

gp24_t400 <- c(homoskedastic = 1.00, x = 0.92, x2 = 0.70, log2 = 0.46, exp = 0.45, step = 0.735)

test_that("GP24 Table 1 (T = 400): efficiency relative to OLS and coverage", {
  skip_if_not_slow()
  for (sk in names(gp24_t400)) {
    r <- gp24_table1(sk, 400, 1000, 2024)
    expect_lt(abs(emse_ratio(r) - gp24_t400[[sk]]), 0.05, label = paste(sk, "eMSE ratio", emse_ratio(r)))
    expect_lt(abs(mean(r$cover) - 0.95), 0.02, label = paste(sk, "coverage", mean(r$cover)))
  }
})

test_that("GP24 Table 1 (T = 100): coverage, no loss, shorter intervals", {
  skip_if_not_slow()
  for (sk in c("homoskedastic", "x2", "log2", "exp", "step")) {
    r <- gp24_table1(sk, 100, 1500, 99)
    expect_lt(abs(mean(r$cover) - 0.95), 0.015, label = paste(sk, "coverage", mean(r$cover)))
    if (sk == "homoskedastic") {
      expect_lte(emse_ratio(r), 1.05, label = "no loss at T = 100")
    } else {
      expect_lt(emse_ratio(r), 0.85, label = paste(sk, "eMSE ratio"))
      expect_lt(mean(r$len), mean(r$len_ols), label = paste(sk, "interval length"))
    }
  }
})

test_that("GP24: no loss under homoskedasticity in larger samples", {
  skip_if_not_slow()
  expect_lte(emse_ratio(gp24_table1("homoskedastic", 800, 500, 7)), 1.02)
})

test_that("GP24 Table 2: interval lengths for v = [log x]^2, T = 100", {
  skip_if_not_slow()
  r <- gp24_table1("log2", 100, 1000, 11)
  expect_lt(abs(mean(r$len) - 0.315), 0.03)
  expect_lt(abs(mean(r$len_ols) - 0.44), 0.03)
})

gp24_extended <- function(case, R = 400, n = 200, k_w = 20) {
  run_mc(R, 99, function() {
    corr <- if (case == "proxy") 0.5 else 0
    sk <- if (case == "homoskedastic") "homoskedastic" else "x2"
    s <- simulate_hetero(n, skedastic = sk, n_extra = k_w, extra_corr = corr)
    W <- as.matrix(s[, paste0("w", seq_len(k_w))])
    blocks <- cbind(W, cos(W), cos(2 * W))
    if (case %in% c("driver", "homoskedastic")) {
      blocks <- cbind(s$x, s$x^2, cos(s$x), cos(2 * s$x), blocks)
    }
    X <- cbind(1, s$x)
    est <- fgls:::estimate_hetero(s$y, X, blocks, sprintf("z%d", seq_len(ncol(blocks))),
      seed = sample.int(.Machine$integer.max, 1)
    )
    list(fgls = est$beta[2] - 1, ols = fgls:::lsq(X, s$y)[2] - 1)
  })
}

test_that("GP24 extended designs: driver included or proxied beats OLS, no noticeable loss otherwise", {
  skip_if_not_slow()
  expect_lt(emse_ratio(gp24_extended("driver")), 0.85)
  expect_lt(emse_ratio(gp24_extended("proxy")), 1.0)
  for (case in c("excluded", "homoskedastic")) {
    expect_lte(emse_ratio(gp24_extended(case)), 1.06, label = case)
  }
})
