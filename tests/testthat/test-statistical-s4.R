# Acceptance tests: PG26 Supplement S-4 (serial correlation + heteroskedasticity),
# Table S.9 designs. y = 1 + x + u, x ~ U[1, 4], u ARMA with innovations
# e_t = sqrt(v(w_t)) eps_t; candidate covariates W = (w, w^2, cos w, cos 2w)
# with w = x (case a) or a Bernoulli(0.75) mixture (case b); T = 200, kmax = 12.

s4_cell <- function(ar = 0.5, ma = numeric(0), skedastic = "x2", driver = "x", gamma = 0,
                    R = 1000, seed = 31) {
  run_mc(R, seed, function() {
    s <- simulate_combined(200, ar = ar, ma = ma, skedastic = skedastic, driver = driver, gamma = gamma)
    w <- s$w
    W <- cbind(w, w^2, cos(w), cos(2 * w))
    X <- cbind(1, s$x)
    sd <- sample.int(.Machine$integer.max, 1)
    est <- fgls:::estimate_combined(s$y, X, W, letters[1:4], x_names = c("c", "x"), kmax = 12, seed = sd)
    fr <- est$serial$filtered
    het <- est$hetero
    k <- est$order
    V_nc <- fgls:::combined_covariance(fr$Xf, drop(fr$yf - fr$Xf %*% het$beta),
      drop(s$y - X %*% het$beta), het$v, k, NULL, fgls:::hat_diag(fr$Xf), het$df)
    lasso_only <- fgls:::estimate_hetero(s$y, X, W, letters[1:4], seed = sd)
    se_fc <- sqrt(est$cov_fc[2, 2])
    se_fch <- sqrt(est$cov_fch[2, 2])
    list(
      fc = est$beta_fc[2] - 1, fch = est$beta[2] - 1, lasso_only = lasso_only$beta[2] - 1,
      ols = fgls:::lsq(X, s$y)[2] - 1,
      cov_fc = covered(est$beta_fc[2], se_fc, 1), cov_fch = covered(est$beta[2], se_fch, 1),
      cov_fch_nofs = covered(est$beta[2], sqrt(V_nc[2, 2]), 1),
      len_fc = se_fc, len_fch = se_fch
    )
  })
}

s4_cells <- list(
  "AR1-x2" = list(ar = 0.5, skedastic = "x2"),
  "AR1-log2" = list(ar = 0.5, skedastic = "log2"),
  "AR1-exp" = list(ar = 0.5, skedastic = "exp"),
  "AR2-x2" = list(ar = c(1.34, -0.42), skedastic = "x2"),
  "MA1-x2" = list(ar = numeric(0), ma = 0.5, skedastic = "x2"),
  "ARMA-x2" = list(ar = 0.8, ma = -0.4, skedastic = "x2"),
  "AR1-x2-mixture" = list(ar = 0.5, skedastic = "x2", driver = "mixture")
)

test_that("Table S.9: F-CH dominates F-C and covers", {
  skip_if_not_slow()
  for (cell in names(s4_cells)) {
    r <- do.call(s4_cell, s4_cells[[cell]])
    expect_lte(mse(r$fch), mse(r$fc), label = paste(cell, "F-CH vs F-C"))
    expect_lte(mse(r$fch), 1.05 * mse(r$lasso_only), label = paste(cell, "F-CH vs Lasso only"))
    expect_lt(mse(r$fch), mse(r$ols), label = paste(cell, "F-CH vs OLS"))
    expect_lt(abs(mean(r$fch)), 0.02, label = paste(cell, "bias"))
    expect_true(mean(r$cov_fch) >= 0.93 && mean(r$cov_fch) <= 0.97,
      label = paste(cell, "F-CH coverage", mean(r$cov_fch)))
    expect_true(mean(r$cov_fc) >= 0.91 && mean(r$cov_fc) <= 0.98,
      label = paste(cell, "F-C coverage", mean(r$cov_fc)))
    expect_lt(mean(r$len_fch), mean(r$len_fc), label = paste(cell, "interval length"))
  }
})

test_that("Table S.9: no loss under homoskedastic errors", {
  skip_if_not_slow()
  r <- s4_cell(ar = 0.5, skedastic = "homoskedastic")
  expect_lte(mse(r$fch), 1.05 * mse(r$fc))
  expect_true(mean(r$cov_fch) >= 0.93 && mean(r$cov_fch) <= 0.97)
})

test_that("open item A: the first-stage correction is needed with non-exogenous regressors", {
  skip_if_not_slow()
  for (gamma in c(0.5, 1)) {
    r <- s4_cell(ar = 0.8, skedastic = "x2", gamma = gamma, R = 800, seed = 77)
    with_fs <- mean(r$cov_fch)
    without_fs <- mean(r$cov_fch_nofs)
    expect_gte(with_fs, without_fs + 0.03, label = paste("gamma", gamma))
    expect_true(with_fs >= 0.91 && with_fs <= 0.97, label = paste("gamma", gamma, with_fs))
  }
})
