# R vs pyfgls distributional agreement (see vignettes/validation.Rmd): the
# data sets of tools/compare_python_r.R are regenerated in R from their seeds and
# refitted; the pyfgls outputs on the same data and folds are stored in
# inst/extdata/cross_language.rds. Slow: FGLS_SLOW_TESTS=true.

xl_designs <- list(
  serial_ar1_nonexog = function(s) simulate_serial(200, ar = 0.8, rho_x = 0.8, gamma = 0.5, burn = 100, seed = s),
  serial_ar2_exog = function(s) simulate_serial(200, ar = c(1.34, -0.42), rho_x = 0.8, gamma = 0, burn = 100, seed = s),
  hetero_log2 = function(s) simulate_hetero(200, skedastic = "log2", seed = s),
  combined_ar1_log2 = function(s) simulate_combined(200, ar = 0.5, skedastic = "log2", seed = s),
  combined_ar1_exp = function(s) simulate_combined(200, ar = 0.5, skedastic = "exp", seed = s)
)

xl_fit <- function(nm, sim, seed) {
  df <- data.frame(y = sim$y, x = sim$x)
  if (startsWith(nm, "serial")) {
    f <- serial_fgls(y ~ x, data = df)
  } else if (startsWith(nm, "hetero")) {
    x <- df$x
    Z <- cbind(x, log(x)^2, x^2, cos(x), cos(2 * x))
    f <- lasso_fgls(y ~ x, data = df, z = Z,
      control = fgls_control(foldid = fgls:::make_folds(200, 10L, seed = seed + 500000L)))
  } else {
    w <- sim$w
    k <- fgls:::estimate_serial(df$y, cbind(1, df$x), c("(Intercept)", "x"))$order
    f <- fgls(y ~ x, data = df, z = cbind(w, w^2, cos(w), cos(2 * w)),
      control = fgls_control(foldid = fgls:::make_folds(200 - k, 10L, seed = seed + 500000L)))
  }
  c(beta = unname(coef(f)[2]), se = sqrt(vcov(f)[2, 2]))
}

test_that("R and pyfgls agree on the same data sets (estimates, standard errors, coverage)", {
  skip_if_not_slow()
  path <- system.file("extdata", "cross_language.rds", package = "fgls")
  skip_if(path == "", "cross_language.rds not installed")
  res <- readRDS(path)
  for (nm in names(xl_designs)) {
    d <- res$designs[[nm]]
    idx <- seq_len(min(200L, length(d$seeds)))
    out <- t(vapply(idx, function(i) {
      suppressWarnings(xl_fit(nm, xl_designs[[nm]](d$seeds[i]), d$seeds[i]))
    }, numeric(2)))
    py <- d$py[idx, ]
    expect_lt(mean(abs(out[, "beta"] - py$beta)), 1e-8, label = paste(nm, "mean |beta_R - beta_Py|"))
    expect_lt(mean(abs(out[, "se"] - py$se)), 1e-4, label = paste(nm, "mean |se_R - se_Py|"))
    cov_r <- mean(abs(out[, "beta"] - 1) <= qnorm(0.975) * out[, "se"])
    cov_p <- mean(abs(py$beta - 1) <= qnorm(0.975) * py$se)
    expect_lt(abs(cov_r - cov_p), 0.01, label = paste(nm, "coverage difference"))
  }
})
