test_that("the pre-determinedness test does not reject for pre-determined regressors", {
  sim <- simulate_serial(400, ar = 0.5, gamma = 0.5, seed = 1)
  expect_warning(t <- test_predetermined(y ~ x, data = sim), class = "fgls_warning_experimental")
  expect_s3_class(t, "fgls_predetermined")
  expect_equal(nrow(t$table), 4L)
  expect_identical(t$table$df, 1:4)
  expect_output(print(t), "Pre-determinedness")
  te <- suppressWarnings(test_predetermined(y ~ x, data = sim, scale = "errors", max_lag = 2, ar_order = 1))
  expect_equal(nrow(te$table), 2L)
  expect_error(suppressWarnings(test_predetermined(y ~ x, data = sim, regressors = "z")), class = "fgls_error_input")
  expect_error(suppressWarnings(test_predetermined(y ~ x, data = sim, max_lag = 0)), class = "fgls_error_input")
})

test_that("the pre-determinedness test rejects when a lag of x is missing", {
  set.seed(2)
  n <- 500
  x <- as.numeric(stats::filter(rnorm(n), 0.5, method = "recursive"))
  u <- fgls:::arma_errors(rnorm(n), 0.5)
  y <- 1 + x + 0.8 * c(0, x[-n]) + u
  df <- data.frame(y = y, x = x)
  t <- suppressWarnings(test_predetermined(y ~ x, data = df))
  expect_true(t$rejected)
  expect_match(t$advice, "Add lag")
  df$xl <- c(0, x[-n])
  t2 <- suppressWarnings(test_predetermined(y ~ x + xl, data = df[-1, ], regressors = "x", max_lag = 1))
  expect_false(t2$rejected)
})

test_that("ADF screening", {
  set.seed(3)
  rw <- cumsum(rnorm(300))
  expect_false(adf_test(rw)$unit_root_rejected)
  expect_true(adf_test(rnorm(300))$unit_root_rejected)
  expect_output(print(adf_test(rnorm(100))), "ADF test")
  expect_equal(adf_test(rnorm(100), max_lags = 0)$lags, 0L)
  df <- data.frame(y = rw + rnorm(300), x = rw)
  expect_warning(
    serial_fgls(y ~ x, data = df, check_stationarity = TRUE),
    class = "fgls_warning_nonstationary"
  )
})

test_that("weight diagnostics, AR roots and the non-exogeneity share", {
  fit <- fgls(y ~ x, data = simulate_combined(200, seed = 4), seed = 1)
  wd <- weight_diagnostics(fit)
  expect_named(wd, c("variance_ratio", "effective_n_share", "top1pct_weight_share"))
  expect_true(wd[["effective_n_share"]] > 0 && wd[["effective_n_share"]] <= 1)
  expect_true(all(Mod(ar_roots(fit)) > 1))
  expect_equal(nonexogeneity_share(fit), fit$nonexogeneity_share)
  s <- serial_fgls(y ~ x, data = simulate_serial(200, seed = 1))
  expect_error(weight_diagnostics(s), class = "fgls_error_input")
  h <- lasso_fgls(y ~ x, data = simulate_hetero(100, seed = 1), seed = 1)
  expect_error(ar_roots(h), class = "fgls_error_input")
  expect_error(nonexogeneity_share(h), class = "fgls_error_input")
})

test_that("formula edge cases", {
  sim <- simulate_serial(200, seed = 6)
  sim$g <- factor(rep(1:2, 100))
  f <- serial_fgls(y ~ x * g + I(x^2), data = sim)
  expect_named(coef(f), c("(Intercept)", "x", "g2", "I(x^2)", "x:g2"))
  fp <- serial_fgls(y ~ poly(x, 2), data = sim)
  expect_length(coef(fp), 3L)
  expect_length(predict(fp, newdata = data.frame(x = c(0, 1))), 2L)
  expect_error(serial_fgls(y ~ x + offset(x), data = sim), "Offsets")
  expect_error(serial_fgls(cbind(y, x) ~ g, data = sim), "single column")
  expect_error(serial_fgls(y ~ x, data = sim[, c("y", "x")], subset = "a"), class = "fgls_error_input")
  expect_error(serial_fgls(y ~ x, data = sim, subset = c(TRUE, FALSE)), class = "fgls_error_input")
  ft <- serial_fgls(y ~ x, data = sim, contrasts = list(g = "contr.sum"))
  expect_equal(coef(ft), coef(serial_fgls(y ~ x, data = sim)))
})
