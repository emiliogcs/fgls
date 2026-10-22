# Behaviour promised by the help pages and the vignettes, so that the text and
# the code do not drift apart.

test_that("test_predetermined() returns the documented components", {
  sim <- simulate_serial(200, seed = 1)
  te <- suppressWarnings(test_predetermined(y ~ x, data = sim, max_lag = 2, ar_order = 1))
  expect_named(te, c("table", "ar_order", "regressors", "alpha", "rejected", "advice", "scale"))
  expect_identical(te$scale, "innovations")
  expect_identical(te$ar_order, 1L)
  expect_named(te$table, c("lags", "LM", "df", "p.value"))
  expect_identical(suppressWarnings(test_predetermined(y ~ x, data = sim, max_lag = 1,
    ar_order = 1, scale = "errors"
  ))$scale, "errors")
})

test_that("the print methods of the diagnostics return their argument invisibly", {
  sim <- simulate_serial(200, seed = 1)
  te <- suppressWarnings(test_predetermined(y ~ x, data = sim, max_lag = 1, ar_order = 1))
  expect_output(v <- withVisible(print(te)), "Pre-determinedness LM test")
  expect_false(v$visible)
  expect_identical(v$value, te)
  a <- adf_test(sim$x)
  expect_output(v <- withVisible(print(a)), "ADF test")
  expect_false(v$visible)
  expect_identical(v$value, a)
})

test_that("weight_diagnostics(), ar_roots() and nonexogeneity_share() return the documented types", {
  sim <- simulate_serial(200, seed = 2)
  f0 <- serial_fgls(y ~ x, data = sim, ar_order = 0)
  expect_identical(ar_roots(f0), complex(0))
  f1 <- serial_fgls(y ~ x, data = sim, ar_order = 1)
  expect_true(is.complex(ar_roots(f1)))
  expect_length(ar_roots(f1), 1L)
  expect_true(is.numeric(nonexogeneity_share(f1)))
  expect_length(nonexogeneity_share(f1), 1L)
  het <- simulate_hetero(100, skedastic = "x2", seed = 3)
  wd <- weight_diagnostics(lasso_fgls(y ~ x, data = het, hetero = "wls_s2"))
  expect_true(is.numeric(wd))
  expect_named(wd, c("variance_ratio", "effective_n_share", "top1pct_weight_share"))
})

test_that("dropped Durbin lag columns are recorded, not signalled as collinearity warnings", {
  sim <- simulate_serial(200, seed = 4)
  sim$trend <- seq_len(200)
  sim$q <- factor(rep_len(1:4, 200))
  expect_no_warning(
    fit <- serial_fgls(y ~ x + trend + q, data = sim, ar_order = 2),
    class = "fgls_warning_collinearity"
  )
  expect_true(all(c("L1.(Intercept)", "L1.trend", "L2.trend") %in% fit$serial$dropped_lags))
  expect_match(summary(fit)$notes, "Collinear lagged columns dropped", all = FALSE)
})

test_that("fgls warnings can be muffled by class with suppressWarnings(classes =)", {
  sim <- simulate_serial(100, seed = 1)
  # no intercept: a specification warning
  expect_warning(serial_fgls(y ~ x - 1, data = sim), class = "fgls_warning_specification")
  expect_no_warning(suppressWarnings(serial_fgls(y ~ x - 1, data = sim), classes = "fgls_warning"))
})

test_that("exog = 'exogenous' only drops the first-stage correction of the combined variance", {
  skip_on_cran()
  sim <- simulate_combined(200, ar = 0.8, skedastic = "x2", gamma = 1, seed = 6)
  pre <- fgls(y ~ x, data = sim, seed = 1)
  exo <- fgls(y ~ x, data = sim, seed = 1, exog = "exogenous")
  nofs <- fgls(y ~ x, data = sim, seed = 1, first_stage_correction = FALSE)
  expect_equal(coef(exo), coef(pre))
  expect_equal(vcov(exo), vcov(nofs))
  expect_false(isTRUE(all.equal(vcov(exo), vcov(pre))))
})

test_that("the homoskedastic example of the combined vignette selects nothing", {
  skip_on_cran()
  hom <- simulate_combined(200, ar = 0.5, skedastic = "homoskedastic", seed = 4)
  fit <- fgls(y ~ x, data = hom, seed = 1)
  expect_length(fit$skedastic$selected, 0L)
  expect_equal(coef(fit), coef(fit$stages$serial))
})
