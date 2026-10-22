# Model frame, data handling and argument validation (pre-submission audit,
# findings 15, 17, 22, 23, 24, 28 and 29).

test_that("subset combined with na.exclude pads residuals and fitted values like lm()", {
  sim <- simulate_serial(200, seed = 1)
  sim$y[51] <- NA
  f <- serial_fgls(y ~ x, data = sim, subset = 51:200, na.action = stats::na.exclude)
  l <- stats::lm(y ~ x, data = sim, subset = 51:200, na.action = stats::na.exclude)
  r <- residuals(f)
  expect_length(r, 150L)
  expect_equal(which(is.na(r)), which(is.na(residuals(l))))
  expect_equal(unname(which(is.na(r))), 1L)
  expect_equal(names(r), names(residuals(l)))
  expect_length(fitted(f), 150L)
  expect_equal(unname(r), unname(sim$y[51:200] - fitted(f)))

  het <- simulate_hetero(100, seed = 1)
  het$y[c(30, 60)] <- NA
  g <- lasso_fgls(y ~ x, data = het, subset = 21:100, na.action = stats::na.exclude, seed = 1)
  expect_equal(unname(which(is.na(residuals(g)))), c(10L, 40L))
  expect_equal(names(residuals(g)), as.character(21:100))
  expect_equal(unname(residuals(g)), unname(het$y[21:100] - fitted(g)))
  expect_equal(unname(which(is.na(predict(g)))), c(10L, 40L))
  # without subset the positions are unchanged
  h <- lasso_fgls(y ~ x, data = het, na.action = stats::na.exclude, seed = 1)
  expect_equal(unname(which(is.na(residuals(h)))), c(30L, 60L))
  expect_equal(as.integer(h$na.action), c(30L, 60L))
  expect_equal(as.integer(g$na.action), c(10L, 40L))
  expect_equal(names(g$na.action), c("30", "60"))
})

test_that("time-series subsets must be increasing blocks of consecutive rows", {
  sim <- simulate_serial(200, seed = 1)
  expect_error(
    serial_fgls(y ~ x, data = sim, subset = seq_len(200) %% 2 == 0),
    "skip rows", class = "fgls_error_time_index"
  )
  expect_error(
    serial_fgls(y ~ x, data = sim, subset = 200:1),
    "out of order", class = "fgls_error_time_index"
  )
  expect_error(
    serial_fgls(y ~ x, data = sim, subset = c(1:50, 1:50)),
    "out of order", class = "fgls_error_time_index"
  )
  expect_error(serial_fgls(y ~ x, data = sim, subset = -100), class = "fgls_error_time_index")
  expect_error(cochrane_orcutt(y ~ x, data = sim, subset = 200:1), class = "fgls_error_time_index")
  # allow_gaps = TRUE accepts skipped rows (positional lags) but not reversed ones
  expect_warning(
    f <- serial_fgls(y ~ x, data = sim, subset = -100, allow_gaps = TRUE),
    "skip rows", class = "fgls_warning_specification"
  )
  expect_equal(coef(f), coef(serial_fgls(y ~ x, data = sim[-100, ])))
  expect_error(
    serial_fgls(y ~ x, data = sim, subset = 200:1, allow_gaps = TRUE),
    class = "fgls_error_time_index"
  )
  # contiguous blocks, given as positions or as a logical vector, are fine
  ref <- coef(serial_fgls(y ~ x, data = sim[51:200, ]))
  expect_equal(coef(serial_fgls(y ~ x, data = sim, subset = 51:200)), ref)
  expect_equal(coef(serial_fgls(y ~ x, data = sim, subset = seq_len(200) > 50)), ref)
  expect_equal(coef(serial_fgls(y ~ x, data = sim, subset = -(1:50))), ref)
  # cross-sections can select rows in any order
  het <- simulate_hetero(100, seed = 2)
  expect_equal(
    coef(lasso_fgls(y ~ x, data = het, subset = seq_len(100) %% 2 == 0, hetero = "none")),
    coef(lasso_fgls(y ~ x, data = het[seq(2, 100, by = 2), ], hetero = "none"))
  )
  expect_s3_class(lasso_fgls(y ~ x, data = het, subset = 100:1, hetero = "none"), "fgls")
  # with a time index the index checks apply instead
  sim$date <- seq(as.Date("2000-01-01"), by = "month", length.out = 200)
  expect_error(
    serial_fgls(y ~ x, data = sim, time = "date", subset = -100),
    "time index has 1 gap", class = "fgls_error_time_index"
  )
})

test_that("missing and infinite values are located by rows of the data", {
  sim <- simulate_serial(200, seed = 1)
  s1 <- sim
  s1$y[60] <- NA
  expect_error(
    serial_fgls(y ~ x, data = s1, subset = 51:200),
    "found at rows 60\\.", class = "fgls_error_missing"
  )
  expect_error(
    serial_fgls(y ~ x, data = s1, subset = 51:200, na.action = stats::na.omit),
    "inside the sample at rows 60\\.", class = "fgls_error_missing"
  )
  s2 <- sim
  s2$y[1] <- NA
  s2$x[50] <- Inf
  expect_error(
    serial_fgls(y ~ x, data = s2, na.action = stats::na.omit),
    "Infinite values found at rows 50\\.", class = "fgls_error_missing"
  )
  het <- simulate_hetero(100, seed = 3)
  het$x[45] <- Inf
  expect_error(
    lasso_fgls(y ~ x, data = het, subset = 41:100, hetero = "none"),
    "rows 45\\.", class = "fgls_error_missing"
  )
  expect_equal(fgls:::format_locations(NULL, c(2, 3), rows = 11:20), "rows 12, 13")
  expect_equal(fgls:::format_locations(NULL, c(2, 3)), "rows 2, 3")
})

test_that("time indexes with missing values or character dates give classed errors", {
  sim <- simulate_serial(120, seed = 1)
  sim$tt <- seq_len(120)
  sim$tt[5] <- NA
  expect_error(
    serial_fgls(y ~ x, data = sim, time = "tt"),
    "missing values at rows 5\\.", class = "fgls_error_time_index"
  )
  expect_equal(
    coef(serial_fgls(y ~ x, data = sim, time = "tt", subset = 6:120)),
    coef(serial_fgls(y ~ x, data = sim[6:120, ]))
  )
  expect_error(fgls:::check_time_index(c(1, NA, 3), rows = 11:13), "rows 12\\.",
    class = "fgls_error_time_index"
  )
  months <- sprintf("%d-%02d", 2000 + (0:119) %/% 12, (0:119) %% 12 + 1)
  sim$tt <- months
  expect_no_warning(expect_error(
    serial_fgls(y ~ x, data = sim, time = "tt"),
    "as.Date", class = "fgls_error_time_index"
  ))
  sim$tt <- factor(months)
  expect_error(serial_fgls(y ~ x, data = sim, time = "tt"), "factor",
    class = "fgls_error_time_index"
  )
  # cross-sections do not check the index
  expect_s3_class(fgls(y ~ x, data = sim, time = "tt", serial = "none", hetero = "none"), "fgls")
  skip_if_not_installed("zoo")
  sim$tt <- zoo::as.yearmon(months)
  expect_equal(
    coef(serial_fgls(y ~ x, data = sim, time = "tt")),
    coef(serial_fgls(y ~ x, data = sim))
  )
})

test_that("known variances follow subset and na.action like the data", {
  het <- simulate_hetero(150, seed = 1)
  v <- het$x^2
  ref <- lasso_fgls(y ~ x, data = het[1:100, ], hetero = "known", variance = v[1:100])
  a <- lasso_fgls(y ~ x, data = het, hetero = "known", variance = v, subset = 1:100)
  expect_equal(coef(a), coef(ref))
  expect_equal(vcov(a), vcov(ref))
  # a vector for the estimation sample is still accepted
  b <- lasso_fgls(y ~ x, data = het, hetero = "known", variance = v[1:100], subset = 1:100)
  expect_equal(coef(b), coef(ref))
  # a function is applied to the model matrix of the sample
  f <- lasso_fgls(y ~ x, data = het, hetero = "known", variance = function(X) X[, 2]^2,
    subset = 1:100
  )
  expect_equal(coef(f), coef(ref))
  # rows reordered by subset in a cross-section: the variances follow the rows
  full <- lasso_fgls(y ~ x, data = het, hetero = "known", variance = v)
  expect_equal(
    coef(lasso_fgls(y ~ x, data = het, hetero = "known", variance = v, subset = 150:1)),
    coef(full)
  )
  hm <- het
  hm$y[c(3, 140)] <- NA
  d <- lasso_fgls(y ~ x, data = hm, hetero = "known", variance = v, na.action = stats::na.omit)
  e <- lasso_fgls(y ~ x, data = het[-c(3, 140), ], hetero = "known", variance = v[-c(3, 140)])
  expect_equal(coef(d), coef(e))
  expect_error(
    lasso_fgls(y ~ x, data = het, hetero = "known", variance = 1:3, subset = 1:100),
    "expected 150 \\(one per row of the data\\) or 100", class = "fgls_error_input"
  )
  expect_error(
    lasso_fgls(y ~ x, data = het, hetero = "known", variance = function(X) 1:3),
    "returned 3 values", class = "fgls_error_input"
  )
  # combined estimator (serial correction + known variances)
  comb <- simulate_combined(200, ar = 0.5, skedastic = "x2", seed = 3)
  vc <- attr(comb, "info")$v
  expect_equal(
    coef(fgls(y ~ x, data = comb, hetero = "known", variance = vc, subset = 21:200)),
    coef(fgls(y ~ x, data = comb[21:200, ], hetero = "known", variance = vc[21:200]))
  )
})
