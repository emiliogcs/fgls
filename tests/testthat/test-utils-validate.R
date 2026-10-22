test_that("prepare_data validates the model matrix", {
  set.seed(3)
  X <- cbind("(Intercept)" = 1, x = rnorm(20))
  d <- fgls:::prepare_data(rnorm(20), X)
  expect_equal(d$const_col, 1L)
  expect_true(d$has_const)
  expect_error(
    fgls:::prepare_data(rnorm(20), cbind(X, z = 2 * X[, 2])),
    class = "fgls_error_collinearity"
  )
  expect_error(fgls:::prepare_data(rnorm(20), cbind(X, z = 0)), class = "fgls_error_collinearity")
  expect_error(fgls:::prepare_data(rnorm(20), cbind(X, c = 2)), class = "fgls_error_collinearity")
  expect_error(fgls:::prepare_data(rnorm(2), X[1:2, ]), class = "fgls_error_insufficient_data")
  expect_error(fgls:::prepare_data(rnorm(19), X), class = "fgls_error_input")
  yi <- rnorm(20)
  yi[3] <- Inf
  expect_error(fgls:::prepare_data(yi, X), class = "fgls_error_missing")
  expect_warning(fgls:::prepare_data(rnorm(20), X[, 2, drop = FALSE]),
    class = "fgls_warning_specification"
  )
  dum <- cbind(a = rep(0:1, 10), b = rep(1:0, 10))
  expect_silent(d2 <- fgls:::prepare_data(rnorm(20), dum))
  expect_true(d2$has_implicit_const)
})

test_that("time index checks", {
  dates <- seq(as.Date("2020-01-01"), by = "month", length.out = 12)
  expect_silent(fgls:::check_time_index(dates))
  expect_error(fgls:::check_time_index(rev(dates)), class = "fgls_error_time_index")
  expect_error(fgls:::check_time_index(dates[c(1, 1:11)]), class = "fgls_error_time_index")
  expect_error(fgls:::check_time_index(dates[-5]), class = "fgls_error_time_index")
  expect_warning(fgls:::check_time_index(dates[-5], allow_gaps = TRUE),
    class = "fgls_warning_specification"
  )
  days <- seq(as.Date("2021-03-01"), by = "day", length.out = 30)
  bdays <- days[as.POSIXlt(days)$wday %in% 1:5]
  expect_silent(fgls:::check_time_index(bdays))
  expect_error(fgls:::check_time_index(bdays[-7]), class = "fgls_error_time_index")
  expect_silent(fgls:::check_time_index(seq(1990, 2000, by = 0.25)))
  expect_error(fgls:::check_time_index(c(1, 2, 3, 5, 6)), class = "fgls_error_time_index")
  expect_silent(fgls:::check_time_index(NULL))
})

test_that("missing-value policies", {
  na_rows <- c(TRUE, FALSE, FALSE, FALSE, TRUE)
  expect_error(fgls:::missing_rows_keep(na_rows, "fail", TRUE, NULL), class = "fgls_error_missing")
  expect_equal(fgls:::missing_rows_keep(na_rows, "omit", TRUE, NULL), !na_rows)
  inner <- c(FALSE, TRUE, FALSE)
  expect_error(fgls:::missing_rows_keep(inner, "omit", TRUE, NULL), class = "fgls_error_missing")
  expect_equal(fgls:::missing_rows_keep(inner, "omit", FALSE, NULL), !inner)
  expect_equal(fgls:::na_policy(stats::na.omit), "omit")
  expect_equal(fgls:::na_policy("na.exclude"), "exclude")
  expect_equal(fgls:::na_policy(NULL), "fail")
  expect_error(fgls:::na_policy(function(x) x), class = "fgls_error_input")
})

test_that("lagged dependent variables are detected", {
  set.seed(4)
  y <- cumsum(rnorm(50))
  X <- cbind(1, c(NA, y[-50]), rnorm(50))
  X[1, 2] <- 0
  colnames(X) <- c("(Intercept)", "ylag", "x")
  found <- fgls:::detect_lagged_dependent(y[2:50], X[2:50, ], colnames(X))
  expect_equal(found[[1]]$name, "ylag")
  expect_equal(found[[1]]$lag, 1L)
  expect_warning(fgls:::warn_lagged_dependent(found), class = "fgls_warning_lagged_dependent")
  expect_equal(fgls:::resolve_columns(c("x", "ylag"), colnames(X), "z"), c(3L, 2L))
  expect_error(fgls:::resolve_columns("nope", colnames(X), "z"), class = "fgls_error_input")
})

test_that("arg_choice errors are classed", {
  expect_equal(fgls:::arg_default(c("a", "b"), c("a", "b"), "opt"), "a")
  expect_equal(fgls:::arg_default("b", c("a", "b"), "opt"), "b")
  expect_error(fgls:::arg_choice("c", c("a", "b"), "opt"), class = "fgls_error_input")
})
