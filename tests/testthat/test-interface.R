sim <- simulate_serial(200, ar = 0.5, gamma = 0.5, burn = 50, seed = 123)

test_that("formula, data and subset behave like lm()", {
  fit <- serial_fgls(y ~ x, data = sim)
  expect_s3_class(fit, "fgls")
  expect_named(coef(fit), c("(Intercept)", "x"))
  expect_equal(formula(fit), y ~ x, ignore_attr = TRUE)
  expect_equal(dim(model.matrix(fit)), c(200L, 2L))
  expect_s3_class(model.frame(fit), "data.frame")
  expect_s3_class(terms(fit), "terms")
  expect_equal(nobs(fit), 200L - fit$ar_order)
  fs <- serial_fgls(y ~ x, data = sim, subset = seq_len(nrow(sim)) <= 150)
  fd <- serial_fgls(y ~ x, data = sim[1:150, ])
  expect_equal(coef(fs), coef(fd))
  expect_equal(serial_fgls(y ~ x, data = sim, subset = 1:150)$coefficients, coef(fd))
  y <- sim$y
  x <- sim$x
  fe <- serial_fgls(y ~ x)
  expect_equal(coef(fe), coef(serial_fgls(y ~ x, data = sim)))
  expect_equal(coef(serial_fgls(log(y + 10) ~ x + I(x^2), data = sim))[3],
    coef(stats::lm(log(y + 10) ~ x + I(x^2), data = sim))[3],
    tolerance = 0.5
  )
})

test_that("factors are expanded and predict() uses the stored levels", {
  df <- sim
  df$g <- factor(rep(c("a", "b", "c", "d"), length.out = 200))
  fit <- serial_fgls(y ~ x + g, data = df)
  expect_named(coef(fit), c("(Intercept)", "x", "gb", "gc", "gd"))
  nd <- data.frame(x = c(0, 1), g = factor(c("b", "d"), levels = c("a", "b", "c", "d")))
  p <- predict(fit, newdata = nd)
  expect_equal(unname(p), unname(c(
    coef(fit)[1] + coef(fit)["gb"],
    coef(fit)[1] + coef(fit)[2] + coef(fit)["gd"]
  )))
})

test_that("missing values follow na.action", {
  df <- sim
  df$y[c(1, 2, 200)] <- NA
  expect_error(serial_fgls(y ~ x, data = df), class = "fgls_error_missing")
  fo <- serial_fgls(y ~ x, data = df, na.action = stats::na.omit)
  expect_equal(fo$nobs, 197L)
  expect_equal(coef(fo), coef(serial_fgls(y ~ x, data = sim[3:199, ])))
  fx <- serial_fgls(y ~ x, data = df, na.action = stats::na.exclude)
  expect_length(residuals(fx), 200L)
  expect_true(all(is.na(residuals(fx)[c(1, 2, 200)])))
  expect_length(fitted(fx), 200L)
  df$y[100] <- NA
  expect_error(serial_fgls(y ~ x, data = df, na.action = stats::na.omit), class = "fgls_error_missing")
  skip("requires the heteroskedasticity engine (Phase 3)")
  fc <- fgls(y ~ x, data = df, serial = "none", hetero = "none", na.action = stats::na.omit)
  expect_equal(fc$nobs, 196L)
})

test_that("arguments are validated with classed errors", {
  expect_error(serial_fgls(y ~ x, data = sim, weights = rep(1, 200)), "not supported")
  expect_error(fgls(y ~ x, data = sim, hetero = "none", vcov_type = "naive"), "Unknown argument")
  expect_error(fgls(y ~ x, data = sim, serial = "bogus"), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, hetero = "none", vcov = "nonrobust"), class = "fgls_error_input")
  expect_error(fgls(~x, data = sim), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, control = list()), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, delta = -1), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, cv_folds = 1), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, psi_grid = -1), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, weight_clip = c(0.9, 0.1)), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, seed = "a"), class = "fgls_error_input")
  expect_error(fgls(y ~ x, data = sim, hetero = "known"), class = "fgls_error_input")
  expect_error(serial_fgls(y ~ x, data = 1:3), class = "fgls_error_input")
  expect_error(serial_fgls(g ~ x, data = transform(sim, g = factor(x > 0))), "numeric")
  expect_error(serial_fgls(y ~ x + z, data = transform(sim, z = 2 * x)), class = "fgls_error_collinearity")
  expect_error(fgls_control(n_boot = 1), class = "fgls_error_input")
  expect_error(fgls_control(transforms = "sqrt"), class = "fgls_error_input")
  expect_error(fgls_control(foldid = "a"), class = "fgls_error_input")
  expect_error(fgls_control(lambda_min_ratio = 2), class = "fgls_error_input")
  expect_error(fgls_control(n_lambda = 1), class = "fgls_error_input")
  expect_error(fgls_control(cores = 0), class = "fgls_error_input")
})

test_that("a model without intercept warns", {
  expect_warning(serial_fgls(y ~ x - 1, data = sim), class = "fgls_warning_specification")
})

test_that("confint() equals coef +/- qnorm * se, and use_t switches to Student-t", {
  fit <- serial_fgls(y ~ x, data = sim)
  se <- sqrt(diag(vcov(fit)))
  ci <- confint(fit, level = 0.9)
  expect_equal(unname(ci[, 1]), unname(coef(fit) - qnorm(0.95) * se))
  expect_equal(unname(ci[, 2]), unname(coef(fit) + qnorm(0.95) * se))
  expect_equal(colnames(ci), c("5 %", "95 %"))
  expect_equal(confint(fit, "x"), confint(fit)["x", , drop = FALSE])
  expect_equal(confint(fit, 2), confint(fit)["x", , drop = FALSE])
  ci_n <- confint(fit, type = "naive")
  expect_equal(unname(ci_n[, 2] - ci_n[, 1]), unname(2 * qnorm(0.975) * sqrt(diag(vcov(fit, "naive")))))
  ft <- serial_fgls(y ~ x, data = sim, use_t = TRUE)
  expect_equal(df.residual(ft), ft$n_eff - 2L)
  expect_equal(df.residual(fit), Inf)
  expect_equal(unname(confint(ft)[, 2]), unname(coef(ft) + qt(0.975, ft$n_eff - 2) * sqrt(diag(vcov(ft)))))
  expect_match(colnames(summary(ft)$coefficients)[3], "t value")
})

test_that("lmtest::coeftest() and coefci() reproduce summary() and confint()", {
  skip_if_not_installed("lmtest")
  fit <- serial_fgls(y ~ x, data = sim)
  ct <- lmtest::coeftest(fit)
  expect_equal(unclass(ct)[, 1:4], summary(fit)$coefficients, ignore_attr = TRUE)
  expect_equal(unclass(lmtest::coefci(fit)), unclass(confint(fit)), ignore_attr = TRUE)
  ft <- serial_fgls(y ~ x, data = sim, use_t = TRUE)
  expect_equal(unclass(lmtest::coeftest(ft))[, 4], summary(ft)$coefficients[, 4], ignore_attr = TRUE)
  expect_equal(unclass(lmtest::coeftest(fit, vcov. = vcov(fit, type = "naive")))[, 2],
    sqrt(diag(vcov(fit, type = "naive"))),
    ignore_attr = TRUE
  )
})

test_that("predict() returns fitted values, standard errors and intervals", {
  fit <- serial_fgls(y ~ x, data = sim)
  expect_equal(predict(fit), fitted(fit))
  nd <- data.frame(x = c(-1, 0, 2))
  p <- predict(fit, nd)
  expect_equal(unname(p), unname(cbind(1, nd$x) %*% coef(fit))[, 1])
  ps <- predict(fit, nd, se.fit = TRUE)
  Xn <- cbind(1, nd$x)
  expect_equal(unname(ps$se.fit), sqrt(rowSums((Xn %*% vcov(fit)) * Xn)))
  pi <- predict(fit, nd, interval = "confidence", level = 0.9)
  expect_equal(colnames(pi), c("fit", "lwr", "upr"))
  expect_equal(unname(pi[, "upr"] - pi[, "fit"]), unname(qnorm(0.95) * ps$se.fit))
  both <- predict(fit, nd, interval = "confidence", se.fit = TRUE)
  expect_named(both, c("fit", "se.fit", "df"))
  expect_error(predict(fit, data.frame(z = 1)))
})

test_that("print and summary output", {
  fit <- serial_fgls(y ~ x, data = sim)
  expect_output(print(fit), "AR order k\\* =")
  s <- summary(fit)
  expect_s3_class(s, "summary.fgls")
  expect_equal(dim(s$coefficients), c(2L, 4L))
  expect_equal(s$coefficients[, 2], sqrt(diag(vcov(fit))))
  expect_output(print(s), "Non-exogeneity share")
  expect_snapshot(print(summary(fit), digits = 4))
  expect_snapshot(print(fit, digits = 4))
})

test_that("ts, zoo and xts data and a time index are accepted", {
  tsd <- stats::ts(as.matrix(sim), start = c(1990, 1), frequency = 4)
  ft <- serial_fgls(y ~ x, data = tsd)
  expect_equal(coef(ft), coef(serial_fgls(y ~ x, data = sim)))
  expect_equal(ft$index[1], 1990)
  y <- tsd[, "y"]
  x <- tsd[, "x"]
  expect_equal(serial_fgls(y ~ x)$index[1:2], c(1990, 1990.25))
  dates <- seq(as.Date("2000-01-01"), by = "month", length.out = 200)
  df <- data.frame(date = dates, sim)
  fd <- serial_fgls(y ~ x, data = df, time = "date")
  expect_s3_class(fd$index, "Date")
  expect_error(serial_fgls(y ~ x, data = df[-50, ], time = "date"), class = "fgls_error_time_index")
  expect_warning(serial_fgls(y ~ x, data = df[-50, ], time = "date", allow_gaps = TRUE),
    class = "fgls_warning_specification"
  )
  expect_error(serial_fgls(y ~ x, data = df[200:1, ], time = "date"), class = "fgls_error_time_index")
  skip("requires the heteroskedasticity engine (Phase 3)")
  expect_silent(fgls(y ~ x, data = df[200:1, ], time = "date", serial = "none", hetero = "none"))
  expect_error(serial_fgls(y ~ x, data = df, time = "when"), class = "fgls_error_input")
  expect_error(serial_fgls(y ~ x, data = df, time = dates[-1]), class = "fgls_error_input")
  skip_if_not_installed("zoo")
  z <- zoo::zoo(as.matrix(sim), order.by = dates)
  fz <- serial_fgls(y ~ x, data = z)
  expect_equal(coef(fz), coef(ft))
  expect_error(serial_fgls(y ~ x, data = z[-30, ]), class = "fgls_error_time_index")
  skip_if_not_installed("xts")
  xx <- xts::xts(as.matrix(sim), order.by = dates)
  expect_equal(coef(serial_fgls(y ~ x, data = xx)), coef(ft))
})
