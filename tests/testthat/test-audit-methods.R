# Regression tests for the pre-submission audit of the S3 methods, the result
# objects, the broom tidiers, plot(), summary() and compare_fgls().

sim <- simulate_serial(200, seed = 1)
comb <- simulate_combined(200, ar = 0.5, skedastic = "x2", seed = 3)

wallis_data <- function(n = 400, seed = 7) {
  set.seed(seed)
  w <- rnorm(n + 1)
  u <- fgls:::arma_errors(rnorm(n + 1), 0.5)
  yy <- numeric(n + 1)
  for (t in 2:(n + 1)) yy[t] <- 0.3 * yy[t - 1] + 1 + w[t] + u[t]
  data.frame(y = yy[-1], ylag = yy[-(n + 1)], w = w[-1])
}

test_that("the F-C stage supports model.frame(), augment() and predict() (audit 13, 21)", {
  fit <- fgls(y ~ x, data = comb, seed = 1)
  fc <- fit$stages$serial
  expect_identical(model.frame(fc), model.frame(fit))
  expect_identical(fc$na.action, fit$na.action)
  expect_equal(predict(fc), fitted(fc))
  expect_length(predict(fc, se.fit = TRUE)$se.fit, 200L)
  skip_if_not_installed("broom")
  au <- broom::augment(fc)
  expect_equal(nrow(au), 200L)
  expect_equal(au$.fitted, unname(fitted(fc)))
  expect_equal(sum(is.na(au$.innovation)), fc$ar_order)
  expect_false(".weight" %in% names(au))
  fm <- fgls(y ~ x, data = comb, seed = 1, model = FALSE)
  expect_error(model.frame(fm$stages$serial), "model = TRUE", class = "fgls_error_input")
})

test_that("vcov(type =) recomputes the F-C stage variance (audit 13, 30)", {
  fit <- fgls(y ~ x, data = comb, seed = 1)
  fc <- fit$stages$serial
  expect_identical(fc$vcov_code, "robust")
  expect_identical(vcov(fc, type = "robust"), vcov(fc))
  for (ct in c("HC0", "HC3")) {
    ref <- fgls(y ~ x, data = comb, seed = 1, vcov = ct)$stages$serial
    expect_equal(vcov(fc, type = ct), vcov(ref), tolerance = 1e-12)
    expect_equal(confint(fc, type = ct), confint(ref), tolerance = 1e-12)
    # and back to the robust variance from a stage fitted with another type
    expect_equal(vcov(ref, type = "robust"), vcov(fc), tolerance = 1e-12)
  }
  # without the first-stage correction (C = NULL)
  fe <- fgls(y ~ x, data = comb, seed = 1, exog = "exogenous")
  fe3 <- fgls(y ~ x, data = comb, seed = 1, exog = "exogenous", vcov = "HC3")
  expect_equal(vcov(fe$stages$serial, type = "HC3"), vcov(fe3$stages$serial), tolerance = 1e-12)
  # the final fit is unchanged
  expect_identical(vcov(fit, type = "robust"), vcov(fit))
  expect_error(vcov(fc, type = "naive"), "\"robust\", \"HC0\", \"HC3\"", class = "fgls_error_input")
  expect_error(vcov(fc, type = "bootstrap"), "F-C stage", class = "fgls_error_input")
})

test_that("benchmark fits return their own variance type (audit 30)", {
  o <- ols_robust(y ~ x, data = sim)
  expect_match(o$vcov_type, "^HAC \\(QS")
  expect_identical(o$vcov_code, "HAC")
  expect_identical(vcov(o, type = "HAC"), vcov(o))
  expect_identical(confint(o, type = "HAC"), confint(o))
  expect_error(vcov(o, type = "HC3"), "fgls\\(\\), serial_fgls\\(\\) or lasso_fgls\\(\\)",
    class = "fgls_error_input"
  )
  expect_error(confint(o, type = "NW"), "type = \"HAC\"", class = "fgls_error_input")
  nw <- ols_robust(y ~ x, data = sim, vcov = "NW")
  expect_identical(vcov(nw, type = "NW"), vcov(nw))
  # the code is the requested type even when the kernel is overridden
  ob <- ols_robust(y ~ x, data = sim, vcov = "HAC", kernel = "bartlett")
  expect_match(ob$vcov_type, "Bartlett")
  expect_identical(vcov(ob, type = "HAC"), vcov(ob))
  expect_error(vcov(ob, type = "NW"), class = "fgls_error_input")
  h1 <- ols_robust(y ~ x, data = sim, vcov = "HC1")
  expect_identical(vcov(h1, type = "HC1"), vcov(h1))
  co <- cochrane_orcutt(y ~ x, data = sim)
  expect_identical(vcov(co, type = "nonrobust"), vcov(co))
  expect_error(vcov(co, type = "robust"), class = "fgls_error_input")
  # objects saved before the code was stored derive it from the label
  o$vcov_code <- NULL
  expect_identical(vcov(o, type = "HAC"), vcov(o))
})

test_that("the variance code is kept apart from the label for bootstrap fits", {
  fit <- serial_fgls(y ~ x, data = sim)
  fb <- serial_fgls(y ~ x, data = sim, vcov = "bootstrap", control = fgls_control(n_boot = 20))
  expect_match(fb$vcov_type, "^bootstrap \\(sieve")
  expect_identical(fb$vcov_code, "bootstrap")
  expect_identical(vcov(fb, type = "bootstrap"), vcov(fb))
  expect_equal(vcov(fb, type = "robust"), vcov(fit), tolerance = 1e-12)
  expect_identical(fit$vcov_code, "robust")
  expect_identical(serial_fgls(y ~ x, data = sim, vcov = "naive")$vcov_code, "naive")
})

test_that("predict() pads standard errors and intervals for na.exclude (audit 20)", {
  d <- sim
  d$y[1:2] <- NA
  f <- serial_fgls(y ~ x, data = d, na.action = stats::na.exclude)
  l <- stats::lm(y ~ x, data = d, na.action = stats::na.exclude)
  p <- predict(f)
  ps <- predict(f, se.fit = TRUE)
  pi <- predict(f, interval = "confidence")
  pb <- predict(f, interval = "confidence", se.fit = TRUE)
  expect_length(p, 200L)
  expect_length(ps$fit, 200L)
  expect_length(ps$se.fit, 200L)
  expect_equal(nrow(pi), 200L)
  expect_equal(nrow(pb$fit), 200L)
  expect_length(pb$se.fit, 200L)
  expect_identical(unname(which(is.na(ps$se.fit))), which(is.na(predict(l, se.fit = TRUE)$se.fit)))
  expect_identical(unname(which(is.na(pi[, "lwr"]))), 1:2)
  expect_equal(ps$fit, p)
  expect_equal(pi[, "fit"], p)
  expect_identical(rownames(pi), rownames(predict(l, interval = "confidence")))
  # na.omit is not padded
  fo <- serial_fgls(y ~ x, data = d, na.action = stats::na.omit)
  expect_length(predict(fo, se.fit = TRUE)$se.fit, 198L)
})

test_that("predict() rejects unsupported types (audit 29)", {
  f <- serial_fgls(y ~ x, data = sim)
  expect_identical(predict(f, type = "response"), predict(f))
  expect_error(predict(f, type = "terms"), "type = \"response\"", class = "fgls_error_input")
  expect_error(predict(f, newdata = data.frame(x = 1), type = "link"), class = "fgls_error_input")
})

test_that("augment() handles na.exclude data, model = FALSE and se_fit/interval (audit 21, 29)", {
  skip_if_not_installed("broom")
  d <- sim
  d$y[c(1, 2, 200)] <- NA
  f <- serial_fgls(y ~ x, data = d, na.action = stats::na.exclude)
  a <- broom::augment(f, data = d)
  expect_equal(nrow(a), 200L)
  expect_identical(which(is.na(a$.fitted)), c(1L, 2L, 200L))
  expect_equal(a$.fitted, unname(fitted(f)))
  expect_equal(a$.resid, unname(residuals(f)))
  expect_equal(a$y, d$y)
  expect_equal(nrow(broom::augment(f)), 197L)
  fo <- serial_fgls(y ~ x, data = d, na.action = stats::na.omit)
  expect_error(broom::augment(fo, data = d), "one row per observation", class = "fgls_error_input")
  fm <- serial_fgls(y ~ x, data = sim, model = FALSE)
  expect_error(broom::augment(fm), "model = FALSE", class = "fgls_error_input")
  expect_equal(nrow(broom::augment(fm, data = sim)), 200L)

  ac <- broom::augment(f, se_fit = TRUE, interval = "confidence", conf.level = 0.9)
  pr <- predict(f, se.fit = TRUE, interval = "confidence", level = 0.9)
  ok <- !is.na(pr$se.fit)
  expect_named(ac, c("y", "x", ".fitted", ".lower", ".upper", ".se.fit", ".resid", ".innovation"))
  expect_equal(ac$.se.fit, unname(pr$se.fit[ok]))
  expect_equal(ac$.lower, unname(pr$fit[ok, "lwr"]))
  expect_equal(ac$.upper, unname(pr$fit[ok, "upr"]))
  ad <- broom::augment(f, data = d, se_fit = TRUE)
  expect_equal(ad$.se.fit, unname(pr$se.fit))
  expect_false(".lower" %in% names(ad))
  nd <- data.frame(x = c(-1, 0, 1))
  an <- broom::augment(f, newdata = nd, se_fit = TRUE, interval = "confidence")
  pn <- predict(f, newdata = nd, se.fit = TRUE, interval = "confidence")
  expect_equal(an$.fitted, unname(pn$fit[, "fit"]))
  expect_equal(an$.lower, unname(pn$fit[, "lwr"]))
  expect_equal(an$.se.fit, unname(pn$se.fit))
  expect_named(broom::augment(f, newdata = nd, se_fit = TRUE), c("x", ".fitted", ".se.fit"))
  expect_error(broom::augment(f, interval = "prediction"), class = "fgls_error_input")
  # heteroskedasticity weights are padded too
  het <- simulate_hetero(100, seed = 4)
  het$y[3] <- NA
  fh <- lasso_fgls(y ~ x, data = het, seed = 1, na.action = stats::na.exclude)
  ah <- broom::augment(fh, data = het)
  expect_equal(nrow(ah), 100L)
  expect_true(is.na(ah$.weight[3]))
  expect_equal(ah$.weight[-3], unname(fh$weights))
})

test_that("plot() accepts main, xlab and ylab (audit 25)", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  f <- serial_fgls(y ~ x, data = sim)
  expect_invisible(plot(f, which = 1, main = "custom"))
  expect_invisible(plot(f, xlab = "foo", ylab = "bar", main = "M", col = 2))
  fb <- fgls(y ~ x, data = comb, seed = 1)
  expect_invisible(plot(fb, main = "both", pch = 20))
  fh <- lasso_fgls(y ~ x, data = simulate_hetero(100, seed = 4), seed = 1)
  expect_invisible(plot(fh, which = 3:4, main = "hetero", ylab = "v"))
})

test_that("the Wallis note names the regressor and is absent when k* = 0 (audit 33)", {
  df <- wallis_data()
  f0 <- serial_fgls(y ~ ylag + w, data = df, lagged_dep_vars = "ylag", ar_order = 0)
  expect_identical(f0$ar_order, 0L)
  expect_false(any(grepl("Wallis", f0$notes)))
  f1 <- serial_fgls(y ~ ylag + w, data = df, lagged_dep_vars = "ylag", ar_order = 1)
  expect_match(f1$notes, "Wallis ratio from regressor\\(s\\) w ", all = FALSE)
})

test_that("compare_fgls() rejects lm fits it cannot convert (audit 19)", {
  f <- serial_fgls(y ~ x, data = sim)
  d <- sim
  d$b <- as.numeric(d$y > stats::median(d$y))
  wls <- stats::lm(y ~ x, data = d, weights = seq(0.1, 5, length.out = 200))
  expect_error(compare_fgls(f, WLS = wls), "`WLS` is a weighted lm fit.*ols_robust",
    class = "fgls_error_input"
  )
  expect_error(compare_fgls(f, stats::lm(y ~ x + offset(0.5 * x), data = d)), "offset",
    class = "fgls_error_input"
  )
  expect_error(compare_fgls(f, stats::lm(y ~ x, offset = 0.5 * x, data = d)), "offset",
    class = "fgls_error_input"
  )
  expect_error(compare_fgls(f, Logit = stats::glm(b ~ x, data = d, family = stats::binomial)),
    "glm", class = "fgls_error_input"
  )
  expect_error(compare_fgls(f, stats::glm(y ~ x, data = d)), "glm", class = "fgls_error_input")
  expect_error(compare_fgls(f, stats::lm(cbind(y, b) ~ x, data = d)), "mlm",
    class = "fgls_error_input"
  )
  # plain lm() and aov() fits are still converted
  tab <- compare_fgls(f, OLS = stats::lm(y ~ x, data = d), AOV = stats::aov(y ~ x, data = d))
  expect_equal(tab$OLS, tab$AOV)
})

test_that("summary() labels the scale of the Lasso penalty (audit 16)", {
  fit <- lasso_fgls(y ~ x, data = simulate_hetero(200, skedastic = "log2", seed = 10), seed = 1)
  h <- fit$skedastic
  expect_equal(h$lambda_paper, 200 * h$lambda)
  val <- summary(fit)$info[["Lasso psi, lambda"]]
  expect_match(val, "glmnet scale", fixed = TRUE)
  expect_match(val, sprintf("%.4g in GP24 eq. 4", h$lambda_paper), fixed = TRUE)
  expect_match(val, sprintf("%g, %.4g", h$psi, h$lambda), fixed = TRUE)
  expect_output(print(summary(fit)), "glmnet scale")
})

test_that("a non-exogeneity share of NA is shown as missing (with audit 26)", {
  fit <- serial_fgls(y ~ x, data = sim)
  fit$serial$nonexogeneity_share <- NA_real_
  fit$nonexogeneity_share <- NA_real_
  expect_identical(summary(fit)$info[["Non-exogeneity share"]], "-")
  expect_output(print(summary(fit)), "Non-exogeneity share: -")
  expect_output(print(fit), "AR order")
  expect_equal(nrow(compare_fgls(fit, ols_robust(y ~ x, data = sim))), 6L)
  info <- fit$serial
  expect_silent(fgls:::serial_notes(info, "robust", "exogenous"))
  skip_if_not_installed("broom")
  expect_true(is.na(broom::glance(fit)$nonexogeneity.share))
  expect_true(is.na(broom::glance(ols_robust(y ~ x, data = sim))$nonexogeneity.share))
})

test_that("print() returns its argument invisibly (audit 3)", {
  fit <- serial_fgls(y ~ x, data = sim)
  out <- NULL
  expect_output(out <- withVisible(print(fit)), "Coefficients")
  expect_false(out$visible)
  expect_identical(out$value, fit)
})
