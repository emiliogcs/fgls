test_that("lag_matrix and ar_filter follow the pyfgls conventions", {
  expect_equal(fgls:::lag_matrix(c(1, 2, 3, 4), 2), matrix(c(2, 3, 1, 2), 2))
  m <- cbind(1:5, 11:15)
  L <- fgls:::lag_matrix(m, 2)
  expect_equal(dim(L), c(3L, 4L))
  expect_equal(L[1, ], c(2, 12, 1, 11))
  expect_equal(dim(fgls:::lag_matrix(1:5, 0)), c(5L, 0L))
  expect_equal(fgls:::ar_filter(c(1, 2, 4, 8), 0.5), c(1.5, 3, 6))
  expect_equal(fgls:::ar_filter(cbind(1, c(1, 2, 4, 8)), 0.5), cbind(0.5, c(1.5, 3, 6)))
  expect_equal(fgls:::ar_filter(1:4, numeric(0)), 1:4)
  expect_error(fgls:::lag_matrix(1:3, 2, start = 1))
})

test_that("AR roots", {
  expect_equal(Mod(fgls:::ar_poly_roots(0.5)), 2)
  r <- fgls:::ar_poly_roots(c(1.34, -0.42))
  expect_equal(sort(Mod(r)), sort(Mod(polyroot(c(1, -1.34, 0.42)))))
  expect_length(fgls:::ar_poly_roots(numeric(0)), 0)
  expect_length(fgls:::ar_poly_roots(c(0.5, 0)), 1)
  expect_length(fgls:::ar_poly_roots(0), 0)
})

test_that("lag-order rules", {
  expect_equal(fgls:::default_kmax(200), 12L)
  expect_equal(fgls:::default_kmax(500), 15L)
  expect_equal(fgls:::theory_kmax(200), floor((200 / log(200))^(1 / 3)))
  expect_equal(fgls:::adjusted_min_lags(200), ceiling(log(200)^1.5))
})

test_that("critical values and p-values", {
  expect_equal(fgls:::critical_value(0.05), qnorm(0.975))
  expect_equal(fgls:::critical_value(0.05, 10), qt(0.975, 10))
  expect_error(fgls:::critical_value(1.5), class = "fgls_error_input")
  expect_equal(fgls:::two_sided_pvalues(1.96), 2 * pnorm(-1.96))
  expect_equal(fgls:::two_sided_pvalues(2, 5), 2 * pt(-2, 5))
})

test_that("with_seed is reproducible and restores the RNG state", {
  set.seed(99)
  before <- .Random.seed
  a <- fgls:::with_seed(1, rnorm(3))
  expect_identical(.Random.seed, before)
  b <- fgls:::with_seed(1, rnorm(3))
  expect_identical(a, b)
  expect_error(fgls:::check_seed("a"), class = "fgls_error_input")
})
