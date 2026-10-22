test_that("package loads", {
  expect_true(is.function(fgls:::.onLoad))
})
