# Monte Carlo helpers for the statistical (paper-replication) tests. These tests
# are slow: they run only when the environment variable FGLS_SLOW_TESTS is
# "true" (nightly workflow) and never on CRAN.
skip_if_not_slow <- function() {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("FGLS_SLOW_TESTS"), "true"), "set FGLS_SLOW_TESTS=true")
}

run_mc <- function(R, seed, one) {
  out <- fgls:::with_seed(seed, lapply(seq_len(R), function(r) suppressWarnings(one())))
  keys <- names(out[[1]])
  stats::setNames(lapply(keys, function(k) vapply(out, function(o) as.numeric(o[[k]]), numeric(1))), keys)
}

rmse <- function(x) sqrt(mean(x^2))
mse <- function(x) mean(x^2)
covered <- function(est, se, true) abs(est - true) <= qnorm(0.975) * se
