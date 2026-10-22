# Monte Carlo helpers for the statistical (paper-replication) tests. These tests
# are slow: they run only when the environment variable FGLS_SLOW_TESTS is
# "true" (.github/workflows/statistical.yaml, on demand or nightly) and never on
# CRAN. FGLS_SLOW_CORES sets the number of forked workers (default 1; the
# workflow uses the runner's core count). Results do not depend on it: every
# replication has its own seed.
skip_if_not_slow <- function() {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("FGLS_SLOW_TESTS"), "true"), "set FGLS_SLOW_TESTS=true")
}

run_mc <- function(R, seed, one) {
  seeds <- fgls:::with_seed(seed, sample.int(.Machine$integer.max, R))
  cores <- as.integer(Sys.getenv("FGLS_SLOW_CORES", "1"))
  f <- function(s) fgls:::with_seed(s, suppressWarnings(one()))
  out <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seeds, f, mc.cores = cores)
  } else {
    lapply(seeds, f)
  }
  keys <- names(out[[1]])
  stats::setNames(lapply(keys, function(k) {
    v <- lapply(out, function(o) as.numeric(o[[k]]))
    if (length(v[[1]]) == 1L) unlist(v) else do.call(rbind, v)
  }), keys)
}

rmse <- function(x) sqrt(mean(x^2))
mse <- function(x) mean(x^2)
covered <- function(est, se, true) abs(est - true) <= qnorm(0.975) * se
