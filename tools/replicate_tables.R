#!/usr/bin/env Rscript
# Monte Carlo replication of the main simulation tables of the papers, for the
# validation vignette (results stored in inst/extdata/replication.rds; never
# recomputed during R CMD check).
#
#   Rscript tools/replicate_tables.R [cores]
#
# Replication counts are smaller than the papers' (10,000), so the Monte Carlo
# standard errors are larger: about 0.007 for a coverage rate near 0.95 with
# 1,000 replications.

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(fgls)
  }
})
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args)) as.integer(args[1]) else 1L

run_mc <- function(R, seed, one) {
  seeds <- fgls:::with_seed(seed, sample.int(.Machine$integer.max, R))
  f <- function(s) fgls:::with_seed(s, suppressWarnings(one()))
  out <- if (cores > 1L) parallel::mclapply(seeds, f, mc.cores = cores) else lapply(seeds, f)
  keys <- names(out[[1]])
  stats::setNames(lapply(keys, function(k) vapply(out, function(o) as.numeric(o[[k]]), numeric(1))), keys)
}
rmse <- function(x) sqrt(mean(x^2))
covered <- function(b, se, true) abs(b - true) <= qnorm(0.975) * se
lsq <- fgls:::lsq
arf <- fgls:::ar_filter
t0 <- Sys.time()

# ------------------------------------------------------------------ PG26 Table 1
message("PG26 Table 1")
r <- run_mc(2000, 20260101, function() {
  s <- simulate_serial(500)
  X <- cbind(1, s$x)
  est <- fgls:::estimate_serial(s$y, X, c("c", "x"), kmax = 12)
  u <- drop(s$y - X %*% lsq(X, s$y))
  rho_co <- sum(u[-1] * u[-500]) / sum(u[-500]^2)
  list(
    ols = lsq(X, s$y)[2] - 1, durbin = est$durbin$beta[2] - 1,
    gls = lsq(arf(X, 0.8), arf(s$y, 0.8))[2] - 1, fgls = est$beta[2] - 1,
    co = lsq(arf(X, rho_co), arf(s$y, rho_co))[2] - 1,
    rho_fgls = sum(est$rho) - 0.8, rho_co = rho_co - 0.8, k = est$order
  )
})
pg26_table1 <- data.frame(
  estimator = c("OLS", "Durbin", "GLS (infeasible)", "FGLS", "Cochrane-Orcutt", "rho (FGLS)", "rho (CO)"),
  rmse = c(rmse(r$ols), rmse(r$durbin), rmse(r$gls), rmse(r$fgls), rmse(r$co), rmse(r$rho_fgls), rmse(r$rho_co)),
  bias = c(mean(r$ols), mean(r$durbin), mean(r$gls), mean(r$fgls), mean(r$co), mean(r$rho_fgls), mean(r$rho_co)),
  paper_rmse = c(0.399, 0.036, 0.025, 0.025, 0.041, 0.034, NA),
  paper_bias = c(0.395, NA, NA, 0.001, 0.031, -0.008, -0.171)
)
pg26_table1_k <- stats::median(r$k)

# ---------------------------------------------------------- PG26 Tables 2-5 cells
message("PG26 Tables 2-5")
cells <- list("AR(1) 0.5" = 0.5, "AR(1) 0.8" = 0.8, "AR(2) (1.34, -0.42)" = c(1.34, -0.42), "AR(2) (0.5, -0.3)" = c(0.5, -0.3))
pg26_general <- do.call(rbind, lapply(names(cells), function(cell) {
  do.call(rbind, lapply(c(0, 0.25, 0.5), function(gamma) {
    ar <- cells[[cell]]
    r <- run_mc(1000, 4242, function() {
      s <- simulate_serial(200, ar = ar, rho_x = 0.8, gamma = gamma, burn = 100)
      X <- cbind(1, s$x)
      est <- fgls:::estimate_serial(s$y, X, c("c", "x"), kmax = 12)
      b <- est$beta[2]
      v <- est$variance
      list(
        ols = lsq(X, s$y)[2] - 1, fgls = b - 1, gls = lsq(arf(X, ar), arf(s$y, ar))[2] - 1,
        cov_r = covered(b, sqrt(est$cov[2, 2]), 1), cov_n = covered(b, sqrt(v$naive_avar[2, 2] / v$n), 1),
        k = est$order
      )
    })
    data.frame(
      errors = cell, gamma = gamma, mse_ols = 100 * mean(r$ols^2), mse_fgls = 100 * mean(r$fgls^2),
      mse_gls = 100 * mean(r$gls^2), coverage_robust = mean(r$cov_r), coverage_naive = mean(r$cov_n),
      median_k = stats::median(r$k)
    )
  }))
}))

# ------------------------------------------------------------------- GP24 Table 1
message("GP24 Table 1")
paper_z <- function(x) cbind(x, log(x)^2, x^2, cos(x), cos(2 * x))
gp24 <- c(homoskedastic = 1.00, x = 0.92, x2 = 0.70, log2 = 0.46, exp = 0.45, step = 0.735)
gp24_table1 <- do.call(rbind, lapply(names(gp24), function(sk) {
  r <- run_mc(1000, 2024, function() {
    s <- simulate_hetero(400, skedastic = sk)
    X <- cbind(1, s$x)
    est <- fgls:::estimate_hetero(s$y, X, paper_z(s$x), letters[1:5], seed = sample.int(.Machine$integer.max, 1))
    list(fgls = est$beta[2] - 1, ols = lsq(X, s$y)[2] - 1, cover = covered(est$beta[2], sqrt(est$cov[2, 2]), 1))
  })
  data.frame(
    skedastic = sk, T = 400, emse_ratio = mean(r$fgls^2) / mean(r$ols^2), paper_ratio = gp24[[sk]],
    coverage = mean(r$cover)
  )
}))

# --------------------------------------------------------------- Supplement S.9
message("Supplement Table S.9 cells")
s9_cells <- list(
  "AR(1) 0.5, x^2" = list(ar = 0.5, skedastic = "x2"),
  "AR(1) 0.5, [log x]^2" = list(ar = 0.5, skedastic = "log2"),
  "AR(1) 0.5, exp" = list(ar = 0.5, skedastic = "exp"),
  "AR(2) (1.34, -0.42), x^2" = list(ar = c(1.34, -0.42), skedastic = "x2"),
  "MA(1) 0.5, x^2" = list(ar = numeric(0), ma = 0.5, skedastic = "x2"),
  "ARMA(1,1) (0.8, -0.4), x^2" = list(ar = 0.8, ma = -0.4, skedastic = "x2"),
  "AR(1) 0.5, x^2, mixture" = list(ar = 0.5, skedastic = "x2", driver = "mixture"),
  "AR(1) 0.5, homoskedastic" = list(ar = 0.5, skedastic = "homoskedastic")
)
s9 <- do.call(rbind, lapply(names(s9_cells), function(cell) {
  a <- s9_cells[[cell]]
  r <- run_mc(1000, 31, function() {
    s <- do.call(simulate_combined, c(list(n = 200), a))
    w <- s$w
    X <- cbind(1, s$x)
    est <- fgls:::estimate_combined(s$y, X, cbind(w, w^2, cos(w), cos(2 * w)), letters[1:4],
      x_names = c("c", "x"), kmax = 12, seed = sample.int(.Machine$integer.max, 1)
    )
    list(
      ols = lsq(X, s$y)[2] - 1, fc = est$beta_fc[2] - 1, fch = est$beta[2] - 1,
      cov_fc = covered(est$beta_fc[2], sqrt(est$cov_fc[2, 2]), 1),
      cov_fch = covered(est$beta[2], sqrt(est$cov_fch[2, 2]), 1),
      se_fc = sqrt(est$cov_fc[2, 2]), se_fch = sqrt(est$cov_fch[2, 2])
    )
  })
  data.frame(
    design = cell, mse_ols = 100 * mean(r$ols^2), mse_fc = 100 * mean(r$fc^2), mse_fch = 100 * mean(r$fch^2),
    coverage_fc = mean(r$cov_fc), coverage_fch = mean(r$cov_fch), se_fc = mean(r$se_fc), se_fch = mean(r$se_fch)
  )
}))

res <- list(
  pg26_table1 = pg26_table1, pg26_table1_median_k = pg26_table1_k, pg26_general = pg26_general,
  gp24_table1 = gp24_table1, s9 = s9,
  replications = c(pg26_table1 = 2000, pg26_general = 1000, gp24_table1 = 1000, s9 = 1000),
  created = format(Sys.time(), "%Y-%m-%d"), fgls_version = as.character(utils::packageVersion("fgls")),
  minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))
)
dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path("inst", "extdata", "replication.rds"), compress = "xz")
print(res[c("pg26_table1", "gp24_table1", "s9")], digits = 3)
message(sprintf("done in %.1f minutes", res$minutes))
