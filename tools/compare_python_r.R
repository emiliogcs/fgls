#!/usr/bin/env Rscript
# Cross-language validation of fgls against pyfgls (see tools/README.md).
#
# Both implementations fit the *same* simulated data sets with the *same*
# cross-validation folds; R generates both (R's random number generator), so a
# data set can be regenerated in R from its seed.
#
#   Rscript tools/compare_python_r.R simulate [n_datasets] [cores]
#       -> validation/cache/<design>.json      (data, folds, seeds)
#   python  tools/make_validation.py validation/cache [workers]
#       -> validation/cache/<design>_py.json   (pyfgls fits)
#   Rscript tools/compare_python_r.R compare [cores]
#       -> validation/cross_language_report.html
#       -> inst/extdata/cross_language.rds      (summary + pyfgls outputs, used by the
#                                                validation vignette and a slow test)
#
# Run from the package root with the package installed (or via devtools::load_all()).

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(fgls)
  }
})

args <- commandArgs(trailingOnly = TRUE)
cmd <- if (length(args)) args[1] else "compare"
cache <- file.path("validation", "cache")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

paper_z <- function(x) {
  Z <- cbind(x, log(x)^2, x^2, cos(x), cos(2 * x))
  colnames(Z) <- c("x", "log(x)^2", "x^2", "cos(x)", "cos(2x)")
  Z
}
paper_w <- function(w) {
  W <- cbind(w, w^2, cos(w), cos(2 * w))
  colnames(W) <- c("w", "w^2", "cos(w)", "cos(2w)")
  W
}

# Design: simulation, model type and seed base. All use T = 200.
designs <- list(
  serial_ar1_nonexog = list(
    kind = "serial", base = 1e6,
    sim = function(s) simulate_serial(200, ar = 0.8, rho_x = 0.8, gamma = 0.5, burn = 100, seed = s)
  ),
  serial_ar2_exog = list(
    kind = "serial", base = 2e6,
    sim = function(s) simulate_serial(200, ar = c(1.34, -0.42), rho_x = 0.8, gamma = 0, burn = 100, seed = s)
  ),
  serial_ma1 = list(
    kind = "ma1", base = 3e6,
    sim = function(s) simulate_serial(200, ar = numeric(0), ma = -0.4, rho_x = 0.8, gamma = 0, burn = 100, seed = s)
  ),
  hetero_log2 = list(
    kind = "hetero", base = 4e6,
    sim = function(s) simulate_hetero(200, skedastic = "log2", seed = s)
  ),
  hetero_exp = list(
    kind = "hetero", base = 5e6,
    sim = function(s) simulate_hetero(200, skedastic = "exp", seed = s)
  ),
  combined_ar1_log2 = list(
    kind = "combined", base = 6e6,
    sim = function(s) simulate_combined(200, ar = 0.5, skedastic = "log2", seed = s)
  ),
  combined_ar1_exp = list(
    kind = "combined", base = 7e6,
    sim = function(s) simulate_combined(200, ar = 0.5, skedastic = "exp", seed = s)
  ),
  combined_ar2_x2 = list(
    kind = "combined", base = 8e6,
    sim = function(s) simulate_combined(200, ar = c(1.34, -0.42), skedastic = "x2", seed = s)
  )
)

# Fold labels for the skedastic fit: n rows (n_eff = T - k* for the combined
# estimator, k* from the deterministic lag selection).
design_folds <- function(d, sim, seed) {
  if (d$kind == "hetero") {
    n <- nrow(sim)
  } else if (d$kind == "combined") {
    X <- cbind(1, sim$x)
    n <- nrow(sim) - fgls:::estimate_serial(sim$y, X, c("(Intercept)", "x"))$order
  } else {
    return(NULL)
  }
  fgls:::make_folds(n, 10L, "random", seed = seed + 500000L)
}

mapper <- function(cores) {
  if (cores > 1L && .Platform$OS.type != "windows") {
    function(X, f) parallel::mclapply(X, f, mc.cores = cores)
  } else {
    lapply
  }
}

simulate_all <- function(n_datasets, cores) {
  map <- mapper(cores)
  for (nm in names(designs)) {
    d <- designs[[nm]]
    seeds <- as.integer(d$base + seq_len(n_datasets))
    sets <- map(seeds, function(s) {
      sim <- d$sim(s)
      out <- list(seed = s, y = sim$y, x = sim$x)
      if (!is.null(sim$w)) out$w <- sim$w
      f <- design_folds(d, sim, s)
      if (!is.null(f)) out$foldid <- f
      out
    })
    jsonlite::write_json(list(design = nm, kind = d$kind, datasets = sets),
      file.path(cache, paste0(nm, ".json")),
      digits = I(17), auto_unbox = TRUE
    )
    message(sprintf("%s: %d data sets", nm, n_datasets))
  }
}

# R fit of one data set (public API), returning the compared quantities.
fit_r <- function(kind, ds) {
  df <- data.frame(y = ds$y, x = ds$x)
  ctl <- if (!is.null(ds$foldid)) fgls_control(foldid = as.integer(ds$foldid)) else fgls_control()
  q <- function(fit) {
    se <- sqrt(diag(vcov(fit)))
    list(beta = unname(coef(fit)[2]), se = unname(se[2]))
  }
  suppressWarnings(switch(kind,
    serial = {
      f <- serial_fgls(y ~ x, data = df)
      c(q(f), list(
        kstar = f$ar_order, rho1 = unname(f$rho[1]),
        se_naive = sqrt(vcov(f, type = "naive")[2, 2]), share = f$nonexogeneity_share
      ))
    },
    ma1 = {
      f <- serial_fgls(y ~ x, data = df)
      g <- serial_fgls(y ~ x, data = df, serial = "ma1")
      c(q(f), list(
        kstar = f$ar_order, rho1 = unname(f$rho[1]), theta = g$serial$extra$theta,
        beta_gma = unname(coef(g)[2]), se_gma = sqrt(vcov(g)[2, 2])
      ))
    },
    hetero = {
      f <- lasso_fgls(y ~ x, data = df, z = paper_z(df$x), control = ctl)
      c(q(f), list(
        psi = f$skedastic$psi, lambda = f$skedastic$lambda, df_hat = f$skedastic$df_hat,
        selected = paste(f$skedastic$selected, collapse = "|")
      ))
    },
    combined = {
      f <- fgls(y ~ x, data = df, z = paper_w(ds$w), control = ctl)
      fc <- f$stages$serial
      c(q(f), list(
        kstar = f$ar_order, beta_fc = unname(coef(fc)[2]), se_fc = sqrt(vcov(fc)[2, 2]),
        psi = f$skedastic$psi, df_hat = f$skedastic$df_hat,
        selected = paste(f$skedastic$selected, collapse = "|")
      ))
    }
  ))
}

as_frame <- function(lst) {
  keys <- names(lst[[1]])
  as.data.frame(lapply(stats::setNames(keys, keys), function(k) sapply(lst, `[[`, k)),
    stringsAsFactors = FALSE
  )
}

compare_all <- function(cores) {
  map <- mapper(cores)
  rows <- list()
  per_design <- list()
  for (nm in names(designs)) {
    f_py <- file.path(cache, paste0(nm, "_py.json"))
    if (!file.exists(f_py)) {
      message("skipping ", nm, ": run tools/make_validation.py first")
      next
    }
    dat <- jsonlite::read_json(file.path(cache, paste0(nm, ".json")), simplifyVector = TRUE,
      simplifyDataFrame = FALSE)
    py <- jsonlite::read_json(f_py, simplifyVector = TRUE)
    kind <- dat$kind
    r_out <- as_frame(map(dat$datasets, function(ds) fit_r(kind, ds)))
    p_out <- as.data.frame(py$results, stringsAsFactors = FALSE)
    p_out <- p_out[, names(r_out), drop = FALSE]
    seeds <- vapply(dat$datasets, function(ds) as.integer(ds$seed), integer(1))
    per_design[[nm]] <- list(kind = kind, seeds = seeds, r = r_out, py = p_out, pyfgls = py$pyfgls)
    for (qty in names(r_out)) {
      a <- r_out[[qty]]
      b <- p_out[[qty]]
      if (is.character(a)) {
        rows[[length(rows) + 1L]] <- data.frame(
          design = nm, quantity = qty, max_abs = NA, median_abs = NA, max_rel = NA,
          share_different = mean(a != b), stringsAsFactors = FALSE
        )
      } else {
        dd <- abs(a - b)
        rel <- dd / pmax(abs(b), 1e-12)
        rows[[length(rows) + 1L]] <- data.frame(
          design = nm, quantity = qty, max_abs = max(dd), median_abs = stats::median(dd),
          max_rel = max(rel), share_different = mean(rel > 1e-6), stringsAsFactors = FALSE
        )
      }
    }
    cov_r <- mean(abs(r_out$beta - 1) <= qnorm(0.975) * r_out$se)
    cov_p <- mean(abs(p_out$beta - 1) <= qnorm(0.975) * p_out$se)
    rows[[length(rows) + 1L]] <- data.frame(
      design = nm, quantity = "coverage_95", max_abs = abs(cov_r - cov_p), median_abs = NA,
      max_rel = NA, share_different = NA, stringsAsFactors = FALSE
    )
    per_design[[nm]]$coverage <- c(R = cov_r, pyfgls = cov_p)
    message(sprintf("%s: compared %d data sets (coverage R %.3f, pyfgls %.3f)", nm, length(seeds), cov_r, cov_p))
  }
  if (!length(rows)) stop("nothing to compare: run tools/make_validation.py first")
  summary <- do.call(rbind, rows)
  res <- list(
    summary = summary, designs = per_design, created = format(Sys.time(), "%Y-%m-%d"),
    fgls_version = as.character(utils::packageVersion("fgls"))
  )
  dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, file.path("inst", "extdata", "cross_language.rds"), compress = "xz")
  if (requireNamespace("rmarkdown", quietly = TRUE)) {
    rds_path <- normalizePath(file.path("inst", "extdata", "cross_language.rds"))
    out_dir <- normalizePath("validation")
    rmarkdown::render(file.path("tools", "cross_language_report.Rmd"),
      output_file = "cross_language_report.html", output_dir = out_dir,
      intermediates_dir = tempdir(), params = list(rds = rds_path), quiet = TRUE
    )
  }
  print(summary, digits = 3, row.names = FALSE)
  invisible(res)
}

if (cmd == "simulate") {
  n <- if (length(args) >= 2) as.integer(args[2]) else 500L
  cores <- if (length(args) >= 3) as.integer(args[3]) else 1L
  simulate_all(n, cores)
} else if (cmd == "compare") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
  compare_all(cores)
} else {
  stop("usage: compare_python_r.R simulate [n] [cores] | compare [cores]")
}
