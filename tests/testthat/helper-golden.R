# Cross-language tolerance policy (CLAUDE_R_PACKAGE.md, Section 3.3).
# Do not loosen any of these without a note in the commit message.
golden_tol <- list(
  exact = 1e-10, # k*, argmin of the BIC path, AR roots, dropped columns
  serial = 1e-8, # rho, beta^D, beta^FGLS given k*; G, Sigma, Q, AVar_*
  gma_theta = 1e-4, # MA(1) MLE (numerical optimiser)
  ridge = 1e-6, # ridge first step given lambda^r
  lasso_objective = 1e-6, # adaptive-Lasso objective at the chosen (lambda, psi)
  lasso_coef_abs = 1e-4, # adaptive-Lasso coefficients (coordinate descent)
  cv_curve_abs = 1e-3, # CV error curves with common folds
  weights = 1e-4, # fitted variances
  given_weights = 1e-8, # WLS beta / eq. (6)-(8) variance / F-CH given pyfgls weights
  fch_beta_e2e = 1e-3, # combined F-CH, R chooses its own weights
  fch_se_e2e = 5e-3
)

golden_dir <- function() test_path("fixtures", "golden")

skip_if_no_golden <- function() {
  skip_if_not_installed("jsonlite")
  skip_if_not(dir.exists(golden_dir()), "golden fixtures not available")
}

read_golden <- function(id) {
  spec <- jsonlite::fromJSON(file.path(golden_dir(), paste0(id, "_golden.json")),
    simplifyVector = TRUE, simplifyMatrix = TRUE
  )
  data <- utils::read.csv(file.path(golden_dir(), paste0(id, "_data.csv")))
  list(spec = spec, data = data, res = spec$results)
}

# Coerce a JSON array to a numeric matrix of the given shape (handles k = 0).
as_mat <- function(x, nrow, ncol) {
  if (nrow * ncol == 0L) {
    return(matrix(0, nrow, ncol))
  }
  m <- as.matrix(x)
  storage.mode(m) <- "double"
  if (!identical(dim(m), c(as.integer(nrow), as.integer(ncol)))) m <- matrix(m, nrow, ncol)
  m
}

# Norm-wise relative difference max|a - b| / max|b|.
rel_diff <- function(a, b) {
  a <- as.numeric(unlist(a))
  b <- as.numeric(unlist(b))
  if (length(a) != length(b)) {
    return(Inf)
  }
  if (!length(b)) {
    return(0)
  }
  max(abs(a - b)) / max(max(abs(b)), .Machine$double.xmin)
}

expect_rel <- function(object, expected, tol, label = NULL) {
  d <- rel_diff(object, expected)
  lab <- if (is.null(label)) deparse(substitute(object)) else label
  expect(d <= tol, sprintf("%s: relative difference %.3g exceeds %.1g", lab, d, tol))
  invisible(d)
}

expect_abs <- function(object, expected, tol, label = NULL) {
  a <- as.numeric(unlist(object))
  b <- as.numeric(unlist(expected))
  d <- if (length(a) == length(b)) max(c(0, abs(a - b))) else Inf
  lab <- if (is.null(label)) deparse(substitute(object)) else label
  expect(d <= tol, sprintf("%s: absolute difference %.3g exceeds %.1g", lab, d, tol))
  invisible(d)
}

# Model matrix of a golden data set (intercept first, as in pyfgls).
golden_X <- function(g) {
  stats::model.matrix(stats::as.formula(g$spec$formula), g$data)
}
