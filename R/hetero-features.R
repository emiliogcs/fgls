# Candidate covariates z_t for the skedastic function, GP24 Section 5 (port of
# pyfgls/hetero/features.py).

.transforms <- c("linear", "poly2", "log", "trig")

# Remove columns that are exact duplicates (up to scale and location) of
# earlier ones. Returns the kept column indices.
dedupe_columns <- function(Z, tol = 1e-10) {
  keep <- integer(0)
  normed <- list()
  for (j in seq_len(ncol(Z))) {
    v <- Z[, j] - mean(Z[, j])
    nv <- sqrt(sum(v^2))
    if (nv == 0) next
    v <- v / nv
    if (any(vapply(normed, function(w) abs(sum(v * w)) > 1 - tol, logical(1)))) next
    keep <- c(keep, j)
    normed[[length(normed) + 1L]] <- v
  }
  keep
}

dictionary_impl <- function(X, names, transforms = .transforms, skip = integer(0)) {
  bad <- setdiff(transforms, .transforms)
  if (length(bad)) {
    fgls_abort(sprintf("Unknown transforms %s.", paste(bad, collapse = ", ")))
  }
  n <- nrow(X)
  cols <- list()
  out_names <- character(0)
  add <- function(v, nm) {
    if (!all(is.finite(v)) || diff(range(v)) <= 1e-12 * max(1, max(abs(v)))) {
      return(invisible(NULL))
    }
    cols[[length(cols) + 1L]] <<- v
    out_names <<- c(out_names, nm)
  }
  for (j in seq_along(names)) {
    if (j %in% skip) next
    x <- X[, j]
    nm <- names[j]
    if (diff(range(x)) == 0) next # constant
    if (length(unique(x)) <= 2L) { # dummy: untransformed
      add(x, nm)
      next
    }
    if ("linear" %in% transforms) add(x, nm)
    if ("poly2" %in% transforms) add(x^2, paste0(nm, "^2"))
    if ("log" %in% transforms && all(x > 0)) {
      lx <- log(x)
      add(lx, paste0("log(", nm, ")"))
      add(lx^2, paste0("log(", nm, ")^2"))
    }
    if ("trig" %in% transforms) {
      add(cos(x), paste0("cos(", nm, ")"))
      add(cos(2 * x), paste0("cos(2", nm, ")"))
    }
  }
  if (!length(cols)) {
    return(list(Z = matrix(0, n, 0L), names = character(0)))
  }
  Z <- do.call(cbind, cols)
  kept <- dedupe_columns(Z)
  Z <- Z[, kept, drop = FALSE]
  colnames(Z) <- out_names[kept]
  if (ncol(Z) > n) {
    fgls_warn(
      sprintf(
        paste0(
          "The candidate dictionary has %d columns for %d observations; the adaptive ",
          "Lasso can handle d > T, but cross-validation will be slow and noisy."
        ),
        ncol(Z), n
      ),
      "sample_size"
    )
  }
  list(Z = Z, names = colnames(Z))
}

#' Candidate covariates for the skedastic function
#'
#' Builds the default dictionary of candidate covariates used by [fgls()] and
#' [lasso_fgls()] with `z = "auto"`: for every non-constant regressor \eqn{x},
#' the columns \eqn{x}, \eqn{x^2}, \eqn{\log x} and \eqn{(\log x)^2} (only when
#' \eqn{x > 0}), \eqn{\cos x} and \eqn{\cos 2x}. This contains the specification
#' \eqn{z = (1, x, [\log x]^2, x^2, \cos x, \cos 2x)} of the simulations of
#' González-Coya and Perron (2024); the intercept of the skedastic regression is
#' always included and never penalised. Binary (dummy) columns are passed through
#' untransformed, undefined or constant transforms are skipped and duplicated
#' columns are removed.
#'
#' @param x A numeric matrix or data frame of regressors (constant columns are
#'   ignored), or a numeric vector.
#' @param transforms Transform families: `"linear"` (\eqn{x}), `"poly2"`
#'   (\eqn{x^2}), `"log"` (\eqn{\log x}, \eqn{(\log x)^2} for positive \eqn{x})
#'   and `"trig"` (\eqn{\cos x}, \eqn{\cos 2x}).
#' @param skip Columns (names or positions) to leave out, e.g. lagged dependent
#'   variables.
#' @return A numeric matrix with one column per candidate covariate and
#'   informative column names.
#' @inherit fgls references
#' @family heteroskedasticity
#' @export
#' @examples
#' X <- cbind(const = 1, x = seq(1, 4, length.out = 10), d = rep(0:1, 5))
#' colnames(make_dictionary(X))
#' colnames(make_dictionary(X, transforms = c("linear", "poly2")))
make_dictionary <- function(x, transforms = c("linear", "poly2", "log", "trig"), skip = NULL) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    nm <- deparse(substitute(x))
    x <- matrix(x, ncol = 1L, dimnames = list(NULL, nm))
  }
  if (!is.numeric(x)) fgls_abort("`x` must be numeric.")
  if (is.null(colnames(x))) colnames(x) <- sprintf("x%d", seq_len(ncol(x)))
  skip_idx <- resolve_columns(skip, colnames(x), "skip")
  dictionary_impl(x, colnames(x), transforms, skip_idx)$Z
}
