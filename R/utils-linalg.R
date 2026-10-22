# Numerically careful linear-algebra helpers (port of pyfgls/_linalg.py).
#
# Rules followed throughout the package:
# * regressions are solved by QR, never through the normal equations;
# * rank decisions use column-normalised matrices so that they do not depend on
#   the units of the regressors;
# * only small (q x q) symmetric positive definite systems are inverted, by
#   Cholesky with a pseudo-inverse fallback.

.eps <- .Machine$double.eps

# Relative tolerance used to declare a column linearly dependent:
# max(n, m) * eps * 100 (about 1e-11 for a few hundred observations).
default_rank_tol <- function(n, m) {
  max(n, m, 1) * .eps * 100
}

# Least-squares coefficients of y (vector or matrix) on X via a QR
# decomposition. Exactly collinear columns get a zero coefficient.
lsq <- function(X, y) {
  p <- ncol(X)
  if (p == 0L) {
    return(if (is.matrix(y)) matrix(0, 0L, ncol(y)) else numeric(0))
  }
  b <- qr.coef(qr(X, tol = 1e-11), y)
  b[is.na(b)] <- 0
  b
}

# OLS coefficients and residuals.
ols_fit <- function(X, y) {
  b <- lsq(X, y)
  list(coef = b, resid = y - X %*% b)
}

# Euclidean norm of each column, zeros replaced by one (safe for scaling).
column_norms <- function(X) {
  nr <- sqrt(colSums(X * X))
  nr[nr == 0] <- 1
  nr
}

# Split the columns of A into linearly independent and dependent ones.
#
# Columns are processed in order by modified Gram-Schmidt with
# re-orthogonalisation on unit-norm columns, so that when a column is a linear
# combination of `base` and of earlier columns of A it is the *later* column that
# is reported as dependent ("the lag of the constant is dropped"). The selection
# for the first j columns does not depend on later columns.
#
# Returns list(keep, drop): 1-based column indices of A.
independent_columns <- function(A, base = NULL, tol = NULL) {
  A <- as.matrix(A)
  n <- nrow(A)
  m <- ncol(A)
  if (is.null(tol)) {
    tol <- default_rank_tol(n, m + if (is.null(base)) 0L else ncol(base))
  }
  if (m == 0L) {
    return(list(keep = integer(0), drop = integer(0)))
  }
  if (!is.null(base) && ncol(base) > 0L) {
    Qb <- qr.Q(qr(sweep(base, 2L, column_norms(base), "/")))
  } else {
    Qb <- matrix(0, n, 0L)
  }
  Q <- matrix(0, n, min(n, m))
  r <- 0L
  keep <- integer(0)
  drop <- integer(0)
  norms <- column_norms(A)
  for (j in seq_len(m)) {
    if (!any(A[, j] != 0)) {
      drop <- c(drop, j)
      next
    }
    v <- A[, j] / norms[j]
    for (pass in 1:2) { # twice is enough (Kahan-Parlett)
      if (ncol(Qb) > 0L) v <- v - Qb %*% crossprod(Qb, v)
      if (r > 0L) {
        Qr <- Q[, seq_len(r), drop = FALSE]
        v <- v - Qr %*% crossprod(Qr, v)
      }
    }
    nv <- sqrt(sum(v * v))
    if (nv <= tol || r >= n) {
      drop <- c(drop, j)
    } else {
      r <- r + 1L
      Q[, r] <- v / nv
      keep <- c(keep, j)
    }
  }
  list(keep = keep, drop = drop)
}

# Numerical rank from a pivoted QR of the column-normalised matrix.
matrix_rank <- function(X, tol = NULL) {
  m <- ncol(X)
  if (m == 0L) {
    return(0L)
  }
  if (is.null(tol)) tol <- default_rank_tol(nrow(X), m)
  d <- abs(diag(qr.R(qr(sweep(X, 2L, column_norms(X), "/"), LAPACK = TRUE))))
  if (length(d) == 0L || d[1L] == 0) {
    return(0L)
  }
  sum(d > tol * d[1L])
}

# Diagonal of the projection matrix X (X'X)^- X' (the leverages), from the
# orthonormal factor of a pivoted QR restricted to the numerical rank, so that
# it is valid for rank-deficient X too.
hat_diag <- function(X) {
  n <- nrow(X)
  m <- ncol(X)
  if (m == 0L) {
    return(numeric(n))
  }
  qx <- qr(X, LAPACK = TRUE)
  d <- abs(diag(qr.R(qx)))
  rank <- if (length(d) && d[1L] > 0) sum(d > default_rank_tol(n, m) * d[1L]) else 0L
  if (rank == 0L) {
    return(numeric(n))
  }
  Q <- qr.Q(qx)[, seq_len(rank), drop = FALSE]
  rowSums(Q * Q)
}

symmetrize <- function(A) {
  0.5 * (A + t(A))
}

# Moore-Penrose pseudo-inverse of a symmetric matrix.
pinv_sym <- function(A) {
  e <- eigen(symmetrize(A), symmetric = TRUE)
  v <- e$values
  cut <- 1e-15 * max(abs(v), 0)
  inv <- ifelse(abs(v) > cut, 1 / v, 0)
  e$vectors %*% (inv * t(e$vectors))
}

# Solve A X = B for a symmetric positive (semi)definite A by Cholesky; if A is
# not numerically positive definite, warn and use the pseudo-inverse.
solve_psd <- function(A, B, name = "matrix") {
  A <- symmetrize(as.matrix(A))
  R <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(R)) {
    return(backsolve(R, forwardsolve(t(R), B)))
  }
  fgls_warn(
    paste0(
      "The ", name, " is not positive definite (it is singular or nearly so); ",
      "a pseudo-inverse is used instead. Check the regressors for collinearity."
    ),
    "collinearity"
  )
  pinv_sym(A) %*% B
}

# Inverse of a small symmetric positive definite matrix (see solve_psd()).
inv_psd <- function(A, name = "matrix") {
  A <- as.matrix(A)
  if (length(A) == 0L) {
    return(A)
  }
  symmetrize(solve_psd(A, diag(nrow(A)), name = name))
}

# Symmetrise A and clip negative eigenvalues to zero. Negative eigenvalues
# smaller than rtol * max|eigenvalue| are rounding noise and are clipped
# silently; larger ones trigger a warning.
make_psd <- function(A, name = "covariance matrix", rtol = 1e-8) {
  A <- symmetrize(as.matrix(A))
  if (length(A) == 0L) {
    return(A)
  }
  e <- eigen(A, symmetric = TRUE)
  w <- e$values
  if (min(w) >= 0) {
    return(A)
  }
  scale <- max(max(abs(w)), .Machine$double.xmin)
  if (min(w) < -rtol * scale) {
    fgls_warn(
      sprintf(
        paste0(
          "The %s was not positive semidefinite (smallest eigenvalue %.3g); ",
          "negative eigenvalues were set to zero."
        ),
        name, min(w)
      ),
      "collinearity"
    )
  }
  w <- pmax(w, 0)
  symmetrize(e$vectors %*% (w * t(e$vectors)))
}

# 2-norm condition number (Inf if singular).
condition_number <- function(A) {
  s <- svd(A, nu = 0L, nv = 0L)$d
  if (length(s) == 0L || s[length(s)] == 0) {
    return(Inf)
  }
  s[1L] / s[length(s)]
}
