# The Durbin (1960) regression underlying the AR filter of PG26 eq. (14)
# (port of pyfgls/serial/durbin.py).
#
# For an AR order k the Durbin regression is
#
#   y_t = sum_{j=1..k} rho_j y_{t-j} + x_t' beta - sum_{j=1..k} x_{t-j}' delta_j + e_kt,
#
# estimated by OLS. Its coefficients on the lagged dependent variable estimate
# the AR coefficients of the errors consistently whether the regressors are
# exogenous or only pre-determined.
#
# The design is laid out lag by lag,
#
#   [x_t | y_{t-1}, x_{t-1} | y_{t-2}, x_{t-2} | ... | y_{t-k}, x_{t-k}],
#
# and lagged columns that are linear combinations of the contemporaneous block
# or of earlier columns (lagged constants, lagged trends, periodic dummies, a lag
# of y that is already a regressor) are dropped *in order*. Hence the design of
# any order k' < k is a prefix of the design of order k on the same sample, which
# lets the lag selection reuse one QR factorization.

# Relative residual-norm threshold for dropping lagged columns: lagged
# deterministic terms are collinear up to rounding (~1e-15).
.lag_block_tol <- 1e-9

# Build the Durbin design of order `order` on observations start + 1, ..., T.
#
# Returns a list with Z (regressors, first q columns = x_t), y, names, n_x,
# ylag_pos (column of each kept y lag, NA if dropped), xlag_pos (order x q
# matrix of columns of the kept x lags), block_end (block_end[k + 1] = number of
# columns of the order-k design), dropped (names), start and order.
durbin_design <- function(y, X, x_names, order, start = order, y_name = "y",
                          tol = .lag_block_tol) {
  n_all <- nrow(X)
  q <- ncol(X)
  if (start < order) stop("start must be at least the AR order")
  rows <- if (start < n_all) (start + 1L):n_all else integer(0)
  base <- X[rows, , drop = FALSE]
  if (order > 0L) {
    parts <- vector("list", 2L * order)
    for (j in seq_len(order)) {
      parts[[2L * j - 1L]] <- y[rows - j]
      parts[[2L * j]] <- X[rows - j, , drop = FALSE]
    }
    A <- do.call(cbind, parts)
    meta_lag <- rep(seq_len(order), each = q + 1L)
    meta_var <- rep(0:q, times = order) # 0 = y, i = regressor i
    ic <- independent_columns(A, base,
      tol = max(tol, default_rank_tol(nrow(A), ncol(A))),
      block_size = q + 1L
    )
    keep <- ic$keep
    drop <- ic$drop
    Qbasis <- ic$Q
  } else {
    A <- matrix(0, length(rows), 0L)
    meta_lag <- meta_var <- integer(0)
    keep <- drop <- integer(0)
    Qbasis <- NULL
  }
  col_name <- function(c) {
    if (meta_var[c] == 0L) {
      paste0("L", meta_lag[c], ".", y_name)
    } else {
      paste0("L", meta_lag[c], ".", x_names[meta_var[c]])
    }
  }
  kept_names <- vapply(keep, col_name, character(1))
  ylag_pos <- rep(NA_integer_, order)
  xlag_pos <- matrix(NA_integer_, order, q)
  counts <- integer(order + 1L)
  for (i in seq_along(keep)) {
    c <- keep[i]
    pos <- q + i
    lag <- meta_lag[c]
    counts[lag + 1L] <- counts[lag + 1L] + 1L
    if (meta_var[c] == 0L) {
      ylag_pos[lag] <- pos
    } else {
      xlag_pos[lag, meta_var[c]] <- pos
    }
  }
  Z <- if (length(keep)) cbind(base, A[, keep, drop = FALSE]) else base
  colnames(Z) <- c(x_names, kept_names)
  list(
    Z = Z, y = y[rows], names = c(x_names, kept_names), n_x = q,
    ylag_pos = ylag_pos, xlag_pos = xlag_pos,
    block_end = as.integer(q + cumsum(counts)),
    dropped = vapply(drop, col_name, character(1)),
    start = as.integer(start), order = as.integer(order), Q = Qbasis
  )
}

# Columns of the order-k design (a prefix of the columns).
design_prefix <- function(design, k) {
  design$Z[, seq_len(design$block_end[k + 1L]), drop = FALSE]
}

# Positions of the kept y lags of order <= k.
design_ylag_columns <- function(design, k = design$order) {
  pos <- design$ylag_pos[seq_len(k)]
  pos[!is.na(pos)]
}

# Diagonal of (Z'Z)^{-1} from a QR factorization (column pivoting undone).
xtx_inv_diag <- function(Z) {
  p <- ncol(Z)
  if (p == 0L) {
    return(numeric(0))
  }
  qz <- qr(Z, tol = 1e-11)
  r <- qz$rank
  R <- qr.R(qz)[seq_len(r), seq_len(r), drop = FALSE]
  Rinv <- backsolve(R, diag(r))
  d <- rep(Inf, p)
  d[qz$pivot[seq_len(r)]] <- rowSums(Rinv * Rinv)
  d
}

# OLS of the Durbin regression with conventional standard errors.
fit_durbin <- function(design) {
  Z <- design$Z
  yv <- design$y
  coef <- lsq(Z, yv)
  resid <- as.numeric(yv - Z %*% coef)
  n <- nrow(Z)
  p <- ncol(Z)
  rss <- sum(resid^2)
  bse <- sqrt(rss / max(n - p, 1L) * xtx_inv_diag(Z))
  names(coef) <- names(bse) <- design$names
  list(
    design = design, coef = coef, resid = resid, sigma2 = rss / n, bse = bse,
    tvalues = coef / bse, beta = coef[seq_len(design$n_x)]
  )
}
