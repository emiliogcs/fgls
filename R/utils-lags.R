# Lags, autoregressive filters, lag-order rules, critical values and RNG helpers
# (port of pyfgls/_utils.py). Sample positions follow the pyfgls convention:
# `start` is the number of initial observations dropped, so the output rows are
# observations start + 1, ..., T.

as_col_matrix <- function(x) {
  if (is.matrix(x)) x else matrix(as.numeric(x), ncol = 1L)
}

# Matrix of lags [x_{t-1}, ..., x_{t-nlags}] for t = start + 1, ..., T, in
# lag-major order (all variables at lag 1, then all at lag 2, ...).
lag_matrix <- function(x, nlags, start = nlags) {
  a <- as_col_matrix(x)
  n <- nrow(a)
  if (nlags < 0L || start < nlags) {
    stop(sprintf("need 0 <= nlags <= start, got nlags=%d, start=%d", nlags, start))
  }
  if (nlags == 0L) {
    return(matrix(0, n - start, 0L))
  }
  rows <- (start + 1L):n
  do.call(cbind, lapply(seq_len(nlags), function(j) a[rows - j, , drop = FALSE]))
}

# Quasi-difference x_t - sum_j rho_j x_{t-j} for t = start + 1, ..., T: the
# filter A(L) = 1 - rho_1 L - ... - rho_k L^k of PG26 eq. (15), applied to every
# column (a constant becomes 1 - sum(rho)). The first k observations are dropped
# (no Prais-Winsten transformation).
ar_filter <- function(x, rho, start = length(rho)) {
  is_mat <- is.matrix(x)
  a <- as_col_matrix(x)
  n <- nrow(a)
  k <- length(rho)
  rows <- if (start < n) (start + 1L):n else integer(0)
  out <- a[rows, , drop = FALSE]
  for (j in seq_len(k)) {
    out <- out - rho[j] * a[rows - j, , drop = FALSE]
  }
  if (is_mat) out else drop(out)
}

# Roots of the lag polynomial 1 - rho_1 z - ... - rho_k z^k. The filter is
# invertible (PG26 Assumption 2) when all roots lie outside the unit circle.
ar_poly_roots <- function(rho) {
  p <- c(1, -as.numeric(rho))
  while (length(p) > 1L && p[length(p)] == 0) p <- p[-length(p)]
  if (length(p) <= 1L) {
    return(complex(0))
  }
  polyroot(p)
}

# Default maximal lag order int(12 * (T/200)^(1/4)) (12 at T = 200).
default_kmax <- function(n) {
  as.integer(floor(12 * (n / 200)^0.25))
}

# Rate-based maximal order floor((T / ln T)^(1/3)) (PG26 Suppl. Lemma A.2).
theory_kmax <- function(n) {
  as.integer(floor((n / log(n))^(1 / 3)))
}

# Lower bound ceil((ln T)^(3/2)) on the AR order (PG26 Suppl. Lemma A.2).
adjusted_min_lags <- function(n) {
  as.integer(ceiling(log(n)^1.5))
}

# Two-sided critical value: standard normal, or Student-t when df is finite.
critical_value <- function(alpha, df = Inf) {
  if (!(alpha > 0 && alpha < 1)) {
    fgls_abort(sprintf("alpha must be in (0, 1), got %s.", format(alpha)))
  }
  if (is.finite(df)) stats::qt(1 - alpha / 2, df) else stats::qnorm(1 - alpha / 2)
}

# Two-sided p-values of t/z statistics.
two_sided_pvalues <- function(stat, df = Inf) {
  a <- abs(stat)
  if (is.finite(df)) 2 * stats::pt(a, df, lower.tail = FALSE) else 2 * stats::pnorm(a, lower.tail = FALSE)
}

# Evaluate `code` with the RNG seeded by `seed`, restoring the caller's RNG
# state (and RNG kind) afterwards. With seed = NULL the current stream is used.
with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(code)
  }
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit(
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    },
    add = TRUE
  )
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  code
}

check_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    fgls_abort("seed must be NULL or a single finite number.")
  }
  invisible(NULL)
}
