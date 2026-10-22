# Romano and Wolf (2017) parametric skedastic functions, GP24 eqs. (2)-(3)
# (port of pyfgls/hetero/wls.py):
#
# * "wls_s1": v(z) = exp(c + g' [log z_1, ..., log z_d]), requires z > 0;
# * "wls_s2": v(z) = exp(c + g' z);
#
# both fitted by OLS of r_t = log(max(u_t^2, delta^2)) on the transformed
# covariates.

wls_skedastic <- function(r, Z, names = colnames(Z), kind = "wls_s1") {
  n <- nrow(Z)
  d <- ncol(Z)
  if (is.null(names)) names <- sprintf("z%d", seq_len(d))
  if (kind == "wls_s1") {
    bad <- names[apply(Z, 2L, function(v) any(v <= 0))]
    if (length(bad)) {
      fgls_abort(sprintf(
        paste0(
          "hetero = \"wls_s1\" takes logs of the covariates, but %s have non-positive ",
          "values. Use \"wls_s2\" or positive covariates."
        ),
        paste(bad, collapse = ", ")
      ))
    }
    Fm <- log(Z)
  } else {
    Fm <- Z
  }
  keep <- independent_columns(Fm, matrix(1, n, 1L))$keep
  D <- cbind(1, Fm[, keep, drop = FALSE])
  b <- lsq(D, r)
  coef <- numeric(d)
  coef[keep] <- b[-1L]
  names(coef) <- names
  list(
    fitted = drop(D %*% b), intercept = b[1L], coef = coef, names = names,
    df = ncol(D), method = kind, selected = names[coef != 0],
    tuning = list(transform = if (kind == "wls_s1") "log" else "level")
  )
}
