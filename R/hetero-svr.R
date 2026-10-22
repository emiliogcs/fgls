# Support-vector-regression estimate of the skedastic function, GP24 Section 3
# (port of pyfgls/hetero/svr.py; libsvm through 'e1071').
#
# epsilon-insensitive SVR with a Gaussian (RBF) kernel of
# r_t = log(max(u_t^2, delta^2)) on the standardized candidate covariates, with
# (epsilon, C) chosen by K-fold cross-validation. GP24 find it dominated by the
# adaptive Lasso. Degrees of freedom for the bias correction: the number of
# observations whose residual lies on or outside the epsilon-tube.

svr_skedastic <- function(r, Z, names = colnames(Z), cost_grid = c(0.1, 1, 10, 100),
                          eps_grid = c(0.05, 0.25, 0.5, 1), cv_folds = 10L,
                          fold_scheme = "random", seed = NULL, foldid = NULL) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    fgls_abort("hetero = \"svr\" needs the 'e1071' package; install it first.")
  }
  n <- length(r)
  d <- ncol(Z)
  if (n > 5000L) {
    fgls_warn(
      sprintf(
        paste(
          "SVR scales quadratically or worse in T (= %d); this may be slow.",
          "Consider hetero = \"lasso\"."
        ),
        n
      ),
      "sample_size"
    )
  }
  Zs <- standardize_cols(Z)$Zs
  sr <- stats::sd(r) * sqrt((n - 1) / n)
  if (!(sr > 0)) sr <- 1
  if (ncol(Zs) == 0L) {
    return(constant_skedastic(r, d, names, "svr", "no non-constant covariates"))
  }
  folds <- if (is.null(foldid)) {
    make_folds(n, min(cv_folds, n), fold_scheme, seed)
  } else {
    check_foldid(foldid, n)
  }
  gamma <- 1 / (ncol(Zs) * mean((Zs - mean(Zs))^2)) # sklearn's gamma = "scale"
  fit_svm <- function(x, y, C, eps) {
    e1071::svm(x, y,
      type = "eps-regression", kernel = "radial", gamma = gamma,
      cost = C, epsilon = eps, scale = FALSE, tolerance = 1e-3, cachesize = 200
    )
  }
  grid <- expand.grid(eps = eps_grid * sr, C = cost_grid * sr) # C varies slowest
  cv <- vapply(seq_len(nrow(grid)), function(g) {
    mean(vapply(sort(unique(folds)), function(k) {
      te <- folds == k
      m <- fit_svm(Zs[!te, , drop = FALSE], r[!te], grid$C[g], grid$eps[g])
      mean((r[te] - stats::predict(m, Zs[te, , drop = FALSE]))^2)
    }, numeric(1)))
  }, numeric(1))
  best <- which.min(cv)
  model <- fit_svm(Zs, r, grid$C[best], grid$eps[best])
  fitted <- as.numeric(stats::predict(model, Zs))
  eps <- grid$eps[best]
  list(
    fitted = fitted, intercept = -model$rho, coef = stats::setNames(numeric(d), names),
    names = names, df = sum(abs(r - fitted) >= eps - 1e-8), method = "svr",
    selected = names,
    tuning = list(
      C = grid$C[best], epsilon = eps, cv_mse = cv[best], n_support = model$tot.nSV,
      foldid = folds
    )
  )
}
