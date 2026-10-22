# Adaptive Lasso estimate of the skedastic function, GP24 eq. (4) (port of
# pyfgls/hetero/adaptive_lasso.py, with glmnet as the Lasso solver).
#
# Given the response r_t = log(max(u_t^2, delta^2)) and candidate covariates
# z_t (d may exceed T), solve
#
#   min_{g0, g}  1/2 sum_t (r_t - g0 - z_t' g)^2 + lambda sum_j w_j |g_j|,
#   w_j = |g_j^ridge|^(-psi),
#
# with an unpenalized intercept (Zou, 2006):
#
# 1. standardize the columns of z (mean 0, variance 1, as glmnet does);
# 2. ridge first step, its penalty chosen by K-fold cross-validation (exact
#    closed form through the SVD, so d > T is fine);
# 3. for each psi, rescale column j by 1/w_j = |g_j^ridge|^psi so that an
#    ordinary Lasso solves the adaptive problem, and back-transform (psi = 0 is
#    the plain Lasso; a zero ridge coefficient means an infinite weight and the
#    column is excluded). Column rescaling is algebraically identical to
#    glmnet's `penalty.factor` and keeps the penalty on one scale for all psi.
#    The Lasso is solved exactly by the LARS homotopy (engine "lars", default,
#    as pyfgls) or by glmnet followed by an exact active-set refinement
#    (engine "glmnet");
# 4. choose (lambda, psi) jointly by K-fold cross-validation of the mean
#    squared prediction error of r, with the *same folds for every psi*. By
#    default (cv_weights = "fold") the ridge step, and hence the weights, are
#    re-estimated on each training fold.
#
# glmnet minimizes (2n)^{-1} RSS + lambda_glmnet * penalty, i.e. the objective
# above divided by n: lambda_paper = n * lambda_glmnet (both reported).
#
# If every slope is zero at the chosen penalty the fitted variance is constant
# and FGLS reduces to OLS (the "no loss under homoskedasticity" property).

.psi_grid <- c(0, 0.25, 0.5, 0.75, 1, 2)

# Fold labels 1..K (random or contiguous blocks).
make_folds <- function(n, K, scheme = "random", seed = NULL) {
  if (K < 2L) fgls_abort(sprintf("cv_folds must be at least 2, got %d.", K))
  if (n < K) fgls_abort(sprintf("cv_folds = %d exceeds the number of observations (%d).", K, n))
  if (scheme == "random") {
    return(with_seed(seed, sample(rep_len(seq_len(K), n))))
  }
  as.integer(((seq_len(n) - 1L) * K) %/% n + 1L)
}

check_foldid <- function(foldid, n) {
  if (length(foldid) != n) {
    fgls_abort(sprintf(
      "`foldid` has %d entries; the skedastic fit uses %d observations.",
      length(foldid), n
    ))
  }
  f <- as.integer(factor(foldid))
  if (max(f) < 2L) fgls_abort("`foldid` must define at least two folds.")
  f
}

# Standardize columns (population sd, as glmnet); constant columns dropped.
standardize_cols <- function(Z) {
  mu <- colMeans(Z)
  Zc <- sweep(Z, 2L, mu, "-")
  sdv <- sqrt(colMeans(Zc^2))
  keep <- which(sdv > 1e-12 * pmax(1, abs(mu)))
  Zs <- sweep(Zc[, keep, drop = FALSE], 2L, sdv[keep], "/")
  list(Zs = Zs, mu = mu, sd = sdv, keep = keep)
}

# Spectral decomposition of a centered design: eigen-decomposition of the
# d x d Gram matrix when n > d (s^2 = eigenvalues), otherwise the SVD. Both give
# ridge solutions V diag(1 / (s^2 + lam)) V' Z'r; the penalties used here are at
# least n * 1e-4, so the Gram route loses no accuracy.
centred_spectrum <- function(Zc) {
  if (nrow(Zc) > ncol(Zc)) {
    e <- eigen(crossprod(Zc), symmetric = TRUE)
    list(s2 = pmax(e$values, 0), v = e$vectors)
  } else {
    s <- svd(Zc, nu = 0L)
    list(s2 = s$d^2, v = s$v)
  }
}

# Ridge slopes argmin ||r - g0 - Zs g||^2 + lam ||g||^2 (intercept unpenalized).
ridge_fit <- function(Zs, r, lam) {
  Zc <- sweep(Zs, 2L, colMeans(Zs), "-")
  sp <- centred_spectrum(Zc)
  vtc <- drop(crossprod(sp$v, crossprod(Zc, r - mean(r))))
  drop(sp$v %*% (vtc / (sp$s2 + lam)))
}

# Ridge with the penalty chosen by K-fold CV over lam = n 10^[-4, 4].
ridge_cv <- function(Zs, r, folds, n_grid = 50L) {
  n <- nrow(Zs)
  grid <- n * 10^seq(-4, 4, length.out = n_grid)
  sse <- numeric(n_grid)
  for (k in sort(unique(folds))) {
    te <- folds == k
    tr <- !te
    Ztr <- Zs[tr, , drop = FALSE]
    zm <- colMeans(Ztr)
    rm <- mean(r[tr])
    Zc <- sweep(Ztr, 2L, zm, "-")
    sp <- centred_spectrum(Zc)
    vtc <- drop(crossprod(sp$v, crossprod(Zc, r[tr] - rm)))
    shrink <- 1 / outer(grid, sp$s2, "+") # (grid, rank)
    coefs <- (shrink * outer(rep(1, n_grid), vtc)) %*% t(sp$v) # (grid, d)
    pred <- sweep(Zs[te, , drop = FALSE], 2L, zm, "-") %*% t(coefs) + rm
    sse <- sse + colSums((r[te] - pred)^2)
  }
  lam <- grid[which.min(sse)]
  list(coef = ridge_fit(Zs, r, lam), lambda = lam, grid = grid, cv_sse = sse)
}

# Column scales 1/w_j = |g_j|^psi (zero where the ridge coefficient is 0).
adaptive_scale <- function(g, psi) {
  a <- abs(g)
  ifelse(a > 0, a^psi, 0)
}

# glmnet path at the given (decreasing) penalties. A path truncated by a
# convergence failure is retried with a looser threshold; remaining penalties
# then reuse the last solution, with a warning.
glmnet_path <- function(x, r, lambdas, thresh) {
  run <- function(th, maxit) {
    withCallingHandlers(
      glmnet::glmnet(x, r,
        family = "gaussian", alpha = 1, lambda = lambdas, standardize = FALSE,
        intercept = TRUE, thresh = th, maxit = maxit
      ),
      warning = function(w) invokeRestart("muffleWarning")
    )
  }
  fit <- run(thresh, 1e5)
  B <- as.matrix(fit$beta)
  if (ncol(B) < length(lambdas)) {
    fit <- run(max(thresh * 1e3, 1e-7), 1e6)
    B <- as.matrix(fit$beta)
  }
  if (ncol(B) < length(lambdas)) {
    fgls_warn(
      "The Lasso path did not converge for the smallest penalties; the last solution is reused.",
      "convergence"
    )
    B <- cbind(B, B[, rep(ncol(B), length(lambdas) - ncol(B)), drop = FALSE])
  }
  B
}

# Exact refinement of Lasso solutions by feature-sign search (Lee, Battle,
# Raina and Ng, 2007), an active-set method that strictly decreases
#
#   f(b) = (2n)^{-1} ||rc - Zc b||^2 + lambda ||b||_1
#
# (written with G = Zc'Zc and c = Zc'rc) and stops at a point satisfying the
# Lasso KKT conditions, i.e. at the exact solution. It is warm-started from
# glmnet's (coordinate-descent) solution, so a few steps suffice; glmnet keeps
# finding the path and the support, and the refinement removes its
# convergence tolerance (pyfgls computes the exact path by LARS).
lasso_kkt_gap <- function(G, cvec, n, lambda, b) {
  g <- drop(G %*% b - cvec) / n
  A <- b != 0
  max(c(abs(g[A] + lambda * sign(b[A])), abs(g[!A]) - lambda, 0))
}

feature_sign <- function(G, cvec, n, lambda, b, tol = 1e-10, max_iter = NULL) {
  d <- length(b)
  if (is.null(max_iter)) max_iter <- 50L + 10L * d
  f <- function(v) (sum(v * (G %*% v)) - 2 * sum(cvec * v)) / (2 * n) + lambda * sum(abs(v))
  theta <- sign(b)
  scale_tol <- tol * max(lambda, max(abs(cvec)) / n)
  for (iter in seq_len(max_iter)) {
    g <- drop(G %*% b - cvec) / n
    A <- which(b != 0)
    # optimality of the non-zero coefficients
    opt_nz <- !length(A) || max(abs(g[A] + lambda * sign(b[A]))) <= scale_tol
    if (opt_nz) {
      inactive <- which(b == 0)
      viol <- inactive[abs(g[inactive]) > lambda + scale_tol]
      if (!length(viol)) {
        return(b)
      }
      i <- viol[which.max(abs(g[viol]))]
      theta[i] <- -sign(g[i])
      A <- sort(c(A, i))
    }
    # feature-sign step: unconstrained minimiser on A with signs theta
    b_active_new <- tryCatch(
      solve(G[A, A, drop = FALSE], cvec[A] - n * lambda * theta[A]),
      error = function(e) NULL
    )
    if (is.null(b_active_new)) {
      return(NULL)
    }
    b_new <- numeric(d)
    b_new[A] <- b_active_new
    # line search over the segment, at its end and at sign changes
    cand <- list(b_new)
    dlt <- b_new - b
    cross <- A[b[A] != 0 & sign(b_new[A]) != sign(b[A])]
    for (j in cross) {
      t <- b[j] / (b[j] - b_new[j])
      bt <- b + t * dlt
      bt[j] <- 0
      cand[[length(cand) + 1L]] <- bt
    }
    vals <- vapply(cand, f, numeric(1))
    best <- cand[[which.min(vals)]]
    if (min(vals) > f(b) + 1e-15 * abs(f(b)) && !opt_nz) {
      return(NULL) # no progress: give up and keep the warm start
    }
    best[abs(best) < 1e-300] <- 0
    b <- best
    theta <- sign(b)
  }
  NULL
}

# Coefficients whose contribution to the fit is below round-off are zero (e.g.
# at lambda = lambda_max, where the solution is the null model).
drop_negligible <- function(B, scale_j, cut) {
  B[abs(B) * scale_j <= cut] <- 0
  B
}

# KKT gap of every column of B (d x L) at the penalties `lambdas`.
kkt_gaps <- function(G, cvec, n, lambdas, B) {
  Gr <- (G %*% B - cvec) / n
  act <- B != 0
  lam <- matrix(lambdas, nrow(B), length(lambdas), byrow = TRUE)
  gap <- ifelse(act, abs(Gr + lam * sign(B)), pmax(abs(Gr) - lam, 0))
  apply(gap, 2L, max)
}

# Exact solution with support and signs given by `sgn` for several
# penalties: solves G_AA b_A = c_A - n lambda s_A. Columns whose signs do not
# agree are NA.
solve_on_support <- function(G, cvec, n, lambdas, sgn) {
  d <- length(sgn)
  out <- matrix(0, d, length(lambdas))
  A <- which(sgn != 0)
  if (!length(A)) {
    return(out)
  }
  rhs <- matrix(cvec[A], length(A), length(lambdas)) - n * outer(sgn[A], lambdas)
  sol <- tryCatch(solve(G[A, A, drop = FALSE], rhs), error = function(e) NULL)
  if (is.null(sol)) {
    out[] <- NA_real_
    return(out)
  }
  bad <- colSums(sign(sol) != sgn[A]) > 0
  out[A, ] <- sol
  out[, bad] <- NA_real_
  out
}

# Refine a whole path. Runs of penalties sharing a sign pattern in glmnet's
# solutions are solved exactly at once; columns that then satisfy the KKT
# conditions are exact Lasso solutions. For the others (support misidentified
# by glmnet) the support of the previous exact solution is tried, and finally
# feature-sign search.
polish_lasso <- function(G, cvec, n, lambdas, B, rss0 = NULL) {
  L <- length(lambdas)
  scale_j <- sqrt(pmax(diag(G), 0) / n)
  cut <- if (is.null(rss0)) 0 else 1e-12 * sqrt(rss0 / n)
  B <- drop_negligible(B, scale_j, cut)
  tol <- 1e-8 * pmax(lambdas, max(abs(cvec)) / n)
  S <- sign(B)
  new_run <- c(TRUE, if (L > 1L) colSums(S[, -1L, drop = FALSE] != S[, -L, drop = FALSE]) > 0)
  run_id <- cumsum(new_run)
  Bx <- matrix(NA_real_, nrow(B), L)
  for (k in unique(run_id)) {
    cols <- which(run_id == k)
    Bx[, cols] <- solve_on_support(G, cvec, n, lambdas[cols], S[, cols[1L]])
  }
  ok <- colSums(is.na(Bx)) == 0
  ok[ok] <- kkt_gaps(G, cvec, n, lambdas[ok], Bx[, ok, drop = FALSE]) <= tol[ok]
  B[, ok] <- Bx[, ok]
  prev <- NULL
  for (l in seq_len(L)) {
    if (ok[l]) {
      prev <- B[, l]
      next
    }
    b <- NULL
    if (!is.null(prev)) {
      cand <- solve_on_support(G, cvec, n, lambdas[l], sign(prev))
      if (!anyNA(cand) && kkt_gaps(G, cvec, n, lambdas[l], cand) <= tol[l]) b <- drop(cand)
    }
    if (is.null(b)) {
      for (b0 in c(list(B[, l]), if (!is.null(prev)) list(prev))) {
        cand <- feature_sign(G, cvec, n, lambdas[l], b0)
        if (!is.null(cand) && kkt_gaps(G, cvec, n, lambdas[l], matrix(cand)) <= tol[l]) {
          b <- cand
          break
        }
      }
    }
    if (!is.null(b)) {
      B[, l] <- drop_negligible(b, scale_j, cut)
      prev <- B[, l]
    }
  }
  B
}

# Lasso coefficients (d x L) of r on Ztil at the penalties `lambdas`
# (decreasing), objective (2n)^{-1} RSS + lambda ||g||_1 with an unpenalized
# intercept, solved by glmnet and refined on the active set. Returns the
# column means and the response mean, so that predictions are
# (Z - zm) %*% coef + rm.
lasso_path <- function(Ztil, r, lambdas, thresh = 1e-10, polish = TRUE) {
  d <- ncol(Ztil)
  n <- nrow(Ztil)
  zm <- colMeans(Ztil)
  rm <- mean(r)
  active <- which(apply(Ztil, 2L, function(v) diff(range(v)) > 0))
  coefs <- matrix(0, d, length(lambdas))
  if (length(active) && stats::sd(r) > 0) {
    x <- Ztil[, active, drop = FALSE]
    xx <- if (ncol(x) == 1L) cbind(x, 0) else x # glmnet needs two columns; a zero one is excluded
    B <- glmnet_path(xx, r, lambdas, thresh)[seq_along(active), , drop = FALSE]
    if (polish) {
      Zc <- sweep(x, 2L, zm[active], "-")
      B <- polish_lasso(crossprod(Zc), drop(crossprod(Zc, r - rm)), n, lambdas, B,
        rss0 = sum((r - rm)^2)
      )
    }
    coefs[active, ] <- B
  }
  list(coef = coefs, zm = zm, rm = rm)
}

# Fold-by-fold mean squared prediction errors: one (K x L) matrix per psi.
# With the exact LARS engine the Gram matrix of each training fold is computed
# once and rescaled for every psi (column rescaling commutes with centering).
cv_errors <- function(Zs, r, folds, psis, grids, g_full, lam_ridge, fold_weights,
                      units = NULL, engine = "lars", thresh = 1e-10) {
  labels <- sort(unique(folds))
  per_fold <- lapply(labels, function(k) {
    te <- folds == k
    tr <- !te
    g <- if (fold_weights) ridge_fit(Zs[tr, , drop = FALSE], r[tr], lam_ridge) else g_full
    if (!is.null(units)) g <- g * units
    if (engine == "lars") {
      zm <- colMeans(Zs[tr, , drop = FALSE])
      rm <- mean(r[tr])
      Zc <- sweep(Zs[tr, , drop = FALSE], 2L, zm, "-")
      G <- crossprod(Zc)
      cv <- drop(crossprod(Zc, r[tr] - rm))
      Zte <- sweep(Zs[te, , drop = FALSE], 2L, zm, "-")
      n_tr <- sum(tr)
    }
    lapply(seq_along(psis), function(i) {
      sc <- adaptive_scale(g, psis[i])
      if (engine == "lars") {
        B <- lars_lasso_coefs(G * outer(sc, sc), cv * sc, n_tr, grids[[i]])
        pred <- Zte %*% (sc * B) + rm
      } else {
        Ztil <- sweep(Zs, 2L, sc, "*")
        p <- lasso_path(Ztil[tr, , drop = FALSE], r[tr], grids[[i]], thresh)
        pred <- sweep(Ztil[te, , drop = FALSE], 2L, p$zm, "-") %*% p$coef + p$rm
      }
      colMeans((r[te] - pred)^2)
    })
  })
  lapply(seq_along(psis), function(i) do.call(rbind, lapply(per_fold, `[[`, i)))
}

# Lasso fit of r on Ztil at the penalties `lambdas` with either engine.
lasso_fit <- function(Ztil, r, lambdas, engine = "lars", thresh = 1e-10) {
  if (engine == "glmnet") {
    return(lasso_path(Ztil, r, lambdas, thresh))
  }
  zm <- colMeans(Ztil)
  rm <- mean(r)
  Zc <- sweep(Ztil, 2L, zm, "-")
  B <- lars_lasso_coefs(crossprod(Zc), drop(crossprod(Zc, r - rm)), nrow(Ztil), lambdas)
  list(coef = B, zm = zm, rm = rm)
}

constant_skedastic <- function(r, d, names, method, reason) {
  rbar <- mean(r)
  list(
    fitted = rep(rbar, length(r)), intercept = rbar, coef = stats::setNames(numeric(d), names),
    names = names, df = 1, method = method, selected = character(0),
    tuning = list(reason = reason)
  )
}

# Adaptive Lasso of r on Z with (lambda, psi) tuned jointly by K-fold CV.
adaptive_lasso <- function(r, Z, names = colnames(Z), psi_grid = .psi_grid, cv_folds = 10L,
                           lambda_rule = "min", fold_scheme = "random", n_lambda = 100L,
                           lambda_min_ratio = NULL, seed = NULL, foldid = NULL,
                           cv_weights = "fold", ridge_scale = "standardized",
                           engine = "lars", thresh = 1e-10) {
  n <- length(r)
  d <- ncol(Z)
  if (is.null(names)) names <- sprintf("z%d", seq_len(d))
  psis <- sort(unique(as.numeric(psi_grid)))
  if (d == 0L) {
    return(constant_skedastic(r, d, names, "lasso", "no candidate covariates"))
  }
  if (is.null(foldid)) {
    K <- min(cv_folds, n)
    folds <- make_folds(n, K, fold_scheme, seed)
  } else {
    folds <- check_foldid(foldid, n)
    K <- max(folds)
  }
  if (n < 2L * K) {
    fgls_warn(
      sprintf(
        "Only %d observations for %d-fold cross-validation; the skedastic fit is noisy.",
        n, K
      ),
      "sample_size"
    )
  }
  st <- standardize_cols(Z)
  Zs <- st$Zs
  keep <- st$keep
  if (ncol(Zs) == 0L) {
    out <- constant_skedastic(r, d, names, "lasso", "all candidate covariates are constant")
    out$tuning$foldid <- folds
    return(out)
  }
  if (stats::sd(r) == 0) {
    out <- constant_skedastic(r, d, names, "lasso", "the response is constant")
    out$tuning$foldid <- folds
    return(out)
  }
  rid <- ridge_cv(Zs, r, folds)
  g_ridge <- rid$coef
  lam_ridge <- rid$lambda
  units <- if (ridge_scale == "original") 1 / st$sd[keep] else NULL
  g_w <- if (is.null(units)) g_ridge else g_ridge * units
  ratio <- if (!is.null(lambda_min_ratio)) lambda_min_ratio else if (n > d) 1e-4 else 1e-2
  rc <- r - mean(r)
  grids <- lapply(psis, function(psi) {
    Ztil <- sweep(Zs, 2L, adaptive_scale(g_w, psi), "*")
    lam_max <- max(abs(crossprod(sweep(Ztil, 2L, colMeans(Ztil), "-"), rc))) / n
    if (!(lam_max > 0)) lam_max <- 1
    lam_max * 10^seq(0, log10(ratio), length.out = n_lambda)
  })
  curves <- cv_errors(Zs, r, folds, psis, grids, g_ridge, lam_ridge,
    fold_weights = cv_weights == "fold", units = units, engine = engine, thresh = thresh
  )
  means <- lapply(curves, colMeans)
  mins <- vapply(means, min, numeric(1))
  i_psi <- which.min(mins)
  i_lam <- which.min(means[[i_psi]])
  cv_min <- mins[i_psi]
  if (lambda_rule == "1se") {
    cm <- curves[[i_psi]]
    se <- apply(cm, 2L, stats::sd) / sqrt(nrow(cm))
    i_lam <- min(which(means[[i_psi]] <= means[[i_psi]][i_lam] + se[i_lam]))
  }
  psi_star <- psis[i_psi]
  alphas <- grids[[i_psi]]
  lam_star <- alphas[i_lam]
  scale <- adaptive_scale(g_w, psi_star)
  p <- lasso_fit(sweep(Zs, 2L, scale, "*"), r, alphas[seq_len(i_lam)], engine, thresh)
  gtil <- p$coef[, i_lam]
  g_std <- gtil * scale # coefficients on standardized z
  coef <- numeric(d)
  coef[keep] <- g_std / st$sd[keep]
  names(coef) <- names
  intercept <- p$rm - sum(p$zm * gtil) - sum(coef[keep] * st$mu[keep])
  fitted <- drop(intercept + Z %*% coef)
  list(
    fitted = fitted, intercept = intercept, coef = coef, names = names,
    df = sum(coef != 0) + 1, method = "lasso", selected = names[coef != 0],
    tuning = list(
      lambda = lam_star, lambda_paper = n * lam_star, psi = psi_star, cv_mse = cv_min,
      lambda_ridge = lam_ridge, ridge_coef = stats::setNames(g_ridge, names[keep]),
      lambda_rule = lambda_rule, cv_folds = K, fold_scheme = fold_scheme,
      cv_weights = cv_weights, engine = engine, foldid = folds,
      cv_path = lapply(seq_along(psis), function(i) {
        list(psi = psis[i], lambda = grids[[i]], cv_mse = means[[i]])
      })
    )
  )
}
