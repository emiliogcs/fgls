# Exact Lasso path by the LARS homotopy with the Lasso modification (Efron,
# Hastie, Johnstone and Tibshirani, 2004), on a Gram matrix.
#
# For centred data the Lasso minimises (2n)^{-1} ||rc - Zc b||^2 + a ||b||_1.
# With G = Zc'Zc, c = Zc'rc and correlations C(b) = c - G b, the solution
# satisfies C_j = n a sign(b_j) on the active set and |C_j| <= n a elsewhere,
# and is piecewise linear in a. Starting from b = 0 at a_max = max|c|/n, the
# homotopy moves along the equiangular direction d_A = G_AA^{-1} s_A until a
# variable joins (an inactive correlation reaches the active level) or a
# coefficient crosses zero (it leaves the active set), and records the knots.
# The coefficients at any penalty are then obtained by linear interpolation
# between knots, exactly. This is what pyfgls does with scikit-learn's
# lars_path_gram, without its absolute stopping threshold at float32 machine
# epsilon (~1.2e-7), which pyfgls hits for very small penalties.

# Knots of the exact path down to a_min: list(alpha (decreasing), coef (d x K)).
lars_lasso_knots <- function(G, cvec, n, a_min = 0) {
  d <- length(cvec)
  b <- numeric(d)
  C <- cvec
  Cmax <- if (d) max(abs(C)) else 0
  knots_a <- Cmax / n
  knots_b <- list(b)
  if (!(Cmax > 0)) {
    return(list(alpha = knots_a, coef = matrix(b, d, 1L)))
  }
  eps <- 1e-12 * Cmax
  active <- integer(0)
  sgn <- numeric(d)
  # a variable that just left the active set has correlation s * Cmax; its
  # same-sign re-entry condition is linear in the step with root 0, so only its
  # opposite-sign root is a genuine event
  just_dropped <- integer(0)
  drop_sign <- 0
  tiny <- .Machine$double.eps * 100
  # variables to add: all those at the maximal correlation
  add <- which(abs(C) >= Cmax - eps)
  for (iter in seq_len(8L * d + 50L)) {
    if (length(add)) {
      active <- sort(c(active, add))
      sgn[add] <- sign(C[add])
      add <- integer(0)
    }
    if (!length(active)) break
    # recompute the correlations exactly at each knot (no drift)
    C <- cvec - drop(G %*% b)
    Cmax <- mean(abs(C[active]))
    GA <- G[active, active, drop = FALSE]
    dA <- tryCatch(solve(GA, sgn[active]), error = function(e) NULL)
    if (is.null(dA) || any(!is.finite(dA))) break # singular active set (e.g. d > n saturation)
    a <- drop(G[, active, drop = FALSE] %*% dA)
    gam_end <- Cmax - n * a_min
    # joins
    inact <- setdiff(seq_len(d), active)
    gam_join <- Inf
    j_join <- NA_integer_
    if (length(inact)) {
      cj <- C[inact]
      aj <- a[inact]
      g1 <- (Cmax - cj) / (1 - aj)
      g2 <- (Cmax + cj) / (1 + aj)
      g1[!is.finite(g1) | g1 <= tiny * Cmax] <- Inf
      g2[!is.finite(g2) | g2 <= tiny * Cmax] <- Inf
      if (length(just_dropped)) {
        jd <- match(just_dropped, inact)
        if (drop_sign > 0) g1[jd] <- Inf else g2[jd] <- Inf
      }
      gj <- pmin(g1, g2)
      if (any(is.finite(gj))) {
        k <- which.min(gj)
        gam_join <- gj[k]
        j_join <- inact[k]
      }
    }
    # drops (coefficient crosses zero)
    gd <- -b[active] / dA
    gd[!(gd > tiny * Cmax)] <- Inf
    gam_drop <- if (length(gd)) min(gd) else Inf
    k_drop <- if (is.finite(gam_drop)) active[which.min(gd)] else NA_integer_
    gam <- min(gam_join, gam_drop, gam_end)
    b[active] <- b[active] + gam * dA
    C <- C - gam * a
    Cmax <- Cmax - gam
    if (gam == gam_end) {
      knots_a <- c(knots_a, Cmax / n)
      knots_b[[length(knots_b) + 1L]] <- b
      break
    }
    if (gam == gam_drop) {
      drop_sign <- sgn[k_drop]
      b[k_drop] <- 0
      sgn[k_drop] <- 0
      active <- setdiff(active, k_drop)
      just_dropped <- k_drop
    } else {
      add <- j_join
      just_dropped <- integer(0)
    }
    knots_a <- c(knots_a, Cmax / n)
    knots_b[[length(knots_b) + 1L]] <- b
    if (Cmax <= n * a_min || Cmax <= eps) break
  }
  list(alpha = knots_a, coef = do.call(cbind, knots_b))
}

# Coefficients (d x L) at the penalties `alphas` by interpolation along the
# exact path (zero above a_max; the last knot below the end of the path).
lars_lasso_coefs <- function(G, cvec, n, alphas) {
  d <- length(cvec)
  if (d == 0L) {
    return(matrix(0, 0L, length(alphas)))
  }
  kn <- lars_lasso_knots(G, cvec, n, a_min = min(alphas))
  ka <- kn$alpha
  KB <- kn$coef
  if (length(ka) == 1L) {
    return(matrix(0, d, length(alphas)))
  }
  out <- matrix(0, d, length(alphas))
  # ka is decreasing; locate each alpha in the segments
  for (l in seq_along(alphas)) {
    al <- alphas[l]
    if (al >= ka[1L]) next
    if (al <= ka[length(ka)]) {
      out[, l] <- KB[, length(ka)]
      next
    }
    i <- max(which(ka >= al)) # segment [ka[i], ka[i + 1]]
    w <- (ka[i] - al) / (ka[i] - ka[i + 1L])
    out[, l] <- KB[, i] + w * (KB[, i + 1L] - KB[, i])
  }
  # safeguard: any solution failing the KKT conditions is refined exactly
  tol <- 1e-8 * pmax(alphas, max(abs(cvec)) / n)
  if (any(kkt_gaps(G, cvec, n, alphas, out) > tol)) out <- polish_lasso(G, cvec, n, alphas, out)
  out
}
