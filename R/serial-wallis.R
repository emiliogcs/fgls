# AR coefficients from the Durbin regression when lagged dependent variables
# are regressors, PG26 Section 4.2 (port of pyfgls/serial/wallis.py).
#
# With lags of y among the regressors the coefficients on the lagged y in the
# Durbin regression no longer estimate rho. Under the common-factor structure
# the coefficient on x_{i,t-j} is -rho_j beta_i for every regressor that is not
# a lagged dependent variable, so rho_j = delta_{j,i} / beta_i (Wallis, 1967).
#
# Choice of i (`wallis_regressor`): "auto" (largest |t| of the contemporaneous
# Durbin coefficient), "median" (median of the ratios, lag by lag) or a
# regressor name. This path is less tested than the standard one.

.weak_t <- 2

# Regressors whose lags 1..k are all in the Durbin design (and not excluded).
eligible_regressors <- function(fit, k, exclude) {
  d <- fit$design
  Filter(function(i) {
    !(i %in% exclude) && all(!is.na(d$xlag_pos[seq_len(k), i]))
  }, seq_len(d$n_x))
}

wallis_single <- function(fit, k, i) {
  b <- fit$coef
  d <- fit$design
  bi <- b[[i]]
  pos <- d$xlag_pos[seq_len(k), i]
  rho <- -b[pos] / bi
  J <- matrix(0, k, length(b))
  for (j in seq_len(k)) {
    J[j, pos[j]] <- -1 / bi
    J[j, i] <- -rho[j] / bi
  }
  list(rho = unname(rho), J = J)
}

# Recover rho_1..rho_k by the Wallis ratio. Returns rho, J (Jacobian with
# respect to the Durbin coefficients), the regressor(s) used and the ratios
# implied by every eligible regressor.
wallis_rho <- function(fit, k, x_names, lagged_dep_idx, choice = "auto") {
  p <- length(fit$coef)
  if (k == 0L) {
    return(list(rho = numeric(0), J = matrix(0, 0L, p), used = character(0), ratios = list()))
  }
  elig <- eligible_regressors(fit, k, lagged_dep_idx)
  if (length(elig) == 0L) {
    fgls_abort(paste0(
      "Cannot recover the AR coefficients with lagged dependent variables: the ",
      "Wallis ratio needs at least one stochastic regressor that is not a lag of y ",
      "(deterministic terms such as the constant are not eligible)."
    ))
  }
  ratios <- stats::setNames(lapply(elig, function(i) wallis_single(fit, k, i)$rho), x_names[elig])
  tvals <- abs(fit$tvalues[elig])
  if (identical(choice, "auto")) {
    pick <- elig[which.max(tvals)]
  } else if (identical(choice, "median")) {
    pick <- elig
  } else {
    i <- match(choice, x_names)
    if (is.na(i)) fgls_abort(sprintf("wallis_regressor '%s' is not a regressor name.", choice))
    if (!(i %in% elig)) {
      fgls_abort(sprintf(
        paste0(
          "wallis_regressor '%s' is not eligible (it is a lagged dependent variable ",
          "or a deterministic term whose lags are collinear)."
        ),
        choice
      ))
    }
    pick <- i
  }
  weak <- x_names[pick][abs(fit$tvalues[pick]) < .weak_t]
  if (length(weak)) {
    fgls_warn(
      sprintf(
        paste0(
          "The Wallis ratio uses regressor(s) %s whose Durbin coefficient has |t| < 2; ",
          "the recovered AR coefficients are weakly identified. Consider ",
          "wallis_regressor = \"median\" or another regressor."
        ),
        paste(weak, collapse = ", ")
      ),
      "specification"
    )
  }
  if (length(pick) == 1L) {
    s <- wallis_single(fit, k, pick)
    return(list(rho = s$rho, J = s$J, used = x_names[pick], ratios = ratios))
  }
  singles <- lapply(pick, function(i) wallis_single(fit, k, i))
  R <- do.call(rbind, lapply(singles, `[[`, "rho")) # (m, k)
  rho <- apply(R, 2L, stats::median)
  J <- matrix(0, k, p)
  m <- nrow(R)
  for (j in seq_len(k)) {
    ord <- order(R[, j], method = "radix")
    mids <- if (m %% 2L == 1L) ord[m %/% 2L + 1L] else ord[c(m %/% 2L, m %/% 2L + 1L)]
    for (s in mids) J[j, ] <- J[j, ] + singles[[s]]$J[j, ] / length(mids)
  }
  list(rho = rho, J = J, used = x_names[pick], ratios = ratios)
}
