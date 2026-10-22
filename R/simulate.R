# Data-generating processes of the simulation studies of the two papers
# (port of pyfgls/simulate.py; R's random number generator is used, so draws
# differ from the Python package while the designs are the same).

# u_t = sum_i ar_i u_{t-i} + e_t + sum_i ma_i e_{t-i}, pre-sample values zero
# (the u_0 = 0 initialisation of PG26 Table 1).
arma_errors <- function(e, ar = numeric(0), ma = numeric(0)) {
  n <- length(e)
  v <- as.numeric(e)
  for (i in seq_along(ma)) {
    if (i < n) v[(i + 1L):n] <- v[(i + 1L):n] + ma[i] * e[seq_len(n - i)]
  }
  if (length(ar)) {
    v <- as.numeric(stats::filter(v, ar, method = "recursive"))
  }
  v
}

skedastic_functions <- list(
  homoskedastic = function(x) rep(1, length(x)),
  x = function(x) x,
  x2 = function(x) x^2,
  log2 = function(x) log(x)^2,
  exp = function(x) exp(0.2 * x + 0.2 * x^2),
  step = function(x) ifelse(x < 2, 1, ifelse(x < 3, 4, 9))
)

resolve_skedastic <- function(v) {
  if (is.function(v)) {
    return(v)
  }
  if (!is.character(v) || length(v) != 1L || !(v %in% names(skedastic_functions))) {
    fgls_abort(sprintf(
      "Unknown skedastic function; choose from %s or pass a function of x.",
      paste0("\"", names(skedastic_functions), "\"", collapse = ", ")
    ))
  }
  skedastic_functions[[v]]
}

new_simulation <- function(df, beta, u, e, info) {
  attr(df, "beta") <- beta
  attr(df, "u") <- u
  attr(df, "e") <- e
  attr(df, "info") <- info
  df
}

#' Simulate the designs of the papers
#'
#' Data-generating processes used in the Monte Carlo studies of Perron and
#' González-Coya (2026) (`simulate_serial()`, `simulate_combined()`,
#' `simulate_predictive()`) and González-Coya and Perron (2024)
#' (`simulate_hetero()`). They are handy for examples, teaching and
#' replication.
#'
#' * `simulate_serial()`: \eqn{y_t = \beta_0 + \beta_1 x_t + u_t} with ARMA
#'   errors driven by standard normal innovations \eqn{e_t} and a regressor
#'   \eqn{x_t = \rho_x x_{t-1} + v_t + \gamma e_{t-1}}, which is pre-determined
#'   and exogenous only when \eqn{\gamma = 0}. The defaults (\eqn{\rho_x = 0},
#'   \eqn{\gamma = 1}, AR(1) coefficient 0.8, no burn-in, \eqn{u_0 = 0}) give the
#'   design of Table 1 of the 2026 paper.
#' * `simulate_hetero()`: \eqn{y_t = \alpha + \beta x_t + \sqrt{v(x_t)}
#'   \varepsilon_t}, \eqn{x_t \sim U[1, 4]} (eq. 9 of the 2024 paper), optionally
#'   with extra candidate covariates \eqn{w_j = c x + (1 - c) q_j},
#'   \eqn{q_j \sim U[1, 4]}.
#' * `simulate_combined()`: ARMA errors whose innovations are heteroskedastic,
#'   \eqn{e_t = \sqrt{v(w_t)} \varepsilon_t}, with \eqn{w_t = x_t} (`driver =
#'   "x"`) or \eqn{w_t = \phi_t x_t + (1 - \phi_t) q_t},
#'   \eqn{\phi_t \sim} Bernoulli(`phi_prob`) (`driver = "mixture"`)
#'   (Supplement S-4).
#' * `simulate_predictive()`: the predictive regression
#'   \eqn{y_{t+2} = \alpha + \beta x_t + u_{t+2}} with MA(1) errors
#'   (Supplement S-3), returned aligned so that row \eqn{t} holds
#'   \eqn{y_{t+2}} and \eqn{x_t}.
#'
#' @param n Number of observations returned.
#' @param ar,ma Autoregressive and moving-average coefficients of the errors.
#' @param rho_x Autoregressive coefficient of the regressor.
#' @param gamma Feedback from past innovations to the regressor; the regressor
#'   is not exogenous when `gamma != 0`.
#' @param beta Intercept and slope.
#' @param burn Number of initial observations discarded.
#' @param skedastic Skedastic function: one of `"homoskedastic"`, `"x"`,
#'   `"x2"`, `"log2"` (\eqn{[\log x]^2}), `"exp"` (\eqn{\exp(0.2x + 0.2x^2)}),
#'   `"step"`, or a function of `x`.
#' @param n_extra Number of extra candidate covariates `w1, w2, ...` (columns of
#'   the result, not regressors).
#' @param extra_corr Weight \eqn{c} of `x` in the extra covariates.
#' @param driver Variance driver of `simulate_combined()`: `"x"` or
#'   `"mixture"`.
#' @param phi_prob Probability that the variance driver equals `x` when
#'   `driver = "mixture"`.
#' @param theta MA(1) coefficient of `simulate_predictive()`.
#' @param seed Optional seed. The random number generator state of the session
#'   is restored afterwards.
#'
#' @return A data frame with columns `y` and `x` (and `w` or `w1, w2, ...` where
#'   relevant), with attributes `beta` (true coefficients), `u` (errors), `e`
#'   (innovations) and `info` (design parameters and, for the heteroskedastic
#'   designs, the true variances `v`).
#' @references
#' Perron, P. and González-Coya, E. (2026). Feasible GLS for time series
#' regression. *Journal of Econometrics*, 258, 106346.
#' \doi{10.1016/j.jeconom.2026.106346}
#'
#' González-Coya, E. and Perron, P. (2024). Estimation in the presence of
#' heteroskedasticity of unknown form: a Lasso-based approach. *Journal of
#' Econometric Methods*, 13(1), 29--48. \doi{10.1515/jem-2023-0007}
#' @family simulation
#' @export
#' @examples
#' sim <- simulate_serial(200, seed = 1)
#' head(sim)
#' attr(sim, "beta")
#'
#' het <- simulate_hetero(100, skedastic = "log2", seed = 2)
#' summary(attr(het, "info")$v)
simulate_serial <- function(n, ar = 0.8, ma = numeric(0), rho_x = 0, gamma = 1,
                            beta = c(1, 1), burn = 0, seed = NULL) {
  check_seed(seed)
  n <- as.integer(n)
  burn <- as.integer(burn)
  with_seed(seed, {
    m <- n + burn
    e_all <- stats::rnorm(m + 1L) # e_0, ..., e_m
    v <- stats::rnorm(m)
  })
  e <- e_all[-1L]
  u <- arma_errors(e, ar, ma)
  x <- as.numeric(stats::filter(v + gamma * e_all[-(m + 1L)], rho_x, method = "recursive"))
  keep <- (burn + 1L):m
  y <- beta[1L] + beta[2L] * x[keep] + u[keep]
  new_simulation(
    data.frame(y = y, x = x[keep]),
    beta = beta, u = u[keep], e = e[keep],
    info = list(ar = ar, ma = ma, rho_x = rho_x, gamma = gamma)
  )
}

#' @rdname simulate_serial
#' @export
simulate_hetero <- function(n, skedastic = "x2", beta = c(1, 1), n_extra = 0L,
                            extra_corr = 0, seed = NULL) {
  check_seed(seed)
  fun <- resolve_skedastic(skedastic)
  n <- as.integer(n)
  with_seed(seed, {
    x <- stats::runif(n, 1, 4)
    eps <- stats::rnorm(n)
    Q <- if (n_extra > 0L) matrix(stats::runif(n * n_extra, 1, 4), n, n_extra)
  })
  v <- fun(x)
  u <- sqrt(v) * eps
  df <- data.frame(y = beta[1L] + beta[2L] * x + u, x = x)
  if (n_extra > 0L) {
    W <- extra_corr * x + (1 - extra_corr) * Q
    colnames(W) <- paste0("w", seq_len(n_extra))
    df <- cbind(df, as.data.frame(W))
  }
  new_simulation(df,
    beta = beta, u = u, e = eps,
    info = list(v = v, skedastic = if (is.character(skedastic)) skedastic else "custom")
  )
}

#' @rdname simulate_serial
#' @export
simulate_combined <- function(n, ar = 0.5, ma = numeric(0), skedastic = "x2",
                              beta = c(1, 1), driver = c("x", "mixture"),
                              phi_prob = 0.75, gamma = 0, burn = 100L, seed = NULL) {
  check_seed(seed)
  driver <- arg_default(driver, c("x", "mixture"), "driver")
  fun <- resolve_skedastic(skedastic)
  n <- as.integer(n)
  burn <- as.integer(burn)
  m <- n + burn
  with_seed(seed, {
    x <- stats::runif(m, 1, 4)
    if (driver == "x") {
      w <- x
    } else {
      phi <- stats::rbinom(m, 1L, phi_prob)
      w <- phi * x + (1 - phi) * stats::runif(m, 1, 4)
    }
    eps <- stats::rnorm(m)
  })
  v <- fun(w)
  e <- sqrt(v) * eps
  u <- arma_errors(e, ar, ma)
  xr <- x
  if (gamma != 0) xr[-1L] <- xr[-1L] + gamma * e[-m]
  keep <- (burn + 1L):m
  new_simulation(
    data.frame(y = beta[1L] + beta[2L] * xr[keep] + u[keep], x = xr[keep], w = w[keep]),
    beta = beta, u = u[keep], e = e[keep],
    info = list(v = v[keep], ar = ar, ma = ma, gamma = gamma, driver = driver)
  )
}

#' @rdname simulate_serial
#' @export
simulate_predictive <- function(n, theta = 0.5, rho_x = 0, gamma = 0, beta = c(1, 1),
                                burn = 100L, seed = NULL) {
  check_seed(seed)
  n <- as.integer(n)
  burn <- as.integer(burn)
  m <- n + burn + 2L
  with_seed(seed, {
    e <- stats::rnorm(m + 1L)
    v <- stats::rnorm(m)
  })
  x <- as.numeric(stats::filter(v + gamma * e[-(m + 1L)], rho_x, method = "recursive"))
  u <- e[-1L] + theta * e[-(m + 1L)] # u_t for t = 1..m
  xs <- x[seq_len(m - 2L)]
  y <- beta[1L] + beta[2L] * xs + u[-(1:2)]
  keep <- (burn + 1L):(m - 2L)
  new_simulation(
    data.frame(y = y[keep], x = xs[keep]),
    beta = beta, u = u[-(1:2)][keep], e = e[-1L][-(1:2)][keep],
    info = list(theta = theta, rho_x = rho_x, gamma = gamma)
  )
}
