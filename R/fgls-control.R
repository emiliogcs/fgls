#' Advanced options for fgls
#'
#' Tuning options of [fgls()] that most users never change. The defaults
#' follow the papers and the Python package 'pyfgls'.
#'
#' @param ic_penalty Penalty of the lag-order criterion: `"order"` (default,
#'   PG26 eq. 14) penalises the AR order; `"params"` penalises every
#'   coefficient a lag adds to the Durbin regression, which avoids
#'   over-selection with many stochastic regressors.
#' @param min_lags Lower bound on the selected AR order: `NULL` (default), a
#'   non-negative integer, or `"adjusted"` for
#'   \eqn{\lceil(\ln T)^{3/2}\rceil} (PG26 Supplement, Lemma A.2).
#' @param fold_scheme Cross-validation folds of the skedastic fit: `"random"`
#'   (default) or `"blocked"` (contiguous blocks).
#' @param foldid Optional integer vector of fold labels (one per observation
#'   entering the skedastic fit), overriding `fold_scheme` and the seed.
#' @param cv_weights `"fold"` (default) re-estimates the ridge first step, and
#'   hence the adaptive weights, on every training fold; `"full"` computes the
#'   weights once on the full sample (the common practice with
#'   `glmnet::cv.glmnet()` and `penalty.factor`).
#' @param ridge_scale Scale on which the ridge coefficients define the adaptive
#'   weights: `"standardized"` (default) or `"original"` units of `z`.
#' @param n_lambda Length of the Lasso penalty grid.
#' @param lambda_min_ratio Smallest penalty relative to the largest; default
#'   `1e-4` when there are more observations than candidate covariates and
#'   `1e-2` otherwise (as in 'glmnet').
#' @param lasso_engine Solver of the Lasso problems: `"lars"` (default) computes
#'   the exact solution path by the LARS homotopy with the Lasso modification,
#'   as 'pyfgls' does; `"glmnet"` uses [glmnet::glmnet()] (coordinate descent)
#'   followed by an exact active-set refinement. Both give the same solutions;
#'   `"lars"` is faster for the moderate numbers of candidate covariates
#'   typical here.
#' @param glmnet_thresh Convergence threshold passed to [glmnet::glmnet()]
#'   when `lasso_engine = "glmnet"`.
#' @param delta_scale Use the scale-aware floor `delta * sigma_OLS` in
#'   \eqn{\log\max(\hat u_t^2, \delta^2)}.
#' @param transforms Transform families used by `z = "auto"`; see
#'   [make_dictionary()].
#' @param svr_cost,svr_epsilon Grids for the SVR skedastic fit, multiplied by
#'   the standard deviation of the response.
#' @param n_boot Number of bootstrap replications for `vcov = "bootstrap"`.
#' @param bootstrap Innovation scheme of the sieve bootstrap: `"wild"`
#'   (default) or `"iid"` resampling.
#' @param wild_dist Wild-bootstrap multipliers: `"rademacher"` (default) or
#'   `"mammen"`.
#' @param boot_seed Seed of the bootstrap; default the `seed` of the fit.
#' @param cores Number of cores for the bootstrap (forked processes via
#'   'parallel'; ignored on Windows). Defaults to 1.
#'
#' @return A list of class `"fgls_control"`.
#' @family fgls
#' @export
#' @examples
#' fgls_control(ic_penalty = "params", n_boot = 199)
fgls_control <- function(ic_penalty = c("order", "params"), min_lags = NULL,
                         fold_scheme = c("random", "blocked"), foldid = NULL,
                         cv_weights = c("fold", "full"),
                         ridge_scale = c("standardized", "original"),
                         n_lambda = 100L, lambda_min_ratio = NULL,
                         lasso_engine = c("lars", "glmnet"), glmnet_thresh = 1e-7,
                         delta_scale = FALSE,
                         transforms = c("linear", "poly2", "log", "trig"),
                         svr_cost = c(0.1, 1, 10, 100), svr_epsilon = c(0.05, 0.25, 0.5, 1),
                         n_boot = 999L, bootstrap = c("wild", "iid"),
                         wild_dist = c("rademacher", "mammen"), boot_seed = NULL,
                         cores = 1L) {
  ic_penalty <- arg_default(ic_penalty, c("order", "params"), "ic_penalty")
  fold_scheme <- arg_default(fold_scheme, c("random", "blocked"), "fold_scheme")
  cv_weights <- arg_default(cv_weights, c("fold", "full"), "cv_weights")
  ridge_scale <- arg_default(ridge_scale, c("standardized", "original"), "ridge_scale")
  lasso_engine <- arg_default(lasso_engine, c("lars", "glmnet"), "lasso_engine")
  bootstrap <- arg_default(bootstrap, c("wild", "iid"), "bootstrap")
  wild_dist <- arg_default(wild_dist, c("rademacher", "mammen"), "wild_dist")
  bad <- setdiff(transforms, c("linear", "poly2", "log", "trig"))
  if (length(bad)) {
    fgls_abort(sprintf(
      "Unknown transforms %s; choose from \"linear\", \"poly2\", \"log\", \"trig\".",
      paste(bad, collapse = ", ")
    ))
  }
  if (!is_count(n_lambda) || n_lambda < 2L) fgls_abort("`n_lambda` must be an integer >= 2.")
  if (!is_count(n_boot) || n_boot < 2L) fgls_abort("`n_boot` must be an integer >= 2.")
  if (!is_count(cores) || cores < 1L) fgls_abort("`cores` must be a positive integer.")
  if (!is.null(lambda_min_ratio) &&
    !(is.numeric(lambda_min_ratio) && lambda_min_ratio > 0 && lambda_min_ratio < 1)) {
    fgls_abort("`lambda_min_ratio` must be NULL or in (0, 1).")
  }
  if (!is.null(foldid) && !(is.numeric(foldid) && all(foldid == round(foldid)))) {
    fgls_abort("`foldid` must be an integer vector of fold labels.")
  }
  check_seed(boot_seed)
  structure(
    list(
      ic_penalty = ic_penalty, min_lags = min_lags, fold_scheme = fold_scheme,
      foldid = if (is.null(foldid)) NULL else as.integer(foldid),
      cv_weights = cv_weights, ridge_scale = ridge_scale, n_lambda = as.integer(n_lambda),
      lambda_min_ratio = lambda_min_ratio, lasso_engine = lasso_engine,
      glmnet_thresh = glmnet_thresh,
      delta_scale = isTRUE(delta_scale), transforms = transforms,
      svr_cost = svr_cost, svr_epsilon = svr_epsilon, n_boot = as.integer(n_boot),
      bootstrap = bootstrap, wild_dist = wild_dist, boot_seed = boot_seed,
      cores = as.integer(cores)
    ),
    class = "fgls_control"
  )
}
