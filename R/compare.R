#' Compare fits side by side
#'
#' A table of coefficients, standard errors and confidence-interval widths for
#' several fits, e.g. FGLS against OLS with HAC standard errors (PG26) or OLS
#' with HC3 standard errors (GP24), or the F-C and F-CH stages of a combined
#' fit.
#'
#' Plain [stats::lm()] fits are accepted and converted to OLS with HAC
#' (Andrews quadratic-spectral) standard errors when any other fit corrects for
#' serial correlation, and with HC3 standard errors otherwise; override with
#' `ols_vcov`. Other models that inherit from `"lm"` (`glm`, multivariate
#' `mlm` fits, weighted fits, fits with an offset, ...) are rejected, because
#' the conversion would replace them by unweighted OLS; use [ols_robust()] for
#' OLS with robust standard errors.
#'
#' @param ... `"fgls"` (or plain `"lm"`) objects, optionally named; the names
#'   label the columns.
#' @param level Confidence level used for the interval widths.
#' @param ols_vcov Covariance type for `lm` fits: `NULL` (automatic), `"HAC"`,
#'   `"NW"` or `"HC0"`--`"HC3"`.
#' @return A data frame with columns `term`, `statistic` (`"estimate"`,
#'   `"std.error"`, `"ci.width"`) and one column per fit.
#' @family fgls
#' @export
#' @examples
#' sim <- simulate_serial(300, seed = 1)
#' compare_fgls(
#'   OLS = stats::lm(y ~ x, data = sim),
#'   FGLS = serial_fgls(y ~ x, data = sim)
#' )
compare_fgls <- function(..., level = 0.95, ols_vcov = NULL) {
  fits <- list(...)
  if (!length(fits)) fgls_abort("compare_fgls() needs at least one fit.")
  labels <- names(fits)
  if (is.null(labels)) labels <- rep("", length(fits))
  any_serial <- any(vapply(fits, function(f) inherits(f, "fgls") && !is.null(f$serial), logical(1)))
  fits <- lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    if (inherits(f, "fgls")) {
      return(f)
    }
    if (inherits(f, "lm")) {
      check_plain_lm(f, if (nzchar(labels[i])) sprintf("`%s`", labels[i]) else sprintf("fit %d", i))
      return(lm_as_fgls(f, if (is.null(ols_vcov)) (if (any_serial) "HAC" else "HC3") else ols_vcov))
    }
    fgls_abort("compare_fgls() accepts fgls and lm fits only.")
  })
  auto <- vapply(fits, function(f) f$method, character(1))
  dup <- auto %in% auto[duplicated(auto)]
  vtypes <- vapply(fits[dup], function(f) f$vcov_type, character(1))
  auto[dup] <- sprintf("%s [%s]", auto[dup], vtypes)
  labels[labels == ""] <- auto[labels == ""]
  labels <- make.unique(labels, sep = " ")
  terms <- unique(unlist(lapply(fits, function(f) names(f$coefficients))))
  stats_ <- c("estimate", "std.error", "ci.width")
  out <- data.frame(
    term = rep(terms, each = 3L), statistic = rep(stats_, length(terms)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(fits)) {
    f <- fits[[i]]
    se <- sqrt(pmax(diag(f$vcov), 0))
    w <- 2 * critical_value(1 - level, fgls_df(f)) * se
    vals <- rbind(f$coefficients, se, w)
    col <- rep(NA_real_, nrow(out))
    idx <- match(out$term, names(f$coefficients))
    ok <- !is.na(idx)
    col[ok] <- vals[cbind(match(out$statistic[ok], stats_), idx[ok])]
    out[[labels[i]]] <- col
  }
  out
}

# compare_fgls() refits an lm object by unweighted OLS from its model matrix and
# response, so it accepts only plain lm() (or aov()) fits: no glm, mlm, weights
# or offset, which the refit would silently drop.
check_plain_lm <- function(f, what) {
  why <- if (inherits(f, "glm")) {
    "a glm fit"
  } else if (inherits(f, "mlm")) {
    "a multivariate lm (mlm) fit"
  } else if (!(class(f)[1L] %in% c("lm", "aov"))) {
    sprintf("a \"%s\" fit", class(f)[1L])
  } else if (!is.null(stats::weights(f))) {
    "a weighted lm fit"
  } else if (!is.null(f$offset)) {
    "an lm fit with an offset"
  } else {
    NULL
  }
  if (!is.null(why)) {
    fgls_abort(sprintf(
      paste0(
        "compare_fgls() converts only plain lm() fits (no weights, offset or glm family) ",
        "to OLS with robust standard errors, but %s is %s. Pass an \"fgls\" fit instead, ",
        "e.g. ols_robust() for OLS with HAC or HC standard errors."
      ),
      what, why
    ))
  }
  invisible(f)
}

# OLS (from an lm fit) with HAC or HC standard errors, as an "fgls"-like list.
lm_as_fgls <- function(fit, vcov_type) {
  arg_choice(vcov_type, c("HAC", "NW", "HC0", "HC1", "HC2", "HC3"), "ols_vcov")
  X <- stats::model.matrix(fit)
  y <- stats::model.response(stats::model.frame(fit))
  const <- which(apply(X, 2L, function(v) all(v == v[1L]) && v[1L] != 0))
  o <- ols_core(X, y, vcov_type, const_col = if (length(const)) const[1L] else NA_integer_)
  names(o$beta) <- colnames(X)
  dimnames(o$cov) <- list(colnames(X), colnames(X))
  structure(list(
    coefficients = o$beta, vcov = o$cov, method = "OLS", vcov_type = o$label,
    vcov_code = vcov_type, use_t = FALSE, df_resid = nrow(X) - ncol(X), serial = NULL
  ), class = "fgls")
}
