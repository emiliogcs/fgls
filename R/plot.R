#' Diagnostic plots for fgls fits
#'
#' Base-graphics diagnostics: (1) residuals against fitted values; (2) the
#' autocorrelation function of the innovations (the residuals of the final
#' filtered and/or weighted regression, which should look like white noise);
#' (3) the lag-order criterion path with the selected order (serial
#' correction); (4) the fitted variances against the fitted values
#' (heteroskedasticity correction), or a normal quantile plot of the
#' innovations otherwise.
#'
#' @param x An `"fgls"` object.
#' @param which Panels to draw (subset of 1:4).
#' @param ask Ask before each new page.
#' @param lag.max Number of lags of the autocorrelation function.
#' @param ... Further graphical arguments passed to [plot()] (or
#'   [stats::qqnorm()] for panel 4). `main`, `xlab` and `ylab` replace the
#'   default titles and axis labels of every panel drawn.
#' @return `x`, invisibly.
#' @family fgls
#' @export
#' @examples
#' fit <- fgls(y ~ x, data = simulate_combined(200, seed = 1), seed = 1)
#' plot(fit)
plot.fgls <- function(x, which = 1:4, ask = FALSE, lag.max = 20L, ...) {
  which <- intersect(which, 1:4)
  if (!length(which)) fgls_abort("`which` must contain values in 1:4.")
  n_pan <- length(which)
  op <- graphics::par(
    mfrow = if (n_pan > 1L) c(ceiling(n_pan / 2), min(n_pan, 2L)) else c(1L, 1L),
    ask = ask && grDevices::dev.interactive()
  )
  on.exit(graphics::par(op), add = TRUE)
  inn <- x$innovations
  # Each panel takes its default titles as formal arguments, so that `main`,
  # `xlab` or `ylab` given in `...` replace them instead of being matched twice.
  panel_resid <- function(xlab = "Fitted values", ylab = "Residuals",
                          main = "Residuals vs fitted", ...) {
    graphics::plot(x$fitted.values, x$residuals, xlab = xlab, ylab = ylab, main = main, ...)
    graphics::abline(h = 0, lty = 3)
  }
  panel_acf <- function(main = "ACF of innovations", ...) {
    graphics::plot(stats::acf(inn, lag.max = lag.max, plot = FALSE), main = main, ...)
  }
  panel_ic <- function(xlab = "AR order k", ylab = toupper(x$serial$method),
                       main = "Lag-order criterion", ...) {
    bp <- x$bic_path
    if (is.null(bp)) {
      graphics::plot.new()
      graphics::title(main = main)
      graphics::text(0.5, 0.5, "no lag-order selection")
    } else {
      graphics::plot(as.numeric(names(bp)), bp,
        type = "b", xlab = xlab, ylab = ylab, main = main, ...
      )
      graphics::abline(v = x$ar_order, lty = 2)
    }
  }
  panel_sked <- function(xlab = "Fitted values", ylab = "Fitted variance (mean 1)",
                         main = "Skedastic fit", ...) {
    v <- x$skedastic$variance
    fv <- utils::tail(x$fitted.values, length(v))
    graphics::plot(fv, v, xlab = xlab, ylab = ylab, main = main, ...)
  }
  panel_qq <- function(main = "Normal Q-Q of innovations", ...) {
    stats::qqnorm(inn, main = main, ...)
    stats::qqline(inn, lty = 3)
  }
  for (w in which) {
    if (w == 1L) {
      panel_resid(...)
    } else if (w == 2L) {
      panel_acf(...)
    } else if (w == 3L) {
      panel_ic(...)
    } else if (!is.null(x$skedastic)) {
      panel_sked(...)
    } else {
      panel_qq(...)
    }
  }
  invisible(x)
}
