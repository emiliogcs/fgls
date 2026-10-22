# Input validation (port of pyfgls/_validation.py).
#
# The formula interface turns user input into a numeric response y and a model
# matrix X; prepare_data() then checks what makes FGLS results silently wrong:
# time order, gaps, missing and infinite values, constant and collinear columns.

.max_locations <- 5L

# Match a scalar option against its choices, with a classed error.
arg_choice <- function(value, choices, name) {
  if (is.null(value) || length(value) != 1L || is.na(value) || !(value %in% choices)) {
    shown <- if (is.null(value)) "NULL" else paste(format(value), collapse = ", ")
    fgls_abort(sprintf(
      "`%s` must be one of %s; got %s.", name,
      paste0("\"", choices, "\"", collapse = ", "), shown
    ))
  }
  value
}

# First element of a default choice vector, or the validated user value.
arg_default <- function(value, choices, name) {
  if (identical(value, choices)) {
    return(choices[[1L]])
  }
  arg_choice(value, choices, name)
}

format_locations <- function(index, positions) {
  pos <- utils::head(positions, .max_locations)
  labels <- if (!is.null(index)) format(index[pos]) else as.character(pos)
  more <- if (length(positions) > length(pos)) {
    sprintf(" (and %d more)", length(positions) - length(pos))
  } else {
    ""
  }
  where <- if (!is.null(index)) "time index" else "rows"
  paste0(where, " ", paste(labels, collapse = ", "), more)
}

# ---------------------------------------------------------------------- time index
modal_value <- function(x) {
  tab <- table(x)
  as.numeric(names(tab)[which.max(tab)])
}

# Positions i such that a gap occurs between observations i and i + 1.
time_index_gaps <- function(index) {
  n <- length(index)
  if (n < 2L) {
    return(integer(0))
  }
  if (inherits(index, "Date")) {
    d <- as.numeric(diff(index))
    if (all(d == d[1L])) {
      return(integer(0))
    }
    lt <- as.POSIXlt(index)
    months <- (lt$year + 1900) * 12 + lt$mon
    dm <- diff(months)
    if (all(dm == dm[1L]) && dm[1L] > 0) {
      return(integer(0)) # monthly, quarterly or annual data
    }
    wday <- lt$wday
    if (all(wday %in% 1:5) && modal_value(d) == 1) {
      # business-day data: a gap is a skipped weekday (holiday or missing day)
      weekdays_between <- vapply(seq_len(n - 1L), function(i) {
        if (d[i] <= 1) {
          return(0L)
        }
        s <- seq(index[i] + 1, index[i + 1L] - 1, by = "day")
        sum(as.POSIXlt(s)$wday %in% 1:5)
      }, integer(1))
      return(which(weekdays_between > 0L))
    }
    return(which(d > 1.5 * modal_value(d)))
  }
  v <- as.numeric(index)
  d <- diff(v)
  if (isTRUE(all.equal(d, rep(d[1L], length(d)), tolerance = 1e-8))) {
    return(integer(0))
  }
  modal <- modal_value(signif(d, 8))
  which(d > 1.5 * modal)
}

# Check that a time index is sorted, unique and (by default) gap-free.
check_time_index <- function(index, allow_gaps = FALSE) {
  if (is.null(index)) {
    return(invisible(NULL))
  }
  if (anyDuplicated(index)) {
    dup <- utils::head(unique(index[duplicated(index)]), .max_locations)
    fgls_abort(
      sprintf(
        "The time index has duplicated values (%s). Each period must appear once.",
        paste(format(dup), collapse = ", ")
      ),
      "time_index"
    )
  }
  v <- as.numeric(index)
  if (is.unsorted(v)) {
    first <- which(diff(v) < 0)[1L] + 1L
    fgls_abort(
      paste0(
        "The data are not sorted in time (first out-of-order value near ",
        format(index[first]), "). FGLS uses lags, so row order matters: sort the ",
        "data by time first."
      ),
      "time_index"
    )
  }
  gaps <- time_index_gaps(index)
  if (length(gaps)) {
    shown <- utils::head(gaps, .max_locations)
    desc <- paste(sprintf("%s -> %s", format(index[shown]), format(index[shown + 1L])),
      collapse = ", "
    )
    more <- if (length(gaps) > .max_locations) {
      sprintf(" and %d more", length(gaps) - .max_locations)
    } else {
      ""
    }
    if (!allow_gaps) {
      fgls_abort(
        sprintf(
          paste0(
            "The time index has %d gap(s) (%s%s). Lags would silently skip the ",
            "missing periods. Fill or interpolate them, or use allow_gaps = TRUE to ",
            "treat lags as positional (e.g. trading days with holidays)."
          ),
          length(gaps), desc, more
        ),
        "time_index"
      )
    }
    fgls_warn(
      sprintf(
        "The time index has %d gap(s) (%s%s); lags are positional (allow_gaps = TRUE).",
        length(gaps), desc, more
      ),
      "specification"
    )
  }
  invisible(NULL)
}

# ------------------------------------------------------------------- missing values
na_policy <- function(na.action) {
  if (is.null(na.action)) {
    return("fail")
  }
  if (is.character(na.action)) {
    nm <- na.action
  } else if (is.function(na.action)) {
    nm <- if (identical(na.action, stats::na.fail)) {
      "na.fail"
    } else if (identical(na.action, stats::na.omit)) {
      "na.omit"
    } else if (identical(na.action, stats::na.exclude)) {
      "na.exclude"
    } else {
      "other"
    }
  } else {
    nm <- "other"
  }
  switch(nm,
    na.fail = "fail",
    na.omit = "omit",
    na.exclude = "exclude",
    fgls_abort("`na.action` must be na.fail (default), na.omit or na.exclude.")
  )
}

# Rows to keep given missing values. For time series only leading and trailing
# rows can be dropped: dropping interior rows would make the lags wrong.
missing_rows_keep <- function(nan_rows, policy, time_series, index) {
  if (!any(nan_rows)) {
    return(rep(TRUE, length(nan_rows)))
  }
  bad <- which(nan_rows)
  if (policy == "fail") {
    fgls_abort(
      sprintf(
        paste0(
          "Missing values (NA) found at %s. Remove or impute them, or use ",
          "na.action = na.omit."
        ),
        format_locations(index, bad)
      ),
      "missing"
    )
  }
  if (!time_series) {
    return(!nan_rows)
  }
  ok <- which(!nan_rows)
  if (length(ok) == 0L) {
    fgls_abort("All observations have missing values.", "missing")
  }
  first <- ok[1L]
  last <- ok[length(ok)]
  interior <- which(nan_rows[first:last]) + first - 1L
  if (length(interior)) {
    fgls_abort(
      sprintf(
        paste0(
          "Missing values inside the sample at %s. With time series only leading ",
          "or trailing missing values can be dropped, because dropping interior ",
          "rows would make the lags wrong. Impute or interpolate them first."
        ),
        format_locations(index, interior)
      ),
      "missing"
    )
  }
  keep <- rep(FALSE, length(nan_rows))
  keep[first:last] <- TRUE
  keep
}

# ------------------------------------------------------------------------ main entry
# Validate y (numeric vector) and X (numeric matrix with column names).
prepare_data <- function(y, X, y_name = "y", index = NULL, time_series = TRUE,
                         allow_gaps = FALSE, warn_no_const = TRUE) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  x_names <- colnames(X)
  if (is.null(x_names)) x_names <- paste0("x", seq_len(ncol(X)))
  colnames(X) <- x_names
  if (nrow(X) != length(y)) {
    fgls_abort(sprintf("y has %d observations but X has %d rows.", length(y), nrow(X)))
  }
  if (ncol(X) == 0L) {
    fgls_abort("The model has no regressors.")
  }
  inf_rows <- is.infinite(y) | rowSums(is.infinite(X)) > 0
  if (any(inf_rows)) {
    fgls_abort(
      sprintf(
        "Infinite values found at %s. Replace them first.",
        format_locations(index, which(inf_rows))
      ),
      "missing"
    )
  }
  if (time_series) check_time_index(index, allow_gaps = allow_gaps)

  n <- nrow(X)
  q <- ncol(X)
  if (n <= q) {
    fgls_abort(
      sprintf(
        "Only %d observations for %d regressors; need more observations than regressors.",
        n, q
      ),
      "insufficient_data"
    )
  }
  is_const <- vapply(seq_len(q), function(j) all(X[, j] == X[1L, j]) && X[1L, j] != 0, logical(1))
  zero <- vapply(seq_len(q), function(j) !any(X[, j] != 0), logical(1))
  if (any(zero)) {
    fgls_abort(
      sprintf(
        "Column(s) %s of the model matrix are identically zero; drop them.",
        paste(x_names[zero], collapse = ", ")
      ),
      "collinearity"
    )
  }
  if (sum(is_const) > 1L) {
    fgls_abort(
      sprintf(
        "The model matrix has several constant columns (%s); keep one.",
        paste(x_names[is_const], collapse = ", ")
      ),
      "collinearity"
    )
  }
  dep <- independent_columns(X)$drop
  if (length(dep)) {
    fgls_abort(
      sprintf(
        paste0(
          "The model matrix is rank deficient: column(s) %s are linear combinations ",
          "of the preceding columns. Drop them (or one level of each factor when an ",
          "intercept is included)."
        ),
        paste(x_names[dep], collapse = ", ")
      ),
      "collinearity"
    )
  }
  const_col <- if (any(is_const)) which(is_const)[1L] else NA_integer_
  implicit <- FALSE
  if (is.na(const_col)) {
    ones <- rep(1, n)
    res <- ones - X %*% lsq(X, ones)
    implicit <- sqrt(sum(res^2)) <= 1e-8 * sqrt(n)
    if (!implicit && warn_no_const) {
      fgls_warn(
        paste0(
          "The model has no intercept. Estimates are valid only if the model truly ",
          "has no intercept; drop `- 1` / `+ 0` from the formula otherwise."
        ),
        "specification"
      )
    }
  }
  list(
    y = y, X = X, y_name = y_name, x_names = x_names, index = index,
    const_col = const_col, has_implicit_const = implicit,
    has_const = !is.na(const_col) || implicit, nobs = n, k_vars = q
  )
}

# Regressors that look like lags of the dependent variable: x_{j,t} equals
# y_{t-p} up to rounding, or their correlation exceeds 0.99999.
detect_lagged_dependent <- function(y, X, names, max_lag = 4L, exclude = character(0)) {
  out <- list()
  n <- length(y)
  for (j in seq_along(names)) {
    if (names[j] %in% exclude) next
    x <- X[, j]
    if (all(x == x[1L])) next
    for (p in seq_len(min(max_lag, n - 3L))) {
      a <- x[(p + 1L):n]
      b <- y[seq_len(n - p)]
      tol_abs <- 1e-12 * (max(abs(b)) + 1)
      if (all(abs(a - b) <= tol_abs + 1e-10 * abs(b))) {
        out[[length(out) + 1L]] <- list(name = names[j], lag = p)
        break
      }
      if (stats::sd(a) > 0 && stats::sd(b) > 0 && stats::cor(a, b) > 0.99999) {
        out[[length(out) + 1L]] <- list(name = names[j], lag = p)
        break
      }
    }
  }
  out
}

warn_lagged_dependent <- function(found) {
  if (length(found) == 0L) {
    return(invisible(NULL))
  }
  desc <- paste(vapply(found, function(f) {
    sprintf("'%s' (looks like y lagged %d)", f$name, f$lag)
  }, character(1)), collapse = ", ")
  nms <- vapply(found, function(f) f$name, character(1))
  fgls_warn(
    sprintf(
      paste0(
        "Regressor(s) %s appear to be lags of the dependent variable. Declare them ",
        "with lagged_dep_vars = c(%s) so the AR filter is estimated consistently ",
        "(PG26 Section 4.2)."
      ),
      desc, paste0("\"", nms, "\"", collapse = ", ")
    ),
    "lagged_dependent"
  )
}

# Resolve a specification of columns (names or 1-based positions) to positions.
resolve_columns <- function(spec, names, what) {
  if (is.null(spec) || length(spec) == 0L) {
    return(integer(0))
  }
  out <- integer(0)
  for (it in as.list(spec)) {
    if (is.character(it)) {
      pos <- match(it, names)
      if (is.na(pos)) {
        fgls_abort(sprintf(
          "%s: '%s' is not a regressor name (%s).", what, it,
          paste(names, collapse = ", ")
        ))
      }
      out <- c(out, pos)
    } else if (is.numeric(it) && it >= 1 && it <= length(names) && it == round(it)) {
      out <- c(out, as.integer(it))
    } else {
      fgls_abort(sprintf("%s: invalid column %s.", what, format(it)))
    }
  }
  unique(out)
}
