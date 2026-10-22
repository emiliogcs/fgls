# Temporary placeholders, replaced in later development phases.
hetero_cov_types <- c("robust", "nonrobust", "HC0", "HC1", "HC2", "HC3")
combined_cov_types <- c("robust", "HC0", "HC3")
not_yet <- function(...) fgls_abort("Not implemented yet.")
estimate_hetero <- not_yet
estimate_combined <- not_yet
build_hetero <- not_yet
build_combined <- not_yet
resolve_candidates <- not_yet
screen_stationarity <- function(d) invisible(NULL)
hetero_covariance_from_estimate <- not_yet
combined_covariance_from_estimate <- not_yet
