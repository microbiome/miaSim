#' @param error_variance Numeric: the variance of measurement error.
#' By default it equals to 0, indicating that the result won't contain any
#' measurement error. This value should be non-negative.
#' (default: \code{error_variance = 0})
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' (default: \code{norm = FALSE})
#' @param t_end Numeric: the end time of the simulationTimes, defining the
#' modeled time length of the community.
#' (default: \code{t_end = 1000})
#' @param ... additional parameters, see \code{\link{utils}} to know more.
