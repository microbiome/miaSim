#' Generate simulation times and the indices of time points to return
#' in simulation functions.
#'
#' @param t.start Numeric scalar indicating the initial time of the simulation.
#' (default: \code{t.start = 0})
#' @param t.end Numeric scalar indicating the final time of the dimulation
#' (default: \code{t.end = 1000})
#' @param t.step Numeric scalar indicating the interval between simulation steps
#' (default: \code{t.step = 0.1})
#' @param t.store Integer scalar indicating the number of evenly distributed
#' time points to keep (default: \code{t.store = 100})
#'
#' @return lists containing simulation times (t.sys) and the indices to keep.
#' @examples
#' Time <- tDyn(t.start = 0, t.end = 100, t.step = 0.5, t.store = 100)
#' DefaultTime <- tDyn(t.start = 0)
#'
#' @docType methods
#' @aliases tDyn-numeric
#' @aliases tDyn,numeric-method
#'
#' @keywords internal
#' @export
setGeneric("tDyn", signature = c("t.start"),
        function(t.start = 0, t.end = 1000, t.step = 0.1, t.store = 1000)
            standardGeneric("tDyn"))

setMethod("tDyn", signature = c(t.start="numeric"),
        function(t.start = 0, t.end = 1000, t.step = 0.1, t.store = 1000){

        t.total <- t.end-t.start

        t.sys <- seq(t.start, t.end, by = t.step)

        t.index <- seq(1, length(t.sys)-1, by=round(length(t.sys)/t.store))

        return(list("t.sys" = t.sys, "t.index" = t.index))
})
