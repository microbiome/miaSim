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
#' Time <- SimulationTimes(t.start = 0, t.end = 100, t.step = 0.5,
#'     t.store = 100)
#' DefaultTime <- SimulationTimes(t.end = 1000)
#'
#' @docType methods
#' @aliases SimulationTimes-numeric
#' @aliases SimulationTimes,numeric-method
#'
#' @keywords internal
#' @export

setGeneric("SimulationTimes", signature = "t.end",
    function(t.start = 0, t.end = 1000, t.step = 0.1, t.store = 1000)
    standardGeneric("SimulationTimes"))

setMethod("SimulationTimes", signature = c(t.end="numeric"),
    function(t.start = 0, t.end = 1000, t.step = 0.1, t.store = 1000){
        t.total <- t.end-t.start
        t.sys <- seq(t.start, t.end, by = t.step)
        t.index <- seq(1, length(t.sys)-1, by=round(length(t.sys)/t.store))
        return(list("t.sys" = t.sys, "t.index" = t.index))
})

setGeneric("isPositiveInteger", signature = "x",
    function(x, tol = .Machine$double.eps^0.5)
        standardGeneric("isPositiveInteger"))
setMethod("isPositiveInteger", signature = c(x="numeric"),
    function(x, tol = .Machine$double.eps^0.5) {
        return(abs(x - round(x)) < tol && x > 0)
    }
)

setGeneric("eventTimes", signature = "t.events",
    function(t.events=c(10,20,30), t.duration=rep(3,3), t.end=1000, ...)
        standardGeneric("eventTimes"))
setMethod("eventTimes", signature = c(t.events="numeric"),
    function(t.events=c(10,20,30), t.duration=rep(3,3), t.end=1000, ...){
        tdyn <- SimulationTimes(t.end = t.end,...)
        t.result = c()
        for (i in seq(length(t.events))){
            p1 <- tdyn$t.sys[(tdyn$t.sys >= t.events[i]) &
                        (tdyn$t.sys <= (t.events[i]+t.duration[i]))]
            t.result <- c(t.result, p1)
        }
        return(t.result)
    })
