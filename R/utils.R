#' Generate simulation times and the indices
#' of time points to return in sumulation functions.
#' @param t.start Numeric scalar. Initial time of the simulation, default 0.
#' @param t.end Numeric scalar. Final time of the simulation, default 1000.
#' @param t.step Numeric scalar. Interval between simulation steps, default 0.1.
#' @param t.keep Integer scalar. Number of evenly distributed time points to keep, default 100.
#' @return list containing simulation times (t.sys) and the indices to keep.\
#' @keywords internal
#' @export
tDyn <- function(t.start = 0, t.end = 1000, t.step = 0.1, t.keep = 100){

  t.total <- t.end-t.start

  t.sys <- seq(t.start, t.end, by = t.step)

  t.index <- seq(1, length(t.sys), by=round(length(t.sys)/t.keep))

  return(list("t.sys" = t.sys, "t.index" = t.index))
}

simulationTimes <- function(t_start = 0, t_end = 1000,
            t_step = 0.1, t_store = 1000){
        t_total <- t_end-t_start
        t_sys <- seq(t_start, t_end, by = t_step)
        t_index <- seq(1, length(t_sys)-1, by=round(length(t_sys)/t_store))
        return(list("t_sys" = t_sys, "t_index" = t_index))
}


isPositiveInteger <- function(x, tol = .Machine$double.eps^0.5) {
        return(abs(x - round(x)) < tol && x > 0)
}

isZeroOrPositiveInteger <- function(x, tol = .Machine$double.eps^0.5) {
    return(abs(x - round(x)) < tol && x >= 0)
}


eventTimes <- function(t_events=c(10,20,30),
                        t_duration=rep(3,3), t_end=1000, ...){
        tdyn <- simulationTimes(t_end = t_end,...)
        t_result = c()
        for (i in seq(length(t_events))){
            p1 <- tdyn$t_sys[(tdyn$t_sys >= t_events[i]) &
                        (tdyn$t_sys <= (t_events[i]+t_duration[i]))]
            t_result <- c(t_result, p1)
        }
        return(t_result)
}