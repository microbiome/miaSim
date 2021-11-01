
simulationTimes <- function(t.start = 0, t.end = 1000,
            t.step = 0.1, t.store = 1000){
        t.total <- t.end-t.start
        t.sys <- seq(t.start, t.end, by = t.step)
        t.index <- seq(1, length(t.sys)-1, by=round(length(t.sys)/t.store))
        return(list("t.sys" = t.sys, "t.index" = t.index))
}


isPositiveInteger <- function(x, tol = .Machine$double.eps^0.5) {
        return(abs(x - round(x)) < tol && x > 0)
}

isZeroOrPositiveInteger <- function(x, tol = .Machine$double.eps^0.5) {
    return(abs(x - round(x)) < tol && x >= 0)
}


eventTimes <- function(t.events=c(10,20,30),
                        t.duration=rep(3,3), t.end=1000, ...){
        tdyn <- simulationTimes(t.end = t.end,...)
        t.result = c()
        for (i in seq(length(t.events))){
            p1 <- tdyn$t.sys[(tdyn$t.sys >= t.events[i]) &
                        (tdyn$t.sys <= (t.events[i]+t.duration[i]))]
            t.result <- c(t.result, p1)
        }
        return(t.result)
    }

