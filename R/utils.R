#' Generate simulation times and the indices
#' of time points to return in sumulation functions.
#' @param t.start Numeric scalar indicating the initial time of the simulation, default is 0. 
#' @param t.end Numeric scalar indicating the final time of the dimulation, default is 1000.
#' @param t.step Numeric scalar indicating the interval between simulation steps.
#' @param t.keep Integer scalar indicating the number of evenly distributed time points to keep
#' @return list containing simulation times (t.sys) and the indices to keep.\
#' @keywords internal
#' @export
tDyn <- function(t.start, t.end, t.step, t.keep){
  
  t.total <- t.end-t.start
  
  
  t.sys <- seq(t.start, t.end, by = t.step)
  
  t.index <- seq(1, length(t.sys), by=round(length(t.sys)/t.keep)-1)
  
  
  
  return(list("t.sys" = t.sys, "t.index" = t.index))
}