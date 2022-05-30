#' Generate simulation times and the indices of time points to return
#' in simulation functions.
#'
#' @param t_start Numeric scalar indicating the initial time of the simulation.
#' (default: \code{t_start = 0})
#' @param t_end Numeric scalar indicating the final time of the dimulation
#' (default: \code{t_end = 1000})
#' @param t_step Numeric scalar indicating the interval between simulation steps
#' (default: \code{t_step = 0.1})
#' @param t_store Integer scalar indicating the number of evenly distributed
#' time points to keep (default: \code{t_store = 100})
#'
#' @return lists containing simulation times (t_sys) and the indices to keep.
#' @examples
#' Time <- simulationTimes(t_start = 0, t_end = 100, t_step = 0.5,
#'     t_store = 100)
#' DefaultTime <- simulationTimes(t_end = 1000)
#'
#' @importFrom stats cov
#' @importFrom stats rbinom
#'
#' @keywords internal
#' @export
simulationTimes <- function(t_start = 0, t_end = 1000,
            t_step = 0.1, t_store = 1000){
        t_total <- t_end-t_start
        t_sys <- seq(t_start, t_end, by = t_step)
    t_index <- seq(1, length(t_sys)-1, by=floor(length(t_sys)/t_store))
    return(list("t_sys" = t_sys, "t_index" = t_index[seq_len(t_store)]))
}

#' check whether a number is a positive integer
#' @param x Numeric number to test
#' @param tol Numeric tolerance of detection
#' @return A logical value: whether the number is a positive integer.
.isPosInt <- function(x, tol = .Machine$double.eps^0.5) {
        return(abs(x - round(x)) < tol && x > 0)
}

#' generate a vector of times when events is happening
#' @param t_events,t_duration Numeric: vector of starting time and durations of
#' the events
#' @param t_end Numeric: end time of the simulation
#' @param ... : additional parameters to pass to simulationTimes, including
#' t_start, t_step, and t_store.
#' @return A vector of time points in the simulation
#' @examples
#' tEvent <-  eventTimes(t_events = c(10, 50, 100),
#'     t_duration = c(1, 2, 3),
#'     t_end = 100,
#'     t_store = 100,
#'     t_step = 1)
#' @export
eventTimes <- function(t_events = NULL, t_duration = NULL,
                       t_end=1000, ...){
        tdyn <- simulationTimes(t_end = t_end,...)
        t_result <- c()
        for (i in seq(length(t_events))){
            p1 <- tdyn$t_sys[(tdyn$t_sys >= t_events[i]) &
            (tdyn$t_sys < (t_events[i]+t_duration[i]))]
            t_result <- c(t_result, p1)
        }
        return(t_result)
}

#' generate pairs of interactions according to interaction types
#' A helper function to be used in combination with .getInteractions()
#' @param I Matrix: defining the interaction between each pair of species
#' @param pair Numeric: a vector with a length of 2, indicating the 2 focusing
#' species in the process of applying the interaction types
#' @param interType Character: one of 'mutualism', 'commensalism', 'parasitism',
#' 'amensalism', or 'competition'. Defining the interaction type
#' @return A matrix of interaction types with one pair changed
.applyInterType <- function(I, pair, interType){
  if (rbinom(1,1,0.5)){
    pair <- rev(pair)
  }
  if (interType=='mutualism'){
    I[pair[1],pair[2]] <- 1
    I[pair[2],pair[1]] <- 1
    return(I)
  }else if (interType=='commensalism'){
    I[pair[1],pair[2]] <- 1
    I[pair[2],pair[1]] <- 0
    return(I)
  }else if (interType=='parasitism'){
    I[pair[1],pair[2]] <- 1
    I[pair[2],pair[1]] <- -1
    return(I)
  }else if (interType=='amensalism'){
    I[pair[1],pair[2]] <- 0
    I[pair[2],pair[1]] <- -1
    return(I)
  }else if (interType=='competition'){
    I[pair[1],pair[2]] <- -1
    I[pair[2],pair[1]] <- -1
    return(I)
  }
}

#' generate interactions according to five types of interactions and
#' their weights
#' @param n_species Integer: defining the dimension of matrix of interaction
#' @param weights Numeric: defining the weights of mutualism, commensalism,
#' parasitism, amensalism, and competition in all interspecies interactions.
#' @param connectance Numeric: defining the density of the interaction network.
#' Ranging from 0 to 1
#' @return A matrix of interactions with all interactions changed according to
#' the weights and connectance.
.getInteractions <- function(n_species, weights, connectance){
    I <- matrix(0, n_species, n_species)
    interactions <- c('mutualism', 'commensalism', 'parasitism',
                      'amensalism', 'competition')
    probs <- abs(weights)/sum(abs(weights))
    combinations <- utils::combn(n_species, 2)
    for (i in sample(seq_along(combinations[1,]),as.integer(ncol(combinations)*connectance), replace=FALSE)){
        I <- .applyInterType(I, combinations[,i], sample(interactions, 1, prob = probs))
    }
    return (I)
}
