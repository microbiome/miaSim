#' Generalized Lotka-Volterra (gLV) simulation
#'
#' Simulates time series with the generalized Lotka-Volterra model.
#'
#' Simulates a community time series using the generalized Lotka-Volterra model,
#' defined as dx/dt = diag(x)(b+Ax), where x is the vector of species abundances
#' ,diag(x) is a diagonal matrix with the diagonal values set to x.
#' A is the interaction matrix and b is the vector of growth rates.
#'
#' @param n_species integer number of species
#' @param A interaction matrix
#' @param x numeric initial abundances
#' @param b numeric growth rates
#' @param sigma_drift numeric degree of drift
#' (turnover of species) in each time step.
#' (default: \code{sigma_drift = 0.01})
#' @param sigma_epoch numeric degree of epoch change of community
#' (default: \code{sigma_epoch = 0.3})
#' @param sigma_external numeric degree of the external events/disturbances
#' (default: \code{sigma_external = 0.3})
#' @param p_epoch numeric value of the probability/frequency of inherit periodic
#' changes of community (default: \code{p_epoch  = 0.01})
#' @param t_external_events numeric value of starting times of external events
#' (default: \code{t_external_events = c(12, 36, 48)})
#' @param t_external_durations numeric durations of external events
#' (default: \code{t_external_durations = c(3, 10, 99)})
#' @param stochastic logical scalar choosing whether the gLV model should be
#' stochastic (default: \code{stochastic = FALSE})
#' @param norm logical scalar returning normalised abundances (proportions
#' in each generation) (default: \code{norm = FALSE})
#' @param t_end numeric value of simulation end time
#' (default: \code{t_end = 1000})
#' @param ... additional parameters including 't_start', 't_step', and 't_store'
#'
#' @seealso
#' \code{\link[miaSim:convertToSE]{convertToSE}}
#'
#' @return
#' \code{simulateGLV} returns an abundance matrix
#'
#' @examples
#' A <- miaSim::powerlawA(4, alpha = 1.01)
#'
#' ExampleGLV <- simulateGLV(n_species = 4, A, t_end = 1000)
#'
#' @importFrom MatrixGenerics colSums2
#' @importFrom deSolve ode
#' @importFrom stats runif
#'
#' @export
simulateGLV <- function(n_species, A,
        x = runif(n_species),
        b = runif(n_species),
        sigma_drift = 0.01,
        sigma_epoch = 0.3,
        sigma_external = 0.3,
        p_epoch  = 0.01,
        t_external_events = c(12, 36, 48),
        t_external_durations = c(3, 10, 99),
        stochastic = FALSE,
        norm = FALSE,
        t_end = 1000, ...){


        # input check
        if(!isPositiveInteger(n_species)){
            stop("n_species must be integer.")}
        if(!all(vapply(list(A,x,b), is.double, logical(1)),
                vapply(list(x,b), length, integer(1)) == n_species)){
            stop("A,x,b must be matrix and the length must be equal to length
                of n_species.")}

        t_dyn <- simulationTimes(t_end = t_end, ...)
        tEvent <- eventTimes(
            t_events = t_external_events,
            t_duration = t_external_durations, t_end = t_end, ...)
        parameters <- list(b=b, A = A, n_species = n_species,
            sigma_drift = sigma_drift, stochastic = stochastic,
            sigma_epoch = sigma_epoch, p_epoch  = p_epoch ,
            sigma_external = sigma_external, tEvent = tEvent)
        out <- ode(
            y = x,
            times = t_dyn$t_sys,
            func = dxdt,
            parms = parameters,
            events = list(func = perturb, time = t_dyn$t_sys))
        counts <- t(out[,2:ncol(out)])
        counts <- counts[,t_dyn$t_index]

        if(norm){
            counts <- t(t(counts)/colSums2(counts))
        }
        return(counts)
    }

dxdt <- function(t, x, parameters){
    b <- parameters$b
    A <- parameters$A
    # rate of change
    dx <- x*(b+A %*% x)
    # return rate of change
    list(dx)
}

perturb <- function(t, y, parameters){
    with(as.list(y),{
        #continuous or episodic perturbation
        epoch_rN <- 0
        external_rN <- 0
        if (rbinom(1,1, parameters$p_epoch )){
            epoch_rN <- rnorm(parameters$n_species, sd=parameters$sigma_epoch)
            epoch_rN <- parameters$stochastic*epoch_rN
        }
        if (t %in% parameters$tEvent){
            external_rN <- rnorm(parameters$n_species,
                                    sd=parameters$sigma_external)
            external_rN <- parameters$stochastic*external_rN
        }
        drift_rN <- rnorm(parameters$n_species, sd=parameters$sigma_drift)
        drift_rN <- parameters$stochastic*drift_rN

        #perturbation is applied to the current population
        y <- y * (1 + drift_rN)*(1 + epoch_rN)*(1 + external_rN)
        return(y)
    })
}