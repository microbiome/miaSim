#' Stochastic Logistic simulation
#'
#' Simulates a community time series using the logistic model.
#' The change rate of the species was defined as
#' `dx/dt = b*x*(1-(x/k))*rN - dr*x`, where
#' b is the vector of growth rates,
#' x is the vector of initial species abundances,
#' k is the vector of maximum carrying capacities,
#' rN is a random number ranged from 0 to 1 which changes in each time step,
#' dr is the vector of constant death rates.
#' Also, the vectors of initial dead species abundances can be set.
#' The number of species will be set to 0 if the dead species abundances
#' surpass the alive species abundances.
#'
#' @param n_species integer number of species
#' @param b numeric growth rates
#' (default: \code{b = runif(n = n_species, min = 0.1, max = 0.2)})
#' @param k numeric value of carrying capacities
#' (default: \code{k = runif(n = n_species, min = 1000, max = 2000)})
#' @param dr numeric value of death rates
#' (default: \code{dr = runif(n = n_species, min = 0.0005, max = 0.0025)})
#' @param x numeric initial abundances
#' (default: \code{x = runif(n = n_species, min = 0.1, max = 10)})
#' @param sigma_drift numeric degree of drift (turnover of species) in each
#' time step.
#' (default: \code{sigma_drift = 0.001})
#' @param sigma_epoch numeric degree of epoch change of community
#' (default: \code{sigma_epoch = 0.1})
#' @param sigma_external numeric degree of external events/disturbances
#' (default: \code{sigma_external = 0.3})
#' @param p_epoch numeric value of probability/frequency of inherit periodic
#' changes of community (default: \code{p_epoch = 0.001})
#' @param t_external_events numeric value of starting times of external events
#' (default: \code{t_external_events = c(0, 240, 480)})
#' @param t_external_durations numeric value of durations of external events
#' (default: \code{t_external_durations = c(0, 1, 1)})
#' @param stochastic logical scalar choosing whether the logistic model should
#' be stochastic (controlled by multiplying the growth rate by a random number)
#' (default: \code{stochastic = TRUE})
#' @param t_end numeric final time of the simulation
#' (default: \code{t_end = 1000})
#' @param ... additional parameters including 't_start', 't_step', and 't_store'
#'
#' @seealso
#' \code{\link[miaSim:convertToSE]{convertToSE}}
#'
#' @examples
#' ## ATTENTION: Don't set a large value to t.step, otherwise the computer won't
#' #give a correct solution to the logistic ODE(ordinary differential equation).
#' #Keeping t_step under 0.05 or 0.01 is a good practice.
#'
#' #while (!exists("ExampleLogistic"))
#' ExampleLogistic <- simulateStochasticLogistic(n_species = 5)
#' #plot the calculated points
#' matplot(ExampleLogistic, type = "l")
#'
#' #calculation by setting initial parameters explicitly
#' ExampleLogistic2 <- simulateStochasticLogistic(n_species = 2,
#' b = c(0.2, 0.1), k = c(1000, 2000), dr = c(0.001, 0.0015), x = c(3, 0.1),
#' sigma_drift = 0.001, sigma_epoch = 0.3, sigma_external = 0.5,p_epoch = 0.001,
#' t_external_events = c(100, 200, 300), t_external_durations = c(1, 2, 3),
#' t_start = 0, t_end = 1500, t_step = 0.01,
#' t_store = 1500, stochastic = TRUE)
#'
#' @return \code{simulateStochasticLogistic} returns an abundance matrix with
#' species abundance as rows and time points as columns
#'
#' @export
simulateStochasticLogistic <- function(n_species,
        b = runif(n = n_species, min = 0.1, max = 0.2),
        k = runif(n = n_species, min = 1000, max = 2000),
        dr = runif(n = n_species, min = 0.0005, max = 0.0025),
        x = runif(n = n_species, min = 0.1, max = 10),
        sigma_drift = 0.001,
        sigma_epoch = 0.1,
        sigma_external = 0.3,
        p_epoch = 0.001,
        t_external_events = c(0, 240, 480),
        t_external_durations = c(0, 1, 1),
        stochastic = TRUE,
        t_end = 1000, ...){

        # define the stochastic logistic model
        stochasticLogisticModel <- function (t, state, params){
            current <- pmax(0,state[names(state) == 'current'])
            live <- state[names(state) == 'live']
            dead <- state[names(state) == 'dead']
            b <- params$b
            k <- params$k
            dr <- params$dr

            dlive <- b*live*(1-(live/(k)))
            ddead <- dr*current
            dcurrent <- (dlive-ddead)
            dxdt <- list(c(dcurrent, dlive, ddead))
            return(dxdt)
        }

        # check the input format
        if(!isPositiveInteger(n_species)){
            stop("n_species must be integer.")}
        if(!all(vapply(list(b,k,dr,x), is.double, logical(1)),
            vapply(list(b,k,dr,x), length, integer(1)) == n_species)){
            stop("b,k,dr,x must be double and of n_species length.")}
        if(!is.logical(stochastic)){
            stop("stochastic should be boolean values.")}

        # select the time points to simulate and to store
        t_dyn <- simulationTimes(t_end = t_end,...)

        #continuous or episodic perturbation
        perturb <- function(t, y, parms){with(as.list(y), {
            epoch_rN <- 0
            external_rN <- 0
            if (rbinom(1,1, parms$p_epoch)){
                epoch_rN <- rnorm(parms$n_species, sd=parms$sigma_epoch)
                epoch_rN <- params$stochastic*epoch_rN
            }
            if (t %in% parms$tEvent){
                external_rN <- rnorm(parms$n_species, sd=parms$sigma_external)
                external_rN <- params$stochastic*external_rN
            }
            drift_rN <- rnorm(parms$n_species, sd=parms$sigma_drift)
            drift_rN <- params$stochastic*drift_rN

            #perturbation is applied to the current population
            current <- pmax(y[names(y)=='current'], 0)
            current <- current*(1+drift_rN)*(1+epoch_rN)*(1+external_rN)
            live <- y[names(y)=='live']
            dead <- y[names(y)=='dead']
            return(c(current, live, dead))})
        }

        tEvent = eventTimes(t_events = t_external_events,
            t_duration = t_external_durations,
            t_end = t_end, ...)

        params <- list(b=b, k=k, dr=dr, n_species = n_species,
            sigma_drift = sigma_drift, stochastic = stochastic,
            sigma_epoch = sigma_epoch, p_epoch = p_epoch,
            sigma_external = sigma_external, tEvent = tEvent)
        yinit <- c(x, x, numeric(n_species))
        names(yinit) <- rep(c("current", "live", "dead"), each = n_species)

        out <- as.data.frame(ode(func=stochasticLogisticModel,
                    y=yinit, times=t_dyn$t_sys, parms=params,
                    events = list(func = perturb, time = t_dyn$t_sys)))

        out_matrix <- out[names(out)=='current']
        names(out_matrix) <- seq_len(n_species)
        out_matrix <- out_matrix[t_dyn$t_index,]

        out_matrix$t <- t_dyn$t_sys[t_dyn$t_index]
        counts <- as.matrix(out_matrix)
        return(counts)

}
