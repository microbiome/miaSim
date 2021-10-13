#' Stochastic Logistic simulation
#'
#' Simulates time series with the (stochastic) logistic model and
#' forms a \linkS4class{SummarizedExperiment} object or a matrix.
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
#' The resulting abundance matrix model is used to construct
#' \linkS4class{SummarizedExperiment} object.
#'
#' @param n.species Integer: number of species
#' @param b Numeric: growth rates
#' (default: \code{b = runif(n = n.species, min = 0.1, max = 0.2)})
#' @param k Numeric: carrying capacities
#' (default: \code{k = runif(n = n.species, min = 1000, max = 2000)})
#' @param dr Numeric: death rates
#' (default: \code{dr = runif(n = n.species, min = 0.0005, max = 0.0025)})
#' @param x Numeric: initial abundances
#' (default: \code{x = runif(n = n.species, min = 0.1, max = 10)})
#' @param sigma.drift Numeric: degree of drift (turnover of species) in each
#' time step.
#' (default: \code{sigma.drift = 0.001})
#' @param sigma.epoch Numeric: degree of epoch change of community
#' (default: \code{sigma.epoch = 0.1})
#' @param sigma.external Numeric: degree of external events/disturbances
#' (default: \code{sigma.external = 0.3})
#' @param p.epoch Numeric: probability/frequency of inherit periodic changes of
#' community
#' (default: \code{p.epoch = 0.001})
#' @param t.external_events Numeric: starting times of external events
#' (default: \code{t.external_events = c(0, 240, 480)})
#' @param t.external_durations Numeric: durations of external events
#' (default: \code{t.external_durations = c(0, 1, 1)})
#' @param t.end Numeric scalar indicating the final time of the dimulation
#' (default: \code{t.end = 2000})
#' @param return.matrix Logical: whether to export only the stored time points
#' in a matrix
#' (default: \code{return.matrix = FALSE})
#' @param stochastic Logical: whether the logistic model should be stochastic
#' (controlled by multiplying the growth rate by a random number)
#' (default: \code{stochastic = TRUE})
#' @param ... additional parameters including 't.start', 't.step', and 't.store'
#' see \code{\link{SimulationTimes}} for more information
#'
#' @docType methods
#' @examples
#' ## ATTENTION: Don't set a large value to t.step, otherwise the computer won't
#' #give a correct solution to the logistic ODE(ordinary differential equation).
#' #Keeping t.step under 0.05 or 0.01 is a good practice.
#'
#' #while (!exists("ExampleLogistic"))
#' ExampleLogistic <- simulateStochasticLogistic(n.species = 5)
#' #plot the calculated points
#' matplot(t(assays(ExampleLogistic)[["counts"]]), type = "l")
#'
#' #calculation by setting initial parameters explicitly
#' ExampleLogistic <- simulateStochasticLogistic(
#' n.species = 2,
#' b = c(0.2, 0.1), k = c(1000, 2000), dr = c(0.001, 0.0015), x = c(3, 0.1),
#' sigma.drift = 0.001, sigma.epoch = 0.3, sigma.external = 0.5,p.epoch = 0.001,
#' t.external_events = c(100, 200, 300), t.external_durations = c(1, 2, 3),
#' t.start = 0, t.end = 1500, t.step = 0.01,
#' t.store = 1500, return.matrix = FALSE, stochastic = TRUE)
#'
#' @return \code{simulateStochasticLogistic} returns a
#' \linkS4class{SummarizedExperiment} class object containing matrix with
#' species abundance as rows and time points as columns
#'
#' @docType methods
#' @aliases simulateStochasticLogistic-numeric
#' @aliases simulateStochasticLogistic,numeric-method
#'
#' @importFrom deSolve ode
#' @importFrom stats runif
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods setGeneric
#'
#' @export

setGeneric("simulateStochasticLogistic",signature = "n.species",
    function(n.species,
        b = runif(n = n.species, min = 0.1, max = 0.2),
        k = runif(n = n.species, min = 1000, max = 2000),
        dr = runif(n = n.species, min = 0.0005, max = 0.0025),
        x = runif(n = n.species, min = 0.1, max = 10),
        sigma.drift = 0.001,
        sigma.epoch = 0.1,
        sigma.external = 0.3,
        p.epoch = 0.001,
        t.external_events = c(0, 240, 480),
        t.external_durations = c(0, 1, 1),
        return.matrix = FALSE,
        stochastic = TRUE,
        t.end = 1000,...)
        standardGeneric("simulateStochasticLogistic"))

setMethod("simulateStochasticLogistic", signature = c(n.species="numeric"),
    function(n.species,
        b = runif(n = n.species, min = 0.1, max = 0.2),
        k = runif(n = n.species, min = 1000, max = 2000),
        dr = runif(n = n.species, min = 0.0005, max = 0.0025),
        x = runif(n = n.species, min = 0.1, max = 10),
        sigma.drift = 0.001,
        sigma.epoch = 0.1,
        sigma.external = 0.3,
        p.epoch = 0.001,
        t.external_events = c(0, 240, 480),
        t.external_durations = c(0, 1, 1),
        return.matrix = FALSE,
        stochastic = TRUE,
        t.end = 1000,...){

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
        if(!isPositiveInteger(n.species)){
            stop("n.species must be integer.")}
        if(!all(vapply(list(b,k,dr,x), is.double, logical(1)),
            vapply(list(b,k,dr,x), length, integer(1)) == n.species)){
            stop("b,k,dr,x must be double and of n.species length.")}
        if(!all(vapply(list(return.matrix,stochastic),is.logical, logical(1)))){
            stop("return.matrix or stochastic should be boolean values.")}

        # select the time points to simulate and to store
        t.dyn <- SimulationTimes(t.end = t.end,...)

        #continuous or episodic perturbation
        perturb <- function(t, y, parms){with(as.list(y), {
            epoch.rN <- 0
            external.rN <- 0
            if (rbinom(1,1, parms$p.epoch)){
                epoch.rN <- rnorm(parms$n.species, sd=parms$sigma.epoch)
                epoch.rN <- params$stochastic*epoch.rN
            }
            if (t %in% parms$tEvent){
                external.rN <- rnorm(parms$n.species, sd=parms$sigma.external)
                external.rN <- params$stochastic*external.rN
            }
            drift.rN <- rnorm(parms$n.species, sd=parms$sigma.drift)
            drift.rN <- params$stochastic*drift.rN

            #perturbation is applied to the current population
            current <- pmax(y[names(y)=='current'], 0)
            current <- current*(1+drift.rN)*(1+epoch.rN)*(1+external.rN)
            live <- y[names(y)=='live']
            dead <- y[names(y)=='dead']
            return(c(current, live, dead))})
        }

        tEvent = eventTimes(t.events = t.external_events,
            t.duration = t.external_durations,
            t.end = t.end, ...)

        params <- list(b=b, k=k, dr=dr, n.species = n.species,
            sigma.drift = sigma.drift, stochastic = stochastic,
            sigma.epoch = sigma.epoch, p.epoch = p.epoch,
            sigma.external = sigma.external, tEvent = tEvent)
        yinit <- c(x, x, numeric(n.species))
        names(yinit) <- rep(c("current", "live", "dead"), each = n.species)

        out <- as.data.frame(ode(func=stochasticLogisticModel,
                    y=yinit, times=t.dyn$t.sys, parms=params,
                    events = list(func = perturb, time = t.dyn$t.sys)))

        out.matrix <- out[names(out)=='current']
        names(out.matrix) <- seq_len(n.species)
        out.matrix <- out.matrix[t.dyn$t.index,]
        if (return.matrix) {
            out.matrix$t <- t.dyn$t.sys[t.dyn$t.index]
            return(as.matrix(out.matrix))
        } else {
            SE <- SummarizedExperiment(assays = list(counts = t(out.matrix)))
            return(SE)
        }
    }
)
