#' Stochastic Logistic simulation
#'
#' Simulates time series with the (stochastic) logistic model
#' 
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
#' @param n_species Integer: number of species
#' @param names_species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n_species))` is used.
#' (default: \code{names_species = NULL})
#' @param growth_rates Numeric: growth rates of simulated species. If NULL,
#'  `runif(n = n_species, min = 0.1, max = 0.2)` is used.
#' (default: \code{growth_rates = NULL})
#' @param carrying_capacities Numeric: The max population of species supported 
#' in the community. If NULL,
#' `runif(n = n_species, min = 1000, max = 2000)` is used.
#' (default: \code{carrying_capacities = NULL})
#' @param death_rates Numeric: death rates of each species. If NULL, 
#' `runif(n = n_species, min = 0.0005, max = 0.0025)` is used.
#' (default: \code{death_rates = NULL})
#' @param x0 Numeric: initial abundances of simulated species. If NULL, 
#' `runif(n = n_species, min = 0.1, max = 10)` is used.
#' (default: \code{x0 = NULL})
#' @param sigma_drift Numeric: standard deviation of a normally distributed 
#' noise applied in each time step (t_step)
#' (default: \code{sigma_drift = 0.001})
#' @param sigma_epoch Numeric: standard deviation of a normally distributed 
#' noise applied to random periods of the community composition with frequency 
#' defined by the epoch_p parameter
#' (default: \code{sigma_epoch = 0.1})
#' @param sigma_external Numeric: standard deviation of a normally distributed 
#' noise applied to user-defined external events/disturbances
#' (default: \code{sigma_external = 0.3})
#' @param sigma_migration Numeric: standard deviation of a normally distributed 
#' variable that defines the intensity of migration at each time step (t_step)
#' (default: \code{sigma_migration = 0.01})
#' @param epoch_p Numeric: the probability/frequency of random periodic 
#' changes introduced to the community composition
#' (default: \code{epoch_p = 0.001})
#' @param t_external_events Numeric: the starting time points of defined 
#' external events that introduce random changes to the community composition
#' (default: \code{t_external_events = c(0, 240, 480)})
#' @param t_external_durations Numeric: respective duration of the external 
#' events that are defined in the 't_external_events' (times) and sigma_external (std).
#' (default: \code{t_external_durations = c(0, 1, 1)})
#' @param migration_p Numeric: the probability/frequency of migration from a 
#' metacommunity.
#' (default: \code{migration_p = 0.01})
#' @param metacommunity_probability Numeric: Normalized probability distribution
#' of the likelihood that species from the metacommunity can enter the community
#' during the simulation. If NULL, 
#' `rdirichlet(1, alpha = rep(1,n_species))` is used.
#' (default: \code{metacommunity_probability = NULL})
#' @param stochastic Logical: whether to introduce noise in the simulation.
#' If False, sigma_drift, sigma_epoch, and sigma_external are ignored.
#' (default: \code{stochastic = TRUE})
#' @param error_variance Numeric: the variance of measurement error.
#' By default it equals to 0, indicating that the result won't contain any 
#' measurement error. This value should be non-negative.
#' (default: \code{error_variance = 0})
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t_end Numeric: the end time of the simulationTimes, defining the 
#' modeled time length of the community. 
#' (default: \code{t_end = 1000})
#' @param ... additional parameters, see \code{\link{utils}} to know more.
#' 
#' @seealso
#' \code{\link[miaSim:convertToSE]{convertToSE}}
#'
#' @importFrom gtools rdirichlet
#' 
#' @examples
#'
#' # Example of logistic model without stochasticity, death rates, or external 
#' # disturbances
#' set.seed(42)
#' ExampleLogistic <- simulateStochasticLogistic(n_species = 5, 
#'     stochastic = FALSE, death_rates=rep(0,5))
#' makePlot(ExampleLogistic$matrix)
#' 
#' # Adding a death rate
#' set.seed(42)
#' ExampleLogistic <- simulateStochasticLogistic(n_species = 5, 
#'     stochastic = FALSE, death_rates=rep(0.01, 5))
#' makePlot(ExampleLogistic$matrix)
#' 
#' # Example of stochastic logistic model
#' set.seed(42)
#' ExampleLogistic <- simulateStochasticLogistic(n_species = 5)
#' makePlot(ExampleLogistic$matrix)
#' 
#' # Example of stochastic logistic model with measurement error
#' set.seed(42)
#' ExampleLogistic <- simulateStochasticLogistic(n_species = 5, 
#'     error_variance = 1000)
#' makePlot(ExampleLogistic$matrix)
#'
#' # example with all the initial parameters defined by the user
#' set.seed(42)
#' ExampleLogistic <- simulateStochasticLogistic(n_species = 2,
#'     names_species = c("species1", "species2"),
#'     growth_rates = c(0.2, 0.1), 
#'     carrying_capacities = c(1000, 2000), 
#'     death_rates = c(0.001, 0.0015), 
#'     x0 = c(3, 0.1),
#'     sigma_drift = 0.001, 
#'     sigma_epoch = 0.3, 
#'     sigma_external = 0.5,
#'     sigma_migration = 0.002,
#'     epoch_p = 0.001,
#'     t_external_events = c(100, 200, 300), 
#'     t_external_durations = c(0.1, 0.2, 0.3),
#'     migration_p = 0.01,
#'     metacommunity_probability = gtools::rdirichlet(1, alpha = rep(1, 2)),
#'     stochastic = TRUE,
#'     error_variance = 0,
#'     norm = TRUE, 
#'     t_end = 400, 
#'     t_start = 0, t_step = 0.01,
#'     t_store = 1500)
#' makePlot(ExampleLogistic$matrix)
#'
#' @return \code{simulateStochasticLogistic} returns a list of community dynamic
#' matrix(species abundance as rows and time points as columns) and its 
#' inputs(including metacommunity_probability and migration_p)
#' 
#' @export
simulateStochasticLogistic <- function(n_species,
    names_species = NULL,
    growth_rates = NULL,
    carrying_capacities = NULL,
    death_rates = NULL,
    x0 = NULL,
    sigma_drift = 0.001,
    sigma_epoch = 0.1,
    sigma_external = 0.3,
    sigma_migration = 0.01,
    epoch_p = 0.001,
    t_external_events = NULL, 
    t_external_durations = NULL,  
    migration_p = 0.01,
    metacommunity_probability = NULL,
    stochastic = TRUE,
    error_variance = 0,
    norm = FALSE,
    t_end = 1000, ...){

    # define the stochastic logistic model
    stochasticLogisticModel <- function (t, state, parameters){
        current <- pmax(0,state[names(state) == 'current'])
        live <- state[names(state) == 'live']
        dead <- state[names(state) == 'dead']
        growth_rates <- parameters$growth_rates
        carrying_capacities <- parameters$carrying_capacities
        death_rates <- parameters$death_rates

        dlive <- growth_rates*live*(1-(live/(carrying_capacities)))
        ddead <- death_rates*current
            dcurrent <- (dlive-ddead)
            dxdt <- list(c(dcurrent, dlive, ddead))
            return(dxdt)
        }

    
    #input check
    if(!isPositiveInteger(n_species)){
      stop("n_species must be integer.")}
    
    # set the default values
    if (is.null(names_species)) {
        names_species <- paste0("sp", seq_len(n_species))
    }
    if (is.null(growth_rates)) {
        growth_rates <- runif(n = n_species, min = 0.1, max = 0.2)
    }
    if (is.null(carrying_capacities)) {
        carrying_capacities <- runif(n = n_species, min = 1000, max = 2000)
    }
    if (is.null(death_rates)) {
        death_rates <- runif(n = n_species, min = 0.0005, max = 0.0025)
    }
    if (is.null(x0)) {
        x0 <- runif(n = n_species, min = 0.1, max = 10)
    }
    if (is.null(metacommunity_probability)){
        metacommunity_probability <- rdirichlet(1, alpha = rep(1,n_species))
    }

        # select the time points to simulate and to store
        t_dyn <- simulationTimes(t_end = t_end,...)

    #continuous or episodic perturbation
        perturb <- function(t, y, parms){with(as.list(y), {
            epoch_rN <- 0
            external_rN <- 0
        migration_rN <- 0
        if (rbinom(1,1, parms$epoch_p)){
                epoch_rN <- rnorm(parms$n_species, sd=parms$sigma_epoch)
            epoch_rN <- parameters$stochastic*epoch_rN
        }
        if (rbinom(1,1, parameters$migration_p)){
            migration_rN <- rmultinom(n = 1, size = 1, prob = parameters$metacommunity_probability)[,]*abs(rnorm(n=1, mean=0, sd = parameters$sigma_migration))
            
            
            }
            if (t %in% parms$tEvent){
                external_rN <- rnorm(parms$n_species, sd=parms$sigma_external)
            external_rN <- parameters$stochastic*external_rN
            }
            drift_rN <- rnorm(parms$n_species, sd=parms$sigma_drift)
        drift_rN <- parameters$stochastic*drift_rN

            #perturbation is applied to the current population
            current <- pmax(y[names(y)=='current'], 0)
        current <- current*(1+drift_rN)*(1+epoch_rN)*(1+external_rN)+ migration_rN
            live <- y[names(y)=='live']
            dead <- y[names(y)=='dead']
            return(c(current, live, dead))})
        }

        tEvent <- eventTimes(t_events = t_external_events,
            t_duration = t_external_durations,
            t_end = t_end, ...)

    parameters <- list(growth_rates=growth_rates, carrying_capacities=carrying_capacities, death_rates=death_rates, n_species = n_species,
            sigma_drift = sigma_drift, stochastic = stochastic,
        sigma_epoch = sigma_epoch, epoch_p = epoch_p,
        sigma_external = sigma_external, tEvent = tEvent,
        migration_p = migration_p, metacommunity_probability = metacommunity_probability,
        sigma_migration = sigma_migration)
    yinit <- c(x0, x0, numeric(n_species))
        names(yinit) <- rep(c("current", "live", "dead"), each = n_species)

        out <- as.data.frame(ode(func=stochasticLogisticModel,
        y=yinit, times=t_dyn$t_sys, parms=parameters,
                    events = list(func = perturb, time = t_dyn$t_sys)))

    out_matrix <- as.matrix(out[,names(out)=='current'])
        out_matrix <- out_matrix[t_dyn$t_index,]

    if(error_variance > 0){
        measurement_error <- rnorm(n = length(t_dyn$t_index)*n_species, 
                                   mean = 0, sd = sqrt(error_variance))
        measurement_error <- matrix(measurement_error, 
                                    nrow = length(t_dyn$t_index))
        out_matrix <- out_matrix + measurement_error
    }
    
    if(norm){
        out_matrix <- out_matrix/rowSums(out_matrix)
    }
    
    colnames(out_matrix) <- names_species
    
    out_matrix <- cbind(out_matrix, time = t_dyn$t_sys[t_dyn$t_index])

    out_list <- list(matrix = out_matrix, 
                     metacommunity_probability = metacommunity_probability,
                     migration_p = migration_p,
                     error_variance = error_variance)
    return(out_list)
}
