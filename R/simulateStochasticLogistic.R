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
#' @template man_spe
#' @param growth_rates \code{Numeric scalar}. Specifies the growth rates 
#' of simulated species. If NULL, `runif(n = n_species, min = 0.1, max = 0.2)` 
#' is used. (Default: \code{NULL})
#' @param carrying_capacities \code{Numeric scalar}. Indicates the max 
#' population of species supported in the community. If NULL,
#' `runif(n = n_species, min = 1000, max = 2000)` is used.
#' (Default: \code{NULL})
#' @param death_rates \code{Numeric scalar}. Indicates the death rates 
#' of each species. If NULL, `runif(n = n_species, min = 0.0005, max = 0.0025)` 
#' is used. (Default: \code{NULL})
#' @param x0 \code{Numeric scalar}. Indicates the initial abundances of simulated 
#' species. If NULL,  `runif(n = n_species, min = 0.1, max = 10)` is used.
#' (Default: \code{NULL})
#' @template man_sto
#' @template man_mig
#' @param stochastic \code{Logical scalar}. Whether to introduce noise in the simulation.
#' If False, sigma_drift, sigma_epoch, and sigma_external are ignored.
#' (Default: \code{TRUE})
#' @template man_mod
#'
#' @examples
#' # Example of logistic model without stochasticity, death rates, or external
#' # disturbances
#' set.seed(42)
#' tse <- simulateStochasticLogistic(
#'     n_species = 5,
#'     stochastic = FALSE, death_rates = rep(0, 5)
#' )
#'
#' # Adding a death rate
#' set.seed(42)
#' tse1 <- simulateStochasticLogistic(
#'     n_species = 5,
#'     stochastic = FALSE, death_rates = rep(0.01, 5)
#' )
#'
#' # Example of stochastic logistic model with measurement error
#' set.seed(42)
#' tse2 <- simulateStochasticLogistic(
#'     n_species = 5,
#'     error_variance = 1000
#' )
#'
#' # example with all the initial parameters defined by the user
#' set.seed(42)
#' tse3 <- simulateStochasticLogistic(
#'     n_species = 2,
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
#'     metacommunity_probability = miaSim::rdirichlet(1, alpha = rep(1, 2)),
#'     stochastic = TRUE,
#'     error_variance = 0,
#'     norm = FALSE, # TRUE,
#'     t_end = 400,
#'     t_start = 0, t_step = 0.01,
#'     t_store = 1500
#' )
#'
#' @return \code{simulateStochasticLogistic} returns a TreeSummarizedExperiment
#' class object
#'
#' @importFrom stats rnorm
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
    t_end = 1000, ...) {

    # define the stochastic logistic model
    stochasticLogisticModel <- function(t, state, parameters) {
        current <- pmax(0, state[names(state) == "current"])
        live <- state[names(state) == "live"]
        dead <- state[names(state) == "dead"]
        growth_rates <- parameters$growth_rates
        carrying_capacities <- parameters$carrying_capacities
        death_rates <- parameters$death_rates

        dlive <- growth_rates * live * (1 - (live / (carrying_capacities)))
        ddead <- death_rates * current
        dcurrent <- (dlive - ddead)
        dxdt <- list(c(dcurrent, dlive, ddead))
        return(dxdt)
    }


    # input check
    if (!.isPosInt(n_species)) {
        stop("n_species must be integer.")
    }

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
    if (is.null(metacommunity_probability)) {
        metacommunity_probability <- rdirichlet(1, alpha = rep(1, n_species))
    }

    # select the time points to simulate and to store
    t_dyn <- .simulationTimes(t_end = t_end, ...)

    # continuous or episodic perturbation
    perturb <- function(t, y, parms) {
        with(as.list(y), {
            epoch_rN <- 0
            external_rN <- 0
            migration_rN <- 0
            if (rbinom(1, 1, parms$epoch_p)) {
                epoch_rN <- rnorm(parms$n_species, sd = parms$sigma_epoch)
                epoch_rN <- parameters$stochastic * epoch_rN
            }
            if (rbinom(1, 1, parameters$migration_p)) {
                migration_rN <- rmultinom(n = 1, size = 1, prob = parameters$metacommunity_probability)[, ] * abs(rnorm(n = 1, mean = 0, sd = parameters$sigma_migration))
            }
            if (t %in% parms$tEvent) {
                external_rN <- rnorm(parms$n_species, sd = parms$sigma_external)
                external_rN <- parameters$stochastic * external_rN
            }
            drift_rN <- rnorm(parms$n_species, sd = parms$sigma_drift)
            drift_rN <- parameters$stochastic * drift_rN

            # perturbation is applied to the current population
            current <- pmax(y[names(y) == "current"], 0)
            current <- current * (1 + drift_rN) * (1 + epoch_rN) * (1 + external_rN) + migration_rN
            live <- y[names(y) == "live"]
            dead <- y[names(y) == "dead"]
            return(c(current, live, dead))
        })
    }

    tEvent <- simulateEventTimes(
        t_events = t_external_events,
        t_duration = t_external_durations,
        t_end = t_end, ...
    )

    parameters <- list(
        growth_rates = growth_rates, carrying_capacities = carrying_capacities, death_rates = death_rates, n_species = n_species,
        sigma_drift = sigma_drift, stochastic = stochastic,
        sigma_epoch = sigma_epoch, epoch_p = epoch_p,
        sigma_external = sigma_external, tEvent = tEvent,
        migration_p = migration_p, metacommunity_probability = metacommunity_probability,
        sigma_migration = sigma_migration
    )
    yinit <- c(x0, x0, numeric(n_species))
    names(yinit) <- rep(c("current", "live", "dead"), each = n_species)

    out <- as.data.frame(ode(
        func = stochasticLogisticModel,
        y = yinit, times = t_dyn$t_sys, parms = parameters,
        events = list(func = perturb, time = t_dyn$t_sys)
    ))

    out_matrix <- as.matrix(out[, names(out) == "current"])
    out_matrix <- out_matrix[t_dyn$t_index, ]

    if (error_variance > 0) {
        measurement_error <- rnorm(
            n = length(t_dyn$t_index) * n_species,
            mean = 0, sd = sqrt(error_variance)
        )
        measurement_error <- matrix(measurement_error,
            nrow = length(t_dyn$t_index)
        )
        out_matrix <- out_matrix + measurement_error
    }

    if (norm) {
        out_matrix <- out_matrix / rowSums(out_matrix)
    }

    colnames(out_matrix) <- names_species

    out_matrix <- cbind(out_matrix, time = t_dyn$t_sys[t_dyn$t_index])

    TreeSE <- TreeSummarizedExperiment(
        assays = list(counts = t(out_matrix[, 1:n_species])),
        colData = DataFrame(time = out_matrix[, "time"]),
        metadata = list(metacommunity_probability = metacommunity_probability,
                        migration_p = migration_p,
                        error_variance = error_variance))

    return(TreeSE)
}
