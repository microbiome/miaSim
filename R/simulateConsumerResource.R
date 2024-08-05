#' Consumer-resource model simulation
#'
#' Simulates time series with the consumer-resource model.
#'
#' @template man_spe
#' @template man_res
#' @param E \code{Matrix}. Defines the efficiency of resource-to-biomass 
#' conversion (positive values) and the relative conversion
#' of metabolic by-products (negative values). If \code{NULL},
#' `randomE(n_species, n_resources)` is used.
#' (Default: \code{NULL})
#' @param x0 \code{Numeric scalar}. Specifies the initial abundances of 
#' simulated species. If \code{NULL}, `runif(n = n_species, min = 0.1, max = 10)` 
#' is used. (Default: \code{NULL})
#' @param resources \code{Numeric scalar}. Specifies the initial 
#' concentrations of resources. If \code{NULL}, `runif(n = n_resources, 
#' min = 1, max = 100)` is used. (Default: \code{NULL})
#' @param resources_dilution \code{Numeric scalar}. Specifies the 
#' concentrations of resources in the continuous inflow 
#' (applicable when inflow_rate > 0). If \code{NULL}, `resources` is used.
#' (Default: \code{NULL})
#' @param growth_rates \code{Numeric vector}. Specifies the maximum 
#' growth rates(mu) of species. If \code{NULL}, `rep(1, n_species)` is used.
#' (Default: \code{NULL})
#' @param monod_constant \code{Matrix}. Specifies the constant of 
#' additive monod growth of n_species consuming n_resources. If \code{NULL},
#' `matrix(rgamma(n = n_species*n_resources, shape = 50*max(resources), 
#' rate = 1), nrow = n_species)` is used. (Default: \code{NULL})
#' @template man_sto
#' @template man_mig
#' @template man_mod
#' @param trophic_priority \code{Matrix}. Defines the orders of resources to
#' be consumed by each species. If \code{NULL}, by default, this feature won't be
#' turned on, and species will consume all resources simultaneously to grow.
#' The dimension should be identical to matrix E.
#' (Default: \code{NULL})
#' @param inflow_rate, outflow_rate \code{Numeric scalar}. The inflow 
#' of a culture process. By default, inflow_rate and is 0, indicating a batch 
#' culture process. When larger than 0, we can simulate a continuous 
#' culture(e.g. chemostat).
#' @param outflow_rate \code{Numeric scalar}. outflow rate of a culture process
#' By default, outflow_rate is 0, indicating a batch culture process. When larger 
#' than 0, we can simulate a continuous culture(e.g. chemostat).
#' @param volume \code{Numeric scalar}. Indicates the volume of the continuous 
#' cultivation. This parameter is important for simulations where inflow_rate 
#' or outflow_rate are not 0. (Default: \code{1000})
#'
#' @examples
#' n_species <- 2
#' n_resources <- 4
#' tse <- simulateConsumerResource(
#'     n_species = n_species,
#'     n_resources = n_resources
#' )
#'
#' \dontrun{
#' # example with user-defined values (names_species, names_resources, E, x0,
#' # resources, growth_rates, error_variance, t_end, t_step)
#'
#' ExampleE <- randomE(
#'     n_species = n_species, n_resources = n_resources,
#'     mean_consumption = 3, mean_production = 1, maintenance = 0.4
#' )
#' ExampleResources <- rep(100, n_resources)
#' tse1 <- simulateConsumerResource(
#'     n_species = n_species,
#'     n_resources = n_resources, names_species = letters[seq_len(n_species)],
#'     names_resources = paste0("res", LETTERS[seq_len(n_resources)]), E = ExampleE,
#'     x0 = rep(0.001, n_species), resources = ExampleResources,
#'     growth_rates = runif(n_species),
#'     error_variance = 0.01,
#'     t_end = 5000,
#'     t_step = 1
#' )
#'
#' # example with trophic levels
#' n_species <- 10
#' n_resources <- 15
#' ExampleEfficiencyMatrix <- randomE(
#'     n_species = 10, n_resources = 15,
#'     trophic_levels = c(6, 3, 1),
#'     trophic_preferences = list(
#'         c(rep(1, 5), rep(-1, 5), rep(0, 5)),
#'         c(rep(0, 5), rep(1, 5), rep(-1, 5)),
#'         c(rep(0, 10), rep(1, 5))
#'     )
#' )
#'
#' ExampleResources <- c(rep(500, 5), rep(200, 5), rep(50, 5))
#' tse2 <- simulateConsumerResource(
#'     n_species = n_species,
#'     n_resources = n_resources,
#'     names_species = letters[1:n_species],
#'     names_resources = paste0(
#'         "res", LETTERS[1:n_resources]
#'     ),
#'     E = ExampleEfficiencyMatrix,
#'     x0 = rep(0.001, n_species),
#'     resources = ExampleResources,
#'     growth_rates = rep(1, n_species),
#'     # error_variance = 0.001,
#'     t_end = 5000, t_step = 1
#' )
#'
#' # example with trophic priority
#' n_species <- 4
#' n_resources <- 6
#' ExampleE <- randomE(
#'     n_species = n_species,
#'     n_resources = n_resources,
#'     mean_consumption = n_resources,
#'     mean_production = 0
#' )
#' ExampleTrophicPriority <- t(apply(
#'     matrix(sample(n_species * n_resources),
#'         nrow = n_species
#'     ),
#'     1, order
#' ))
#' # make sure that for non-consumables resources for each species,
#' # the priority is 0 (smaller than any given priority)
#' ExampleTrophicPriority <- (ExampleE > 0) * ExampleTrophicPriority
#' tse3 <- simulateConsumerResource(
#'     n_species = n_species,
#'     n_resources = n_resources,
#'     E = ExampleE,
#'     trophic_priority = ExampleTrophicPriority,
#'     t_end = 2000
#' )
#' }
#'
#' @return an TreeSummarizedExperiment class object
#'
#' @export
simulateConsumerResource <- function(n_species, n_resources,
    names_species = NULL,
    names_resources = NULL,
    E = NULL,
    x0 = NULL,
    resources = NULL,
    resources_dilution = NULL,
    growth_rates = NULL,
    monod_constant = NULL,
    sigma_drift = 0.001,
    sigma_epoch = 0.1,
    sigma_external = 0.3,
    sigma_migration = 0.01,
    epoch_p = 0.001,
    t_external_events = NULL,
    t_external_durations = NULL,
    stochastic = FALSE,
    migration_p = 0.01,
    metacommunity_probability = NULL,
    error_variance = 0,
    norm = FALSE,
    t_end = 1000,
    trophic_priority = NULL,
    inflow_rate = 0,
    outflow_rate = 0,
    volume = 1000,
    ...) {
    t_dyn <- .simulationTimes(t_end = t_end, ...)

    # calculate the time points influenced by the disturbances
    tEvent <- simulateEventTimes(
        t_events = t_external_events,
        t_duration = t_external_durations,
        t_end = t_end,
        ... = ...
    )

    if (!is.null(trophic_priority)) {
        if (!identical(dim(E), dim(trophic_priority))) {
            stop("The dimension of trophic priority is not correct.")
        }
    }

    # set the default values
    if (is.null(names_species)) {
        names_species <- paste0("sp", seq_len(n_species))
    }
    if (is.null(names_resources)) {
        names_resources <- paste0("res", seq_len(n_resources))
    }
    if (is.null(E)) {
        E <- randomE(n_species, n_resources)
    }
    if (is.null(x0)) {
        x0 <- runif(n = n_species, min = 0.1, max = 10)
    }
    if (is.null(resources)) {
        resources <- runif(n = n_resources, min = 1, max = 100)
    }
    if (is.null(resources_dilution)) {
        resources_dilution <- resources
    }
    if (is.null(growth_rates)) {
        growth_rates <- rep(1, n_species)
    }
    if (is.null(monod_constant)) {
        monod_constant <- matrix(rgamma(
            n = n_species * n_resources,
            shape = 50 * max(resources), rate = 1
        ), nrow = n_species)
    }
    if (is.null(metacommunity_probability)) {
        metacommunity_probability <- rdirichlet(1, alpha = rep(1, n_species))
    }

    # normalize metacommunity_probability
    metacommunity_probability <- metacommunity_probability /
        sum(metacommunity_probability)

    # define the perturbation event
    perturb <- function(t, y, parameters) {
        with(as.list(y), {

            # continuous or episodic perturbation
            epoch_rN <- 0
            external_rN <- 0
            migration_rN <- 0
            if (rbinom(1, 1, parameters$epoch_p)) {
                epoch_rN <- rnorm(parameters$n_species, sd = parameters$sigma_epoch)
                epoch_rN <- parameters$stochastic * epoch_rN
            }

            if (rbinom(1, 1, parameters$migration_p)) {
                migration_rN <- rmultinom(
                    n = 1, size = 1,
                    prob = parameters$metacommunity_probability
                )[, ] * abs(rnorm(
                    n = 1,
                    mean = 0, sd = parameters$sigma_migration
                ))
            }

            if (t %in% parameters$tEvent) {
                external_rN <- rnorm(parameters$n_species,
                    sd = parameters$sigma_external
                )
                external_rN <- parameters$stochastic * external_rN
            }

            drift_rN <- rnorm(parameters$n_species, sd = parameters$sigma_drift)
            drift_rN <- parameters$stochastic * drift_rN

            # perturbation is applied to the current population
            consumer <- pmax(y[grepl("consumer", names(y))], 0)

            consumer <- consumer * (1 + drift_rN) * (1 + epoch_rN) * (1 + external_rN) + migration_rN
            resource <- y[grepl("resource", names(y))]

            volume <- y[grepl("volume", names(y))]
            return(c(consumer, resource, volume))
        })
    }


    state_init <- c(x0, resources, volume)
    names(state_init) <- c(
        paste0("consumer", seq(length(x0))),
        paste0("resource", seq(length(resources))),
        "volume"
    )

    parameters <- list(
        growth_rates = growth_rates, E = E, monod_constant = monod_constant,
        n_species = n_species, sigma_drift = sigma_drift, stochastic = stochastic,
        sigma_epoch = sigma_epoch, epoch_p = epoch_p,
        sigma_external = sigma_external, tEvent = tEvent,
        migration_p = migration_p, metacommunity_probability = metacommunity_probability,
        sigma_migration = sigma_migration,
        resources_dilution = resources_dilution,
        trophic_priority = trophic_priority,
        inflow_rate = inflow_rate,
        outflow_rate = outflow_rate,
        volume = volume, threshold = 0.01
    )

    out <- as.data.frame(ode(
        y = state_init, times = t_dyn$t_sys,
        func = consumerResourceModel, parms = parameters,
        events = list(func = perturb, time = t_dyn$t_sys)
    ))

    # end_row and index_to_return helps to solve returning early of ode solver.
    end_row <- max(t_dyn$t_index)
    if (nrow(out) <= end_row) {
        warning("returned output are earlier finished than t_end.")
        index_to_return <- t_dyn$t_index[t_dyn$t_index <= nrow(out)]
    } else {
        index_to_return <- t_dyn$t_index
    }
    species_index <- grepl("consumer", names(out))
    resource_index <- grepl("resource", names(out))
    volume_index <- grepl("volume", names(out))

    out_species_matrix <- as.matrix(out[, species_index])
    out_species_matrix <- as.matrix(out_species_matrix[index_to_return, ])

    out_resource_matrix <- as.matrix(out[, resource_index])
    out_resource_matrix <- as.matrix(out_resource_matrix[index_to_return, ])


    out_volume_matrix <- as.matrix(out[, volume_index])
    out_volume_matrix <- as.matrix(out_volume_matrix[index_to_return, ])

    if (error_variance > 0) {
        measurement_error <- rgamma(length(out_species_matrix), 1 / error_variance, 1 / error_variance)
        out_species_matrix <- out_species_matrix * matrix(measurement_error, ncol = n_species)
    }

    if (norm) {
        out_species_matrix <- out_species_matrix / rowSums(out_species_matrix)
        out_resource_matrix <- out_resource_matrix / rowSums(out_resource_matrix)
    }

    colnames(out_species_matrix) <- names_species
    colnames(out_resource_matrix) <- names_resources
    colnames(out_volume_matrix) <- c("volume")

    out_species_matrix <- cbind(
        out_species_matrix,
        time = t_dyn$t_sys[index_to_return]
    )
    out_resource_matrix <- cbind(
        out_resource_matrix,
        time = t_dyn$t_sys[index_to_return]
    )
    out_volume_matrix <- cbind(
        out_volume_matrix,
        time = t_dyn$t_sys[index_to_return]
    )

    TreeSE <- TreeSummarizedExperiment(
        assays = list(counts = t(out_species_matrix[, 1:n_species])),
        colData = DataFrame(time = out_species_matrix[, "time"]),
        metadata = list(resources = out_resource_matrix,
                        volume = out_volume_matrix,
                        E = E,
                        monod_constant = monod_constant,
                        error_variance = error_variance))

    return(TreeSE)
}

# define the consumer-resource model
consumerResourceModel <- function(t, state, params) {
    with(as.list(c(state, params)), {
        volume <- state[startsWith(names(state), "volume")]
        volume <- volume * (volume > 10^(-10))


        x0 <- (volume > 0) * pmax(0, state[startsWith(names(state), "consumer")])

        resources <- (volume > 0) * pmax(0, state[startsWith(names(state), "resource")])



        growth_rates <- params[["growth_rates"]]
        E <- params[["E"]]
        monod_constant <- params[["monod_constant"]]

        resources_dilution <- params[["resources_dilution"]]
        inflow_rate <- params[["inflow_rate"]]
        outflow_rate <- params[["outflow_rate"]]

        dilution_rate <- inflow_rate / volume


        trophic_priority <- params[["trophic_priority"]]

        threshold <- params[["threshold"]]

        if (!is.null(trophic_priority)) {

            # modify E in each step
            Emod <- trophic_priority
            ## resources <= a relatively small number(instead of 0) ######
            Emod[, resources <= threshold] <- 0
            Emod <- t(apply(Emod, 1, function(x) {
                ifelse(x == max(x), x, 0)
            })) > 0
            ## 1. E*(Emod): consumption of prior resource
            ## 2. E*(E<0): production
            ## 3. E*apply(!Emod, 1, all): if all resources are run out, then
            ## select all of them (instead of none of them), to avoid the
            ## production without consumption
            Emod <- E * (Emod) + E * (E < 0) + E * apply(!Emod, 1, all)
            E <- Emod
        }



        B <- matrix(rep(resources, length(x0)),
            ncol = length(resources), byrow = TRUE
        ) + monod_constant
        growth <- ((E * (E > 0) / B) %*% resources) * x0


        consumption <- -resources * (t(1 / B * (E > 0)) %*% x0)

        production <- -t(E * (E < 0)) %*% growth

        dVolume <- max(-volume, inflow_rate - outflow_rate)
        dResources <- consumption + production - (dilution_rate * (resources - resources_dilution))
        dConsumers <- growth_rates * growth - dilution_rate * x0

        dxdt <- list(c(dConsumers, dResources, dVolume))
        return(dxdt)
    })
}
