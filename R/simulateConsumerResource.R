#' Consumer-resource model simulation
#'
#' Simulates time series with the consumer-resource model
#'
#' @param n_species Integer: number of species
#' @param n_resources Interger: number of resources
#' @param names_species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n_species))` is used.
#' (default: \code{names_species = NULL})
#' @param names_resources Character: names of resources. If NULL,
#' `paste0("res", seq_len(n_resources))` is used.
#' @param E matrix: matrix of efficiency. A matrix defining the efficiency of
#' resource-to-biomass conversion (positive values) and the relative conversion
#' of metabolic by-products (negative values). If NULL, 
#' `randomE(n_species, n_resources)` is used.
#' (default: \code{E = NULL})
#' @param x0 Numeric: initial abundances of simulated species. If NULL,
#' `runif(n = n_species, min = 0.1, max = 10)` is used.
#' (default: \code{x0 = NULL})
#' @param resources Numeric: initial concentrations of resources. If NULL,
#' `runif(n = n_resources, min = 1, max = 100)` is used.
#' (default: \code{resources = NULL})
#' @param resources_dilution Numeric: concentrations of resources in the 
#' continuous inflow (applicable when inflow_rate > 0). If NULL,
#' `resources` is used. 
#' (default: \code{resources_dilution = NULL}) 
#' @param growth_rates Numeric: vector of maximum growth rates(mu) of species.
#' If NULL, `rep(1, n_species)` is used.
#' (default: \code{growth_rates = NULL})
#' @param monod_constant matrix: the constant of additive monod growth of
#' n_species consuming n_resources. If NULL,
#' `matrix(rgamma(n = n_species*n_resources, shape = 50*max(resources), rate = 1), nrow = n_species)`
#' is used.
#' (default: \code{monod_constant = NULL})
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
#' (default: \code{t_external_events = NULL})
#' @param t_external_durations Numeric: respective duration of the external 
#' events that are defined in the 't_external_events' (times) and 
#' sigma_external (std).
#' (default: \code{t_external_durations = NULL})
#' @param stochastic Logical: whether to introduce noise in the simulation.
#' If False, sigma_drift, sigma_epoch, and sigma_external are ignored.
#' (default: \code{stochastic = FALSE})
#' @param migration_p Numeric: the probability/frequency of migration from a 
#' metacommunity.
#' (default: \code{migration_p = 0.01})
#' @param metacommunity_probability Numeric: Normalized probability distribution
#' of the likelihood that species from the metacommunity can enter the community
#' during the simulation. If NULL, `rdirichlet(1, alpha = rep(1,n_species))` is 
#' used.
#' (default: \code{metacommunity_probability = NULL})
#' @param error_variance Numeric: the variance of measurement error.
#' By default it equals to 0, indicating that the result won't contain any 
#' measurement error. This value should be non-negative.
#' (default: \code{error_variance = 0})
#' @param norm Logical scalar: returns normalized abundances (proportions
#' in each generation) 
#' (default: \code{norm = FALSE})
#' @param t_end Numeric scalar indicating the final time of the simulation
#' (default: \code{t_end = 1000})
#' @param trophic_priority Matrix: a matrix defining the orders of resources to 
#' be consumed by each species. If NULL, by default, this feature won't be 
#' turned on, and species will consume all resources simultaneously to grow. 
#' The dimension should be identical to matrix E.
#' (default: \code{trophic_priority = NULL})
#' @param inflow_rate,outflow_rate Numeric: the inflow and outflow rate of a 
#' culture process. By default, inflow_rate and outflow_rate are 0, indicating a
#' batch culture process. By setting them equally larger than 0, we can simulate
#' a continuous culture(e.g. chemostat).
#' @param volume Numeric: the volume of the continuous cultivation. This 
#' parameter is important for simulations where inflow_rate or outflow_rate are
#' not 0. 
#' (default: \code{volume = 1000})
#' @param ... additional parameters, see \code{\link{utils}} to know more.
#' @seealso
#' \code{\link[miaSim:convertToSE]{convertToSE}}
#'
#' @examples
#' ExampleCR <- simulateConsumerResource(n_species = 2, 
#'     n_resources = 4)
#' makePlot(ExampleCR$matrix)
#' makePlot(ExampleCR$resources)
#' 
#' # example to get relative abundance and relative proportion of resources
#' ExampleCR <- simulateConsumerResource(n_species = 2, 
#'     n_resources = 4, norm = TRUE)
#' makePlot(ExampleCR$matrix)
#' makePlot(ExampleCR$resources)
#'
#' # example with user-defined values (names_species, names_resources, E, x0, 
#' # resources, growth_rates, error_variance, t_end, t_step)
#' ExampleE <- randomE(n_species = 4, n_resources = 6, 
#'     mean_consumption = 3, mean_production = 1, maintenance = 0.4)
#' ExampleResources <- rep(100, 6)
#' ExampleCR <- simulateConsumerResource(n_species = 4, 
#'     n_resources = 6, names_species = letters[1:4], 
#'     names_resources = paste0("res",LETTERS[1:6]), E = ExampleE, 
#'     x0 = rep(0.001, 4), resources = ExampleResources, 
#'     growth_rates = runif(4),
#'     error_variance = 0.01,
#'     t_end = 5000,
#'     t_step = 1)
#' makePlot(ExampleCR$matrix)
#' makePlot(ExampleCR$resources)
#' 
#' # example with trophic levels
#' n_species <- 10
#' n_resources <- 15
#' 
#' ExampleEfficiencyMatrix <- randomE(n_species = 10, n_resources = 15,
#'                                    trophic_levels = c(6,3,1),
#'                                    trophic_preferences = list(
#'                                        c(rep(1,5), rep(-1, 5), rep(0, 5)), 
#'                                        c(rep(0,5), rep(1, 5), rep(-1, 5)), 
#'                                        c(rep(0,10), rep(1, 5))))
#' makeHeatmap(ExampleEfficiencyMatrix)
#' 
#' # ExampleResources <- rep(100, n_resources)
#' ExampleResources <- c(rep(500, 5), rep(200, 5), rep(50, 5))
#' ExampleCR <- simulateConsumerResource(n_species = n_species, 
#'                                       n_resources = n_resources, 
#'                                       names_species = letters[1:n_species], 
#'                                       names_resources = paste0(
#'                                           "res",LETTERS[1:n_resources]), 
#'                                       E = ExampleEfficiencyMatrix, 
#'                                       x0 = rep(0.001, n_species), 
#'                                       resources = ExampleResources, 
#'                                       growth_rates = rep(1, n_species),
#'                                       # error_variance = 0.001,
#'                                       t_end = 5000, t_step = 1)
#' makePlot(ExampleCR$matrix)
#' makePlotRes(ExampleCR$resources)
#' 
#' 
#' # example with trophic priority
#' n_species <- 4
#' n_resources <- 6
#' ExampleE <- randomE(n_species = n_species, 
#'                     n_resources = n_resources, 
#'                     mean_consumption = n_resources, 
#'                     mean_production = 0)
#' ExampleTrophicPriority <- t(apply(matrix(sample(n_species * n_resources), 
#'                                          nrow = n_species), 
#'                                   1, order))
#' # make sure that for non-consumables resources for each species, 
#' # the priority is 0 (smaller than any given priority) 
#' ExampleTrophicPriority <- (ExampleE > 0) * ExampleTrophicPriority
#' ExampleCR <- simulateConsumerResource(n_species = n_species, 
#'                                       n_resources = n_resources, 
#'                                       E = ExampleE, 
#'                                       trophic_priority = ExampleTrophicPriority,
#'                                       t_end = 2000
#'                                       )
#' makePlot(ExampleCR$matrix)
#' makePlotRes(ExampleCR$resources)
#'
#' @return an abundance matrix with species and resources abundance as rows and
#' time points as columns
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
    ...){
    
    t_dyn <- simulationTimes(t_end = t_end,...)
    
    # calculate the time points influenced by the disturbances
    tEvent <- eventTimes(
        t_events = t_external_events,
        t_duration = t_external_durations,
        t_end = t_end,
        ... = ...)
    
    if(!is.null(trophic_priority)){
        if (!identical(dim(E), dim(trophic_priority))){
            stop("The dimension of trophic priority is not correct.")}}

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
        monod_constant <- matrix(rgamma(n = n_species*n_resources,
                                        shape = 50*max(resources), rate = 1), nrow = n_species)
    }
    if (is.null(metacommunity_probability)) {
        metacommunity_probability <- rdirichlet(1, alpha = rep(1,n_species))
    }
    
    #normalize metacommunity_probability
    metacommunity_probability <- metacommunity_probability/
        sum(metacommunity_probability)
    
    # define the perturbation event
    perturb <- function(t, y, parameters){
        with(as.list(y),{
            
            #continuous or episodic perturbation
            epoch_rN <- 0
            external_rN <- 0
            migration_rN <- 0
            if (rbinom(1,1, parameters$epoch_p)){
                epoch_rN <- rnorm(parameters$n_species, sd=parameters$sigma_epoch)
                epoch_rN <- parameters$stochastic*epoch_rN
            }
            
            if (rbinom(1,1, parameters$migration_p)){
                migration_rN <- rmultinom(n = 1, size = 1, 
                                          prob = parameters$metacommunity_probability)[,]*abs(rnorm(n=1, 
                                                                                                    mean=0, sd = parameters$sigma_migration))
                
            }
            
            if (t %in% parameters$tEvent){
                external_rN <- rnorm(parameters$n_species,
                                     sd=parameters$sigma_external)
                external_rN <- parameters$stochastic*external_rN
            }
            
            drift_rN <- rnorm(parameters$n_species, sd=parameters$sigma_drift)
            drift_rN <- parameters$stochastic*drift_rN
            
            #perturbation is applied to the current population
            consumer <- pmax(y[grepl('consumer', names(y))], 0)
            
            consumer <- consumer*(1+drift_rN)*(1+epoch_rN)*(1+external_rN)+ migration_rN
            resource <- y[grepl('resource', names(y))]
            
            volume <- y[grepl('volume', names(y))]
            return(c(consumer, resource, volume))})
    }
    
    
    state_init <- c(x0, resources, volume)
    names(state_init) <- c(paste0("consumer", seq(length(x0))),
                           paste0("resource", seq(length(resources))),
                           "volume")
    
    parameters <- list(growth_rates = growth_rates, E = E, monod_constant = monod_constant,  
                       n_species = n_species, sigma_drift = sigma_drift, stochastic = stochastic,
                       sigma_epoch = sigma_epoch, epoch_p = epoch_p,
                       sigma_external = sigma_external, tEvent = tEvent,
                       migration_p = migration_p, metacommunity_probability = metacommunity_probability,
                       sigma_migration = sigma_migration, 
                       resources_dilution = resources_dilution,
                       trophic_priority = trophic_priority, 
                       inflow_rate = inflow_rate,
                       outflow_rate = outflow_rate,
                       volume = volume, threshold=0.01
                       )
    
    out <- as.data.frame(ode(y = state_init, times = t_dyn$t_sys,
                             func = consumerResourceModel, parms = parameters, 
                             events = list(func = perturb, time = t_dyn$t_sys)))
    
    
    species_index <- grepl('consumer', names(out))
    resource_index <- grepl('resource', names(out))
    volume_index <- grepl('volume', names(out))
    
    out_species_matrix <- as.matrix(out[,species_index])
    out_species_matrix <- as.matrix(out_species_matrix[t_dyn$t_index,])
    
    out_resource_matrix <- as.matrix(out[,resource_index])
    out_resource_matrix <- as.matrix(out_resource_matrix[t_dyn$t_index,])
    
    
    out_volume_matrix <- as.matrix(out[,volume_index])
    out_volume_matrix <- as.matrix(out_volume_matrix[t_dyn$t_index,])
    
    if(error_variance > 0){
        measurement_error <- rgamma(length(out_species_matrix), 1/error_variance, 1/error_variance)
        out_species_matrix <- out_species_matrix * matrix(measurement_error, ncol = n_species)
    }
    
    if(norm){
        out_species_matrix <- out_species_matrix/rowSums(out_species_matrix)
        out_resource_matrix <- out_resource_matrix/rowSums(out_resource_matrix)
    }
    
    colnames(out_species_matrix) <- names_species
    colnames(out_resource_matrix) <- names_resources
    colnames(out_volume_matrix) <- c("volume")

    out_species_matrix <- cbind(
        out_species_matrix, 
        time = t_dyn$t_sys[t_dyn$t_index])
    out_resource_matrix <- cbind(
        out_resource_matrix, 
        time = t_dyn$t_sys[t_dyn$t_index])
    out_volume_matrix <- cbind(
        out_volume_matrix,
        time = t_dyn$t_sys[t_dyn$t_index])

    out_list <- list(matrix = out_species_matrix, 
        resources = out_resource_matrix,
        volume = out_volume_matrix,
        E = E,
        monod_constant = monod_constant,
        error_variance = error_variance)
    return(out_list)
}

# define the consumer-resource model
consumerResourceModel <- function(t, state, params){
  with(as.list(c(state, params)),{
    volume <- state[startsWith(names(state), "volume")]
    volume <- volume * (volume>10^(-10))
    
    
    x0 <- (volume>0) * pmax(0, state[startsWith(names(state), "consumer")])
    
    resources <- (volume>0) * pmax(0, state[startsWith(names(state), "resource")])
    
    
    
    growth_rates <- params[['growth_rates']]
    E <- params[['E']]
    monod_constant <- params[['monod_constant']]
    
    resources_dilution <- params[['resources_dilution']]
    inflow_rate <- params[['inflow_rate']]
    outflow_rate <- params[['outflow_rate']]
    
    dilution_rate <- inflow_rate/volume
    
    
    trophic_priority <- params[['trophic_priority']]
    
    threshold <- params[['threshold']]
    
    if(!is.null(trophic_priority)){
      
      # modify E in each step
      Emod <- trophic_priority
      ## resources <= a relatively small number(instead of 0) ######
      Emod[, resources <= threshold] <- 0 
      Emod <- t(apply(Emod, 1, getRowMax))>0
      ## 1. E*(Emod): consumption of prior resource 
      ## 2. E*(E<0): production 
      ## 3. E*apply(!Emod, 1, all): if all resources are run out, then
      ## select all of them (instead of none of them), to avoid the 
      ## production without consumption 
      Emod <- E*(Emod) + E*(E<0) + E*apply(!Emod, 1, all)
      E <- Emod
    }
    
    
    
    B <- matrix(rep(resources, length(x0)),
                ncol = length(resources), byrow = TRUE) + monod_constant
    growth <- ((E*(E>0)/B) %*% resources)*x0
    
    
    consumption <- -resources*(t(1/B*(E>0))%*%x0)
    
    production <- -t(E*(E<0)) %*% growth
    
    dVolume <- max(-volume, inflow_rate - outflow_rate)
    dResources <- consumption + production - (dilution_rate*(resources - resources_dilution))
    dConsumers <- growth_rates*growth - dilution_rate*x0
    
    dxdt <- list(c(dConsumers, dResources, dVolume))
    return(dxdt)
  })
}

