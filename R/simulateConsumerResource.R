#' Consumer-resource model simulation
#'
#' Simulates a community time series using the consumer-resource model.
#' The change rate of each species was defined as
#' `dx/dt = mumax*sum(monod)*X`, where
#' mumax is the vector of maximum growth rates for the species,
#' monod is the monod growth rate, S/(Ks+S), where S is the concentration of the
#' limiting resource, and Ks is the half-velocity constant for species X and S.
#' X is the vector of abundances of species.
#' The concentrations of resource will be set to 0 if they were calculated
#' less than 0 during the iteration.
#'
#' @param n_species integer number of species
#' @param n_resources interger number of resources
#' @param eff a matrix of efficiency. How efficient are resources
#' converted into biomass, negative values represent excreted resources
#' (default: \code{eff = randomE(n_species, n_resources)})
#' @param consumers numeric vector of species
#' (default: \code{consumers = runif(n = n_species, min = 0.1, max = 10)})
#' @param resources numeric vector of resources
#' (default: \code{resources = runif(n = n_resources, min = 1, max = 100)})
#' @param mumax numeric vector of maximum mu of species
#' (default: \code{mumax = rep(1, n_species)})
#' @param k_table a matrix of K values in monod model
#' (default: \code{k_table = matrix(rgamma(n=n_species*n_resources,
#' shape = 50, rate = 0.25), nrow = n_species)})
#' @param t_end numeric scalar indicating the final time of the simulation
#' (default: \code{t_end = 1000})
#' @param ... additional parameters including 't_start', 't_step', and 't_store'
#'
#' @seealso
#' \code{\link[miaSim:convertToSE]{convertToSE}}
#'
#' @examples
#' x <- simulateConsumerResource(n_species = 2, n_resources = 4)
#'
#' @return an abundance matrix with species and resources abundance as rows and
#' time points as columns
#'
#' @export
simulateConsumerResource <- function(n_species, n_resources,
    eff = randomE(n_species, n_resources),
    consumers = runif(n = n_species, min = 0.1, max = 10),
    resources = runif(n = n_resources, min = 1, max = 100),
    mumax = rep(1, n_species),
    k_table = matrix(rgamma(n=n_species*n_resources, shape = 50,rate = 0.25),
                        nrow = n_species), t_end = 1000, ...){

    # define the consumer-resource model
    consumerResourceModel <- function(t, state, params){
        with(as.list(c(state, params)),{
            consumers <- pmax(0, state[startsWith(names(state), "consumer")])
            resources <- pmax(0, state[startsWith(names(state), "resource")])
            mumax <- params[['mumax']]
            eff <- params[['eff']]
            k_table <- params[['k_table']]
            # B = R + K, sums of resources and K,
            # this matrix is created to facilitate the monod growth calculation
            B <- matrix(rep(resources, length(consumers)),
                ncol = length(resources), byrow = TRUE) + k_table
            growth <- ((eff*(eff>0)/B) %*% resources)*consumers
            consumption <- (t(growth) %*% ((eff>0)/B))*resources
            production <- -(t(growth) %*% (eff*(eff<0)/B))*resources
            dResources <- - consumption + production
            dConsumers <- mumax*growth
            dxdt <- list(c(dConsumers, dResources))
            return(dxdt)
        })
    }
    state_init <- c(consumers, resources)
    names(state_init) <- c(paste0("consumer", seq(length(consumers))),
        paste0("resource", seq(length(resources))))
    parameters <- list(mumax = mumax, eff = eff, k_table = k_table)

    t_dyn <- simulationTimes(t_end = t_end, ...)
    out_matrix <- as.data.frame(ode(y = state_init, times = t_dyn$t_sys,
        func = consumerResourceModel, parms = parameters))
    out_matrix <- out_matrix[t_dyn$t_index,]
    counts <- as.matrix(out_matrix[,2:ncol(out_matrix)])

    return(counts)

}
