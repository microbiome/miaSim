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
#' @param n.species Integer: number of species
#' @param n.resources Interger: number of resources
#' @param eff matrix: matrix of efficiency. How efficient are resources
#' converted into biomass, negative values represent excreted resources
#' (default: \code{eff = randomE(n.species, n.resources)})
#' @param consumers Numeric: vector of species
#' (default: \code{consumers = runif(n = n.species, min = 0.1, max = 10)})
#' @param resources Numeric: vector of resources
#' (default: \code{resources = runif(n = n.resources, min = 1, max = 100)})
#' @param mumax Numeric: vector of maximum mu of species
#' (default: \code{mumax = rep(1, n.species)})
#' @param k.table matrix: matrix of K values in monod model
#' (default: \code{k.table = matrix(rgamma(n=n.species*n.resources,
#' shape = 50, rate = 0.25), nrow = n.species)})
#' @param t.end Numeric scalar indicating the final time of the dimulation
#' (default: \code{t.end = 1000})
#' @param ... additional parameters including 't.start', 't.step', and 't.store'
#'
#' @examples
#' # example1 users provide least parameters.
#' ExampleConsumerResource <- simulateConsumerResource(n.species = 2,
#' n.resources = 4)
#' # visualize the dynamics of the model
#' matplot(ExampleConsumerResource, type = "l")
#'
#' @return an abundance matrix with species and resources abundance as rows and
#' time points as columns
#'
#' @export
simulateConsumerResource <- function(n.species, n.resources,
    eff = randomE(n.species, n.resources),
    consumers = runif(n = n.species, min = 0.1, max = 10),
    resources = runif(n = n.resources, min = 1, max = 100),
    mumax = rep(1, n.species),
    k.table = matrix(rgamma(n=n.species*n.resources,
        shape = 50,rate = 0.25), nrow = n.species),
    t.end = 1000, ...){

    # define the consumer-resource model
    consumerResourceModel <- function(t, state, params){
        with(as.list(c(state, params)),{
            consumers <- pmax(0, state[startsWith(names(state), "consumer")])
            resources <- pmax(0, state[startsWith(names(state), "resource")])
            mumax <- params[['mumax']]
            eff <- params[['eff']]
            k.table <- params[['k.table']]
            # B = R + K, sums of resources and K,
            # this matrix is created to facilitate the monod growth calculation
            B <- matrix(rep(resources, length(consumers)),
                ncol = length(resources), byrow = TRUE) + k.table
            growth <- ((eff*(eff>0)/B) %*% resources)*consumers
            consumption <- (t(growth) %*% ((eff>0)/B))*resources
            production <- -(t(growth) %*% (eff*(eff<0)/B))*resources
            dResources <- - consumption + production
            dConsumers <- mumax*growth
            dxdt <- list(c(dConsumers, dResources))
            return(dxdt)
        })
    }
    state.init <- c(consumers, resources)
    names(state.init) <- c(paste0("consumer", seq(length(consumers))),
        paste0("resource", seq(length(resources))))
    parameters <- list(mumax = mumax, eff = eff, k.table = k.table)

    t.dyn <- simulationTimes(t.end = t.end, ...)
    out.matrix <- as.data.frame(ode(y = state.init, times = t.dyn$t.sys,
        func = consumerResourceModel, parms = parameters))
    out.matrix <- out.matrix[t.dyn$t.index,]

    return(as.matrix(out.matrix[,2:ncol(out.matrix)]))

}
