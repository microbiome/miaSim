#' Hubbell's neutral model simulation applied to time series
#'
#' Neutral species abundances simulation according to the Hubbell model.This
#' model shows that losses in society can be replaced either by the birth of
#' individuals or by immigration depending on their probabilities.
#' The specific time between the events of birth or migration is calculated and
#' time effect is considered to determine the next event.
#'
#' @param community.initial Numeric: a vector of integers, containing species
#' counts greater or equal to zero.
#' @param migration.p Numeric: immigration possibility. It defines the
#' probability of replacement in the community. The value can be equal
#' be between 0 and 1. The sum of probability of migration and the probability
#' birth has to be 1.
#' @param metacommunity.p Numeric: the probability to find a species in the
#' metacommunity is the weighted average of probabilities in communities.
#' @param k.events Integer:
#' @param growth.rates Numeric: the rate of the change in the community size
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t.end Numeric: simulation end time (default: \code{t.end = 1000})
#' @param ... additional parameters including 't.start', 't.step', and 't.store'
#'
#' @examples
#' com_1 <- c(0, 5, 10)
#' com_2 <- c(20, 40, 60)
#' community.initial <- t(array(c(com_1, com_2)))
#' HubbellRates <- simulateHubbellRates(community.initial, migration.p = 0.01,
#'              metacommunity.p = NULL, k.events = 1, growth.rates = NULL,
#'              norm = FALSE, t.end=1000)
#'
#' @return a list object containing the abundance matrix of each species,
#' number of individuals of each species before and after the simulation and the
#' metacommunity probability
#'
#' @importFrom gtools rdirichlet
#'
#' @references Rosindell, James et al. "The unified neutral theory of
#' biodiversity and biogeography at age ten." Trends in ecology & evolution
#' vol. 26,7 (2011).
#
#' @export
simulateHubbellRates <- function(community.initial,
                                 migration.p = 0.01,
                                 metacommunity.p = NULL,
                                 k.events = 1,
                                 growth.rates = NULL,
                                 norm = FALSE,
                                 t.end=1000, ...){

    t.dyn <- simulationTimes(t.end = t.end, ...)

    n.species <- length(community.initial)

    birth.p <- 1 - migration.p

    community <- community.initial

    if (is.null(metacommunity.p)){
        metacommunity.p <- rdirichlet(1, alpha = rep(1,n.species))
    }
    metacommunity.p <- metacommunity.p/sum(metacommunity.p)

    if (is.null(growth.rates)){
        growth.rates <- rep(1,n.species)
    }

    out.matrix <- matrix(0, nrow=length(t.dyn$t.index), ncol = n.species)

    current_t <- t.dyn$t.sys[1]

    current.sample.index <- t.dyn$t.index[1]

    next.sample.index <- t.dyn$t.index[2]

    while(current_t <= t.end){
        tau.events <- min(min(community[community>0]),k.events)
        propensities <- growth.rates*community
        probabilities <- propensities/(sum(propensities))
        tau <- rgamma(n = 1, shape = tau.events, scale = 1/(sum(propensities)))
        current_t <- current_t + tau

    if (current_t >= t.dyn$t.sys[next.sample.index]) {
        limit.sample.index <-
            max(seq(t.dyn$t.index)[t.dyn$t.sys[t.dyn$t.index]<=current_t])
        out.matrix[current.sample.index:limit.sample.index,] = community
        current.sample.index <- limit.sample.index
        next.sample.index <- limit.sample.index + 1
        }
    else if (current_t >= t.end){
        break
        }
    }

    #k deaths
    community <-
        community - t(rmultinom(n = 1, size = k.events, prob = probabilities))
    n.births <- sum(rbinom(n=tau.events, size=1, p = birth.p))
    n.migration <- tau.events-n.births

    community <-
        community + t(rmultinom(n = 1, size = n.births, prob = probabilities)) +
        t(rmultinom(n = 1, size = n.migration, prob = metacommunity.p))

    if(norm){
        out.matrix <- out.matrix/rowSums(out.matrix)
    }

    colnames(out.matrix) <- seq_len(n.species)

    out.matrix <- cbind(out.matrix, time = t.dyn$t.sys[t.dyn$t.index])

    #out.matrix$t <- t.dyn$t.sys[t.dyn$t.index]
    #SE <- SummarizedExperiment(assays = list(counts = out.matrix))
    out.list <- list(matrix = out.matrix,
                   community = community,
                   community.initial = community.initial,
                   metacommunity.p = metacommunity.p)
    return(out.list)
}
