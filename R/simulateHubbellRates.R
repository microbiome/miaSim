#' Hubbell's neutral model simulation applied to time series
#'
#' Neutral species abundances simulation according to the Hubbell model.This
#' model shows that losses in society can be replaced either by the birth of
#' individuals or by immigration depending on their probabilities.
#' The specific time between the events of birth or migration is calculated and
#' time effect is considered to determine the next event.
#'
#' @param community_initial Numeric: a vector of integers, containing species
#' counts greater or equal to zero.
#' @param migration_p Numeric: immigration possibility. It defines the
#' probability of replacement in the community. The value can be equal
#' be between 0 and 1. The sum of probability of migration and the probability
#' birth has to be 1.
#' @param metacommunity_p Numeric: the probability to find a species in the
#' metacommunity is the weighted average of probabilities in communities.
#' @param k_events Integer:
#' @param growth_rates Numeric: the rate of the change in the community size
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t.end Numeric: simulation end time (default: \code{t.end = 1000})
#' @param ... additional parameters including 't.start', 't.step', and 't.store'
#'
#' @examples
#' HubbellRates <- simulateHubbellRates(community_initial = c(0,5,10),
#'              migration_p = 0.01, metacommunity_p = NULL, k_events = 1,
#'              growth_rates = NULL, norm = FALSE, t.end=1000)
#'
#' @return a \linkS4class{SummarizedExperiment} object containing the community
#' abundance matrix and generated values: metacommunity probabilities, growth
#' rates and time points
#'
#' @importFrom gtools rdirichlet
#'
#' @references Rosindell, James et al. "The unified neutral theory of
#' biodiversity and biogeography at age ten." Trends in ecology & evolution
#' vol. 26,7 (2011).
#
#' @export
simulateHubbellRates <- function(community_initial,
                                migration_p = 0.01,
                                metacommunity_p = NULL,
                                k_events = 1,
                                growth_rates = NULL,
                                norm = FALSE,
                                t.end=1000,...){

    t.dyn <- simulationTimes(t.end = t.end, ...)

    t.store <- length(t.dyn$t.index)

    n.species <- length(community_initial)

    birth_p <- 1 - migration_p

    community <- community_initial

    if (is.null(metacommunity_p)){
        metacommunity_p <- rdirichlet(1, alpha = rep(1,n.species))
    }

    metacommunity_p <- metacommunity_p/sum(metacommunity_p)

    if (is.null(growth_rates)){
        growth_rates <- rep(1,n.species)
    }

    out_matrix <- matrix(0, nrow=length(t.dyn$t.index), ncol = n.species)

    out_matrix[1,] = community_initial

    current_t <- t.dyn$t.sys[2]

    current_sample_index <- t.dyn$t.index[2]

    next_sample_index <- t.dyn$t.index[3]

    while(current_t <= t.end){

        tau_events <- min(min(community[community>0]),k_events)

        propensities <- growth_rates*community

        probabilities <- propensities/(sum(propensities))

        tau <- rgamma(n = 1, shape = tau_events, scale = 1/(sum(propensities)))

        current_t <- current_t + tau



        # if reached end of simulation:
        if(current_t >= t.end){
            break
        }

        if (current_t >= t.dyn$t.sys[next_sample_index]) {

            limit_sample_index <-
                max(t.dyn$t.index[t.dyn$t.sys[t.dyn$t.index]<=current_t])

            out_matrix[seq(t.store)[t.dyn$t.index==current_sample_index]:seq(t.store)[t.dyn$t.index==limit_sample_index],] = community
            current_sample_index <- limit_sample_index
            next_sample_index <- limit_sample_index + 1
        }



        #k deaths
        community <- community -
            t(rmultinom(n = 1, size = k_events, prob = probabilities))

        n_births <- sum(rbinom(n=tau_events, size=1, prob = birth_p))
        n_migration <- tau_events-n_births

        community <- community +
            t(rmultinom(n = 1, size = n_births, prob = probabilities)) +
            t(rmultinom(n = 1, size = n_migration, prob = metacommunity_p))

    }

    if(norm){
        out_matrix <- out_matrix/rowSums(out_matrix)
    }

    colnames(out_matrix) <- seq_len(n.species)
    colnames(out_matrix) <- paste("s", colnames(out_matrix), sep = "_")
    rownames(out_matrix) <- t(t.dyn$t.sys[t.dyn$t.index])

    col_data <- DataFrame(t(out_matrix))

    SE <- SummarizedExperiment(assays = list(counts = out_matrix),
                               colData = col_data,
                metadata = list(metacommunity_p = metacommunity_p,
                                growth_rates = growth_rates))
    return(SE)
}
