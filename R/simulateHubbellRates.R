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
#' probability of migration when replacement is needed in the community.
#' The value can be between 0 and 1. The sum of the probability of
#' migration and the probability birth must be 1.
#' @param metacommunity_p Numeric: the probability of a species being found in
#' the metacommunity.
#' @param k_events Integer: number of times the replacement happening
#' @param growth_rates Numeric: the rate of change in community size.
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

    out_matrix <- matrix(0, nrow = n.species, ncol = length(t.dyn$t.index))

    out_matrix[,1] = community_initial

    stored_time <- t.dyn$t.sys[t.dyn$t.index]
    current_t <- stored_time[1]
    last_stored_t <- stored_time[1]

    while(current_t <= t.end){

        tau_events <- min(min(community[community>0]),k_events)

        propensities <- sum(community)*(c(migration_p, 1-migration_p))
        event_probabilities <- propensities/(sum(propensities))

        tau <- rgamma(n = 1, shape = tau_events, scale = 1/(sum(propensities)))

        current_t <- current_t + tau

        composition_probabilities <- community/sum(community)

        #k deaths

        community <- community -
            (rmultinom(n = 1, size = tau_events, prob = composition_probabilities))

        n_births <- rbinom(n=1, size=tau_events, prob = event_probabilities[2])

        n_migration <- tau_events-n_births

        community <- community +
            (rmultinom(n = 1, size = n_births, prob = composition_probabilities)) +
            (rmultinom(n = 1, size = n_migration, prob = metacommunity_p))

        index <- ((current_t >= stored_time )  & (last_stored_t < stored_time))

        if (sum(index)>0) {

            out_matrix[,index] <-
                t(matrix(rep(t(community),sum(index)), ncol = sum(index)))
            last_stored_t <- stored_time[max(seq(t.store)[index])]

        }
    }
    if(norm){
        out_matrix <- out_matrix/rowSums(out_matrix)
    }

    rownames(out_matrix) <- seq_len(n.species)
    rownames(out_matrix) <- paste("s", rownames(out_matrix), sep = "_")
    colnames(out_matrix) <- (t.dyn$t.sys[t.dyn$t.index])


    timepoints <- c(t.dyn$t.sys[t.dyn$t.index])
    time_int <- diff(timepoints)
    time_int[t.end] <- NA

    col_data <- DataFrame(time = timepoints,
                          time_interval = time_int)

    SE <- SummarizedExperiment(assays = list(counts = out_matrix),
                              colData = col_data,
                metadata = list(metacommunity_p = metacommunity_p,
                                growth_rates = growth_rates))
    return(SE)
}
