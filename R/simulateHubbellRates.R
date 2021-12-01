#' Hubbell's neutral model simulation applied to time series
#'
#' Neutral species abundances simulation according to the Hubbell model. This
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
#' @param k_events Integer: the number of steps performed at a time point.
#' It can be equal or more than 1. Bigger k_events increases speed while
#' decreasing precision.
#' @param growth_rates Numeric: the rate of change in community size.
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t.end Numeric: simulation end time (default: \code{t.end = 1000})
#' @param list Logical : decides whether output is a list object or not
#' (default: \code{norm = TRUE})
#' @param ... additional parameters including 't.start', 't.step', and 't.store'
#'
#' @examples
#' ExampleHubbellRates <- simulateHubbellRates(community_initial = c(0,5,10),
#'              migration_p = 0.01, metacommunity_p = NULL, k_events = 1,
#'              growth_rates = NULL, norm = FALSE, t.end=1000)
#'
#' @return a community abundance matrix or a list object that contains
#' growth rates, time points and metacommunity probabilities
#'
#' @importFrom gtools rdirichlet
#' @importFrom stats rgamma
#' @importFrom S4Vectors DataFrame
#'
#' @references Rosindell, James et al. "The unified neutral theory of
#' biodiversity and biogeography at age ten." Trends in ecology & evolution
#' vol. 26,7 (2011).
#
#' @export
simulateHubbellRates <- function(community_initial,
                                migration_p = 0.1,
                                metacommunity_p = NULL,
                                k_events = 1,
                                growth_rates = NULL,
                                norm = FALSE,
                                t.end = 1000,
                                list = TRUE,
                                ...){

    #input check
    i <- seq_len(length(community_initial))
    if(!all(vapply(community_initial[i],
            FUN = isZeroOrPositiveInteger, logical(1)))){
        stop("community size for each species must be equal to zero or greater
            than zero.")}

    if(!all(vapply(list(migration_p, k_events),
            FUN = is.numeric, logical(1)))){
        stop("migration possibility and k_events parameter must be positive
            numeric values.")}

    if(!is.logical(norm)){
        stop("'norm' must be TRUE or FALSE.", call. = FALSE)
    }

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
            (rmultinom(n = 1, size = tau_events,
                prob = composition_probabilities))

        n_births <- rbinom(n=1, size=tau_events, prob = event_probabilities[2])

        n_migration <- tau_events-n_births

        community <- community +
            (rmultinom(n = 1, size = n_births,
                prob = composition_probabilities)) +
            (rmultinom(n = 1, size = n_migration,
                prob = metacommunity_p))

        index <- ((current_t >= stored_time )  & (last_stored_t < stored_time))

        if (sum(index)>0) {

            out_matrix[,index] <-
                t(matrix(rep(t(community),sum(index)), ncol = sum(index)))
            last_stored_t <- stored_time[max(seq(t.store)[index])]

        }
    }

    rownames(out_matrix) <- seq_len(n.species)
    rownames(out_matrix) <- paste("s", rownames(out_matrix), sep = "_")
    colnames(out_matrix) <- (t.dyn$t.sys[t.dyn$t.index])

    if(norm){
        output <- out_matrix/rowSums(out_matrix)
    }
    if(list){
        timepoints <- c(t.dyn$t.sys[t.dyn$t.index])
        time_int <- diff(timepoints)
        time_int[t.end] <- NA

        col_data <- DataFrame(time = timepoints, time_interval = time_int)
        meta_data <- list(metacommunity_p = metacommunity_p,
                            growth_rates = growth_rates)
        out_matrix<- list(abundancematrix = out_matrix, time = col_data,
                            metadata = meta_data)
    }
    return(out_matrix)
}
