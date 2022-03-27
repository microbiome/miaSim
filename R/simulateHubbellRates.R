#' Hubbell's neutral model simulation applied to time series
#'
#' Neutral species abundances simulation according to the Hubbell model. This
#' model shows that losses in society can be replaced either by the birth of
#' individuals or by immigration depending on their probabilities.
#' The specific time between the events of birth or migration is calculated and
#' time effect is considered to determine the next event.
#'
#' @param community_initial numeric value a vector of integers,
#' containing species counts greater or equal to zero.
#' @param migration_p numeric immigration possibility. It defines the
#' probability of migration when replacement is needed in the community.
#' The value can be between 0 and 1. The sum of the probability of
#' migration and the probability birth must be 1.
#' @param metacommunity_p numeric value the probability of a species being
#' found in the metacommunity.
#' @param k_events integer number of steps performed at a time point.
#' It can be equal or more than 1. Bigger k_events increases speed while
#' decreasing precision.
#' @param growth_rates numeric rate of change in community size.
#' @param norm logical scalar choosing whether the time series should be
#' returned with the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t_end numeric value of simulation end time
#' (default: \code{t_end = 1000})
#' @param list logical scalar deciding whether output is a list object or not
#' (default: \code{norm = TRUE})
#' @param ... additional parameters including 't_start', 't_step', and 't_store'
#'
#' @seealso
#' \code{\link[miaSim:convertToSE]{convertToSE}}
#'
#' @examples
#' x <- simulateHubbellRates(community_initial = c(0,5,10),
#'              migration_p = 0.01, metacommunity_p = NULL, k_events = 1,
#'              growth_rates = NULL, norm = FALSE, t_end=1000)
#'
#' @return a community abundance matrix or a list object that contains
#' growth rates, time points and metacommunity probabilities
#'
#' @importFrom MatrixGenerics rowSums2
#' @importFrom gtools rdirichlet
#' @importFrom stats rgamma
#' @importFrom S4Vectors DataFrame
#'
#' @references Rosindell J, Hubbell SP, Etienne RS. The unified neutral theory
#' of biodiversity and biogeography at age ten. Trends Ecol Evol.
#' 2011 Jul;26(7):340-8. doi: 10.1016/j.tree.2011.03.024.
#' Epub 2011 May 10. PMID: 21561679.
#'
#' @export
simulateHubbellRates <- function(community_initial,
                                migration_p = 0.1,
                                metacommunity_p = NULL,
                                k_events = 1,
                                growth_rates = NULL,
                                norm = FALSE,
                                t_end = 1000,
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

    t_dyn <- simulationTimes(t_end = t_end, ...)

    t_store <- length(t_dyn$t_index)

    n_species <- length(community_initial)

    birth_p <- 1 - migration_p

    community <- community_initial

    if (is.null(metacommunity_p)){
        metacommunity_p <- rdirichlet(1, alpha = rep(1,n_species))
    }

    metacommunity_p <- metacommunity_p/sum(metacommunity_p)

    if (is.null(growth_rates)){
        growth_rates <- rep(1,n_species)
    }

    counts <- matrix(0, nrow = n_species, ncol = length(t_dyn$t_index))

    counts[,1] = community_initial

    stored_time <- t_dyn$t_sys[t_dyn$t_index]
    current_t <- stored_time[1]
    last_stored_t <- stored_time[1]

    while(current_t <= t_end){

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

            counts[,index] <-
                t(matrix(rep(t(community),sum(index)), ncol = sum(index)))
            last_stored_t <- stored_time[max(seq(t_store)[index])]

        }
    }

    rownames(counts) <- seq_len(n_species)
    rownames(counts) <- paste("s", rownames(counts), sep = "_")
    colnames(counts) <- (t_dyn$t_sys[t_dyn$t_index])

    if(norm){
        output <- counts/rowSums2(counts)
    }
    if(list){
        timepoints <- c(t_dyn$t_sys[t_dyn$t_index])
        time_int <- diff(timepoints)
        time_int[t_end] <- NA

        col_data <- DataFrame(time = timepoints, time_interval = time_int)
        meta_data <- list(metacommunity_p = metacommunity_p,
                            growth_rates = growth_rates)
        counts<- list(counts = counts, time = col_data,
                            metadata = meta_data)
    }
    return(counts)
}
