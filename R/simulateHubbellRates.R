#' Hubbell's neutral model simulation applied to time series
#'
#' Neutral species abundances simulation according to the Hubbell model. This
#' model shows that losses in society can be replaced either by the birth of
#' individuals or by immigration depending on their probabilities.
#' The specific time between the events of birth or migration is calculated and
#' time effect is considered to determine the next event.
#'
#' @template man_spe
#' @param x0 \code{Numeric scalar}. Indicates the initial species 
#' composition. If \code{NULL}, `rep(100, n_species)` is used.
#' @template man_mig
#' @param k_events \code{Integer scalar}. Indicates the number of 
#' events to simulate before updating the sampling distributions.
#' (Default: \code{1})
#' @param growth_rates \code{Numeric scalar}. Indicates the maximum 
#' growth rates(mu) of species. If \code{NULL}, `rep(1, n_species)` is used.
#' (Default: \code{NULL})
#' @template man_mod
#'
#' @examples
#' set.seed(42)
#' tse <- simulateHubbellRates(n_species = 5)
#'
#' miaViz::plotSeries(tse, x = "time")
#'
#' # no migration, all stochastic birth and death
#' set.seed(42)
#' tse1 <- simulateHubbellRates(n_species = 5, migration_p = 0)
#'
#' # all migration, no stochastic birth and death
#' set.seed(42)
#' tse2 <- simulateHubbellRates(
#'     n_species = 5,
#'     migration_p = 1,
#'     metacommunity_probability = c(0.1, 0.15, 0.2, 0.25, 0.3),
#'     t_end = 20,
#'     t_store = 200
#' )
#'
#' # all migration, no stochastic birth and death, but with measurement errors
#' set.seed(42)
#' tse3 <- simulateHubbellRates(
#'     n_species = 5,
#'     migration_p = 1,
#'     metacommunity_probability = c(0.1, 0.15, 0.2, 0.25, 0.3),
#'     t_end = 20,
#'     t_store = 200,
#'     error_variance = 100
#' )
#'
#' # model with specified inputs
#' set.seed(42)
#' tse4 <- simulateHubbellRates(
#'     n_species = 5,
#'     migration_p = 0.1,
#'     metacommunity_probability = c(0.1, 0.15, 0.2, 0.25, 0.3),
#'     t_end = 200,
#'     t_store = 1000,
#'     k_events = 5,
#'     growth_rates = c(1.1, 1.05, 1, 0.95, 0.9)
#' )
#'
#' @return \code{simulateHubbellRates} returns a TreeSummarizedExperiment class
#' object
#'
#' @docType methods
#' @aliases simulateHubbellRates-numeric
#' @aliases simulateHubbellRates,numeric-method
#'
#' @importFrom stats rbinom rgamma rmultinom rnorm
#'
#' @references Rosindell, James et al. "The unified neutral theory of
#' biodiversity and biogeography at age ten." Trends in ecology & evolution
#' vol. 26,7 (2011).
#
#' @export
simulateHubbellRates <- function(n_species = NULL,
    x0 = NULL,
    names_species = NULL,
    migration_p = 0.01,
    metacommunity_probability = NULL,
    k_events = 1,
    growth_rates = NULL,
    error_variance = 0,
    norm = FALSE,
    t_end = 1000, ...) {
    # set the default values
    if (is.null(n_species) & !is.null(x0)) {
        n_species <- length(x0)
    } else if (is.null(x0) & !is.null(n_species)) {
        x0 <- rep(100, n_species)
    } else if (is.null(x0) & is.null(n_species)) {
        stop("At least one of x0 and n_species shall be given.")
    } else {
        if (n_species != length(x0)) stop("n_species and length of x0 is not identical.")
    }

    # input check
    if (!.isPosInt(n_species)) {
        stop("n_species must be positive integer.")
    }
    if (!is.logical(norm)) stop('norm" should be a logical variable.')
    if (!all(x0 >= 0)) stop("x0 should be non negative.")

    if (is.null(names_species)) {
        names_species <- paste0("sp", seq_len(n_species))
    }

    if (is.null(metacommunity_probability)) {
        metacommunity_probability <- rdirichlet(1, alpha = rep(1, n_species))
    }

    # normalize metacommunity_probability
    metacommunity_probability <- metacommunity_probability /
        sum(metacommunity_probability)

    if (is.null(growth_rates)) {
        growth_rates <- rep(1, n_species)
    }

    t_dyn <- .simulationTimes(t_end = t_end, ...)
    t_store <- length(t_dyn$t_index)

    birth_p <- 1 - migration_p
    community <- x0

    propensities <- sum(community) * (c(migration_p, birth_p))

    out_matrix <- matrix(0, nrow = length(t_dyn$t_index), ncol = n_species)
    out_matrix[1, ] <- x0

    stored_time <- t_dyn$t_sys[t_dyn$t_index]
    current_t <- stored_time[1]
    last_stored_t <- stored_time[1]

    while (current_t <= t_end) {
        tau_events <- min(min(community[community > 0]), k_events)
        tau <- rgamma(n = 1, shape = tau_events, scale = 1 / (sum(propensities)))

        current_t <- current_t + tau

        composition_propensities <- community * growth_rates

        composition_probabilities <- composition_propensities / sum(composition_propensities)

        # k deaths
        community <- community -
            (rmultinom(n = 1, size = tau_events, prob = community))

        n_births <- rbinom(n = 1, size = tau_events, prob = birth_p)
        n_migration <- tau_events - n_births

        community <- community +
            (rmultinom(n = 1, size = n_births, prob = composition_probabilities)) +
            (rmultinom(n = 1, size = n_migration, prob = metacommunity_probability))

        index <- ((current_t >= stored_time) & (last_stored_t < stored_time))

        if (sum(index) > 0) {
            out_matrix[index, ] <- t(matrix(rep(t(community), sum(index)),
                ncol = sum(index)
            ))
            last_stored_t <- stored_time[max(seq(t_store)[index])]
        }
    }

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
        metadata = list(x0 = x0,
                        metacommunity_probability = metacommunity_probability,
                        error_variance = error_variance,
                        growth_rates = growth_rates))

    return(TreeSE)
}
