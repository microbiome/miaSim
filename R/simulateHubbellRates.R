#' Hubbell's neutral model simulation applied to time series
#'
#' Neutral species abundances simulation according to the Hubbell model. This
#' model shows that losses in society can be replaced either by the birth of
#' individuals or by immigration depending on their probabilities.
#' The specific time between the events of birth or migration is calculated and
#' time effect is considered to determine the next event.
#'
#' @param n_species Integer: number of species
#' @param x0 Numeric: initial species composition. If NULL, 
#' `rep(100, n_species)` is used.
#' @param names_species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n_species))` is used.
#' (default: \code{names_species = NULL})
#' @param migration_p Numeric: the probability/frequency of migration from a 
#' metacommunity
#' (default: \code{migration_p = 0.01})
#' @param metacommunity_probability Numeric: normalized probability distribution
#' of the likelihood that species from the metacommunity can enter the community
#' during the simulation. If NULL, `rdirichlet(1, alpha = rep(1,n_species))` is 
#' used.
#' (default: \code{metacommunity_probability = NULL})
#' @param k_events Integer: number of events to simulate before updating the 
#' sampling distributions.
#' (default: \code{k_events = 1})
#' @param growth_rates Numeric: maximum growth rates(mu) of species.
#' If NULL, `rep(1, n_species)` is used.
#' (default: \code{growth_rates = NULL})
#' @param error_variance Numeric: the variance of measurement error.
#' By default it equals to 0, indicating that the result won't contain any 
#' measurement error. This value should be non-negative.
#' (default: \code{error_variance = 0})
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t_end Numeric: simulation end time (default: \code{t_end = 1000})
#' @param ... additional parameters including 't_start', 't_step', and 't_store'
#' see \code{\link{utils}} for more information.
#'
#' @seealso
#' \code{\link[miaSim:convertToSE]{convertToSE}}
#'
#' @examples
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(n_species = 5)
#' makePlot(ExampleHubbellRates$matrix)
#' 
#' # no migration, all stochastic birth and death
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(n_species = 5, migration_p = 0)
#' makePlot(ExampleHubbellRates$matrix)
#' 
#' # all migration, no stochastic birth and death
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(
#'     n_species = 5, 
#'     migration_p = 1, 
#'     metacommunity_probability = c(0.1, 0.15, 0.2, 0.25, 0.3), 
#'     t_end = 20, 
#'     t_store = 200)
#' makePlot(ExampleHubbellRates$matrix)
#' 
#' # all migration, no stochastic birth and death, but with measurement errors
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(
#'     n_species = 5, 
#'     migration_p = 1, 
#'     metacommunity_probability = c(0.1, 0.15, 0.2, 0.25, 0.3), 
#'     t_end = 20, 
#'     t_store = 200, 
#'     error_variance = 100)
#' makePlot(ExampleHubbellRates$matrix) 
#' 
#' # model with specified inputs
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(
#'     n_species = 5,
#'     migration_p = 0.1,
#'     metacommunity_probability = c(0.1, 0.15, 0.2, 0.25, 0.3), 
#'     t_end = 200, 
#'     t_store = 1000, 
#'     k_events = 5,
#'     growth_rates = c(1.1, 1.05, 1, 0.95, 0.9))
#' makePlot(ExampleHubbellRates$matrix)
#' 
#' @return \code{simulateHubbellRates} returns a list of initial states, 
#' parameters of the model, including a matrix with species abundance as rows 
#' and time points as columns.
#' 
#' @docType methods
#' @aliases simulateHubbellRates-numeric
#' @aliases simulateHubbellRates,numeric-method
#'
#' @importFrom MatrixGenerics rowSums2
#' @importFrom gtools rdirichlet
#' @importFrom stats rbinom
#' @importFrom stats rgamma
#' @importFrom stats rmultinom
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
    t_end=1000,...){
    # set the default values
    if (is.null(n_species) & !is.null(x0)) {
        n_species <- length(x0)
    } else if (is.null(x0) & !is.null(n_species)) {
        x0 <- rep(100, n_species)
    } else if (is.null(x0) & is.null(n_species)){
        stop("At least one of x0 and n_species shall be given.")
    } else {
        if (n_species != length(x0)) stop("n_species and length of x0 is not identical.")
    }

    #input check
    if(!isPositiveInteger(n_species)){
        stop("n_species must be positive integer.")}
    if (!is.logical(norm)) stop('norm" should be a logical variable.')
    if (!all(x0 >= 0)) stop("x0 should be non negative.")

    if (is.null(names_species)) {
        names_species <- paste0("sp", seq_len(n_species))
    }

    if (is.null(metacommunity_probability)){
        metacommunity_probability <- rdirichlet(1, alpha = rep(1,n_species))
    }

    # normalize metacommunity_probability
    metacommunity_probability <- metacommunity_probability/
        sum(metacommunity_probability)

    if (is.null(growth_rates)){
        growth_rates <- rep(1,n_species)
    }

    t_dyn <- simulationTimes(t_end = t_end,...)
    t_store <- length(t_dyn$t_index)
    
    birth_p <- 1 - migration_p
    community <- x0

    propensities <- sum(community)*(c(migration_p, birth_p))

    out_matrix <- matrix(0, nrow=length(t_dyn$t_index), ncol = n_species)
    out_matrix[1,] <- x0
    
    stored_time <- t_dyn$t_sys[t_dyn$t_index]
    current_t <- stored_time[1]
    last_stored_t <- stored_time[1]

    while(current_t <= t_end){

        tau_events <- min(min(community[community>0]),k_events)
        tau <- rgamma(n = 1, shape = tau_events, scale = 1/(sum(propensities)))

        current_t <- current_t + tau

        composition_propensities <- community*growth_rates

        composition_probabilities <- composition_propensities/sum(composition_propensities)

        #k deaths
        community <- community -
            (rmultinom(n=1, size=tau_events, prob=community))

        n_births <- rbinom(n=1, size=tau_events, prob=birth_p)
        n_migration <- tau_events-n_births

        community <- community +
          (rmultinom(n=1, size=n_births, prob=composition_probabilities)) +
          (rmultinom(n=1, size=n_migration, prob=metacommunity_probability))

        index <- ((current_t >= stored_time )  & (last_stored_t < stored_time))

        if (sum(index)>0) {
            out_matrix[index,] <- t(matrix(rep(t(community),sum(index)), 
                ncol = sum(index)))
            last_stored_t <- stored_time[max(seq(t_store)[index])]

        }
    }

    if(error_variance > 0){
        measurement_error <- rnorm(n = length(t_dyn$t_index)*n_species, 
                                   mean = 0, sd = sqrt(error_variance))
        measurement_error <- matrix(measurement_error, 
                                    nrow = length(t_dyn$t_index))
        out_matrix <- out_matrix + measurement_error
    }

    if(norm){
        out_matrix <- out_matrix/rowSums(out_matrix)
    }

    colnames(out_matrix) <- names_species
    
    out_matrix <- cbind(out_matrix, time = t_dyn$t_sys[t_dyn$t_index])
    
    #out_matrix$t <- t_dyn$t_sys[t_dyn$t_index]
    #SE <- SummarizedExperiment(assays = list(counts = out_matrix))
    out_list <- list(matrix = out_matrix, 
        community = community, 
        x0 = x0,
        metacommunity_probability = metacommunity_probability,
        error_variance = error_variance,
                            growth_rates = growth_rates)
    return(out_list)
}
