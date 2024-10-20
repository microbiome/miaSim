#' Generalized Lotka-Volterra (gLV) simulation
#'
#' Simulates time series with the generalized Lotka-Volterra model.
#'
#' Simulates a community time series using the generalized Lotka-Volterra model,
#' defined as dx/dt = x(b+Ax), where x is the vector of species abundances,
#' diag(x) is a diagonal matrix with the diagonal values set to x.
#' A is the interaction matrix and b is the vector of growth rates.
#'
#' @template man_spe
#' @param A \code{Matrix}. Interaction matrix defining the positive and negative
#' interactions between n_species. If \code{NULL}, `randomA(n_species)` is used.
#' (Default: \code{NULL})
#' @param x0 \code{Numeric scalar}. Indicates the initial abundances of 
#' simulated species. If \code{NULL}, `runif(n = n_species, min = 0, max = 1)` 
#' is used. (Default: \code{NULL})
#' @param growth_rates \code{Numeric scalar}. Indicates the growth rates 
#' of simulated species. If \code{NULL}, `runif(n = n_species, min = 0, max = 1)` 
#' is used. (Default: \code{NULL})
#' @template man_sto
#' @template man_mig
#' @template man_mod
#'
#' @return
#' \code{simulateGLV} returns a TreeSummarizedExperiment class object
#'
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @examples
#'
#' # generate a random interaction matrix
#' ExampleA <- randomA(n_species = 4, diagonal = -1)
#'
#' # run the model with default values (only stochastic migration considered)
#' tse <- simulateGLV(n_species = 4, A = ExampleA)
#'
#' # run the model with two external disturbances at time points 240 and 480
#' # with durations equal to 1 (10 time steps when t_step by default is 0.1).
#' ExampleGLV <- simulateGLV(
#'     n_species = 4, A = ExampleA,
#'     t_external_events = c(0, 240, 480), t_external_durations = c(0, 1, 1)
#' )
#'
#' # run the model with no perturbation nor migration
#' set.seed(42)
#' tse1 <- simulateGLV(
#'     n_species = 4, A = ExampleA, stochastic = FALSE,
#'     sigma_migration = 0
#' )
#'
#' # run the model with no perturbation nor migration but with measurement error
#' set.seed(42)
#' tse2 <- simulateGLV(
#'     n_species = 4, A = ExampleA, stochastic = FALSE,
#'     error_variance = 0.001, sigma_migration = 0
#' )
#'
#' @importFrom deSolve ode
#' @importFrom stats runif rnorm
#'
#' @export
simulateGLV <- function(n_species,
    names_species = NULL,
    A = NULL,
    x0 = NULL,
    growth_rates = NULL,
    sigma_drift = 0.001,
    sigma_epoch = 0.1,
    sigma_external = 0.3,
    sigma_migration = 0.01,
    epoch_p = 0.001,
    t_external_events = NULL,
    t_external_durations = NULL,
    stochastic = TRUE,
    migration_p = 0.01,
    metacommunity_probability = NULL,
    error_variance = 0,
    norm = FALSE,
    t_end = 1000, ...) {

    # set the default values
    if (is.null(names_species)) {
        names_species <- paste0("sp", seq_len(n_species))
    }
    if (is.null(A)) {
        A <- randomA(n_species)
    }
    if (is.null(x0)) {
        x0 <- runif(n = n_species, min = 0, max = 1)
    }
    if (is.null(growth_rates)) {
        growth_rates <- runif(n = n_species, min = 0, max = 1)
    }
    if (is.null(metacommunity_probability)) {
        metacommunity_probability <- rdirichlet(1, alpha = rep(1, n_species))
    }

    # normalize metacommunity_probability
    metacommunity_probability <- metacommunity_probability /
        sum(metacommunity_probability)

    # select the time points to simulate and to store
    t_dyn <- .simulationTimes(t_end = t_end, ...)

    # calculate the time points influenced by the disturbances
    tEvent <- simulateEventTimes(
        t_events = t_external_events,
        t_duration = t_external_durations,
        t_end = t_end,
        ... = ...
    )

    parameters <- list(
        growth_rates = growth_rates, A = A, n_species = n_species,
        sigma_drift = sigma_drift, stochastic = stochastic,
        sigma_epoch = sigma_epoch, epoch_p = epoch_p,
        sigma_external = sigma_external, tEvent = tEvent,
        migration_p = migration_p,
        metacommunity_probability = metacommunity_probability,
        sigma_migration = sigma_migration
    )

    out <- ode(
        y = x0,
        times = t_dyn$t_sys,
        func = glvModel,
        parms = parameters,
        events = list(func = perturb, time = t_dyn$t_sys),
        maxsteps = 10^9, method = "ode45"
    )


    out_matrix <- out[, 2:ncol(out)]

    out_matrix <- out_matrix[t_dyn$t_index, ]

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
        metadata = list(migration_p = migration_p,
                        error_variance = error_variance))


    return(TreeSE)
}

# define the GLV Model
glvModel <- function(t, x0, parameters) {
    x0[x0 < 10^-8] <- 0
    growth_rates <- parameters$growth_rates
    A <- parameters$A
    # rate of change
    dx <- x0 * (growth_rates + A %*% x0)
    # return rate of change
    list(dx)
}

# define the perturbation event
perturb <- function(t, y, parameters) {
    with(as.list(y), {
        # continuous or episodic perturbation
        epoch_rN <- 0
        external_rN <- 0
        migration_rN <- 0
        if (rbinom(1, 1, parameters$epoch_p)) {
            epoch_rN <- rnorm(parameters$n_species, sd = parameters$sigma_epoch)
            epoch_rN <- parameters$stochastic * epoch_rN
        }

        if (rbinom(1, 1, parameters$migration_p)) {
            migration_rN <- rmultinom(
                n = 1, size = 1,
                prob = parameters$metacommunity_probability
            )[, ] * abs(rnorm(
                n = 1,
                mean = 0, sd = parameters$sigma_migration
            ))
            # TODO: is migration also stochastic? if so, add the following:####
            # migration_rN <- parameters$stochastic*migration_rN
        }

        if (t %in% parameters$tEvent) {
            external_rN <- rnorm(parameters$n_species,
                sd = parameters$sigma_external
            )
            external_rN <- parameters$stochastic * external_rN
        }
        drift_rN <- rnorm(parameters$n_species, sd = parameters$sigma_drift)
        drift_rN <- parameters$stochastic * drift_rN

        # perturbation is applied to the current population
        y <- y * (1 + drift_rN) * (1 + epoch_rN) * (1 + external_rN) + migration_rN
        return(y * (y > 0))
    })
}
