#' Self-Organised Instability model (SOI) simulation
#'
#' Generate time-series with The Self-Organised Instability (SOI) model.
#' Implements a K-leap method for accelerating stochastic simulation.
#'
#' @template man_spe
#' @param x0 \code{Numeric scalar}. Specifies initial community abundances
#' If \code{NULL}, based on migration rates. (Default: \code{NULL})
#' @param carrying_capacity \code{Integer scalar}. Indicates community size, 
#' number of available sites (individuals). (Default: \code{1000})
#' @param A \code{Matrix}. Defines the positive and negative
#' interactions between n_species. If \code{NULL}, `powerlawA(n_species)` is used.
#' (Default: \code{NULL})
#' @param k_events \code{Integer scalar}. Indicates the number of transition 
#' events that are allowed to take place during one leap. Higher values reduce 
#' runtime, but also accuracy of the simulation. (Default: \code{5}).
#' @param t_end \code{Numeric scalar}. Specifies the end time of the simulation, 
#' defining the modeled time length of the community. (Default: \code{1000})
#' @param metacommunity_probability \code{Numeric scalar}: Indicates the 
#' normalized probability distribution of the likelihood that species from 
#' the metacommunity can enter the community during the simulation. 
#' (Default: \code{runif(n_species, min = 0.1, max = 0.8)})
#' @param death_rates \code{Numeric scalar}. Indicates the death rates of each species.
#' (Default: \code{runif(n_species, min = 0.01, max = 0.08)})
#' @param norm \code{Logical scalar}. Whether the time series should be
#' returned with the abundances as proportions (\code{norm = TRUE})
#' or the raw counts. (Default: \code{FALSE})
#'
#' @return
#' \code{simulateSOI} returns a TreeSummarizedExperiment class object
#'
#' @examples
#' # Generate interaction matrix
#' A <- miaSim::powerlawA(10, alpha = 1.2)
#' # Simulate data from the SOI model
#' tse <- simulateSOI(
#'     n_species = 10, carrying_capacity = 1000, A = A,
#'     k_events = 5, x0 = NULL, t_end = 150, norm = TRUE
#' )
#'
#' @importFrom stats rgamma rmultinom
#' @importFrom MatrixGenerics colSums2
#'
#' @export
simulateSOI <- function(n_species,
    x0 = NULL,
    names_species = NULL,
    carrying_capacity = 1000,
    A = NULL,
    k_events = 5,
    t_end = 1000,
    metacommunity_probability = runif(n_species, min = 0.1, max = 0.8),
    death_rates = runif(n_species, min = 0.01, max = 0.08),
    norm = FALSE) {

    # input check
    if (!all(vapply(
        list(
            n_species,
            carrying_capacity,
            k_events
        ),
        .isPosInt, logical(1)
    ))) {
        stop("n_species,carrying_capacity and k_events values must be integer.")
    }

    if (is.null(x0)) {
        counts <- c(metacommunity_probability, runif(1))
        counts <- round((counts / sum(counts)) * carrying_capacity)
    } else if (length(x0) == n_species + 1) {
        counts <- x0
    } else if (length(x0) == n_species) {
        empty_count <- carrying_capacity - sum(x0)
        counts <- c(x0, empty_count)
    } else {
        stop("x0 vector argument can only be either of length N, or N+1
             in case the number of empty sites are included")
    }

    if (is.null(names_species)) {
        names_species <- paste0("sp", seq_len(n_species))
    }
    if (is.null(A)) {
        A <- powerlawA(n_species = n_species)
    }

    # making sure counts vector adds up to carrying_capacity with the empty sites
    # (if not, the difference is added or subtracted from the empty
    # sites to get a total of carrying_capacity)
    if (sum(counts) != carrying_capacity) {
        counts[n_species + 1] <- counts[n_species + 1] + (carrying_capacity - sum(counts))
    }
    # TRANSITION MATRIX & INTERACTION RATES
    # Initialise trans_mat transition matrix representing the actual
    # reactions
    # First n_species columns = migration: species_i counts go up with 1
    # empty sites go down with -1
    # Next n_species columns = death: species_i counts go down with -1,
    # empty sites go up with +1
    # Last columns: the competition / negative interaction jumps:
    trans_mat <- cbind(
        rbind(
            diag(1, ncol = n_species, nrow = n_species),
            rep(-1, times = n_species)
        ),
        rbind(
            diag(-1, ncol = n_species, nrow = n_species),
            rep(1, times = n_species)
        )
    )
    # Interaction:
    # happens when sp_stronger interacts stronger with sp_weaker,
    # than sp_weaker with sp_stronger
    pos_inter_rates <- matrix(0, nrow = n_species, ncol = n_species)
    neg_inter_rates <- pos_inter_rates # copy to initialise
    # pos interaction and immigration
    # AND sp_stronger interacts positively with sp_weaker
    for (stronger in seq_len(n_species)) {
        weaker <- which(A[stronger, ] > A[, stronger] & A[, stronger] >= 0)
        if (length(weaker) > 0) {
            rate <- A[stronger, weaker] + A[weaker, stronger]
            pos_inter_rates[stronger, weaker] <- rate
        }
    }
    # neg interaction: competition
    # AND sp_weaker interacts negatively with sp_stronger
    for (stronger in seq_len(n_species)) {
        weaker_vec <- which(A[stronger, ] > A[, stronger] & A[, stronger] < 0)
        if (length(weaker_vec) > 0) {
            rate <- A[stronger, weaker_vec] - A[weaker_vec, stronger]
            neg_inter_rates[stronger, weaker_vec] <- rate
            for (weaker in weaker_vec) {
                jump <- rep(0, times = n_species + 1)
                jump[stronger] <- 1
                jump[weaker] <- -1
                trans_mat <- cbind(trans_mat, jump)
            }
        }
    }

    # INITIAL PROPENSITIES
    propensities <- updatePropensities(
        carrying_capacity, counts, death_rates,
        metacommunity_probability, pos_inter_rates, neg_inter_rates
    )
    # K-LEAPS SIMULATION
    sys_t <- seq(from = 0, to = t_end, length.out = t_end)
    # time variables
    current_t <- 0 # current time
    sample_t <- sys_t[1] # sampled time
    series_t <- 1 # column to be filled with counts in series matrix
    # (next generation)
    series <- matrix(0, nrow = n_species, ncol = t_end)
    colnames(series) <- paste0("t", seq_len(t_end))
    rownames(series) <- names_species
    while (series_t <= t_end) {
        c <- sum(propensities)
        p <- propensities / c
        tau <- rgamma(n = 1, shape = k_events, scale = 1 / c)
        current_t <- current_t + tau
        # if reached end of simulation:
        if (current_t >= t_end) {
            series[, t_end] <- series
            break
        }
        # if current_t exceeds sample_t, update series with current series
        if (current_t > sample_t) {
            series[, series_t] <- counts[seq_len(n_species)]
            series_t <- series_t + 1
            index <- which(sys_t >= sample_t & sys_t < current_t)
            sample_t <- sys_t[index][length(index)]
        }
        # which reactions allowed?
        reactionsToFire <- as.vector(rmultinom(n = 1, size = k_events, prob = p))
        # update counts with allowed transitions (reactions)
        counts <- vapply(counts + trans_mat %*% reactionsToFire,
            FUN = function(x) {
                x <- max(0, x)
            }, numeric(1)
        )
        # update propensities with updated counts
        propensities <- updatePropensities(
            carrying_capacity, counts, death_rates,
            metacommunity_probability, pos_inter_rates, neg_inter_rates
        )
    }
    if (norm) {
        series <- t(t(series) / colSums2(series))
    }
    counts <- series
    out_matrix <- cbind(t(counts), time = seq_len(ncol(counts)))

    TreeSE <- TreeSummarizedExperiment(
        assays = list(counts = t(out_matrix[, 1:n_species])),
        colData = DataFrame(time = out_matrix[, "time"]),
        metadata = list(x0 = x0,
                        A = A,
                        carrying_capacity = carrying_capacity,
                        k_events = k_events))

    return(TreeSE)
}

updatePropensities <- function(carrying_capacity, # total nr of sites
    counts, # species counts vector incl the empty sites
    death_rates, # species-specific death rates
    metacommunity_probability, # species-specific migration rates
    pos_inter_rates, # positive interaction rates matrix
    neg_inter_rates # negative interaction rates matrix
) {
    N <- length(counts) - 1
    # migration AND positive interaction
    m_propensities <- counts[N + 1] * metacommunity_probability / N
    for (stronger in seq_len(N)) {
        if (sum(pos_inter_rates[stronger, ] != 0) > 0) {
            weaker <- which(pos_inter_rates[stronger, ] != 0)
            inter_rate <- pos_inter_rates[stronger, weaker]
            m_propensities[stronger] <- m_propensities[stronger] +
                sum(counts[weaker] *
                    inter_rate) / carrying_capacity * counts[stronger] / carrying_capacity * counts[N + 1]
        }
    }
    # death / extinction
    d_propensities <- (counts[seq_len(N)] * death_rates)
    # competition / negative interaction
    j_propensities <- c()
    for (stronger in seq_len(N)) {
        if (sum(neg_inter_rates[stronger, ] != 0) > 0) {
            weaker <- which(neg_inter_rates[stronger, ] != 0)
            inter_rate <- neg_inter_rates[stronger, weaker]
            j_propensities <- c(
                j_propensities,
                counts[stronger] * counts[weaker] * inter_rate / carrying_capacity
            )
        }
    }
    propensities <- c(m_propensities, d_propensities, j_propensities)
    return(propensities)
}
