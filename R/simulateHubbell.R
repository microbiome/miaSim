#' Hubbell's neutral model simulation
#'
#' Neutral species abundances simulation according to the Hubbell model.
#'
#' @param n_species \code{Integer scalar}. Specifies the amount of 
#' different species initially in the local community.
#' @param M \code{Integer scalar}. Specifies the amount of different 
#' species in the metacommunity, including those of the local community
#' @param carrying_capacity \code{Integer scalar}. Indicates the fixed 
#' amount of individuals in the local community. (Default: \code{1000})
#' @param k_events \code{Inteer scalar}. Indicates the fixed amount of 
#' deaths of local community individuals in each generation (Default: \code{10})
#' @param migration_p \code{Numeric scalar}. The immigration rate; specifies 
#' the probability that a death in the local community is replaced by a 
#' migrant of the metacommunity rather than by the birth of a local 
#' community member (Default: \code{0.02})
#' @param t_skip \code{Integer scalar}. Indicates the number of generations 
#' that should not be included in the outputted species abundance matrix. 
#' (Default: \code{0})
#' @param t_end \code{Integer scalar}. Indicates the number of simulations 
#' to be simulated
#' @param norm \code{Logical scalar}. Whether the time series should be
#' returned with the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts. (Default: \code{FALSE})
#'
#' @examples
#' tse <- simulateHubbell(
#'     n_species = 8, M = 10, carrying_capacity = 1000, k_events = 50,
#'     migration_p = 0.02, t_end = 100
#' )
#'
#' @return \code{simulateHubbell} returns a TreeSummarizedExperiment class
#' object
#'
#' @importFrom stats rbinom rmultinom
#' @importFrom MatrixGenerics colSums2
#'
#' @references Rosindell, James et al. "The unified neutral theory of
#' biodiversity and biogeography at age ten." Trends in ecology & evolution
#' vol. 26,7 (2011).
#
#' @export
simulateHubbell <- function(n_species, M, carrying_capacity = 1000,
    k_events = 10, migration_p = 0.02, t_skip = 0,
    t_end, norm = FALSE) {

    # input check
    if (!.isPosInt(n_species)) stop("n_species must be positive integer.")
    if (!.isPosInt(M)) stop("M must be positive integer.")
    if (!.isPosInt(carrying_capacity)) stop("carrying_capacity must be positive integer.")
    if (!.isPosInt(k_events)) stop("k_events must be positive integer.")

    pbirth <- runif(n_species, min = 0, max = 1)
    pmigr <- runif(M, min = 0, max = 1)
    pbirth <- c(pbirth, rep(0, times = (M - n_species)))
    pbirth <- pbirth / sum(pbirth)
    pmigr <- pmigr / sum(pmigr)
    com <- ceiling(carrying_capacity * pbirth)
    if (sum(com) > carrying_capacity) {
        ind <- sample(seq_len(M), size = sum(com) - carrying_capacity, prob = 1 - pbirth)
        com[ind] <- com[ind] - 1
    }

    tseries <- matrix(0, nrow = M, ncol = t_end)
    colnames(tseries) <- paste0("t", seq_len(t_end))
    rownames(tseries) <- seq_len(M)
    com[which(com < 0)] <- 0
    tseries[, 1] <- com
    for (t in seq(2, t_end, 1)) {
        pbirth <- com / sum(com)
        pbirth[which(pbirth < 0)] <- 0
        deaths <- rmultinom(n = 1, size = k_events, prob = pbirth)
        while (sum(com - deaths < 0) > 0) { # species with count 0 have probability
            # 0 and species not present in the community can also not die
            neg_sp <- which(com - deaths < 0)
            pbirth[neg_sp] <- 0
            deaths <- rmultinom(n = 1, size = k_events, prob = pbirth)
        }
        event <- rbinom(k_events, 1, prob = migration_p) # immigration rate migration_p: probability
        # death replaced by immigrant; immigration 1, birth 0
        n_migrants <- sum(event)
        n_births <- length(event) - n_migrants
        births <- rmultinom(1, n_births, prob = pbirth)
        migr <- rmultinom(1, n_migrants, prob = pmigr)
        com <- com - deaths + births + migr
        com[which(com < 0)] <- 0
        tseries[, t] <- com
    }
    if (norm) {
        tseries <- t(t(tseries) / colSums2(tseries))
    }
    counts <- tseries[, seq((t_skip + 1), t_end, 1)]
    # return(counts)
    matrix_out <- cbind(t(counts), time = seq_len(t_end))

    TreeSE <- TreeSummarizedExperiment(
      assays = list(counts = t(matrix_out[, 1:n_species])),
      colData = DataFrame(time = matrix_out[, "time"]),
      metadata = list(carrying_capacity = carrying_capacity,
                      k_events = k_events,
                      migration_p = migration_p,
                      t_skip = t_skip,
                      norm = norm))

    return(TreeSE)
}
