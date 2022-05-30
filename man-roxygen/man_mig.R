#' @param migration_p Numeric: the probability/frequency of migration from a
#' metacommunity.
#' (default: \code{migration_p = 0.01})
#' @param metacommunity_probability Numeric: Normalized probability distribution
#' of the likelihood that species from the metacommunity can enter the community
#' during the simulation. If NULL, `rdirichlet(1, alpha = rep(1,n_species))` is
#' used.
#' (default: \code{metacommunity_probability = NULL})
