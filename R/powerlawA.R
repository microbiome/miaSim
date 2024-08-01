#' Interaction matrix with Power-Law network adjacency matrix
#'
#' N is the an Interspecific Interaction matrix with values drawn from
#' a normal distribution H the interaction strength heterogeneity drawn from
#' a power-law distribution with the parameter alpha, and G the adjacency matrix
#' of with out-degree that reflects the heterogeneity of the powerlaw.
#' A scaling factor s may be used to constrain the values of the interaction
#' matrix to be within a desired range. Diagonal elements of A are defined
#' by the parameter d.
#'
#' @param n_species \code{Integer scalar}. Indicates the number of species.
#' @param alpha \code{Numeric scalar}. Specifies the power-law distribution. 
#' Should be > 1. Larger values will give lower interaction
#' strength heterogeneity, whereas values closer to 1 give strong heterogeneity
#' in interaction strengths between the species. In other words, values of alpha
#' close to 1 will give Strongly Interacting Species (SIS). (Default: \code{3.0}) 
#' @param stdev \code{Numeric scalar}. Specifies the standard deviation of the normal
#' distribution with mean 0 from which the elements of the nominal interspecific
#' interaction matrix N are drawn. (Default: \code{1})
#' @param s \code{Numeric scalar}. Specifies the scaling with which the final global
#' interaction matrix A is multiplied. (Default: \code{0.1})
#' @param d \code{Numeric scalar}. Diagonal values, indicating self-interactions (use
#' negative values for stability). (Default: \code{1.0})
#' @param symmetric \code{Logical scalar}. Whether a symmetric interaction matrix
#' is returned. (Default: \code{FALSE})
#'
#' @return The interaction matrix A with dimensions (n_species x n_species)
#'
#' @references Gibson TE, Bashan A, Cao HT, Weiss ST, Liu YY (2016)
#' On the Origins and Control of Community Types in the Human Microbiome.
#' PLOS Computational Biology 12(2): e1004688.
#' https://doi.org/10.1371/journal.pcbi.1004688
#'
#' @importFrom stats rnorm
#' @importFrom poweRlaw rplcon
#'
#' @examples
#' # Low interaction heterogeneity
#' A_low <- powerlawA(n_species = 10, alpha = 3)
#' # Strong interaction heterogeneity
#' A_strong <- powerlawA(n_species = 10, alpha = 1.01)
#'
#' @export
powerlawA <- function(n_species, alpha = 3.0, stdev = 1, s = 0.1, d = -1,
    symmetric = FALSE) {

    # input check
    if (!.isPosInt(n_species)) {
        stop("n_species must be integer.")
    }
    if (!all(vapply(
        list(alpha, stdev, s, d),
        is.numeric, logical(1)
    ))) {
        stop("alpha, stdev, s and d values must be numeric.")
    }

    # Nominal Interspecific Interaction matrix N
    N <- matrix(
        data = rnorm(n_species^2, mean = 0, sd = stdev),
        nrow = n_species,
        ncol = n_species
    )

    # power law sample
    pl <- rplcon(n = n_species, xmin = 1, alpha = alpha)
    pl[is.infinite(pl)] <- 10^308

    # Interaction strength heterogeneity
    H <- diag(1 + (pl - min(pl)) / (max(pl) - min(pl)))

    # Adjacency matrix G of power-law out-degree digraph ecological
    # network
    deg <- 0.1 * n_species

    h <- pmin(ceiling(deg * pl / mean(pl)), n_species)

    G <- matrix(0, nrow = n_species, ncol = n_species)
    for (i in seq_len(n_species)) {
        index <- sample(x = seq_len(n_species), size = h[i])
        G[index, i] <- 1
    }

    # G[t(G) == 1] <- 1
    A <- N %*% H * G
    A <- A * s / max(A)

    if (symmetric) {
        A[lower.tri(A)] <- t(A)[lower.tri(A)]
    }

    diag(A) <- d
    colnames(A) <- seq_len(n_species)
    rownames(A) <- seq_len(n_species)
    return(A)
}
