#' Generate random uniform interaction matrix
#'
#' Generate random simplified interaction matrix from a uniform distribution.
#'
#' @param n.species integer number of species
#' @param d numeric diagonal values (should be negative)
#' (default: \code{d = -0.5})
#' @param min.strength numeric minimal off-diagonal interaction strength
#' (default: \code{min.strength = -0.5})
#' @param max.strength numeric maximal off-diagonal interaction strength
#' (default: \code{max.strength = 0.5})
#' @param connectance numeric (interaction probability)
#' (default: \code{connectance = 0.02})
#' @param symmetric logical return a symmetric interaction matrix
#' (default: \code{symmetric = FALSE})
#'
#' @examples
#' high_inter_A <- randomA(n.species = 10, d = -0.4, min.strength = -0.8,
#'                                     max.strength = 0.8, connectance = 0.5)
#'
#' low_inter_A <- randomA(n.species = 10, connectance = 0.01)
#'
#' @return
#' \code{randomA} returns a matrix A with dimensions (n.species x n.species)
#'
#' @export
randomA <- function(n.species, d = -0.5, min.strength = -0.5,
                    max.strength = 0.5, connectance = 0.02, symmetric = FALSE){
            A = runif(n.species^2, min = min.strength, max = max.strength)

            #an efficient approximation to reach the desired connectance
            setZeros = n.species^2*(1-connectance)
            A[sample(seq_len(length(A)), setZeros, replace=0)] = 0

            A = matrix(A,
                        nrow = n.species,
                        ncol = n.species
            )

            if(symmetric){
                A[lower.tri(A)] <- t(A)[lower.tri(A)]
            }

            diag(A) <- d
            return(A)
}
