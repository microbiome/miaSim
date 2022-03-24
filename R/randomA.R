#' Generate random uniform interaction matrix
#'
#' Generate random simplified interaction matrix from a uniform distribution.
#'
#' @param n_species integer number of species
#' @param d numeric diagonal values (should be negative)
#' (default: \code{d = -0.5})
#' @param min_strength numeric value of minimal off-diagonal interaction
#' strength (default: \code{min_strength = -0.5})
#' @param max_strength numeric value of maximal off-diagonal interaction
#' strength (default: \code{max_strength = 0.5})
#' @param connectance numeric interaction probability
#' (default: \code{connectance = 0.02})
#' @param symmetric logical scalar returning a symmetric interaction matrix
#' (default: \code{symmetric=FALSE})
#'
#' @examples
#' high_inter_A <- randomA(n_species = 10, d = -0.4, min_strength = -0.8,
#'                                     max_strength = 0.8, connectance = 0.5)
#'
#' low_inter_A <- randomA(n_species = 10, connectance = 0.01)
#'
#' @return
#' \code{randomA} returns a matrix A with dimensions (n_species x n_species)
#'
#' @export
randomA <- function(n_species, d = -0.5, min_strength = -0.5,
                    max_strength = 0.5,connectance = 0.02, symmetric = FALSE){
            A = runif(n_species^2, min = min_strength, max = max_strength)

            #input check
            if(!isPositiveInteger(n_species)){
                stop("n_species must be integer.")}
            if(!all(vapply(list(d, min_strength, max_strength, connectance),
                    is.numeric, logical(1)))){
                stop("d, min_strength, max_strength and connectance values
                    must be numeric.")}

            #an efficient approximation to reach the desired connectance
            setZeros = n_species^2*(1-connectance)
            A[sample(seq_len(length(A)), setZeros, replace=0)] = 0

            A = matrix(A,
                        nrow = n_species,
                        ncol = n_species
            )

            if(symmetric){
                A[lower.tri(A)] <- t(A)[lower.tri(A)]
            }

            diag(A) <- d
            return(A)
}
