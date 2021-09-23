#' @title Generate random uniform interaction matrix
#' @description Generate random simplified interaction matrix from a uniform distribution.
#'
#' @param n.species number of species
#' @param d diagonal values, indicating self-interactions (use negative values for stability)
#' @param min.strength minimal off-diagonal interaction strength
#' @param max.strength maximal off-diagonal interaction strength
#' @param connectance (interaction probability)
#' @examples
#' high_inter_A <- randomA(10, d = -0.4, min.strength = -0.8, max.strength = 0.8, c = 0.5)
#' low_inter_A <- randomA(10, c = 0.01)
#' @export
#'

randomA <- function(
  n.species = 100,
  d = -0.5,
  min.strength = -0.5,
  max.strength = 0.5,
  connectance = 0.02 
){
  A = runif(n.species^2, min = min.strength, max = max.strength)
  setZeros = n.species^2*(1-connectance) #an efficient approximation to reach the desired connectance 
  A[sample(1:length(A), setZeros, replace=0)] = 0
  
  A = matrix(A,
             nrow = n.species,
             ncol = n.species
  )
  diag(A) <- d
  
  return(A)
}
