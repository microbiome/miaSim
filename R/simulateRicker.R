#' Generate time series with the Ricker model
#'
#' The Ricker model is a discrete version of the generalized
#' Lotka-Volterra model and is implemented here as proposed by Fisher and Mehta
#' in PLoS ONE 2014.
#'
#' @param n_species integer number of species
#' @param A interaction matrix
#' @param K numeric carrying capacities
#' @param x numeric initial abundances
#' @param sigma numeric value of noise level, if set to a non-positive value,
#' no noise is added (default: \code{sigma = 0.05})
#' @param explosion_bound numeric value of boundary for explosion
#' (default: \code{explosion_bound = 10^8})
#' @param tskip integer number of generations that should not be included in the
#' outputted species abundance matrix (default: \code{tskip = 0})
#' @param tend integer number of simulations to be simulated
#' @param norm logical scalar returning normalised abundances (proportions
#' in each generation) (default: \code{norm = FALSE})
#'
#' @seealso
#' \code{\link[miaSim:convertToSE]{convertToSE}}
#'
#' @return
#' \code{simulateRicker} returns an abundance matrix with species abundance
#' as rows and time points as columns
#'
#' @references Fisher & Mehta (2014). Identifying Keystone Species in the Human
#' Gut Microbiome from Metagenomic Timeseries using Sparse Linear Regression.
#' PLoS One 9:e102451
#'
#' @examples
#' A <- powerlawA(10, alpha = 1.01)
#' ExampleRicker <- simulateRicker(n_species=10,A,tend=100)
#'
#' @importFrom stats rlnorm
#'
#' @export
simulateRicker <- function(n_species, A, x = runif(n_species),
                        K = runif(n_species),sigma=0.05,
                        explosion_bound=10^8, tskip = 0, tend, norm = FALSE){
        if(length(x) != n_species){
            stop("x needs to have n_species entries.")
        }
        if(nrow(A)!=n_species || ncol(A)!=n_species){
            stop("A needs to have n_species rows and n_species columns.")
        }
        if(length(K)!=n_species){
            stop("K needs to have n_species entries.")
        }
        counts <- matrix(nrow=n_species, ncol=tend-tskip)
        counts[,1] <- x
        # simulate difference equation
        for(t in seq(from=2, to=tend)){
            if(sigma > 0){
                b=rlnorm(n_species,meanlog=0,sdlog=sigma)
            }else{
                b=rep(1,n_species)
            }

            x <- b*x*exp(A%*%(x-K))

            if(max(x) > explosion_bound){
                # report which species explodes
                stop("Explosion for species ", which(x==max(x)))
            }
            if(length(x[x<0]) > 0){
                stop("Species below 0!")
            }
            if(t > tskip){
                counts[,t-tskip] <- x
            }
        }
        if(norm){
            counts <- t(t(counts)/colSums2(counts))
        }

        return(counts)
}
