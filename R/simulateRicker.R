#' @title Generate time series with the Ricker model
#'
#' @description The Ricker model is a discrete version of the generalized
#' Lotka-Volterra model and is implemented here as proposed by Fisher and Mehta
#' in PLoS ONE 2014. The function returns a
#' \linkS4class{SummarizedExperiment} object.
#'
#' @param n.species Integer: number of species
#' @param A interaction matrix
#' @param K Numeric: carrying capacities
#' @param x Numeric: initial abundances
#' @param sigma noise level, if set to a non-positive value, no noise is added
#' (default: \code{sigma = 0.05})
#' @param explosion.bound boundary for explosion
#' (default: \code{explosion.bound = 10^8})
#' @param norm Logical scalar: returns normalised abundances (proportions
#' in each generation) (default: \code{norm = FALSE})
#' @param tskip  Integer: number of generations that should
#' not be included in the outputted species abundance matrix.
#' (default: \code{tskip = 0})
#' @param tend Integer: number of simulations to be simulated
#' @return
#' \code{simulateRicker} returns a \linkS4class{SummarizedExperiment}
#' class object containing matrix with species abundance as rows and time
#' points as columns
#' @references Fisher & Mehta (2014). Identifying Keystone Species in the Human
#' Gut Microbiome from Metagenomic Timeseries using Sparse Linear Regression.
#' PLoS One 9:e102451
#' @examples
#' A <- powerlawA(10, alpha = 1.01)
#' ExampleRicker <- simulateRicker(n.species=10,A,tend=100)
#'
#' @importFrom stats rlnorm
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
simulateRicker <- function(n.species, A, x = runif(n.species),
                        K = runif(n.species),sigma=0.05,
                        explosion.bound=10^8, tskip = 0, tend, norm = FALSE){
        if(length(x) != n.species){
            stop("x needs to have n.species entries.")
        }
        if(nrow(A)!=n.species || ncol(A)!=n.species){
            stop("A needs to have n.species rows and n.species columns.")
        }
        if(length(K)!=n.species){
            stop("K needs to have n.species entries.")
        }
        tseries <- matrix(nrow=n.species, ncol=tend-tskip)
        tseries[,1] <- x
        # simulate difference equation
        for(t in seq(2:tend)){
            if(sigma > 0){
                b=rlnorm(n.species,meanlog=0,sdlog=sigma)
            }else{
                b=rep(1,n.species)
            }

            x <- b*x*exp(A%*%(x-K))

            if(max(x) > explosion.bound){
                # report which species explodes
                stop("Explosion for species ", which(x==max(x)))
            }
            if(length(x[x<0]) > 0){
                stop("Species below 0!")
            }
            if(t > tskip){
                tseries[,t-tskip] <- x
            }
        }
        if(norm){
            tseries <- t(t(tseries)/colSums(tseries))
        }
        SE <- SummarizedExperiment(assays = list(counts = tseries))
        return(SE)
}
