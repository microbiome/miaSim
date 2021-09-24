#' Generalized Lotka-Volterra (gLV) simulation
#'
#' Simulates time series with the generalized Lotka-Volterra model and forms a
#' \linkS4class{SummarizedExperiment} object.
#'
#' Simulates a community time series using the generalized Lotka-Volterra model,
#' defined as dx/dt = x(b+Ax), where x is the vector of species abundances,
#' A is the interaction matrix and growth_rates the vector of growth rates.
#'
#' The resulting abundance matrix model is used to construct
#'\linkS4class{SummarizedExperiment} object.
#'
#' @param n.species Integer scalar specifying the number of species
#' @param A  a interaction matrix
#' @param b Numeric scalar indicating growth rates
#' @param x Numeric scalar indicating initial abundances
#' @param tend Integer scalar indicating timepoints
#'      (default: \code{tend = 1000})
#' @param norm Logical scalar \code{TRUE} or \code{FALSE},
#'  return normalised abundances (proportions in each generation)
#'  (default: \code{norm = FALSE})
#'
#' @return
#' \code{simulateGLV} returns a \linkS4class{SummarizedExperiment} object
#'
#' @examples
#' row_data <- data.frame(Kingdom = "Animalia",
#'                 Phylum = rep(c("Chordata", "Echinodermata"), c(500, 500)),
#'                 Class = rep(c("Mammalia", "Asteroidea"), each = 500),
#'                 ASV = paste0("X", seq_len(1000)),
#'                 row.names = rownames(paste0("species", seq_len(1000))),
#'                 stringsAsFactors = FALSE)
#'
#' row_data <- t(row_data)
#'
#' col_data <- DataFrame(sampleID = seq_len(1000),
#'                     time = as.Date(1000, origin = "2000-01-01"),
#'                     row.names = colnames(paste0("sample", seq_len(1000))))
#'
#' SEobject <- simulateGLV(n.species = 4,
#'                         A = powerlawA(n.species = 4, alpha = 2),tend = 1000)
#' rowData(SEobject) <- row_data
#' colData(SEobject) <- col_data
#'
#' @docType methods
#' @aliases simulateGLV-numeric
#' @aliases simulateGLV,numeric-method
#'
#' @importFrom deSolve ode
#' @importFrom stats runif
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("simulateGLV", signature = "n.species",
    function(n.species, A, x = runif(n.species), b = runif(n.species),
            tend = 1000, norm = FALSE)
        standardGeneric("simulateGLV"))

dxdt <- function(t, x, parameters){
    b <- parameters[,1]
    A <- parameters[,2:ncol(parameters)]
    # rate of change
    dx <- x*(b+A %*% x)
    # return rate of change
    list(dx)
}

setMethod("simulateGLV", signature = c(n.species="numeric"),
    function(n.species, A, x = runif(n.species), b = runif(n.species),
            tend = 1000, norm = FALSE){
        parameters <- cbind(b, A)
        times <- seq(0, tend, by = 0.01)

        out <- ode(
                y = x,
                times = times,
                func = dxdt,
                parms = parameters
            )
        spab <- t(out[,2:ncol(out)])
        spab <- spab[,round(seq(1, tend*100, length.out = tend))]
        if(norm){
        spab <- t(t(spab)/colSums(spab))
            }
        spab
        SE <- SummarizedExperiment(assays = list(counts=spab))
        return(SE)
    }
)
