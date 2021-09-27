#' Generalized Lotka-Volterra (gLV) simulation
#'
#' Simulates time series with the generalized Lotka-Volterra model and forms a
#' \linkS4class{SummarizedExperiment} object.
#'
#' Simulates a community time series using the generalized Lotka-Volterra model,
#' defined as dx/dt = diag(x)(b+Ax), where x is the vector of species abundances
#' ,diag(x) is a diagonal matrix with the diagonal values set to x.
#' A is the interaction matrix and b is the vector of growth rates.
#'
#' The resulting abundance matrix model is used to construct
#'\linkS4class{SummarizedExperiment} object.
#'
#' @param n.species Integer: number of species
#' @param A interaction matrix
#' @param b Numeric: growth rates
#' @param x Numeric: initial abundances
#' @param norm Logical scalar: returns normalised abundances (proportions
#' in each generation) (default: \code{norm = FALSE})
#' @param ... additional arguments that can be called from miaSim::tDyn
#'
#' @return
#' \code{simulateGLV} returns a \linkS4class{SummarizedExperiment} object
#' containing abundance matrix
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
#' A <- miaSim::powerlawA(4, alpha = 1.01)
#'
#' SEobject <- simulateGLV(n.species = 4, A, t.start = 0, t.store = 1000)
#' rowData(SEobject) <- row_data
#' colData(SEobject) <- col_data
#'
#' @docType methods
#' @aliases simulateGLV-numeric
#' @aliases simulateGLV,numeric-method
#'
#' @importFrom utils str
#' @importFrom deSolve ode
#' @importFrom stats runif
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("simulateGLV", signature = "n.species",
    function(n.species, A, x = runif(n.species),
            b = runif(n.species), t.start = 0, t.store, norm = FALSE, ...)
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
    function(n.species, A, x = runif(n.species),
                b = runif(n.species), t.start = 0, t.store, norm = FALSE, ...){
        parameters <- cbind(b, A)
        t.dyn <- tDyn(t.start, ..., t.store)

        out <- ode(
                y = x,
                times = t.dyn$t.sys,
                func = dxdt,
                parms = parameters
        )
        spab <- t(out[,2:ncol(out)])
        spab <- spab[,t.dyn$t.index]
        print(str(spab))
        if(norm){
        spab <- t(t(spab)/colSums(spab))
            }
        spab
        SE <- SummarizedExperiment(assays = list(counts=spab))
        return(SE)
    }
)
