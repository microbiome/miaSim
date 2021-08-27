#' simulateGLV
#'
#' Simulates time series with the generalized Lotka-Volterra model and forms a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]
#' {SummarizedExperiment}} object.
#'
#' Simulates a community time series using the generalized Lotka-Volterra model,
#' defined as dx/dt = x(b+Ax), where x is the vector of species abundances,
#' A is the interaction matrix and growth_rates the vector of growth rates.
#'
#' The resulting abundance matrix model is used to construct
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#'
#' @param N Integer scalar specifying the number of species
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
#' \code{simulateGLV} returns a \code{\link{SummarizedExperiment}} object
#'
#' @examples
#' result <- simulateGLV(N = 4, A = powerlawA(n = 4, alpha = 2), tend = 1000)
#'
#' @importFrom deSolve ode
#' @importFrom stats runif
#' @importFrom microsim methods
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export

setGeneric("simulateGLV",signature = "N",
    function(N, A, b = runif(N), x = runif(N), tend = 1000, norm = FALSE)
    standardGeneric("simulateGLV"))

.dxdt <- function(t, x, parameters){
        b <- parameters[,1]
        A <- parameters[,2:ncol(parameters)]

        dx <- x*(b+A %*% x)
        list(dx)
}

row_data <- data.frame(Kingdom = "A",
        Phylum = rep(c("B1", "B2"), c(500, 500)),
        Class = rep(c("C1", "C2"), each = 500),
        ASV = paste0("D", 1:1000),
        row.names = rownames(paste0("species", seq_len(1000))),
        stringsAsFactors = FALSE)

row_data <- t(row_data)

col_data <- data.frame(sampleID = c(1:1000),
        time = as.Date(sample( as.numeric(as.Date('2000-01-01')):
        as.numeric(as.Date('2014-01-01')), 1000,
        replace = T),
        origin = '1970-01-01'),
        row.names = colnames(paste0("sample", 1:1000)),
        stringsAsFactors = FALSE)

setMethod("simulateGLV", signature = c(N="numeric"),
    function(N, A, b = runif(N), x = runif(N), tend = 1000, norm = FALSE){
        parameters <- cbind(b, A)
        times <- seq(0, tend, by = 0.01)

        .dxdt(t, x, parameters)

        out <- ode(
                y = x,
                times = times,
                func = .dxdt,
                parms = parameters
            )
        spab <- t(out[,2:ncol(out)])
        spab <- spab[,round(seq(1, tend*100, length.out = tend))]
        if(norm){
        spab <- t(t(spab)/colSums(spab))
            }
        spab
        SE <- SummarizedExperiment(assays = spab,
                            rowData = row_data,
                            colData = col_data)

        return(SE)
    }
)
