#' simulateGLV
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
#' \code{simulateGLV} returns a \linkS4class{SummarizedExperiment} object
#'
#' @examples
#' row_data <- data.frame(Kingdom = "A",
#'                     Phylum = rep(c("B1", "B2"), c(500, 500)),
#'                     Class = rep(c("C1", "C2"), each = 500),
#'                     ASV = paste0("D", seq_len(1000)),
#'                     row.names = rownames(paste0("species", seq_len(1000))),
#'                     stringsAsFactors = FALSE)
#'                     row_data <- t(row_data)
#'
#' col_data <- data.frame(sampleID = seq_len(1000),
#'                     time = as.Date(sample(as.numeric(as.Date('2000-01-01')):
#'                                 as.numeric(as.Date('2014-01-01')), 1000,
#'                                             replace = TRUE),
#'                                                     origin = '1970-01-01'),
#'                     row.names = colnames(paste0("sample", seq_len(1000))),
#'                     stringsAsFactors = FALSE)
#'
#' simulateGLV(N = 4, A = powerlawA(n = 4, alpha = 2), tend = 1000)
#'
#' @importFrom deSolve ode
#' @importFrom stats runif
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom poweRlaw rplcon
#'
#' @export


setGeneric("simulateGLV",signature = "N",
    function(N, A, b = runif(N), x = runif(N), tend = 1000, norm = FALSE)
    standardGeneric("simulateGLV"))

#' @export

powerlawA <- function(
    n, # number of species
    alpha, # power-law distribution parameter
    stdev = 1, # sd normal distribution
    s = 0.1 # scaling parameter, default: 0.1/max(A)
){
    # Nominal Interspecific Interaction matrix N
    N <- matrix(
        data = rnorm(n^2, mean = 0, sd = stdev),
        nrow = n,
        ncol = n
    )
    # power law sample
    pl <- rplcon(n = n, xmin = 1, alpha = alpha)
    # Interaction strength heterogeneity H
    H <- sapply(seq_len(n), FUN = function(i){
        1 + ((pl[i]-min(pl))/(max(pl)-min(pl)))
    })
    H <- diag(H)
    # Adjacency matrix G of power-law out-degree digraph ecological network
    d <- 0.1*n
    h <- sapply(seq_len(n), FUN = function(i){
        min(
            ceiling(d*pl[i]/mean(pl)),
            n
        )
    })
    G <- matrix(0, nrow = n, ncol = n)
    for(i in seq_len(n)){
        index <- sample(x = seq_len(n), size = h[i])
        G[index, i] <- 1
    }
    A <- N %*% H * G
    A <- A*s/max(A)
    diag(A) <- -1
    colnames(A) <- seq_len(n)
    rownames(A) <- seq_len(n)
    return(A)
}

#' @export

dxdt <- function(t, x, parameters){
        b <- parameters[,1]
        A <- parameters[,2:ncol(parameters)]

        dx <- x*(b+A %*% x)
        list(dx)
}

setMethod("simulateGLV", signature = c(N="numeric"),
    function(N, A, b = runif(N), x = runif(N), tend = 1000, norm = FALSE){
        parameters <- cbind(b, A)
        times <- seq(0, tend, by = 0.01)

        dxdt(t, x, parameters)
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
        SE <- SummarizedExperiment(assays = list(counts=spab),
                                            colData=col_data,
                                            rowData=row_data)
        return(SE)
    }
)
