#' Lowest-RMSE sample from time series
#'
#' Computes the Root Mean Square Error between subsequent generations (columns)
#' in given species abundances matrix (time series) that is in
#' \linkS4class{SummarizedExperiment} object as a measure of convergence.
#'
#' The timepoint with the lowest RMSE-difference with respect to the previous
#' timepoint is then extracted and a normalised sample is
#' returned (proportions).
#'
#' @param SE a \linkS4class{SummarizedExperiment} object that has species
#' abundances matrix with OTUs in the rows and the timepoints as columns.
#' @param warn logical: (default: \code{warn = TRUE}), set to FALSE
#' to suppress convergence warning
#' @param norm logical: (default: \code{warn = TRUE}) compute RMSE on
#' compositional time series (abundances per time point sum to 1)
#' @param cutoff numeric: (default: \code{cutoff = 1e-04}) The value the minimum
#' RMSE cannot exceed in order to declare convergence.
#'
#' @return A sample vector with the length equal to the number of rows of
#' species abundances matrix.
#'
#' @examples
#' A <- miaSim::powerlawA(n.species = 10, alpha = 1.2)
#' SE <- miaSim::simulateGLV(n.species = 10, A, tend = 100)
#' sample <- rmse_sample(SE)
#'
#' @export

setGeneric("rmse_sample",signature = "SE",
           function(SE, warn = TRUE, norm = TRUE, cutoff = 1e-04)
             standardGeneric("rmse_sample"))

setMethod("rmse_sample", signature = c(SE="SummarizedExperiment"),
            function(SE, warn = TRUE, norm = TRUE, cutoff = 1e-04){
            spab <- assay(SE)
            if(norm){
                rmse_vec <- rmse_t(spab = t(t(spab)/colSums(spab)))
            } else {
                rmse_vec <- rmse_t(spab)
            }
            if(warn){
                if(min(rmse_vec) >= (cutoff)){
                warning(paste0("WARNING -  simulation did not converge.",
                    "Run the simulation longer or change parameters."))
                }
            }
                sample <- spab[,which.min(rmse_vec)]
                return(sample)
            })

            # Helper function to compute rmse vector over all timepoints of
            # a given species abundances x time points matrix
            rmse_t <- function(spab){
                N <- nrow(spab)
                K <- ncol(spab)

                rmse_vec <- rep.int(0, times = (K-1)) # initialize
                rmse_vec <- sapply(X = 1:(K-1), FUN = function(j){
                rmse_vec[j] = rmse(spab[,j], spab[,(j+1)])
            })
            return(rmse_vec)
            }
)
