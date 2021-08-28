#' Microbiome Phyloseq wrapper for community samples
#'
#' Takes in samples in vector or matrix format
#' and first, merges them to a microbiome phyloseq class object.The returning
#' phyloseq object is used to construct
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#'
#' @param S vector or matrix : Each parameter given to the function is a sample vector or matrix to be merged
#' into one phyloseq object of samples.
#' The vectors/matrices must all contain the same number of species in order to be merged.
#' A phyloseq object requires at least a matrix of min. 2 columns,
#' therefore an error is thrown if only a vector is given without additional vector/matrix objects.
#' @param ... vector(s)/matrice(s): see description of S
#' @return A \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object with species/OTUs as rows and samples as columns.
#' @details Requires phyloseq package (Bioconductor). Function will install package when not installed yet.
#' @examples
#'
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' \dontrun{
#' A = powerlawA(n = 80, alpha = 1.6),
#' glv_timeseries <- glv(N = 80, A = A, tend = 100)
#' soi_sample <- rmse_sample(soi(N = 80, I = 1000, A = A, tend = 1000))
#' physeq <- asPhyloseq(glv_timeseries, soi_sample)
#' }
#' @rdname simulatePhyloseq
#' @export
#'
setGeneric("simulatePhyloseq",signature = c(V = "vector", M = "matrix"),
           function(S, ...)
             standardGeneric("simulatePhyloseq"))

setMethod("simulatePhyloseq", signature = c(V = "vector", M = "matrix"),
          function(S, ...){
            param=c(as.list(environment()), list(...))
            names(param)=1:length(param)
            merge = as.matrix(param[[1]])
            if(length(param)>1){
              for(i in 2:length(param)){
                param[[i]] <- as.matrix(param[[i]])
                if(nrow(merge) != nrow(param[[i]])){
                  stop("Samples must contain the same number of species.")
                }
                if(ncol(param[[i]]) == sum(colSums(param[[i]]))){
                  stop("Samples should consist of absolute count data (no relative abundances)")
                }
                merge = cbind(merge, param[[i]])
              }
            } else {
              if(ncol(param[[1]])==1){
                stop("Need more samples.")
              }
            }
            rownames(merge) <- paste0("OTU", 1:nrow(merge))
            colnames(merge) <- paste0("Sample", 1:ncol(merge))

            physeq = phyloseq(otu_table(merge, taxa_are_rows = TRUE))
            tree <- phy_tree(physeq)
            TreeSE <- TreeSummarizedExperiment(assays = list(counts = M),
                                               colData = ...,
                                               rowData = ...,
                                               rowTree = tree
                                               )
          }
)
