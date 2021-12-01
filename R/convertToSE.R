#' It helps to store any matrix in \linkS4class{SummarizedExperiment}
#' object class data container(`SE`). The resulting abundance matrix from the
#' simulation functions used in `miaSim` can be used.
#'
#' @param matrix : A matrix-like or list of matrix-like objects.
#' Rows: genes, genomic coordinates, etc. Columns: samples, cells, etc.
#' @param tse : Logical, decides if the object is stored in
#' \linkS4class{TreeSummarizedExperiment} or \linkS4class{SummarizedExperiment}
#' @param ... : additional objects that can be implemented in
#' the `SE` object.
#'
#' @examples
#' ExampleHubbellRates <- miaSim::simulateHubbellRates(
#'     community_initial = c(0,5,10), migration_p = 0.01,
#'     metacommunity_p = NULL, k_events = 1, growth_rates = NULL, norm = FALSE,
#'     t.end=1000)
#'
#' HubbellSE <- convertToSE(matrix = ExampleHubbellRates$abundancematrix,
#'                         colData = ExampleHubbellRates$time,
#'                         metadata = ExampleHubbellRates$metadata)
#'
#' @return \linkS4class{SummarizedExperiment} an object containing abundance
#' matrix
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
#' @export
convertToSE<- function(matrix, tse = FALSE, ...){

    if(tse){
        objectContainer <- TreeSummarizedExperiment(
                                    assays = list(counts = matrix),
                                    ...)
    }

    objectContainer <- SummarizedExperiment(assays = list(counts = matrix),
                                            ...)
    return(objectContainer)
}
