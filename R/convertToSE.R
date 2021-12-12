#' \linkS4class{SummarizedExperiment}(`SE`) or
#' `TreeSE` construction function
#'
#' Storing the data in \linkS4class{SummarizedExperiment} enables access to
#' various tools for further analysis of data. A large number of Bioconductor
#' packages contain extension of \linkS4class{SummarizedExperiment} class.
#' \linkS4class{SummarizedExperiment} class offers data and metadata
#' synchronization, while still accommodating specialized data structures for
#' particular scientific applications.
#'
#' Further examples for `SE` object manipulation and analysis can be found at
#' https://microbiome.github.io/OMA/
#'
#' The resulting abundance matrix from the simulation functions used
#' in `miaSim` can be easily converted to \linkS4class{SummarizedExperiment}
#' class object.
#'
#' @param matrix is a matrix-like or list of matrix-like object.
#' Rows refer to taxa and columns refer to samples.
#' @param output character value for storing the matrix in either
#' \linkS4class{TreeSummarizedExperiment} (\code{output = TSE}) or
#' \linkS4class{SummarizedExperiment} (default: \code{output = SE})
#'
#' @param ... : additional parameters that can be implemented in
#' the `SE` object.
#'
#' @examples
#' ExampleHubbellRates <- simulateHubbellRates(
#'     community_initial = c(0,5,10), migration_p = 0.01,
#'     metacommunity_p = NULL, k_events = 1, growth_rates = NULL, norm = FALSE,
#'     t.end=1000)
#'
#' HubbellSE <- convertToSE(matrix = ExampleHubbellRates$counts,
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
convertToSE<- function(matrix, output, ...){

    if(missing(output)) {
        objectContainer <- SummarizedExperiment(assays = list(counts = matrix),
                                                ...)
    } else {
        objectContainer <- TreeSummarizedExperiment(
                                    assays = list(counts = matrix),
                                    ...)
    }

    return(objectContainer)
}
