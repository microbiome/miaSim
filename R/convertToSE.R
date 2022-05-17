#' \linkS4class{SummarizedExperiment}(`SE`) or
#' `TreeSE` construction function
#'
#' The abundance matrix from the simulation functions
#' in `miaSim` can be converted to \linkS4class{SummarizedExperiment}
#' class object.
#'
#' Storing the data in \linkS4class{SummarizedExperiment} enables access to
#' various Bioconductor packages and tools that extend the
#' \linkS4class{SummarizedExperiment} class.
#' This offers data and metadata
#' synchronization, while accommodating specialized data structures for
#' particular scientific applications.
#'
#' Further examples for `SE` object manipulation and analysis can be found at
#' \url{https://microbiome.github.io/OMA}
#'
#' @param assay is a matrix-like or list of matrix-like object.
#' Rows refer to taxa and columns refer to samples.
#' @param output character value for storing the matrix in either
#' \linkS4class{TreeSummarizedExperiment} (\code{output = TSE}) or
#' \linkS4class{SummarizedExperiment} (default: \code{output = SE})
#'
#' @param ... : additional parameters to pass for
#' the `SE` object.
#'
#' @examples
#' x <- simulateHubbellRates(
#'     community_initial = c(0,5,10), migration_p = 0.01,
#'     metacommunity_p = NULL, k_events = 1, growth_rates = NULL, norm = FALSE,
#'     t_end=1000)
#'
#' HubbellSE <- convertToSE(matrix = x$counts,
#'                          colData = x$time,
#'                          metadata = x$metadata)
#'
#' @return \linkS4class{SummarizedExperiment} an object containing abundance
#' matrix
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
#' @export
convertToSE<- function(assay, output, ...){

  if(missing(output)) {
    objectContainer <- SummarizedExperiment(assays = list(counts = assay),
                                            ...)
  } else {
    objectContainer <- TreeSummarizedExperiment(
      assays = list(counts = assay),
      ...)
  }

  return(objectContainer)
}