#' Convert simulated microbiome data to \linkS4class{TreeSummarizedExperiment}(`TreeSE`) data container
#'
#' This function converts simulated microbial community data to
#' \linkS4class{TreeSummarizedExperiment}(`TreeSummarizedExperiment`) format.
#'
#' The abundance matrix from the simulation functions
#' in `miaSim` can be converted to \linkS4class{TreeSummarizedExperiment}
#' class object.
#'
#' Storing the data in \linkS4class{TreeSummarizedExperiment} enables access to
#' various Bioconductor packages and tools that extend the
#' this class. This offers data and metadata
#' synchronization, while accommodating specialized data structures for
#' particular scientific applications.
#'
#' Examples for `TreeSE` object manipulation and analysis can be found at
#' \url{https://microbiome.github.io/OMA}
#'
#' @param assay is a matrix-like or list of matrix-like object.
#' Rows refer to taxa and columns refer to samples.
#' @param output character value for storing the matrix in
#' \linkS4class{TreeSummarizedExperiment} (\code{output = TreeSE}).
#' @param ... : additional parameters to pass
#'
#' @examples
#' # Simulate time series data
#' x <- simulateHubbellRates(
#'     n_species = 3,
#'     migration_p = 0.01,
#'     metacommunity_probability = NULL,
#'     k_events = 1,
#'     growth_rates = NULL,
#'     norm = FALSE,
#'     t_end = 1000
#' )
#' # Convert into TreeSE format.
#' HubbellSE <- TreeSummarizedExperiment(
#'     assays = list(counts = t(x$matrix[, 1:3])),
#'     colData = DataFrame(time = x$matrix[, "time"]),
#'     metadata = x[-which(names(x) == "matrix")]
#' )
#' miaViz::plotSeries(HubbellSE, x = "time")
#' @return \linkS4class{TreeSummarizedExperiment} data object containing abundance matrix
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @export
convertToTreeSE <- function(assay, output, ...) {
    .Deprecated(msg = "The convertToTreeSE is replaced with TreeSummarizedExperiment.")

    TreeSummarizedExperiment(
        assays = list(counts = assay), ...
    )
}

#' Old conversion function (to be deprecated)
#' @param assay is a matrix-like or list of matrix-like object.
#' Rows refer to taxa and columns refer to samples.
#' @param output character value for storing the matrix in
#' \linkS4class{TreeSummarizedExperiment} (\code{output = TreeSE}).
#' @param ... : additional parameters to pass
#' @return SummarizedExperiment data object containing abundance matrix
#' @export
convertToSE <- function(assay, output, ...) {
    .Deprecated(msg = "The convertToSE function is replaced with TreeSummarizedExperiment.")
}
