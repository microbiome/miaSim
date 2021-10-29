#' Generate random efficiency matrix
#'
#' Generate random efficiency matrix for consumer resource model from a normal
#' distribution. Positive efficiencies indicate the consumption of resources,
#' whilst negatives indicate that the species would produce the resource.
#' Efficiencies larger than 1 would be replaced by 1.
#'
#' @param n.species Integer: number of species
#' @param n.resources Integer: number of resources
#' @param min.con Integer: minimum number of resources consumed by each species
#' @param max.con Integer: maximum number of resources consumed by each species
#' @param min.prod Integer: minimum number of resources produced by each species
#' @param max.prod Integer: maximum number of resources produced by each species
#' @param sd Numeric: standard deviation of the normal distribution
#'
#' @examples
#' # example with specific parameters
#' ExampleEfficiencyMatrix <- randomE(n.species = 3, n.resources = 6,
#' min.con = 3, max.con = 4, min.prod = 1, max.prod = 1, sd = 0.4)
#' # example with minimum parameters
#' ExampleEfficiencyMatrix2 <- randomE(n.species = 5, n.resources = 12)
#'
#' @return
#' \code{randomE} returns a matrix E with dimensions (n.species x n.resources),
#' and each row represents a species.
#'
#' @export
randomE <- function(n.species,
                    n.resources,
                    min.con = round(n.resources/4),
                    max.con = round(n.resources/3),
                    min.prod = round(n.resources/6),
                    max.prod = round(n.resources/4),
                    sd = 0.5){
    if(min.con > max.con) {
        warning("min.con surpassed max.con, modified to max.con")
        min.con <- max.con
    }
    if(min.prod > max.prod) {
        warning("min.prod surpassed max.prod, modified to max.prod")
        min.prod <- max.prod
    }
    if(min.con + min.prod > n.resources) {
        stop("min.con + min.prod surpassed n.resources.")
    }
    efficiency.matrix <- matrix(0, nrow = n.species, ncol = n.resources)
    for (i in seq(n.species)){
        irow <- abs(rnorm(n.resources,sd = sd))
        irow[irow >1] <- 1
        # a larger sd might cause exponential growth in consumer-resource model
        # efficiencies larger than 1 is not meaningful.
        index.zeros <- sample(seq(n.resources),
                              size = sample(seq(n.resources - max.con - max.prod,
                                                n.resources-min.con-min.prod)))
        irow[index.zeros] <- 0
        index.neg <- sample(seq(n.resources)[irow>0],
                            size = sample(seq(min.prod, max.prod)))
        irow[index.neg] <- -irow[index.neg]
        efficiency.matrix[i,] <- irow
    }
    return(efficiency.matrix)
}
