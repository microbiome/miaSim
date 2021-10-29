#' Generate random efficiency matrix
#'
#' Generate random efficiency matrix for consumer resource model from dirichlet
#' distribution. Positive efficiencies indicate the consumption of resources,
#' whilst negatives indicate that the species would produce the resource.
#'
#' @param n.species Integer: number of species
#' @param n.resources Integer: number of resources
#' @param min.con Integer: minimum number of resources consumed by each species
#' @param max.con Integer: maximum number of resources consumed by each species
#' @param min.prod Integer: minimum number of resources produced by each species
#' @param max.prod Integer: maximum number of resources produced by each species
#' @param maintenance Numeric: 0~1 the proportion of resources used to maintain
#' the living of microorganisms. 0 means all the resources will be used for the
#' reproduction of microorganisms, and 1 means all the resources would be used
#' to maintain the living of organisms and no resources would be left for
#' their growth(reproduction).
#'
#' @examples
#' # example with specific parameters
#' ExampleEfficiencyMatrix <- randomE(n.species = 3, n.resources = 6,
#' min.con = 3, max.con = 4, min.prod = 1, max.prod = 1, maintenance = 0.4)
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
                    maintenance = 0.5){
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
        irow <- efficiency.matrix[i,]
        consumption <- irow
        production <- irow
        index.consumption <- sample(seq(n.resources),
            size = sample(min.con:max.con))
        consumption[index.consumption] <- 1
        irow <- rdirichlet(1, consumption)[,]
        index.production <- sample(seq(n.resources)[irow==0], size=sample(min.prod:max.prod))
        production[index.production] <- 1
        prod <- (-1)*(1-maintenance)* rdirichlet(1, production)[,]
        irow[index.production] <- prod[index.production]
        efficiency.matrix[i,] <- irow
    }

    return(efficiency.matrix)
}
