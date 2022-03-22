#' Generate random efficiency matrix
#'
#' Generate random efficiency matrix for consumer resource model from Dirichlet
#' distribution. Positive efficiencies indicate the consumption of resources,
#' whilst negatives indicate that the species would produce the resource.
#'
#' @param n_species integer number of species
#' @param n_resources integer number of resources
#' @param min_con integer minimum number of resources consumed by each species
#' @param max_con integer maximum number of resources consumed by each species
#' @param min_prod integer minimum number of resources produced by each species
#' @param max_prod integer maximum number of resources produced by each species
#' @param maintenance numeric value between 0~1 the proportion of resources used
#' to maintain the living of microorganisms. 0 means all the resources will be
#' used for the reproduction of microorganisms, and 1 means all the resources
#' would be used to maintain the living of organisms and no resources would be
#' left for their growth(reproduction).
#'
#' @examples
#' # example with minimum parameters
#' ExampleEfficiencyMatrix2 <- randomE(n_species = 5, n_resources = 12)
#'
#' @return
#' \code{randomE} returns a matrix E with dimensions (n_species x n_resources),
#' and each row represents a species.
#'
#' @export
randomE <- function(n_species,
                    n_resources,
                    min_con = round(n_resources/4),
                    max_con = round(n_resources/3),
                    min_prod = round(n_resources/6),
                    max_prod = round(n_resources/4),
                    maintenance = 0.5){

    if(!all(vapply(list(n_species, n_resources), isPositiveInteger,
            logical(1)))){
        stop("n_species and/or n_resources must be integer.")}

    if(min_con > max_con) {
        warning("min_con surpassed max_con, modified to max_con")
        min_con <- max_con
    }
    if(min_prod > max_prod) {
        warning("min_prod surpassed max_prod, modified to max_prod")
        min_prod <- max_prod
    }
    if(min_con + min_prod > n_resources) {
        stop("min_con + min_prod surpassed n_resources.")
    }
    efficiency_matrix <- matrix(0, nrow = n_species, ncol = n_resources)

    for (i in seq(n_species)){
        irow <- efficiency_matrix[i,]
        consumption <- irow
        production <- irow
        index_consumption <- sample(seq(n_resources),
            size = min_con)
        consumption[index_consumption] <- 1
        irow <- rdirichlet(1, consumption)[,]
        index_production <-
            sample(seq(n_resources)[irow==0], size= min_prod)
        production[index_production] <- 1
        prod <- (-1)*(1-maintenance)* rdirichlet(1, production)[,]
        irow[index_production] <- prod[index_production]
        efficiency_matrix[i,] <- irow
        E <- efficiency_matrix
    }

    return(E)
}
