#' Generate random efficiency matrix
#'
#' Generate random efficiency matrix for consumer resource model from Dirichlet
#' distribution, where positive efficiencies indicate the consumption of resources,
#' whilst negatives indicate that the species would produce the resource.
#'
#' @template man_spe
#' @template man_res
#' @param mean_consumption \code{Numeric scalar}. Specifies the mean number 
#' of resources consumed by each species drawn from a poisson distribution
#' (Default: \code{n_resources/4})
#' @param mean_production \code{Numeric scalar}. Specifies the mean number 
#' of resources produced by each species drawn from a poisson distribution
#' (Default: \code{n_resources/6})
#' @param maintenance \code{Numeric scalar}. Specifies the proportion of 
#' resources that cannot be converted into products between 0~1 the 
#' proportion of resources used to maintain the living of microorganisms. 
#' 0 means all the resources will be used for the reproduction of 
#' microorganisms, and 1 means all the resources would be used to maintain 
#' the living of organisms and no resources would be left for their 
#' growth(reproduction). (Default: \code{0.5})
#' @param trophic_levels \code{Integer scalar}. Indicates the number of 
#' species in microbial trophic levels. If NULL, by default, microbial 
#' trophic levels would not be considered. (Default: \code{NULL})
#' @param trophic_preferences \code{List}. Indicates the preferred resources 
#' and productions of each trophic level. Positive values indicate the 
#' consumption of resources, whilst negatives indicate that the species 
#' would produce the resource. (Default: \code{NULL})
#' @param exact \code{Logical scalar}. Whether to set the number of 
#' consumption/production to be exact as mean_consumption/mean_production 
#' or to set them using a Poisson distribution. (Default: \code{FALSE})
#' If `length(trophic_preferences)` is smaller than `length(trophic_levels)`,
#' then NULL values would be appended to lower trophic levels.
#' If NULL, by default, the consumption preference will be defined randomly.
#' (Default: \code{trophic_preferences = NULL})
#'
#' @examples
#' # example with minimum parameters
#' ExampleEfficiencyMatrix <- randomE(n_species = 5, n_resources = 12)
#'
#' # examples with specific parameters
#' ExampleEfficiencyMatrix <- randomE(
#'     n_species = 3, n_resources = 6,
#'     names_species = letters[1:3],
#'     names_resources = paste0("res", LETTERS[1:6]),
#'     mean_consumption = 3, mean_production = 1
#' )
#' ExampleEfficiencyMatrix <- randomE(
#'     n_species = 3, n_resources = 6,
#'     maintenance = 0.4
#' )
#' ExampleEfficiencyMatrix <- randomE(
#'     n_species = 3, n_resources = 6,
#'     mean_consumption = 3, mean_production = 1, maintenance = 0.4
#' )
#'
#' # examples with microbial trophic levels
#' ExampleEfficiencyMatrix <- randomE(
#'     n_species = 10, n_resources = 15,
#'     trophic_levels = c(6, 3, 1),
#'     trophic_preferences = list(
#'         c(rep(1, 5), rep(-1, 5), rep(0, 5)),
#'         c(rep(0, 5), rep(1, 5), rep(-1, 5)),
#'         c(rep(0, 10), rep(1, 5))
#'     )
#' )
#' ExampleEfficiencyMatrix <- randomE(
#'     n_species = 10, n_resources = 15,
#'     trophic_levels = c(6, 3, 1),
#'     trophic_preferences = list(c(rep(1, 5), rep(-1, 5), rep(0, 5)), NULL, NULL)
#' )
#' ExampleEfficiencyMatrix <- randomE(
#'     n_species = 10, n_resources = 15,
#'     trophic_levels = c(6, 3, 1)
#' )
#'
#' @return
#' \code{randomE} returns a matrix E with dimensions (n_species x n_resources),
#' and each row represents a species.
#'
#' @importFrom stats rpois
#' @export
randomE <- function(n_species,
    n_resources,
    names_species = NULL,
    names_resources = NULL,
    mean_consumption = n_resources / 4,
    mean_production = n_resources / 6,
    maintenance = 0.5,
    trophic_levels = NULL,
    trophic_preferences = NULL,
    exact = FALSE) {
    if (!all(vapply(
        list(n_species, n_resources), .isPosInt,
        logical(1)
    ))) {
        stop("n_species and/or n_resources must be integer.")
    }

    # set the default values
    if (is.null(names_species)) {
        names_species <- paste0("sp", seq_len(n_species))
    }
    if (is.null(names_resources)) {
        names_resources <- paste0("res", seq_len(n_resources))
    }
    if (is.null(trophic_levels)) {
        trophic_levels <- n_species
    }
    if (sum(trophic_levels) != n_species) {
        stop("Sum of 'trophic_levels' should equal to 'n_species'.")
    }
    if (!is.null(trophic_preferences)) {
        if (!is.list(trophic_preferences) && length(trophic_preferences) == n_resources) {
            trophic_preferences <- list(trophic_preferences)
        }
        while (length(trophic_preferences) < length(trophic_levels)) {
            warning("Autofilling 'trophic_preferences' with NULL")
            trophic_preferences <- c(trophic_preferences, list(NULL))
        }
    }
    efficiency_matrix <- matrix(0,
        nrow = n_species, ncol = n_resources,
        dimnames = list(names_species, names_resources)
    )

    list_auto_trophic_preference <- list(NULL)
    for (j in seq_len(length(trophic_levels))) {
        n_species_this_level <- trophic_levels[j]

        for (i in seq(n_species_this_level)) {
            irow <- efficiency_matrix[i + sum(trophic_levels[0:(j - 1)]), ]
            consumption <- irow
            production <- irow
            # calculate consumption
            consumption_pref <- trophic_preferences[[j]] * (trophic_preferences[[j]] > 0)
            if (length(consumption_pref) == 0 && is.null(list_auto_trophic_preference[[j]])) {
                # no consumption preference nor auto_trophic_preference
                # consumption_pref <- NULL
                consumption_pref <- rep(1, n_resources)
                if (exact) {
                    index_consumption <- sample(seq(n_resources),
                        size = min(max(1, round(mean_consumption)), n_resources)
                    )
                } else {
                    index_consumption <- sample(seq(n_resources),
                        size = min(max(1, rpois(1, mean_consumption)), n_resources)
                    )
                }
            } else {
                # with consumption preference
                if (length(consumption_pref) == 0) {
                    consumption_pref <- list_auto_trophic_preference[[j]]
                }
                if (exact) {
                    index_consumption <- sample(seq(n_resources),
                        size = min(
                            sum(consumption_pref > 0),
                            max(1, round(mean_consumption))
                        ),
                        replace = FALSE,
                        prob = consumption_pref
                    )
                } else {
                    index_consumption <- sample(seq(n_resources),
                        size = min(
                            sum(consumption_pref > 0),
                            max(1, rpois(1, mean_consumption))
                        ),
                        replace = FALSE,
                        prob = consumption_pref
                    )
                }
            }
            consumption[index_consumption] <- 1
            irow <- rdirichlet(1, consumption * consumption_pref * 100)

            # calculate production
            production_pref <- trophic_preferences[[j]] * (trophic_preferences[[j]] < 0)
            if (sum(production_pref) == 0) { # no production preference
                production_pref <- NULL
                setprod <- setdiff(seq(n_resources), index_consumption)
                if (length(setprod) > 0) {
                    if (exact) {
                        index_production <- unique(
                            sample(setprod,
                                size = round(mean_production),
                                replace = TRUE
                            )
                        )
                    } else {
                        index_production <- unique(
                            sample(setprod,
                                size = rpois(1, mean_production),
                                replace = TRUE
                            )
                        )
                    }
                    index_production <- setdiff(index_production, index_consumption)
                } else {
                    index_production <- c()
                }
            } else { # with production preference
                if (exact) {
                    index_production <- sample(seq(n_resources),
                        size = min(
                            sum(production_pref < 0),
                            round(mean_production)
                        ),
                        replace = FALSE,
                        prob = abs(production_pref)
                    )
                } else {
                    index_production <- sample(seq(n_resources),
                        size = min(
                            sum(production_pref < 0),
                            rpois(1, mean_production)
                        ),
                        replace = FALSE,
                        prob = abs(production_pref)
                    )
                }
            }

            production[index_production] <- 1
            prod <- (-1) * (1 - maintenance) * rdirichlet(1, production)
            irow[index_production] <- prod[index_production]


            efficiency_matrix[i + sum(trophic_levels[0:(j - 1)]), ] <- irow
        }

        # automatically generate consumption of next level according to
        # the production of this level
        if (j < length(trophic_levels)) {
            if (j + 1 > length(list_auto_trophic_preference) || is.null(trophic_preferences[[j + 1]])) {
                eff_mat <- efficiency_matrix[seq_len(n_species_this_level) + sum(trophic_levels[seq(0, j - 1, 1)]), ]
                eff_mat[eff_mat > 0] <- 0
                eff_mat <- -eff_mat
                list_auto_trophic_preference[[j + 1]] <- colSums(eff_mat)
            }
        }
    }
    return(efficiency_matrix)
}
