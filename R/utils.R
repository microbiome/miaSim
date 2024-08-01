#' Generate simulation times and the indices of time points to return
#' in simulation functions.
#'
#' @param t_start \code{Numeric scalar}. Indicates the initial time of the simulation.
#' (Default: \code{0})
#' @param t_end \code{Numeric scalar}. Indicates the final time of the simulation
#' (Default: \code{1000})
#' @param t_step \code{Numeric scalar}. Indicates the interval between simulation steps
#' (Default: \code{0.1})
#' @param t_store \code{Integer scalar}. Indicates the number of evenly distributed
#' time points to keep (Default: \code{100})
#'
#' @return lists containing simulation times (t_sys) and the indices to keep.
#' @examples
#' Time <- .simulationTimes(
#'     t_start = 0, t_end = 100, t_step = 0.5,
#'     t_store = 100
#' )
#' DefaultTime <- .simulationTimes(t_end = 1000)
#'
#' @importFrom stats cov
#' @importFrom stats rbinom
#'
#' @keywords internal
#' @export
.simulationTimes <- function(t_start = 0, t_end = 1000,
    t_step = 0.1, t_store = 1000) {
    t_total <- t_end - t_start
    t_sys <- seq(t_start, t_end, by = t_step)
    t_index <- seq(1, length(t_sys) - 1, by = floor(length(t_sys) / t_store))
    return(list("t_sys" = t_sys, "t_index" = t_index[seq_len(t_store)]))
}

#' Check whether a number is a positive integer
#' @param x Numeric number to test
#' @param tol Numeric tolerance of detection
#' @return A logical value: whether the number is a positive integer.
.isPosInt <- function(x, tol = .Machine$double.eps^0.5) {
    return(abs(x - round(x)) < tol && x > 0)
}

#' Generate dirichlet random deviates
#'
#' @param n Number of random vectors to generate.
#' @param alpha Vector containing shape parameters.
#'
#' @importFrom stats rgamma
#'
#' @return a vector containing the Dirichlet density
#'
#' @examples
#' dirichletExample <- rdirichlet(1, c(1, 2, 3))
#'
#' @export
rdirichlet <- function(n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x / as.vector(sm)
}

#' Generate a vector of event times
#'
#' @param t_events Numeric vector; starting time of the events
#' @param t_duration Numeric vector; duration of the events
#' @param t_end Numeric: end time of the simulation
#' @param ... : additional parameters to pass to simulationTimes, including
#' t_start, t_step, and t_store.
#' @return A vector of time points in the simulation
#' @examples
#' tEvent <- simulateEventTimes(
#'     t_events = c(10, 50, 100),
#'     t_duration = c(1, 2, 3),
#'     t_end = 100,
#'     t_store = 100,
#'     t_step = 1
#' )
#' @export
simulateEventTimes <- function(t_events = NULL, t_duration = NULL,
    t_end = 1000, ...) {
    tdyn <- .simulationTimes(t_end = t_end, ...)
    t_result <- c()
    for (i in seq(length(t_events))) {
        p1 <- tdyn$t_sys[(tdyn$t_sys >= t_events[i]) &
            (tdyn$t_sys < (t_events[i] + t_duration[i]))]
        t_result <- c(t_result, p1)
    }
    return(t_result)
}

#' Generate pairs of interactions according to interaction types
#'
#' A helper function to be used in combination with .getInteractions()
#' @param I Matrix: defining the interaction between each pair of species
#' @param pair Numeric: a vector with a length of 2, indicating the 2 focusing
#' species in the process of applying the interaction types
#' @param interType Character: one of 'mutualism', 'commensalism', 'parasitism',
#' 'amensalism', or 'competition'. Defining the interaction type
#' @return A matrix of interaction types with one pair changed
.applyInterType <- function(I, pair, interType) {
    if (rbinom(1, 1, 0.5)) {
        pair <- rev(pair)
    }
    if (interType == "mutualism") {
        I[pair[1], pair[2]] <- 1
        I[pair[2], pair[1]] <- 1
        return(I)
    } else if (interType == "commensalism") {
        I[pair[1], pair[2]] <- 1
        I[pair[2], pair[1]] <- 0
        return(I)
    } else if (interType == "parasitism") {
        I[pair[1], pair[2]] <- 1
        I[pair[2], pair[1]] <- -1
        return(I)
    } else if (interType == "amensalism") {
        I[pair[1], pair[2]] <- 0
        I[pair[2], pair[1]] <- -1
        return(I)
    } else if (interType == "competition") {
        I[pair[1], pair[2]] <- -1
        I[pair[2], pair[1]] <- -1
        return(I)
    }
}

#' Generate interactions according to five types of interactions and their
#' weights
#' @param n_species Integer: defining the dimension of matrix of interaction
#' @param weights Numeric: defining the weights of mutualism, commensalism,
#' parasitism, amensalism, and competition in all interspecies interactions.
#' @param connectance Numeric: defining the density of the interaction network.
#' Ranging from 0 to 1
#' @return A matrix of interactions with all interactions changed according to
#' the weights and connectance.
.getInteractions <- function(n_species, weights, connectance) {
    I <- matrix(0, n_species, n_species)
    interactions <- c(
        "mutualism", "commensalism", "parasitism",
        "amensalism", "competition"
    )
    probs <- abs(weights) / sum(abs(weights))
    combinations <- utils::combn(n_species, 2)
    for (i in sample(seq_along(combinations[1, ]), as.integer(ncol(combinations) * connectance), replace = FALSE)) {
        I <- .applyInterType(I, combinations[, i], sample(interactions, 1, prob = probs))
    }
    return(I)
}

#' Replace one element with zero in a list
#'
#' If the list contains m elements, then lengths of each element must be m, too.
#' This function is intended to generate a list of x0 (the initial community)
#' with one missing species, to prepare the parameter simulations_compare in
#' `estimateAFromSimulations`.
#' @param input_list A list containing m elements, and lengths of each element
#' must be m, too.
#' @return A list of same dimension as input_list, but with 0 at specific
#' positions in the elements of the list.
.replaceByZero <- function(input_list) { # params_iter$x0 as input_list
    if (!all(length(input_list) == unlist(unique(lapply(input_list, length))))) {
        stop("Length of input_list doesn't match length of element in it.")
    }
    for (i in seq_along(input_list)) {
        input_list[[i]][[i]] <- 0
    }
    return(input_list)
}

#' Get the interspecies interaction matrix A using leave-one-out method
#'
#' generate matrix A from the comparisons between simulations with one absent
#' species and a simulation with complete species (leave-one-out)
#'
#' @param simulations A list of simulation(s) with complete species
#' @param simulations2 A list of simulation(s), each with one absent
#' species
#' @param n_instances Integer: number of instances to generate
#' (default: \code{n_instances = 1})
#' @param t_end Numeric: end time of the simulation. If not identical with t_end
#' in params_list, then it will overwrite t_end in each simulation
#' (default: \code{t_end = 1000})
#' @param scale_off_diagonal Numeric: scale of the off-diagonal elements
#' compared to the diagonal. Same to the parameter in function `randomA`.
#' (default: \code{scale_off_diagonal = 0.1})
#' @param diagonal Values defining the strength of self-interactions. Input can
#' be a number (will be applied to all species) or a vector of length n_species.
#' Positive self-interaction values lead to exponential growth. Same to the
#' parameter in function `randomA`.
#' (default: \code{diagonal = -0.5})
#' @param connectance Numeric frequency of inter-species interactions.
#' i.e. proportion of non-zero off-diagonal terms. Should be in the interval
#'  0 <= connectance <= 1. Same to the parameter in function `randomA`.
#' (default: \code{connectance = 0.2})

#' @return a matrix A with dimensions (n_species x n_species) where n_species
#' equals to the number of elements in simulations2
#' @importFrom SummarizedExperiment assay
.estimateAFromSimulations <- function(simulations,
    simulations2,
    n_instances = 1,
    t_end = NULL,
    scale_off_diagonal = 0.1,
    diagonal = -0.5,
    connectance = 0.2) {

    # Use last time point if t_end is not given
    if (is.null(t_end)) {t_end <- which.max(simulations[[1]]$time)}

    # simulations_means should be a vector of n_species length, indicating the original average simulation result
    # simulations should be generated using generateSimulations with param_iter = NULL
    simulations_means <- t(as.matrix(rowMeans(sapply(simulations, function (x) {assay(x, "counts")[, which.max(x$time)]}))))

    # simulations_compare_means should be a dataframe with a dimension of n_species X n_species
    # and each row indicates the result of simulation with one species absent
    simulations2_matrix <- t(sapply(simulations2, function (x) {assay(x, "counts")[, which.max(x$time)]}))

    simulations_compare_means <- matrix(
        data = NA,
        nrow =  nrow(simulations2_matrix) / n_instances,
        ncol =  ncol(simulations2_matrix)
    )

    colnames(simulations_compare_means) <- colnames(simulations2_matrix)
    for (i in seq_len(nrow(simulations_compare_means))) {
        simulations_compare_means[i, ] <- simulations2_matrix[(i - 1) * n_instances + seq_len(n_instances), ]
    }
    if (nrow(simulations_compare_means) != length(simulations_means)) {
        stop("Number of simulations to compare should be equal to the numbers of species!")
    }

    # normalize and substitute
    simulations_compare_means_norm <- simulations_compare_means / rowSums(simulations_compare_means)
    simulations_effect <- simulations_compare_means_norm
    for (i in seq_len(nrow(simulations_compare_means_norm))) {
        simulations_means_i <- simulations_means
        simulations_means_i[1, i] <- 0
        simulations_means_norm_i <- simulations_means_i / rowSums(simulations_means_i)
        simulations_effect[i, ] <- simulations_compare_means_norm[i, ] - simulations_means_norm_i
    }

    matrixA <- as.matrix(simulations_effect) * scale_off_diagonal
    value_threshold <- sort(abs(matrixA))[max(1, length(as.vector(matrixA)) * (1 - connectance))]
    matrixA <- matrixA * (abs(matrixA) >= value_threshold)
    diag(matrixA) <- diagonal
    return(matrixA)
}


