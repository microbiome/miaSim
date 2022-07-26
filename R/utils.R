#' Generate simulation times and the indices of time points to return
#' in simulation functions.
#'
#' @param t_start Numeric scalar indicating the initial time of the simulation.
#' (default: \code{t_start = 0})
#' @param t_end Numeric scalar indicating the final time of the dimulation
#' (default: \code{t_end = 1000})
#' @param t_step Numeric scalar indicating the interval between simulation steps
#' (default: \code{t_step = 0.1})
#' @param t_store Integer scalar indicating the number of evenly distributed
#' time points to keep (default: \code{t_store = 100})
#'
#' @return lists containing simulation times (t_sys) and the indices to keep.
#' @examples
#' Time <- simulationTimes(t_start = 0, t_end = 100, t_step = 0.5,
#'     t_store = 100)
#' DefaultTime <- simulationTimes(t_end = 1000)
#'
#' @importFrom stats cov
#' @importFrom stats rbinom
#'
#' @keywords internal
#' @export
simulationTimes <- function(t_start = 0, t_end = 1000,
            t_step = 0.1, t_store = 1000){
        t_total <- t_end-t_start
        t_sys <- seq(t_start, t_end, by = t_step)
    t_index <- seq(1, length(t_sys)-1, by=floor(length(t_sys)/t_store))
    return(list("t_sys" = t_sys, "t_index" = t_index[seq_len(t_store)]))
}

#' check whether a number is a positive integer
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
#' dirichletExample <- rdirichlet(1, c(1,2,3))
#'
#' @export
rdirichlet <- function(n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x / as.vector(sm)
}

#' generate a vector of times when events is happening
#' @param t_events,t_duration Numeric: vector of starting time and durations of
#' the events
#' @param t_end Numeric: end time of the simulation
#' @param ... : additional parameters to pass to simulationTimes, including
#' t_start, t_step, and t_store.
#' @return A vector of time points in the simulation
#' @examples
#' tEvent <-  eventTimes(t_events = c(10, 50, 100),
#'     t_duration = c(1, 2, 3),
#'     t_end = 100,
#'     t_store = 100,
#'     t_step = 1)
#' @export
eventTimes <- function(t_events = NULL, t_duration = NULL,
                       t_end=1000, ...){
        tdyn <- simulationTimes(t_end = t_end,...)
        t_result <- c()
        for (i in seq(length(t_events))){
            p1 <- tdyn$t_sys[(tdyn$t_sys >= t_events[i]) &
            (tdyn$t_sys < (t_events[i]+t_duration[i]))]
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
.applyInterType <- function(I, pair, interType){
  if (rbinom(1,1,0.5)){
    pair <- rev(pair)
  }
  if (interType=='mutualism'){
    I[pair[1],pair[2]] <- 1
    I[pair[2],pair[1]] <- 1
    return(I)
  }else if (interType=='commensalism'){
    I[pair[1],pair[2]] <- 1
    I[pair[2],pair[1]] <- 0
    return(I)
  }else if (interType=='parasitism'){
    I[pair[1],pair[2]] <- 1
    I[pair[2],pair[1]] <- -1
    return(I)
  }else if (interType=='amensalism'){
    I[pair[1],pair[2]] <- 0
    I[pair[2],pair[1]] <- -1
    return(I)
  }else if (interType=='competition'){
    I[pair[1],pair[2]] <- -1
    I[pair[2],pair[1]] <- -1
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
.getInteractions <- function(n_species, weights, connectance){
    I <- matrix(0, n_species, n_species)
    interactions <- c('mutualism', 'commensalism', 'parasitism',
                      'amensalism', 'competition')
    probs <- abs(weights)/sum(abs(weights))
    combinations <- utils::combn(n_species, 2)
    for (i in sample(seq_along(combinations[1,]),as.integer(ncol(combinations)*connectance), replace=FALSE)){
        I <- .applyInterType(I, combinations[,i], sample(interactions, 1, prob = probs))
    }
    return (I)
}


#' Generate multiple simulations and store them in a list
#'
#' This function is useful when generating simulations with different parameters
#' or simulating multiple instances to evaluate stochasticity therein.
#'
#' @param model Character: name of the model to use, must be one of
#' "simulateConsumerResource", "simulateGLV", "simulateHubbellRates", and
#' "simulateStochasticLogistic"
#' @param params_list List: a list containing all parameters used in the
#' aforementionned called model
#' @param param_iter List of parameters to overwrite in different simulations.
#' If NULL, then simulate different instances using the same set of parameters.
#' (default: \code{param_iter = NULL})
#' @param n_instances Integer: number of instances to generate
#' (default: \code{n_instances = 1})
#' @param t_end Numeric: end time of the simulation. If not identical with t_end
#' in params_list, then it will overwrite t_end in each simulation
#' (default: \code{t_end = 1000})
#' @return A list containing multiple simulating results
#' @examples
#'
#' crm_params <- list(n_species = 10,
#'                    n_resources = 5,
#'                    E = randomE(n_species = 10, n_resources = 5,
#'                        mean_consumption = 1, mean_production = 3),
#'                    x0 = rep(0.001, 10),
#'                    resources = rep(1000,5),
#'                    monod_constant = matrix(rbeta(10*5, 10,10),nrow=10, ncol=5),
#'                    inflow_rate = .5,
#'                    outflow_rate = .5,
#'                    migration_p = 0,
#'                    stochastic = TRUE,
#'                    t_start = 0,
#'                    t_end = 20,
#'                    t_store = 100,
#'                    growth_rates = runif(10),
#'                    norm=FALSE)
#' CRMSimus <- generateSimulations(model = "simulateConsumerResource",
#'     params_list = crm_params, param_iter = NULL, n_instances = 1, t_end = 20)
#'
#' @export
generateSimulations <- function(model,
                                params_list,
                                param_iter = NULL,
                                n_instances = 1,
                                t_end = 1000){
    simulations <- list()
    if (params_list$t_end != t_end){
        warning("t_end in the variable list not equals to t_end")
        params_list$t_end <- t_end
    }
    if (is.null(param_iter)){
        for (i in seq_len(n_instances)){
            # print(paste(i, "of", n_instances,"instances in one set of params."))
            simulation <- do.call(model, params_list)
            # simulation$matrix[,colnames(simulation$model)!= "time"]
            simulations[[length(simulations)+1]] <- simulation
        }
        return(simulations)
    } else {
        if (!all(names(param_iter) %in% names(params_list))){
            stop("not recognized parameter(s) :",
                 names(param_iter)[!names(param_iter) %in% names(params_list)]
            )
        }
        if (length(unique(unlist(lapply(param_iter, length)))) > 1) {
            stop("parameters in list param_iter are not of the same length.")
        }

        simulations_full <- list()
        for (i in seq_along(param_iter[[1]])){
            # print(paste(i, "of", length(param_iter[[1]]), "sets of params."))
            params_list_local <- utils::modifyList(params_list,
                                                   lapply(param_iter, "[[", i))
            simulations_local <- generateSimulations(
                model,
                params_list_local,
                param_iter = NULL,
                n_instances = n_instances,
                t_end = t_end)
            simulations_full <- append(simulations_full, simulations_local)
        }
        return(simulations_full)
    }
}

#' Get the community from one simulation
#'
#' Get the community composition at specific time point in a simulation
#' (default: last time point)
#'
#' @param simulation A list generated from simulateXXX function
#' @param t_end Numeric: end time of the simulation. If NULL, the last t_end in
#' the simulation will be used to extract the final state.
#' (default: \code{t_end = NULL})
#' @return a subset of dataframe simulation$matrix, containing the composition
#' at t_end.
#' @examples
#' # single use: extract community composition at t_end = 200 from a simulation
#' getCommunity(simulateGLV(n_species = 3, t_end = 1000), t_end = 200)
#' # apply this function to previously generated list of simulations
#' GLVSimus <- generateSimulations("simulateGLV", list(n_species = 3, t_end = 1000), n_instances = 5)
#' GLVSimusCom <- as.data.frame(do.call(rbind, lapply(GLVSimus, getCommunity)))
#'
#' @export
getCommunity <- function(simulation, t_end = NULL){
    if (is.null(t_end)){
        t_end <- max(simulation$matrix[,"time"])
    }
    return(simulation$matrix[simulation$matrix[,"time"] == t_end, -ncol(simulation$matrix)])
}
#' Replace one element with zero in a list.
#'
#' If the list contains m elements, then lengths of each element must be m, too.
#' This function is intended to generate a list of x0 (the initial community)
#' with one missing species, to prepare the parameter simulations_compare in
#' `estimateAFromSimulations`.
#' @param input_list A list containing m elements, and lengths of each element
#' must be m, too.
#' @return A list of same dimension as input_list, but with 0 at specific
#' positions in the elements of the list.
#'
.replaceByZero <- function(input_list){ # crm_params_iter$x0 as input_list
    if(!all(length(input_list) == unlist(unique(lapply(input_list, length))))){
        stop("Length of input_list doesn't match length of element in it.")
    }
    for (i in seq_along(input_list)){
        input_list[[i]][[i]] <- 0
    }
    return(input_list)
}
#' Create a list of parameter.
#'
#' This function is intended to generate a list of parameters such as initial
#' community, to prepare the parameter param_iter in `generateSimulations`.
#' @param input_param Vector: the input parameter for a simulation.
#' @param n_repeat Integer: the number to repeat the parameter.
#' @param replace_by_zero Boolean: whether to replace certain elements in result with
#' 0 using internal function `.replaceByZero`.
#' @return A list of parameters which is ready as input of param_iter for
#' `generateSimulations`.
#' @examples
#' paramx0 <- createParamList(input_param = rep(99, 10), n_repeat = 10, replace_by_zero = TRUE)
#' @export
createParamList <- function(input_param, n_repeat, replace_by_zero = FALSE){
    res_list <- vector(mode = "list", length = n_repeat)
    for (i in seq_len(n_repeat)){
        res_list[[i]] <- input_param
    }
    if (replace_by_zero) res_list <- .replaceByZero(res_list)
    return(res_list)
}

#' Get the interspecies interaction matrix A using leave-one-out method
#'
#' generate matrix A from the comparisons between simulations with one absent
#' species and a simulation with complete species (leave-one-out)
#'
#' @param simulations A list of simulation(s) with complete species
#' @param simulations_compare A list of simulation(s), each with one absent
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
#' @examples
#' # example of generateSimulations
#' crm_params <- list(n_species = 10,
#'                    n_resources = 5,
#'                    E = randomE(n_species = 10, n_resources = 5,
#'                        mean_consumption = 1, mean_production = 3),
#'                    x0 = rep(0.001, 10),
#'                    resources = rep(1000,5),
#'                    monod_constant = matrix(rbeta(10*5, 10,10),nrow=10, ncol=5),
#'                    inflow_rate = .5,
#'                    outflow_rate = .5,
#'                    migration_p = 0,
#'                    stochastic = TRUE,
#'                    t_start = 0,
#'                    t_end = 20,
#'                    t_store = 100,
#'                    growth_rates = runif(10),
#'                    norm=FALSE)
#' CRMSimus <- generateSimulations(model = "simulateConsumerResource",
#'     params_list = crm_params, param_iter = NULL, n_instances = 1, t_end = 20)
#' CRMSimus_SE <- convertToSE(assay = t(CRMSimus[[1]]$matrix[,1:10]),
#'                            colData = DataFrame(time = CRMSimus[[1]]$matrix[,"time"]))
#' miaViz::plotSeries(CRMSimus_SE, x = "time")
#'
#' # get average of all instances
#' CRMSimusCom <- as.data.frame(do.call(rbind, lapply(CRMSimus, getCommunity)))
#' CRMSimusMeans <- colMeans(CRMSimusCom)
#' CRMSimusVariance <- apply(CRMSimusCom, 2, var)
#'
#'
#' # test leave-one-out in CRM
#'
#' .replaceByZero <- function(input_list){ # crm_params_iter$x0 as input_list
#'     if(!all(length(input_list) == unlist(unique(lapply(input_list, length))))){
#'         stop("Length of input_list doesn't match length of element in it.")
#'     }
#'     for (i in seq_along(input_list)){
#'         input_list[[i]][[i]] <- 0
#'     }
#'     return(input_list)
#' }
#' createParamList <- function(input_param, n_repeat, replace_by_zero = FALSE){
#'     res_list <- vector(mode = "list", length = n_repeat)
#'     for (i in seq_len(n_repeat)){
#'         res_list[[i]] <- input_param
#'     }
#'     if (replace_by_zero) res_list <- .replaceByZero(res_list)
#'     return(res_list)
#' }
#'
#' paramx0 <- createParamList(input_param = rep(0.001, 10), n_repeat = 10, replace_by_zero = TRUE)
#' paramresources <- createParamList(input_param = rep(1000,5), n_repeat = 10)
#' # test overwrite params
#' crm_params_iter <- list(x0 = paramx0, resources = paramresources)
#' CRMSimus2 <- generateSimulations(model = "simulateConsumerResource",
#'     params_list = crm_params, param_iter = crm_params_iter, n_instances = 1, t_end = 20)
#' # get average of all instances
#' CRMSimusCom2 <- as.data.frame(do.call(rbind, lapply(CRMSimus2, getCommunity)))
#'
#' # get only one average for 1 set of params
#' n_instances <- 1
#' CRMSimusCom2simp <- matrix(NA,
#'                            nrow = nrow(CRMSimusCom2)/n_instances,
#'                            ncol = ncol(CRMSimusCom2))
#' colnames(CRMSimusCom2simp) <- colnames(CRMSimusCom2)
#' for (i in seq_len(nrow(CRMSimusCom2simp))){
#'     CRMSimusCom2simp[i,] <- colMeans(CRMSimusCom2[(i-1)*n_instances + (1:n_instances),])
#' }
#' # View(CRMSimusCom2simp)
#' estimatedA <- estimateAFromSimulations(CRMSimus, CRMSimus2, n_instances = 1,
#'     scale_off_diagonal = 1, diagonal = -0.5, connectance = 0.2)/1000
#'
#' estimatedGLVmodel <- simulateGLV(n_species = 10, x0 = crm_params$x0,
#'     A = estimatedA, growth_rates = crm_params$growth_rates, t_end = 20, t_store = 100)
#' estimatedGLVmodel_SE <- convertToSE(assay = t(estimatedGLVmodel$matrix[,1:10]),
#'                                     colData = DataFrame(time = estimatedGLVmodel$matrix[, "time"]))
#' miaViz::plotSeries(estimatedGLVmodel_SE, x = "time")
#' @return a matrix A with dimensions (n_species x n_species) where n_species
#' equals to the number of elements in simulations_compare
#' @export
estimateAFromSimulations <- function(
        simulations,
        simulations_compare,
        n_instances = 1,
        t_end = NULL,
        scale_off_diagonal = 0.1,
        diagonal = -0.5,
        connectance = 0.2
){
    # simulations_means should be a vector of n_species length, indicating the original average simulation result
    # simulations should be generated using generateSimulations with param_iter = NULL
    simulations_means <- colMeans(as.data.frame(do.call(rbind, lapply(simulations, getCommunity, t_end = t_end))))
    simulations_means <- t(as.data.frame(simulations_means))
    # simulations_compare_means should be a dataframe with a dimension of n_species X n_species
    # and each row indicates the result of simulation with one species absent
    simulations_compare_df <- as.data.frame(do.call(rbind, lapply(simulations_compare, getCommunity, t_end = t_end)))
    simulations_compare_means <- matrix(data = NA,
                                        nrow = nrow(simulations_compare_df)/n_instances,
                                        ncol = ncol(simulations_compare_df))
    colnames(simulations_compare_means) <- colnames(simulations_compare_df)
    for (i in seq_len(nrow(simulations_compare_means))){
        simulations_compare_means[i,] <- colMeans(simulations_compare_df[(i-1)*n_instances + seq_len(n_instances),])
    }
    if (nrow(simulations_compare_means) != length(simulations_means)){
        stop("numbers of simulations to compare not equals to the numbers of species")
    }
    simulations_compare_means <- as.data.frame(simulations_compare_means)

    # normalize and substrate
    simulations_compare_means_norm <- simulations_compare_means/rowSums(simulations_compare_means)
    simulations_effect <- simulations_compare_means_norm
    for (i in seq_len(nrow(simulations_compare_means_norm))){
        simulations_means_i <- simulations_means
        simulations_means_i[1,i] <- 0
        simulations_means_norm_i <- simulations_means_i/rowSums(simulations_means_i)
        simulations_effect[i,] <- simulations_compare_means_norm[i,] - simulations_means_norm_i
    }

    matrixA <- as.matrix(simulations_effect) * scale_off_diagonal
    value_threshold <- sort(abs(matrixA))[max(1, length(as.vector(matrixA)) * (1 - connectance))]
    matrixA <- matrixA * (abs(matrixA)>=value_threshold)
    diag(matrixA) <- diagonal
    return(matrixA)
}
