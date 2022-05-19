#' @title Generate simulation times and the indices of time points to return
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
#' @importFrom utils combn
#' @importFrom stats cov
#' @importFrom stats rbinom
#' @importFrom reshape2 melt
#' @import ggplot2
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

#' @title check whether a number is a positive integer
#' @param x Numeric number to test
#' @param tol Numeric tolerance of detection
#' @return A logical value: whether the number is a positive integer.
#' @export
isPositiveInteger <- function(x, tol = .Machine$double.eps^0.5) {
        return(abs(x - round(x)) < tol && x > 0)
}

#' @title check whether a number is zero or positive integer
#' @param x Numeric number to test
#' @param tol Numeric tolerance of detection
#' @return A logical value: whether the number is zero or positive integer.
#' @export
isZeroOrPositiveInteger <- function(x, tol = .Machine$double.eps^0.5) {
    return(abs(x - round(x)) < tol && x >= 0)
}

#' @title generate a vector of times when events is happening
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

#' @title generate pairs of interactions according to interaction types
#' @description A helper function to be used in combination with 
#' getInteractions()
#' @param I Matrix: defining the interaction between each pair of species
#' @param pair Numeric: a vector with a length of 2, indicating the 2 focusing 
#' species in the process of applying the interaction types
#' @param interType Character: one of 'mutualism', 'commensalism', 'parasitism',
#' 'amensalism', or 'competition'. Defining the interaction type
#' @examples
#' applyInteractionType(matrix(0, nrow = 3, ncol = 3), 
#'                      pair = c(1,2), 
#'                      interType = "parasitism")
#' @return A matrix of interaction types with one pair changed
#' @export
applyInteractionType <- function(I, pair, interType){
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

#' @title generate interactions according to five types of interactions and
#' their weights
#' @param n_species Integer: defining the dimension of matrix of interaction
#' @param weights Numeric: defining the weights of mutualism, commensalism, 
#' parasitism, amensalism, and competition in all interspecies interactions.
#' @param connectance Numeric: defining the density of the interaction network.
#' Ranging from 0 to 1
#' @examples
#' I <- getInteractions(n_species = 10, 
#'                      weights = c(1,1,1,1,1), 
#'                      connectance = 0.5)
#' @return A matrix of interactions with all interactions changed according to
#' the weights and connectance.
#' @export
getInteractions <- function(n_species, weights, connectance){
    I <- matrix(0, n_species, n_species)
    interactions <- c('mutualism', 'commensalism', 'parasitism', 
                      'amensalism', 'competition')
    probs <- abs(weights)/sum(abs(weights))
    combinations <- combn(n_species, 2)
    for (i in sample(seq_along(combinations[1,]),as.integer(ncol(combinations)*connectance), replace=FALSE)){
        I <- applyInteractionType(I, combinations[,i], sample(interactions, 1, prob = probs))
    }
    return (I)
}

#' @title only get the maximum value from a row 
#' @description A helper function to be used in combination with apply()
#' @param aRow Vector: a row in a matrix/dataframe
#' @return A vector of numbers with only the first maximun value kept, and other
#' numbers changed to 0
#' @export
getRowMax <- function(aRow){
    maxPos <- which.max(aRow)
    newRow <- rep(FALSE, length(aRow))
    newRow[maxPos] <- TRUE
    return(aRow * newRow)
}

# plotting functions ####
#' @title simple line plot of species
#' @param out_matrix a result matrix. Rows are timepoints and columns are 
#' species. Last column named "time" indicates the timepoint of each row in the
#' simulation.
#' @param title,obj,y.label title, x, and y axis of the plot
#' @examples
#' ExampleCR <- simulateConsumerResource(n_species = 2, 
#'     n_resources = 4)
#' makePlot(ExampleCR$matrix)
#' @return a ggplot type of plot.
#' @export
makePlot <- function(out_matrix, title = "abundance of species by time", obj = "species", y.label = "x.t"){
    df <- as.data.frame(out_matrix)
    dft <-  melt(df, id="time")
    names(dft)[2] <- obj
    names(dft)[3] <- y.label
    lgd <- ncol(df)<= 20
    ggplot(dft, aes_string(names(dft)[1], names(dft)[3], col = names(dft)[2])) +
        geom_line(show.legend = lgd, lwd=0.5) +
        ggtitle(title) + 
        theme_linedraw() +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
}

#' @title simple line plot of resources
#' @param out_matrix a result matrix of resources. Rows are timepoints and 
#' columns are resources. Last column named "time" indicates the timepoint of 
#' each row in the simulation.
#' @param title title of the plot
#' @return a ggplot type of plot.
#' @examples 
#' ExampleCR <- simulateConsumerResource(n_species = 2, 
#'     n_resources = 4)
#' makePlotRes(ExampleCR$resources)
#' 
#' @return a ggplot type of plot.
#' @export
makePlotRes <- function(out_matrix, title = "quantity of compounds by time"){
    time <- S.t <- resources <- NULL
    
    df <- as.data.frame(out_matrix)
    dft <-  melt(df, id="time")
    names(dft)[2] <- "resources"
    names(dft)[3] <- "S.t"
    lgd <- ncol(df)<= 20
    ggplot(dft, aes(time, S.t, col = resources)) + 
        geom_line(show.legend = lgd, lwd=0.5) + 
        ggtitle(title) + 
        theme_linedraw() + 
        theme(plot.title = element_text(hjust = 0.5, size = 14))
}

#' @title simple heatmap plot of matrix
#' @param matrix.A a matrix. Can be a matrix generated with randomA or randomE
#' @param title,x.label,y.label,midpoint_color,lowColor,midColor,highColor 
#' detailed controls of the heatmap
#' @examples 
#' dense_A <- randomA(n_species = 10, 
#'     scale_off_diagonal = 1, 
#'     diagonal = -1.0, 
#'     connectance = 0.9)
#' makeHeatmap(dense_A, lowColor = 'blue', highColor = 'red')
#' @return a ggplot type of plot.
#' @export
makeHeatmap <-function(matrix.A, 
                       title = "Consumption/production matrix",
                       y.label = 'resources',
                       x.label = 'species',
                       midpoint_color = NULL, 
                       lowColor = "red", 
                       midColor = "white", 
                       highColor = "blue"){
    x <- y <- strength <- NULL
    df <- melt(t(matrix.A))
    if (is.null(midpoint_color)) {
        midpoint_color <- 0
    }
    names(df)<- c("x", "y", "strength")
    df$y <- factor(df$y, levels=rev(unique(sort(df$y))))
    fig <- ggplot(df, aes(x,y,fill=strength)) + geom_tile() + coord_equal() +
        theme(axis.title = element_blank()) + 
        scale_fill_gradient2('strength', low = lowColor, mid = midColor, high = highColor, midpoint = midpoint_color)+
        theme_void() + ggtitle(title)
    
    if (ncol(matrix.A)<=10 & nrow(matrix.A)<=10){
        fig <- fig + geom_text(aes(label = round(strength, 2))) 
    } else if (ncol(matrix.A)<=15 & nrow(matrix.A)<=15){
        fig <- fig + geom_text(aes(label = round(strength, 1)))
    } else {
        fig <- fig
    }
    
    fig <- fig + labs(x = x.label, y = y.label)+
        theme_linedraw() + 
        theme(plot.title = element_text(hjust = 0.5, size = 14), axis.text.x = element_text(
            angle = 90))
    
    if (nrow(matrix.A) >= 20){
        # too many species 
        fig <- fig + theme(
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
        )
    }
    if (ncol(matrix.A) >= 20){
        # too many resources
        fig <- fig + theme(
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()
        )
    }
    fig
}
