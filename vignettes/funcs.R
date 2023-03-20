# TODO:
# To be replaced by the TreeSummarizedExperiment equivalents.

makePlot <- function(out_matrix, title = "abundance of species by time", obj = "species", y.label = "x.t"){
    df <- as.data.frame(out_matrix)
    dft <-  melt(df, id="time")
    names(dft)[2] = obj
    names(dft)[3] = y.label
    lgd = ncol(df)<= 20
    ggplot(dft, aes_string(names(dft)[1], names(dft)[3], col = names(dft)[2])) +
        geom_line(show.legend = lgd, lwd=0.5) +
        ggtitle(title) +
        theme_linedraw() +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
}



makePlotRes <- function(out_matrix, title = "quantity of compounds by time"){
    df <- as.data.frame(out_matrix)
    dft <-  melt(df, id="time")
    names(dft)[2] = "resources"
    names(dft)[3] = "S.t"
    lgd = ncol(df)<= 20
    ggplot(dft, aes(time, S.t, col = resources)) +
        geom_line(show.legend = lgd, lwd=0.5) +
        ggtitle(title) +
        theme_linedraw() +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
}







makeHeatmap <-function(matrix.A,
                       title = "Consumption/production matrix",
                       y.label = 'resources',
                       x.label = 'species',
                       midpoint_color = NULL,
                       lowColor = "red",
                       midColor = "white",
                       highColor = "blue"){
    df <- melt(t(matrix.A))
    if (is.null(midpoint_color)) {
        midpoint_color <- 0
    }
    names(df)<- c("x", "y", "strength")
    df$y <- factor(df$y, levels=rev(unique(sort(df$y))))
    fig <- ggplot(df, aes(x,y,fill=strength)) + geom_tile() + coord_equal() +
        theme(axis.title = element_blank()) +
        scale_fill_gradient2('strength', low = lowColor,
	    mid = midColor, high = highColor, midpoint = midpoint_color)+
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
        theme(plot.title = element_text(hjust = 0.5, size = 14),
	      axis.text.x = element_text(
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

## Construct a function taking umap_CRM_coor as df and return the mean distance
average_distance <- function(df, res_conc_type, com_type, method = "euclidean"){
    sub_df <- df[df$concentration == res_conc_type & df$community == com_type,]
    combines <- combn(sub_df$medium, 2)
    distances <- NULL
    for (i in seq_len(ncol(combines))) {
        distances[i] <- dist(sub_df[combines[,i], c(1, 2)])
    }
    return(mean(distances))
}



# This function generates a data frame, where each row is arranged in an
# increasing dissimilarity to the first row.
gradient.df.generator <- function(n_row, n_col, density_row, max_gradient, error_interval){
    list_initial <- list()
    dissimilarity.gradient <- seq(from = 0, to = max_gradient, length.out = n_row)
    for (i in seq_len(n_row)){
        print(i)
        if (i == 1){
            row_temp <- rbeta(n_col, 1, 1/n_col)
            col_to_remove <- sample(x = seq_len(n_col), size = n_col-n_col*density_row)
            row_temp[col_to_remove] <- 0
            list_initial[[i]] <- row_temp
        } else {
            while (length(list_initial) < i) {
                row_temp <- rbeta(n_col, 1, 1/n_col)
                col_to_remove <- sample(x = seq_len(n_col), size = n_col-n_col*density_row)
                row_temp[col_to_remove] <- 0
                diff_temp <- abs(vegdist(rbind(list_initial[[1]], row_temp), method = "bray") - dissimilarity.gradient[i])
                if (diff_temp < error_interval) {
                    list_initial[[i]] <- row_temp
                }
            }
        }
    }
    dataframe_to_return <- as.data.frame(t(matrix(unlist(list_initial), ncol = n_row)))
    return(dataframe_to_return)
}