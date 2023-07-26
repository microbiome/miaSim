## ---- echo=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(
  cache = FALSE,
  fig.width = 9,
  message = FALSE,
  warning = FALSE)


## -----------------------------------------------------------------------------
library(ggplot2)
library(miaSim)
library(vegan)
library(ggplot2)
library(reshape2)


## -----------------------------------------------------------------------------
set.seed(42)


## -----------------------------------------------------------------------------
n_rep <- 50


## -----------------------------------------------------------------------------
result_df <- data.frame(
    n_species = integer(),
    theta = numeric(),
    i = integer(),
    n_resources = integer(),
    value = numeric()
)

result_df2 <- data.frame(
  matrix(NA, nrow = 0, ncol = 13, dimnames = list(NULL, paste0("sp", seq_len(13))))
)

sorensen_df <- data.frame(
    n_species = integer(),
    theta = numeric(),
    rho_mean = numeric(),
    rho_sd = numeric()
)



## -----------------------------------------------------------------------------
gradient_df_generator <- function(n_row, n_col, density_row, max_gradient, error_interval){
  list_initial <- list()
  dissimilarity.gradient <- seq(from = 0, to = max_gradient, length.out = n_row)
  for (i in seq_len(n_row)){
    print(i)
    if (i == 1){
      row_temp <- rbeta(n_col, 1, n_col)
      col_to_remove <- sample(x = seq_len(n_col), size = n_col-n_col*density_row)
      row_temp[col_to_remove] <- 0
      list_initial[[i]] <- row_temp
    } else {
      while (length(list_initial) < i) {
        row_temp <- rbeta(n_col, 1, n_col)
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


## -----------------------------------------------------------------------------
n_species_types <- c(13, 3, 4)
theta_types <- c(1, 0.75, 0.5, 0.25, 0.1, 0.05)
n_resources_types <- c(1,2,4,8,16,32)

community.initial.df <- as.list(
  lapply(n_species_types, 
         gradient_df_generator, 
         n_row = n_rep, 
         density_row = 1,
         max_gradient = 0.7,
         error_interval = 0.15)
)


## ---- eval=FALSE--------------------------------------------------------------
## for (n_species in n_species_types){
##     for (theta in theta_types) {
##         sorensen <- c()
##         for (i in seq_len(n_rep)){
##             for (n_resources in n_resources_types) {
##                 ### generate E ####
##                 Etest <- randomE(n_species = n_species,
##                                  n_resources = n_resources,
##                                  mean_consumption = theta*n_resources,
##                                  exact = TRUE)
## 
##                 ### calculate rho ####
##                 if (n_resources == max(n_resources_types)){
##                     Etest_pos <- Etest
##                     Etest_pos[Etest_pos<0] <- 0
##                     for (j in seq_len(n_species - 1)){
##                         for (k in 2:n_species){
##                             sorensen <- c(sorensen,
##                                           sum(apply(Etest_pos[c(j,k),], 2, min)))
##                         }
##                     }
##                 }
## 
##                 if (n_resources > 1){
##                     Priority <- t(apply(matrix(sample(n_species * n_resources), nrow = n_species), 1, order))
##                     Priority <- (Etest > 0) * Priority
##                 } else {
##                     Priority <- NULL
##                 }
## 
##                 print(paste0("n_species=",n_species, " theta=",theta, " i=", i, " n_resources=", n_resources))
##                 x0temp <- as.numeric(community.initial.df[[match(n_species, n_species_types)]][i,])
##                 x0temp <- 10*x0temp/sum(x0temp)
##                 CRMtest <- simulateConsumerResource(n_species = n_species,
##                                                     n_resources = n_resources,
##                                                     x0 = x0temp, #rep(10, n_species),
##                                                     resources = rep(100, n_resources),
##                                                     E = Etest,
##                                                     # trophic_priority = Priority,
##                                                     stochastic = TRUE,
##                                                     t_end = 1000,
##                                                     t_step = 1,
##                                                     t_store = 1000)
##                 CRMspecies <- getCommunity(CRMtest)
##                 CRMspeciesTotal <- sum(CRMspecies)
##                 result_df[nrow(result_df)+1,] <- c(n_species, theta, i, n_resources, CRMspeciesTotal)
##                 result_df2[nrow(result_df2)+1,] <- c(CRMspecies, rep(NA, 13-length(CRMspecies)))
##                 # makePlotRes(CRMtest$resources)
##                 # makePlot(CRMtest$matrix)
##             }
##         }
##         rho_mean <- mean(sorensen)
##         rho_sd <- var(sorensen)
##         sorensen_df[nrow(sorensen_df)+1, ] <- c(n_species, theta, rho_mean, rho_sd)
##     }
## }


## ---- eval=FALSE--------------------------------------------------------------
## p_fig2_result_df <- ggplot(result_df, aes(x = n_resources, y = value, group = n_resources)) +
##     geom_boxplot(outlier.shape = NA) +
##     geom_jitter(alpha = 0.2, width = 0.2) +
##     facet_grid(. ~ factor(n_species, levels = n_species_types)) +
##     theme_bw() +
##     scale_x_continuous(trans = "log2", breaks = n_resources_types) +
##     xlab("number of resources") +
##     ylab("growth yield (number of individuals)")
## p_fig2_result_df

