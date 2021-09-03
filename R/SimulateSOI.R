#' Self-Organised Instability Model Community Simulation
#'
#' Generate time-series with The Self-Organised Instability (SOI) model.
#' Implements a K-leap method for accelerating stochastic simulation.
#'
#' @param N Integer: number of species
#' @param I Integer: community size, number of available sites (individuals)
#' @param A a interaction matrix of dimension NxN
#' @param com a vector inditicating initial community abundances.
#' If NULL (default), based on migration rates.
#' @param tend Integer: number of timepoints to be returned in the time series
#' (nr of generations)
#' @param k Integer: the number of transition events that are allowed to take
#' place during one leap. By default set to 5. Higher values reduce runtime,
#' but also accuracy of the simulation.
#' @param norm logical to indicate whether the time series should be returned
#' with the abundances as proportions (norm = TRUE)
#' or the raw counts (norm = FALSE, default)
#' @return \linkS4class{SummarizedExperiment} object with abundance matrix
#' consisting of species abundance as rows and time points as columns
#' @examples simulateSOI(N = 10, I = 1000, A = powerlawA(n = 10, alpha = 1.2),
#'                     k=5, com = NULL, tend = 150, norm = TRUE)
#' @export

setGeneric("simulateSOI",signature = "N",
           function(N, I, A, k, com, tend, norm = FALSE)
             standardGeneric("simulateSOI"))

#' @export

powerlawA <- function(
  n, # number of species
  alpha, # power-law distribution parameter
  stdev = 1, # sd normal distribution
  s = 0.1 # scaling parameter, default: 0.1/max(A)
){
  # Nominal Interspecific Interaction matrix N
  N <- matrix(
    data = rnorm(n^2, mean = 0, sd = stdev),
    nrow = n,
    ncol = n
  )
  # power law sample
  pl <- rplcon(n = n, xmin = 1, alpha = alpha)
  # Interaction strength heterogeneity H
  H <- sapply(seq_len(n), FUN = function(i){
    1 + ((pl[i]-min(pl))/(max(pl)-min(pl)))
  })
  H <- diag(H)
  # Adjacency matrix G of power-law out-degree digraph ecological network
  d <- 0.1*n
  h <- sapply(seq_len(n), FUN = function(i){
    min(
      ceiling(d*pl[i]/mean(pl)),
      n
    )
  })
  G <- matrix(0, nrow = n, ncol = n)
  for(i in seq_len(n)){
    index <- sample(x = seq_len(n), size = h[i])
    G[index, i] <- 1
  }
  A <- N %*% H * G
  A <- A*s/max(A)
  diag(A) <- -1
  colnames(A) <- seq_len(n)
  rownames(A) <- seq_len(n)
  return(A)
}

#' @export

update_propensities <- function(
  I, #total nr of sites
  counts, # species counts vector incl the empty sites
  death_rates, # species-specific death rates
  migr_rates, # species-specific migration rates
  pos_inter_rates, # positive interaction rates matrix
  neg_inter_rates # negative interaction rates matrix
){
  N <- length(counts)-1
  # migration AND positive interaction
  m_propensities <- counts[N+1]*migr_rates/N
  for(stronger in seq_len(N)){
    if(sum(pos_inter_rates[stronger,]!=0) >0){
      weaker <- which(pos_inter_rates[stronger,]!=0)
      inter_rate <- pos_inter_rates[stronger,weaker]
      m_propensities[stronger] <- m_propensities[stronger] +
        sum(counts[weaker] * inter_rate)/I*counts[stronger]/I*counts[N+1]
    }
  }
  # death / extinction
  d_propensities <- (counts[seq_len(N)]*death_rates)
  # competition / negative interaction
  j_propensities <- c()
  for(stronger in seq_len(N)){
    if(sum(neg_inter_rates[stronger,]!=0) >0){
      weaker <- which(neg_inter_rates[stronger,]!=0)
      inter_rate <- neg_inter_rates[stronger, weaker]
      j_propensities <- c(j_propensities,
                          counts[stronger]*counts[weaker]*inter_rate/I)
    }
  }
  propensities <- c(m_propensities, d_propensities, j_propensities)
  return(propensities)
}

setMethod("simulateSOI", signature = c(N="numeric"),
          function(N, I, A, k, com, tend, norm = FALSE){
            #species-specific immigration probabilities
            migr_rates <- runif(N, min = 0.1, max = 0.8)
            #species-specific extinction probabilities
            death_rates <- runif(N, min = 0.01, max = 0.08)
            # COUNTS VECTOR
            # initial counts vector based on migration rates of the species
            # additional draw from the uniform distribution for the empty sites
            if(is.null(com)){
              counts <- c(migr_rates, runif(1))
              counts <- round((counts/sum(counts))*I)
            } else if(length(com) == N+1){
              counts <- com
            } else if(length(com) == N){
              empty_count <- I - sum(com)
              counts <- c(com, empty_count)
            } else {
              stop("com vector argument can only be either of length N, or N+1
                    in case the number of empty sites are included")
            }
            # making sure counts vector adds up to I with the empty sites
            # (if not, the difference is added or subtracted from the empty
            # sites to get a total of I)
            if(sum(counts)!=I){
              counts[N+1] <- counts[N+1] + (I-sum(counts))
            }
            # TRANSITION MATRIX & INTERACTION RATES
            # Initialise trans_mat transition matrix reprenting the actual
            # reactions
            # First N columns = migration: species_i counts go up with 1,
            # empty sites go down with -1
            # Next N columns = death: species_i counts go down with -1,
            # empty sites go up with +1
            # Last columns: the competition / negative interaction jumps:
            trans_mat <- cbind(
              rbind(diag(1, ncol = N, nrow = N), rep(-1, times = N)),
              rbind(diag(-1, ncol = N, nrow = N), rep(1, times = N))
            )
            # Interaction:
            # happens when happens when sp_stronger interacts more strongle
            # with sp_weaker, than sp_weaker with sp_stronger
            pos_inter_rates <- matrix(0, nrow = N, ncol = N)
            neg_inter_rates <- pos_inter_rates # copy to initialise
            # pos interaction and immigration
            # AND sp_stronger interacts positively with sp_weaker
            for(stronger in seq_len(N)){
              weaker <- which(A[stronger,] > A[,stronger] & A[,stronger] >= 0)
              if(length(weaker) > 0){
                rate <- A[stronger, weaker] + A[weaker, stronger]
                pos_inter_rates[stronger,weaker] <- rate
              }
            }
            # neg interaction: competition
            # AND sp_weaker interacts negatively with sp_stronger
            for(stronger in seq_len(N)){
              weaker_vec <- which(A[stronger,] > A[,stronger] & A[,stronger] < 0)
              if(length(weaker_vec) > 0){
                rate <- A[stronger, weaker_vec] - A[weaker_vec, stronger]
                neg_inter_rates[stronger, weaker_vec] <- rate
                for(weaker in weaker_vec){
                  jump <- rep(0, times = N+1)
                  jump[stronger] <- 1
                  jump[weaker] <- -1
                  trans_mat <- cbind(trans_mat, jump)
                }
              }
            }
            # INITIAL PROPENSITIES
            propensities <- update_propensities(I, counts, death_rates,
                                                migr_rates,pos_inter_rates, neg_inter_rates)

            # K-LEAPS SIMULATION
            sys_t <- seq(from = 0, to = tend, length.out = tend)
            # time variables
            current_t <- 0 # current time
            sample_t <- sys_t[1]  # sampled time
            series_t <- 1 # column to be filled with counts in series matrix
            #(next generation)
            series <- matrix(0, nrow = N, ncol = tend)
            colnames(series) <- paste0("t", seq_len(tend))
            rownames(series) <- paste0("sp_", seq_len(N))
            while(series_t <= tend){
              c <- sum(propensities)
              p <- propensities/c
              tau <- rgamma(n = 1, shape = k, scale = 1/c)
              current_t <- current_t + tau
              # if reached end of simulation:
              if(current_t >= tend){
                series[,tend] <- counts
                break
              }
              # if current_t exceeds sample_t, update series with current counts
              if(current_t > sample_t){
                series[,series_t] <- counts[seq_len(N)]
                series_t <- series_t +1
                index <- which(sys_t >= sample_t & sys_t < current_t)
                sample_t <- sys_t[index][length(index)]
              }
              # which reactions allowed?
              reactionsToFire <- as.vector(rmultinom(n= 1, size= k, prob = p))
              # update counts with allowed transitions (reactions)
              counts <- sapply(counts + trans_mat%*%reactionsToFire,
                               FUN = function(x){x = max(0,x)})
              # update propensities with updated counts
              propensities <- update_propensities(I, counts, death_rates,
                                                  migr_rates, pos_inter_rates, neg_inter_rates)
            }
            if(norm){
              series <- t(t(series)/colSums(series))
            }
            series
            SOI <- SummarizedExperiment(assays = list(counts=series))
            return(SOI)
          })
