#' Self-Organised Instability model (SOI) simulation
#'
#' Generate time-series with The Self-Organised Instability (SOI) model.
#' Implements a K-leap method for accelerating stochastic simulation.
#'
#' @param n.species Integer: number of species
#' @param I Integer: community size, number of available sites (individuals)
#' @param A interaction matrix
#' @param com a vector: initial community abundances
#' If (default: \code{com = NULL}), based on migration rates
#' @param tend Integer: number of timepoints to be returned in the time series
#' (number of generations)
#' @param k Integer: the number of transition events that are allowed to take
#' place during one leap. (default: \code{k = 5}). Higher values reduce runtime,
#' but also accuracy of the simulation.
#' @param norm Logical: indicates whether the time series should be returned
#' with the abundances as proportions (\code{norm = TRUE})
#' or the raw counts (default: \code{norm = FALSE})
#'
#' @return \linkS4class{SummarizedExperiment} object containing abundance matrix
#' consisting of species abundance as rows and time points as columns
#'
#' @examples
#' A <- miaSim::powerlawA(10, alpha = 1.2)
#'
#' ExampleSOI <- simulateSOI(n.species = 10, I = 1000, A, k=5, com = NULL,
#'                                             tend = 150, norm = TRUE)
#' @docType methods
#' @aliases simulateSOI-numeric
#' @aliases simulateSOI,numeric-method
#'
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#'
#' @export
simulateSOI <- function(n.species, I, A, k = 5, com = NULL, tend, norm = FALSE){
            #species-specific immigration probabilities
            migr_rates <- runif(n.species, min = 0.1, max = 0.8)
            #species-specific extinction probabilities
            death_rates <- runif(n.species, min = 0.01, max = 0.08)
            # COUNTS VECTOR
            # initial counts vector based on migration rates of the species
            # additional draw from the uniform distribution for the empty sites

            #input check
            if(!all(vapply(list(n.species,I,k), isPositiveInteger,
                    logical(1)))){
                stop("n.species,I and k values must be integer.")}

            if(is.null(com)){
                    counts <- c(migr_rates, runif(1))
                    counts <- round((counts/sum(counts))*I)
    }       else if(length(com) == n.species+1){
                counts <- com
    }       else if(length(com) == n.species){
                empty_count <- I - sum(com)
                counts <- c(com, empty_count)
    }       else {
                stop("com vector argument can only be either of length N, or N+1
                    in case the number of empty sites are included")
    }
            # making sure counts vector adds up to I with the empty sites
            # (if not, the difference is added or subtracted from the empty
            # sites to get a total of I)
            if(sum(counts)!=I){
                counts[n.species+1] <- counts[n.species+1] + (I-sum(counts))
    }
            # TRANSITION MATRIX & INTERACTION RATES
            # Initialise trans_mat transition matrix representing the actual
            # reactions
            # First n.species columns = migration: species_i counts go up with 1
            # empty sites go down with -1
            # Next n.species columns = death: species_i counts go down with -1,
            # empty sites go up with +1
            # Last columns: the competition / negative interaction jumps:
            trans_mat <- cbind(
                rbind(diag(1, ncol = n.species, nrow = n.species),
                        rep(-1, times = n.species)),
                rbind(diag(-1, ncol = n.species, nrow = n.species),
                        rep(1, times = n.species))
    )
            # Interaction:
            # happens when sp_stronger interacts stronger with sp_weaker,
            # than sp_weaker with sp_stronger
            pos_inter_rates <- matrix(0, nrow = n.species, ncol = n.species)
            neg_inter_rates <- pos_inter_rates # copy to initialise
            # pos interaction and immigration
            # AND sp_stronger interacts positively with sp_weaker
            for(stronger in seq_len(n.species)){
                weaker <- which(A[stronger,] > A[,stronger] & A[,stronger] >= 0)
                    if(length(weaker) > 0){
                    rate <- A[stronger, weaker] + A[weaker, stronger]
                    pos_inter_rates[stronger,weaker] <- rate
                    }
                }
            # neg interaction: competition
            # AND sp_weaker interacts negatively with sp_stronger
            for(stronger in seq_len(n.species)){
            weaker_vec <- which(A[stronger,] > A[,stronger] & A[,stronger] < 0)
                if(length(weaker_vec) > 0){
                rate <- A[stronger, weaker_vec] - A[weaker_vec, stronger]
                neg_inter_rates[stronger, weaker_vec] <- rate
                for(weaker in weaker_vec){
                    jump <- rep(0, times = n.species+1)
                    jump[stronger] <- 1
                    jump[weaker] <- -1
                    trans_mat <- cbind(trans_mat, jump)
                }
                }
            }

            # INITIAL PROPENSITIES
            propensities <- updatePropensities(I, counts, death_rates,
                                migr_rates,pos_inter_rates, neg_inter_rates)
            # K-LEAPS SIMULATION
            sys_t <- seq(from = 0, to = tend, length.out = tend)
            # time variables
            current_t <- 0 # current time
            sample_t <- sys_t[1]  # sampled time
            series_t <- 1 # column to be filled with counts in series matrix
            #(next generation)
            series <- matrix(0, nrow = n.species, ncol = tend)
            colnames(series) <- paste0("t", seq_len(tend))
            rownames(series) <- paste0("sp_", seq_len(n.species))
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
                series[,series_t] <- counts[seq_len(n.species)]
                series_t <- series_t +1
                index <- which(sys_t >= sample_t & sys_t < current_t)
                sample_t <- sys_t[index][length(index)]
    }
            # which reactions allowed?
            reactionsToFire <- as.vector(rmultinom(n= 1, size= k, prob = p))
            # update counts with allowed transitions (reactions)
            counts <- vapply(counts + trans_mat%*%reactionsToFire,
                                FUN = function(x){x = max(0,x)}, numeric(1))
            # update propensities with updated counts
            propensities <- updatePropensities(I, counts, death_rates,
                                migr_rates, pos_inter_rates, neg_inter_rates)
    }
            if(norm){
                series <- t(t(series)/colSums(series))
    }
            series
            SOI <- SummarizedExperiment(assays = list(counts=series))
            return(SOI)
}

updatePropensities <- function(
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
                sum(counts[weaker] *
                        inter_rate)/I*counts[stronger]/I*counts[N+1]
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