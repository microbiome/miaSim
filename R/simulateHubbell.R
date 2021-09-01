#' Hubbell community simulation
#'
#' Neutral species abundances simulation turned into
#' \linkS4class{SummarizedExperiment} according to the Hubbell model.
#'
#' @param N Integer: the amount of different species initially
#' in the local community
#' @param M Integer: amount of different species in the
#' metacommunity,including those of the local community
#' @param I Integer: fixed amount of individuals in the
#' local community (default: \code{I = 1000})
#' @param d Integer: fixed amount of deaths of local community
#' individuals in each generation (default: \code{d = 10})
#' @param m Numeric: immigration rate: the probability that a
#' death in the local community is replaced by a migrant of the metacommunity
#' rather than by the birth of a local community member
#' (default: \code{m = 0.02})
#' @param tskip  Integer: number of generations that should
#' not be included in the outputted species abundance matrix.
#' (default: \code{tskip = 0})
#' @param tend Integer: number of simulations
#' to be simulated
#' @param norm Logical: whether the time series should
#' be returned with the abundances as proportions (norm = TRUE) or the raw
#' counts (norm = FALSE, default)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @aliases simulateNeutral
#'
#' @examples
#' col_data <- data.frame(sampleID = c(1:100),
#'                         stringsAsFactors = FALSE)
#'
#' row_data <- data.frame(Kingdom = "A",
#'                         Phylum = rep(c("B1", "B2"), c(50, 50)),
#'                         Class = rep(c("C1", "C2"), each = 50),
#'                         ASV1 = paste0("D", 1:100),
#'                         ASV2 = paste0("E", 1:100),
#'                         ASV3 = paste0("F", 1:100),
#'                         ASV4 = paste0("G", 1:100),
#'                         ASV5 = paste0("H", 1:100),
#'                         ASV6 = paste0("J", 1:100),
#'                         ASV7 = paste0("K", 1:100),
#'                         row.names = rownames(paste0("species", seq_len(10))),
#'                         stringsAsFactors = FALSE)
#' row_data <- t(row_data)
#'
#' simulateHubbell(N = 8, M = 10, I = 1000, d = 50, m = 0.02, tend = 100)
#'
#' @return \code{simulateHubbell} returns a \code{\link{SummarizedExperiment}}
#' object
#'
#' @references Rosindell, James et al. "The unified neutral theory of
#' biodiversity and biogeography at age ten." Trends in ecology & evolution
#' vol. 26,7 (2011).
#' @export

setGeneric("simulateHubbell",signature = "N",
    function(N, M, I, d, m, tskip = 0, tend, norm = FALSE)
        standardGeneric("simulateHubbell"))


setMethod("simulateHubbell", signature = c(N="numeric"),
    function(N, M, I, d, m, tskip = 0 , tend, norm = FALSE){
            # First setting the function arguments right
            pbirth <- runif(N, min = 0, max = 1)
            pmigr <- runif(M, min = 0, max = 1)
            if(length(pbirth)!=N | length(pmigr)!=M){
                stop("Either length of pbirth vector does not match with N or
                    length of pmigr vector does not match with M")
            }
            pbirth <- c(pbirth, rep(0, times = (M-N)))
            pbirth <- pbirth/sum(pbirth)
            pmigr <- pmigr/sum(pmigr)
            com <- ceiling(I*pbirth)
            if(sum(com)<I){
                ind <- sample(1:M, size = I-sum(com), prob = pbirth)
                com[ind] <- com[ind] +1
            } else if(sum(com)>I){
                ind <- sample(1:M, size = sum(com)-I, prob = 1-pbirth)
                com[ind] <- com[ind] -1
            }
            # The simulation
            tseries <- matrix(0, nrow = M, ncol = tend)
            colnames(tseries) <- paste0("t", 1:tend)
            rownames(tseries) <- 1:M

            com[which(com < 0)] <- 0
            tseries[,1] <- com
            for (t in 2:tend){
            # Each iteration the probability of births is updated by the counts
                pbirth <- com/sum(com)
                pbirth[which(pbirth < 0)] <- 0
            # Probability of births is used to pick the species that will die
            # because species with count 0 will have probability 0 and species
            #not present in the community can also not die
                deaths <- rmultinom(n = 1, size = d, prob = pbirth)
            while(sum(com-deaths <0) >0){
                neg_sp <- which(com-deaths <0)
                pbirth[neg_sp] <- 0
                deaths <- rmultinom(n = 1, size = d, prob = pbirth)
            }

            event <- rbinom(d, 1, prob = m) # immigration rate m: probability
            #death replaced by immigrant
            #immigration 1, birth 0
            n_migrants <- sum(event)
            n_births <- length(event) - n_migrants

            births <- rmultinom(1, n_births, prob = pbirth)
            migr <- rmultinom(1, n_migrants, prob = pmigr)
            com <- com - deaths + births + migr
            com[which(com < 0)] <- 0
            tseries[,t] <- com
            }
            if(norm){
                tseries <- t(t(tseries)/colSums(tseries))
            }

            AbundanceM <- tseries[, (tskip +1):tend]
            SE <- SummarizedExperiment(assays = list(counts = AbundanceM))

            return(SE)
})
