#' Hubbell's neutral model simulation
#'
#' Neutral species abundances simulation according to the Hubbell model.
#'
#' @param n.species Integer: the amount of different species initially
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
#' be returned with the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#'
#' @docType methods
#' @aliases simulateHubbell-numeric
#' @aliases simulateHubbell,numeric-method
#' @aliases simulateNeutral
#'
#' @examples
#' colData <- DataFrame(sampleID = c(seq_len(100)),
#'                         time = as.Date(100, origin = "2000-01-01"))
#'
#' rowData <- data.frame(Kingdom = "Animalia",
#'                 Phylum = rep(c("Platyhelminthes", "Mollusca"), c(50, 50)),
#'                 Class = rep(c("Turbellaria", "Polyplacophora"), each = 50),
#'                 ASV1 = paste0("D", seq_len(100)),
#'                 ASV2 = paste0("E", seq_len(100)),
#'                 ASV3 = paste0("F", seq_len(100)),
#'                 ASV4 = paste0("G", seq_len(100)),
#'                 ASV5 = paste0("H", seq_len(100)),
#'                 ASV6 = paste0("J", seq_len(100)),
#'                 ASV7 = paste0("K", seq_len(100)),
#'                 row.names = rownames(paste0("species", seq_len(10))),
#'                 stringsAsFactors = FALSE)
#'
#' rowData <- t(rowData)
#'
#' ExampleHubbell <- simulateHubbell(n.species = 8, M = 10, I = 1000, d = 50,
#'                                                         m = 0.02, tend = 100)
#' rowData(ExampleHubbell) <- rowData
#' colData(ExampleHubbell) <- colData
#'
#' @return \code{simulateHubbell} returns a \linkS4class{SummarizedExperiment}
#' object containing matrix with species abundance as rows and
#' time points as columns
#'
#' @importFrom stats rbinom
#' @importFrom stats rmultinom
#'
#' @references Rosindell, James et al. "The unified neutral theory of
#' biodiversity and biogeography at age ten." Trends in ecology & evolution
#' vol. 26,7 (2011).
#
#' @export

setGeneric("simulateHubbell",signature = "n.species",
    function(n.species, M, I = 1000, d = 10, m = 0.02, tskip = 0,
            tend, norm = FALSE)
        standardGeneric("simulateHubbell"))

setMethod("simulateHubbell", signature = c(n.species="numeric"),
    function(n.species, M, I = 1000, d = 10, m = 0.02, tskip = 0,
            tend, norm = FALSE){
            pbirth <- runif(n.species, min = 0, max = 1)
            pmigr <- runif(M, min = 0, max = 1)
            pbirth <- c(pbirth, rep(0, times = (M-n.species)))
            pbirth <- pbirth/sum(pbirth)
            pmigr <- pmigr/sum(pmigr)
            com <- ceiling(I*pbirth)
            if(sum(com)>I){
                ind <- sample(seq_len(M), size = sum(com)-I, prob = 1-pbirth)
                com[ind] <- com[ind] -1
            }

            tseries <- matrix(0, nrow = M, ncol = tend)
            colnames(tseries) <- paste0("t", seq_len(tend))
            rownames(tseries) <- seq_len(M)
            com[which(com < 0)] <- 0
            tseries[,1] <- com
            for (t in 2:tend){
                pbirth <- com/sum(com)
                pbirth[which(pbirth < 0)] <- 0
                deaths <- rmultinom(n = 1, size = d, prob = pbirth)
            while(sum(com-deaths <0) >0){ #species with count 0 have probability
            # 0 and species not present in the community can also not die
                neg_sp <- which(com-deaths <0)
                pbirth[neg_sp] <- 0
                deaths <- rmultinom(n = 1, size = d, prob = pbirth)
            }
            event <- rbinom(d, 1, prob = m) # immigration rate m: probability
            # death replaced by immigrant; immigration 1, birth 0
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
