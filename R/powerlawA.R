#' Interaction matrix with Power-Law network adjacency matrix
#'
#' Generate an interaction matrix A that can be decomposed as NH.Gs .
#' Where N is the Normal Interspecific Interaction matrix, H the interaction
#' strength heterogeneity drawn from a power-law distribution with given
#' parameter alpha, and G the adjacency matrix of power-law out-degree digraph
#' ecological network, and s a scaling factor. Diagonal elements of A are
#' subsequently set to -1.
#'
#' @references Gibson TE, Bashan A, Cao HT, Weiss ST, Liu YY (2016)
#' On the Origins and Control of Community Types in the Human Microbiome.
#' PLOS Computational Biology 12(2): e1004688.
#' https://doi.org/10.1371/journal.pcbi.1004688
#'
#' @param n the number of species
#' @param alpha the power-law distribution parameter. Should be > 1.
#' Larger values will give lower interaction strength heterogeneity,
#' whereas values closer to 1 give strong heterogeneity in interaction strengths
#' between the species. In other words, values of alpha close to 1 will give
#' Strongly Interacting Species (SIS).
#' @param stdev the standard deviation parameter of the normal distribution
#' with mean 0 from which the elements of the nominal interspecific
#' interaction matrix N are drawn. (default: \code{stdev = 1})
#' @param s scaling parameter with which the final global interaction matrix A
#' is multiplied. Default set to NULL where s is set to 0.1*max(A) after
#' constructing the matrix A = NH*G. (default: \code{s = 0.1})
#'
#' @return The global interaction matrix A with n rows and n columns.
#'
#' @docType methods
#' @aliases powerlawA-numeric
#' @aliases powerlawA,numeric-method
#'
#' @importFrom poweRlaw rplcon
#' @examples
#' # Low interaction heterogeneity
#' A_low <- powerlawA(n = 10, alpha = 3)
#' # Strong interaction heterogeneity
#' A_strong <- powerlawA(n = 10, alpha = 1.01)
#' @export

setGeneric("powerlawA",signature = "n",
            function(n, alpha, stdev, s)
                standardGeneric("powerlawA"))

setMethod("powerlawA", signature = c(n="numeric"),
            function(n, alpha, stdev, s){
                if (missing(stdev) || missing(s)) {
                stdev <- 1
                s <- 0.1
            }
            # Nominal Interspecific Interaction matrix N
            N <- matrix(
                data = rnorm(n^2, mean = 0, sd = stdev),
                nrow = n,
                ncol = n
            )
            #diag(N) <- 0

            # power law sample
            pl <- rplcon(n = n, xmin = 1, alpha = alpha)
            # Interaction strength heterogeneity H
            H <- vapply(seq_len(n), FUN = function(i){
                1 + ((pl[i]-min(pl))/(max(pl)-min(pl)))
            }, numeric(1))
            H <- diag(H)

            # Adjacency matrix G of power-law out-degree digraph ecological
            #network
            d <- 0.1*n
            h <- vapply(seq_len(n), FUN = function(i){
                min(
                    ceiling(d*pl[i]/mean(pl)),
                n
            )
            },numeric(1))
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
)
