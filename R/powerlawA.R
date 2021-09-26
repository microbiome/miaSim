#' Interaction matrix with Power-Law network adjacency matrix
#'
#' Where N is the an Interspecific Interaction matrix with values drawn from
#' a normal distribution H the interaction strength heterogeneity drawn from
#' a power-law distribution with the parameter alpha, and G the adjacency matrix
#' of with out-degree that reflects the heterogeneity of the powerlaw.
#' A scaling factor s may be used to constrain the values of the interaction
#' matrix to be within a desired range. Diagonal elements of A are defined
#' by the parameter d.
#'
#' @param n.species integer: the number of species
#' @param alpha numeric: the power-law distribution parameter. Should be > 1.
#' (default: \code{alpha = 3.0}) Larger values will give lower interaction
#' strength heterogeneity, whereas values closer to 1 give strong heterogeneity
#' in interaction strengths between the species. In other words, values of alpha
#' close to 1 will give Strongly Interacting Species (SIS).
#' @param stdev numeric: the standard deviation parameter of the normal
#' distribution with mean 0 from which the elements of the nominal interspecific
#' interaction matrix N are drawn. (default: \code{stdev = 1})
#' @param s numeric: scaling parameter with which the final global
#' interaction matrix A is multiplied. (default: \code{s = 0.1})
#' @param d numeric: diagonal values, indicating self-interactions (use
#' negative values for stability). (default: \code{s = 1.0})
#'
#' @return The interaction matrix A with n rows and n columns.
#'
#' @references Gibson TE, Bashan A, Cao HT, Weiss ST, Liu YY (2016)
#' On the Origins and Control of Community Types in the Human Microbiome.
#' PLOS Computational Biology 12(2): e1004688.
#' https://doi.org/10.1371/journal.pcbi.1004688
#'
#' @docType methods
#' @aliases powerlawA-numeric
#' @aliases powerlawA,numeric-method
#'
#' @importFrom poweRlaw rplcon
#' @examples
#' # Low interaction heterogeneity
#' A_low <- powerlawA(n.species = 10, alpha = 3)
#' # Strong interaction heterogeneity
#' A_strong <- powerlawA(n.species = 10, alpha = 1.01)
#' @export

setGeneric("powerlawA",signature = "n.species",
            function(n.species, alpha, stdev = 1, s = 0.1, d = -1)
                standardGeneric("powerlawA"))

setMethod("powerlawA", signature = c(n.species = "numeric"),
            function(n.species, alpha, stdev = 1, s = 0.1, d = -1){
            # Nominal Interspecific Interaction matrix N
            N <- matrix(
                data = rnorm(n.species^2, mean = 0, sd = stdev),
                nrow = n.species,
                ncol = n.species
            )

            # power law sample
            pl <- rplcon(n = n.species, xmin = 1, alpha = alpha)
            pl[is.infinite(pl)] = 10^308
<<<<<<< HEAD
=======
            
            
>>>>>>> e220ab5a825005bf8dcd8f7ccf505d4817fc2e61

            # Interaction strength heterogeneity
            H <- diag(1 + (pl-min(pl))/(max(pl)-min(pl)))

            # Adjacency matrix G of power-law out-degree digraph ecological
            #network
            deg <- 0.1*n.species

            h <- pmin(ceiling(deg*pl/mean(pl)), n.species)

            G <- matrix(0, nrow = n.species, ncol = n.species)
            for(i in seq_len(n.species)){
                index <- sample(x = seq_len(n.species), size = h[i])
                G[index, i] <- 1
            }
            #G[t(G) == 1] <- 1
            A <- N %*% H * G
            A <- A*s/max(A)
            diag(A) <- d
            colnames(A) <- seq_len(n.species)
            rownames(A) <- seq_len(n.species)
            return(A)
})

