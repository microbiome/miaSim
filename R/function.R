

library(SummarizedExperiment)
library(microsim)
library(deSolve)

setGeneric("glv",signature = "x",
           function(N, A, b = runif(N), x = runif(N), tend = 1000, norm = FALSE)
             standardGeneric("glv"))

setGeneric("conversionSE",signature = "x",
           function(x)
             standardGeneric("conversionSE"))

.dxdt <- function(t, x, parameters){
  b <- parameters[,1]
  A <- parameters[,2:ncol(parameters)]

  dx <- x*(b+A %*% x)
  list(dx)
}

setMethod("glv", signature = c(x="matrix"),
          function(N, A, b = runif(N), x = runif(N), tend = 1000, norm = FALSE){
            parameters <- cbind(b, A)
            times <- seq(0, tend, by = 0.01)

            .dxdt(t, x, parameters)

            out <- ode(
              y = x,
              times = times,
              func = dxdt,
              parms = parameters
            )
            spab <- t(out[,2:ncol(out)])
            spab <- spab[,round(seq(1, tend*100, length.out = tend))]
            if(norm){
              spab <- t(t(spab)/colSums(spab))
            }
            return(spab)
          }
)



row_data <- data.frame(Kingdom = "A",
                       Phylum = rep(c("B1", "B2"), c(500, 500)),
                       Class = rep(c("C1", "C2"), each = 500),
                       ASV = paste0("D", 1:1000),
                       row.names = rownames(paste0("species", seq_len(1000))),
                       stringsAsFactors = FALSE)

row_data <- t(row_data)

col_data <- data.frame(sampleID = c(1:1000),
                       time = as.Date(sample( as.numeric(as.Date('2000-01-01')): as.numeric(as.Date('2014-01-01')), 1000,
                                              replace = T),
                                      origin = '1970-01-01'),
                       row.names = colnames(paste0("sample", 1:1000)),
                       stringsAsFactors = FALSE)


setMethod("conversionSE", signature = c(x="SummarizedExperiment"),
          function(x){
           SummarizedExperiment(assays = x,
                                rowData = row_data,
                                colData = col_data)
          }
)

result <- glv(N = 4, A = powerlawA(n = 4, alpha = 2), tend = 1000)
SE <- conversionSE(result)

