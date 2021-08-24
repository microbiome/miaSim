library(SummarizedExperiment)
library(microsim)
library(deSolve)

setGeneric("glv",signature = c("x"),
           function(N, A, b = runif(N), x = runif(N), tend = 1000, norm = FALSE)
             standardGeneric("glv"))

setGeneric("conversionSE",signature = c("x"),
           function(x, ...)
             standardGeneric("conversionSE"))

setMethod("glv", signature = c(x="matrix"),
          function(N, A, b = runif(N), x = runif(N), tend = 1000, norm = FALSE){
            parameters <- cbind(b, A)
            times <- seq(0, tend, by = 0.01)
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

setMethod("conversionSE", signature = c(x="matrix"),
          function(x){
           SummarizedExperiment(assays = x)
          }
)

result <- glv(N = 4, A = powerlawA(n = 4, alpha = 2), tend = 1000)

SE <- conversionSE(result)
