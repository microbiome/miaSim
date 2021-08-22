# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

library(miaSim)
glv <- function(
  N,
  A,
  b = runif(N),
  x = runif(N),
  tend = 1000,
  norm = FALSE
){
  # 1 Model specification
  # Model parameters
  parameters <- cbind(b, A)

  # 2 Model application
  # Time specification
  times <- seq(0, tend, by = 0.01)
  # Model integration
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

x <- miaSim::glv(N = 4, A = powerlawA(n = 4, alpha = 2), tend = 1000)


