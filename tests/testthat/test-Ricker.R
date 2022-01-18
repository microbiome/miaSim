test_that("simulateRicker", {
    #check simulateRicker with 10 species interaction matrix (alpha = 1.01)
    A <- powerlawA(10, alpha = 1.01)
    ExampleRicker <- simulateRicker(n_species=10,A,tend=100)
    expect_type(ExampleRicker, "double")
    expect_equal(dim(ExampleRicker), c(10, 100))

    #check simulateRicker with custom inputs
    ExampleRicker2<-simulateRicker(n_species=10,A,sigma=-0.05,tend=100,
        explosion_bound=10^4,norm=TRUE)
    expect_type(ExampleRicker2, "double")
    expect_equal(dim(ExampleRicker), c(10, 100))

    #check simulateRicker with errors in inputs
    expect_error(ErrorRicker1<-simulateRicker(10,A,x=runif(9),tend=100))
    expect_error(ErrorRicker2<-simulateRicker(9,A,tend=100))
    expect_error(ErrorRicker3<-simulateRicker(10,A,K=runif(9),tend=100))

    expect_error(ErrorLogistic2 <- simulateStochasticLogistic(
        n_species = 3, b = c(0.2, 0.1), k = c(1000, 2000),
        dr = c(0.001, 0.0015), x = c(3, 1)))
    expect_error(ErrorLogistic3 <- simulateStochasticLogistic(
        n_species = 4, partial = FALSE, stochastic = 1))

})
