test_that("simulateRicker", {
    # check simulateRicker with 10 species interaction matrix (alpha = 1.01)
    A <- powerlawA(10, alpha = 1.01)
    ExampleRicker <- simulateRicker(n_species = 10, A, t_end = 100, t_store = 100)

    expect_s4_class(ExampleRicker, "TreeSummarizedExperiment")
    expect_type(ExampleRicker@assays@data@listData[["counts"]], "double")
    expect_equal(dim(ExampleRicker@assays@data@listData[["counts"]]), c(10, 100))


    # check simulateRicker with custom inputs
    ExampleRicker2 <- simulateRicker(
        n_species = 10, A, error_variance = -0.05, t_end = 100,
        explosion_bound = 10^4, norm = TRUE, t_store = 100
    )
    expect_s4_class(ExampleRicker2, "TreeSummarizedExperiment")
    expect_type(ExampleRicker2@assays@data@listData[["counts"]], "double")
    expect_equal(dim(ExampleRicker@assays@data@listData[["counts"]]), c(10, 100))

    # check simulateRicker with errors in inputs
    expect_error(ErrorRicker1 <- simulateRicker(10, A, x0 = runif(9), t_end = 100))
    expect_error(ErrorRicker2 <- simulateRicker(9, A, t_end = 100))
    expect_error(ErrorRicker3 <- simulateRicker(10, A, carrying_capacities = runif(9), t_end = 100))

    expect_error(ErrorLogistic2 <- simulateStochasticLogistic(
        n_species = 3, b = c(0.2, 0.1), carrying_capacities = c(1000, 2000),
        dr = c(0.001, 0.0015), x0 = c(3, 1)
    ))
    expect_error(ErrorLogistic3 <- simulateStochasticLogistic(
        n_species = 4, partial = FALSE, stochastic = 1
    ))
})
