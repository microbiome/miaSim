test_that("simulateStochasticLogistic", {
    #check simulateStochasticLogistic with 5 species and default inputs
    ExampleLogistic <- simulateStochasticLogistic(n.species = 5)
    expect_type(assay(ExampleLogistic), "double")
    expect_equal(dim(assay(ExampleLogistic)), c(5,1000))
    expect_s4_class(ExampleLogistic, "SummarizedExperiment")

    #check simulateStochasticLogistic with custom inputs
    ExampleLogistic2 <- simulateStochasticLogistic(
        n.species = 2, b = c(0.2, 0.1), k = c(1000, 2000),
        dr = c(0.001, 0.0015), x = c(3, 0.1),
        t.start = 0, t.end = 1500, t.step = 0.01,
        t.store = 1200, partial = FALSE, stochastic = FALSE)
    expect_type(assay(ExampleLogistic2), "double")
    expect_equal(dim(assay(ExampleLogistic2)), c(2,150001))
    expect_s4_class(ExampleLogistic2, "SummarizedExperiment")

    #check simulateStochasticLogistic with errors in inputs
    expect_error(ErrorLogistic1 <- simulateStochasticLogistic(n.species = 4.1))
    expect_error(ErrorLogistic2 <- simulateStochasticLogistic(
        n.species = 3, b = c(0.2, 0.1), k = c(1000, 2000),
        dr = c(0.001, 0.0015), x = c(3, 1)))
    expect_error(ErrorLogistic3 <- simulateStochasticLogistic(
        n.species = 4, partial = FALSE, stochastic = 1))
})
