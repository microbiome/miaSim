test_that("simulateStochasticLogistic", {
    # check simulateStochasticLogistic with 5 species and default inputs
    tse <- simulateStochasticLogistic(n_species = 5)
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(dim(tse), c(5, 1000))    

    # check simulateStochasticLogistic with custom inputs
    tse2 <- simulateStochasticLogistic(
        n_species = 2,
	growth_rates = c(0.2, 0.1),
	carrying_capacities = c(100, 2000),
        death_rates = c(0.001, 0.0015),
	x0 = c(3, 0.1),
        t_start = 0,
	t_end = 500,
	t_step = 0.01,
        t_store = 200,
	stochastic = TRUE)
	
    expect_s4_class(tse2, "TreeSummarizedExperiment")
    expect_equal(dim(tse2), c(2, 200))    

    # check simulateStochasticLogistic with errors in inputs
    expect_error(ErrorLogistic1 <- simulateStochasticLogistic(n_species = 4.1))
    expect_error(ErrorLogistic2 <- simulateStochasticLogistic(
        n_species = 3, b = c(0.2, 0.1), k = c(100, 2000),
        dr = c(0.001, 0.0015), x = c(3, 1)
    ))
    expect_error(ErrorLogistic3 <- simulateStochasticLogistic(
        n_species = 4, partial = FALSE, stochastic = 1
    ))
})
