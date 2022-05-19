test_that("convertToSE", {
    ExampleHubbellRates <- simulateHubbellRates(x0 = c(0,5,10),
                                              migration_p = 0.01,
                                              metacommunity_probability = NULL,
                                              k_events = 1,
                                              growth_rates = NULL,
                                              norm = FALSE,
                                              t_end=1000)

    HubbellSE <- convertToSE(assay = t(ExampleHubbellRates$matrix[,1:3]))

    expect_s4_class(HubbellSE, "SummarizedExperiment")

    HubbellTSE <- convertToSE(assay = t(ExampleHubbellRates$matrix[,1:3]), output = TSE)

    expect_s4_class(HubbellTSE, "TreeSummarizedExperiment")

})