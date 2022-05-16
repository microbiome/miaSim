test_that("convertToSE", {
    ExampleHubbellRates <- simulateHubbellRates(community_initial = c(0,5,10),
                                              migration_p = 0.01,
                                              metacommunity_p = NULL,
                                              k_events = 1,
                                              growth_rates = NULL,
                                              norm = FALSE,
                                              t_end=1000)

    HubbellSE <- convertToSE(assay = ExampleHubbellRates$counts)

    expect_s4_class(HubbellSE, "SummarizedExperiment")

    HubbellTSE <- convertToSE(assay = ExampleHubbellRates$counts, output = TSE)

    expect_s4_class(HubbellTSE, "TreeSummarizedExperiment")

})