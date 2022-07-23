test_that("convertToTreeSE", {
    ExampleHubbellRates <- simulateHubbellRates(community_initial = c(0,5,10),
                                              migration_p = 0.01,
                                              metacommunity_p = NULL,
                                              k_events = 1,
                                              growth_rates = NULL,
                                              norm = FALSE,
                                              t_end=1000)

    HubbellSE <- convertToTreeSE(assay = ExampleHubbellRates$counts)
    expect_s4_class(HubbellSE, "TreeSummarizedExperiment")

})