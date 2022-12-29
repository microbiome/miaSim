test_that("simulateConsumerResource", {
    tse <- simulateConsumerResource(
        n_species = 2,
        n_resources = 4
    )
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(dim(tse), c(2, 1000))    
})
