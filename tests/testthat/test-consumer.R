test_that("simulateConsumerResource", {
    ExampleCR <- simulateConsumerResource(
        n_species = 2,
        n_resources = 4
    )
    expect_s4_class(ExampleCR, "TreeSummarizedExperiment")
    expect_equal(dim(ExampleCR@assays@data@listData[["counts"]]), c(2, 1000))
})
