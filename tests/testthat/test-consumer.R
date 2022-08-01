test_that("simulateConsumerResource", {
    ExampleCR <- simulateConsumerResource(
        n_species = 2,
        n_resources = 4
    )
    expect_type(ExampleCR, "list")
    expect_equal(dim(ExampleCR$matrix), c(1000, 3))
    expect_equal(dim(ExampleCR$resources), c(1000, 5))
})
