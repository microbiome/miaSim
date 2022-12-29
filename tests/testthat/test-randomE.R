test_that("randomE", {
    x <- randomE(n_species = 5, n_resources = 12)
    expect_type(x, "double")
    expect_equal(dim(x), c(5, 12))
})
