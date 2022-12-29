test_that("randomA", {
    expect_true(is.numeric(miaSim:::randomA(5, connectance = 0.02)))
    expect_type(miaSim:::randomA(5, connectance = 0.02), "double")
    expect_error(Error1 <- randomA(n_species = 0.5))
    expect_error(Error2 <- randomA(n_species = 5, d = "5.8", min_strength = "nonnumeric"))
    low_inter_A <- randomA(10, connectance = 0.01, symmetric = TRUE)
    expect_type(low_inter_A, "double")
})
