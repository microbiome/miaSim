test_that("randomA", {

    expect_true(is.numeric(miaSim:::randomA(5, connectance = 0.02)))
    expect_type(miaSim:::randomA(5, connectance = 0.02), "double")

    low_inter_A <- randomA(10, connectance = 0.01, symmetric = TRUE)
    expect_type(low_inter_A, "double")
})
