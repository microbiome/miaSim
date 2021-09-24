test_that("randomA", {

    expect_true(is.numeric(miaSim:::randomA(5, 0.02)))
    expect_type(miaSim:::randomA(5, 0.02), "double")
})
