test_that("powerlawA", {
    expect_true(is.numeric(miaSim:::powerlawA(6, 4)))
    expect_type(miaSim:::powerlawA(6, 4), "double")

    expect_error(Error1 <- powerlawA(n_species = 3.2))
    expect_error(Error1 <- powerlawA(n_species = 3, alpha = "nonnumeric"))

    A_strong <- powerlawA(n_species = 10, alpha = 1.01, symmetric = TRUE)
    expect_type(A_strong, "double")
})
