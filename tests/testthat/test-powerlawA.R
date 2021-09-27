test_that("powerlawA", {

expect_true(is.numeric(miaSim:::powerlawA(6,4)))
expect_type(miaSim:::powerlawA(6,4), "double")

A_strong <- powerlawA(n.species = 10, alpha = 1.01, symmetric = TRUE)
expect_type(A_strong, "double")

})
