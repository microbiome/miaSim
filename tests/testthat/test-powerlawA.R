test_that("powerlawA", {
expect_error(miaSim:::powerlawA(n.species = 5),
             'argument "alpha" is missing')
expect_true(is.numeric(miaSim:::powerlawA(6,4)))
expect_type(miaSim:::powerlawA(6,4), "double")

})
