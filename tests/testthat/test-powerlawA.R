test_that("powerlawA", {

expect_true(is.numeric(miaSim:::powerlawA(6,4)))
expect_type(miaSim:::powerlawA(6,4), "double")
})
