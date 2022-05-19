test_that("simulateGLV", {

    #check simulateGLV
    set.seed(533)
    A <- miaSim::powerlawA(4, alpha = 1.01)
    simulateGLV(n_species = 4, A = A, t_start = 0, t_store = 1000)
    SEobject <- simulateGLV(n_species = 4, A = A, t_start = 0,
                                     t_store = 1000)

    expect_type(SEobject, "list")
    expect_equal(dim(t(SEobject$matrix)), c(5, 1000))

    expect_error(Error1 <- simulateGLV(n_species = 0.5))
    expect_error(Error2 <- simulateGLV(n_species = 5, A = 2, x = 2, b = 3))


    #check norm = TRUE
    set.seed(53)
    SEobject2 <- miaSim:::simulateGLV(n_species = 3,
                A = powerlawA(n_species = 3, alpha = 2), t_start = 0,
                t_store = 1000, norm = TRUE)
    expect_type(SEobject2$matrix, "double")
    expect_equal(dim(t(SEobject2$matrix)), c(4,1000))

})
