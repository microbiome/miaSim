test_that("simulateGLV", {
    #check dxdt
    expect_error(miaSim:::dxdt(x),
                'argument "parameters" is missing')

    #check simulateGLV
    A <- miaSim::powerlawA(4, alpha = 1.01)
    simulateGLV(n_species = 4, A, t_start = 0, t_store = 1000)
    SEobject <- simulateGLV(n_species = 4, A, t_start = 0,
                                     t_store = 1000)

    expect_type(SEobject, "double")
    expect_equal(dim(SEobject), c(4,1000))

    expect_error(Error1 <- simulateGLV(n_species = 0.5))
    expect_error(Error2 <- simulateGLV(n_species = 5, A = 2, x = 2, b = 3))


    #check norm = TRUE
    SEobject2 <- miaSim:::simulateGLV(n_species = 4,
                A = powerlawA(n_species = 4, alpha = 2), t_start = 0,
                t_store = 1000, norm = TRUE)
    expect_type(SEobject2, "double")
    expect_equal(dim(SEobject2), c(4,1000))
})
