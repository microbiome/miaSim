test_that("simulateGLV", {
    #check dxdt
    expect_error(miaSim:::dxdt(x),
                'argument "parameters" is missing')

    #check simulateGLV
    A <- miaSim::powerlawA(4, alpha = 1.01)
    simulateGLV(n.species = 4, A, t.start = 0, t.store = 1000)
    SEobject <- simulateGLV(n.species = 4, A, t.start = 0,
                                     t.store = 1000)

    expect_type(SEobject, "double")
    expect_equal(dim(SEobject), c(4,1000))

    expect_error(Error1 <- simulateGLV(n.species = 0.5))
    expect_error(Error2 <- simulateGLV(n.species = 5, A = 2, x = 2, b = 3))


    #check norm = TRUE
    SEobject2 <- miaSim:::simulateGLV(n.species = 4,
                A = powerlawA(n.species = 4, alpha = 2), t.start = 0,
                t.store = 1000, norm = TRUE)
    expect_type(SEobject2, "double")
    expect_equal(dim(SEobject2), c(4,1000))
})
