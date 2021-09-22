test_that("simulateGLV", {
    #check dxdt
    expect_error(miaSim:::dxdt(x),
                'argument "parameters" is missing')

    #check simulateGLV
    SEobject <- miaSim:::simulateGLV(N = 4, A = powerlawA(n = 4, alpha = 2), tend = 1000)
    InterMatx <- assay(SEobject)
    expect_true(is.matrix(InterMatx))
    expect_equal(dim(InterMatx), c(4,1000))
    expect_s4_class(SEobject, "SummarizedExperiment")

    #check norm = TRUE
    SEobject2 <- miaSim:::simulateGLV(N = 4, A = powerlawA(n = 4, alpha = 2), tend = 1000, norm = TRUE)
    InterMatx2 <- assay(SEobject)
    expect_true(is.matrix(InterMatx2))
    expect_equal(dim(InterMatx2), c(4,1000))
})
