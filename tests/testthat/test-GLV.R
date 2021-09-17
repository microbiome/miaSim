#usethis::use_r("simulateGLV")
test_that("simulateGLV", {
    # check powerlawA
    expect_error(miaSim:::powerlawA(),
               'argument "n" is missing')
    expect_error(miaSim:::powerlawA(n = 5),
               'argument "alpha" is missing')
    expect_true(is.numeric(miaSim:::powerlawA(6,4)))

    #check dxdt
    expect_error(miaSim:::dxdt(x),
                'argument "parameters" is missing')

    #check simulateGLV
    SEobject <- miaSim:::simulateGLV(N = 4, A = powerlawA(n = 4, alpha = 2), tend = 1000)
    InterMatx <- assay(SEobject)
    expect_true(is.matrix(InterMatx))
    expect_equal(dim(InterMatx), c(4,1000))
    expect_s4_class(SEobject, "SummarizedExperiment")
})
