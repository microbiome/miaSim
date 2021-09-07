usethis::use_r("simulateSOI")
test_that("simulateSOI", {
  SOI <- simulateSOI(N = 10, I = 1000, A = powerlawA(n = 10, alpha = 1.2),
                     k=5, com = NULL, tend = 150, norm = TRUE)
  # check powerlawA
  expect_error(miaSim:::powerlawA(),
               'argument "n" is missing')
  expect_error(miaSim:::powerlawA(n = 10),
               'argument "alpha" is missing')
  expect_true(is.numeric(miaSim:::powerlawA(10, 1.2)))
  expect_false(is.integer(miaSim:::powerlawA(10, 1.2)),
               'arguments must be integers')

  SOI_matrix <- assay(SOI)
  expect_true(is.matrix(SOI_matrix))
  expect_equal(dim(SOI_matrix), c(10,150))
  expect_s4_class(SOI, "SummarizedExperiment")
})
