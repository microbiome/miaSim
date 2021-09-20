test_that("simulateSOI", {
  SOI <- simulateSOI(N = 10, I = 1000, A = powerlawA(n = 10, alpha = 1.2),
                     k=5, com = NULL, tend = 150, norm = TRUE)

  SOI_matrix <- assay(SOI)
  expect_true(is.matrix(SOI_matrix))
  expect_equal(dim(SOI_matrix), c(10,150))
  expect_s4_class(SOI, "SummarizedExperiment")
})
