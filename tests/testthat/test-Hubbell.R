test_that("simulateHubbell", {
  ExampleHubbell <- simulateHubbell(N = 8, M = 10, I = 1000, d = 50,
                                        m = 0.02, tend = 100)

  AbundanceM <- assay(ExampleHubbell)
  expect_true(is.matrix(AbundanceM))
  expect_equal(dim(AbundanceM), c(10,100))
  expect_s4_class(ExampleHubbell, "SummarizedExperiment")
})