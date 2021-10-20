test_that("simulateHubbell", {
  ExampleHubbell <- simulateHubbell(n.species = 8, M = 10, I = 1000, d = 50,
                                     m = 0.02, tend = 100)
  AbundanceM <- assay(ExampleHubbell)
  expect_type(AbundanceM, "double")
  expect_equal(dim(AbundanceM), c(10,100))
  expect_s4_class(ExampleHubbell, "SummarizedExperiment")

  #check norm = TRUE
  ExampleHubbell_2 <- simulateHubbell(n.species = 8, M = 10, I = 1000, d = 50,
                                      m = 0.02, tend = 100, norm = TRUE)
})
