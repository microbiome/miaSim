usethis::use_r("simulateHubbell")
test_that("simulateHubbell", {
  ExampleHubbell <- simulateHubbell(N = 8, M = 10, I = 1000, d = 50,
                                        m = 0.02, tend = 100)

  AbundanceM <- assay(ExampleHubbell)
  expect_true(is.matrix(AbundanceM))
  expect_equal(dim(AbundanceM), c(10,100))
<<<<<<< HEAD
  expect_s4_class(ExampleHubbell, "SummarizedExperiment")


=======
  expect_true(class(ExampleHubbell)== "SummarizedExperiment")
>>>>>>> aaf6f49c85378edc57fe0174245804b61a0dcdae
})
