test_that("simulateConsumerResource", {

  ExampleConsumerResource <- simulateConsumerResource(n.species = 2,
                                                      n.resources = 4)
  A_Matrix <- assay(ExampleConsumerResource)
  expect_type(A_Matrix, "double")
  expect_equal(dim(A_Matrix), c(6,1000))
  expect_s4_class(ExampleConsumerResource, "SummarizedExperiment")

  #check return.matrix = TRUE
  ExampleConsumerResource2 <- simulateConsumerResource(n.species = 2,
                                                      n.resources = 4,
                                                      return.matrix = TRUE)
})
