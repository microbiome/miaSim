test_that("simulateConsumerResource", {

  ExampleConsumerResource <- simulateConsumerResource(n.species = 2,
                                                      n.resources = 4)
  expect_type(ExampleConsumerResource, "double")
  expect_equal(dim(ExampleConsumerResource), c(1000,6))

})
