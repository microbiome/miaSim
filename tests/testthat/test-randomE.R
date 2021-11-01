test_that("randomE", {

  ExampleEfficiencyMatrix2 <- randomE(n.species = 5, n.resources = 12)

  expect_type(ExampleEfficiencyMatrix2, "double")
  expect_equal(dim(ExampleEfficiencyMatrix2), c(5,12))

  expect_error(E3 <- randomE(n.species = 12, n.resources = -1))

})
