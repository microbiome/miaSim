test_that("randomE", {

  ExampleEfficiencyMatrix2 <- randomE(n_species = 5, n_resources = 12)

  expect_type(ExampleEfficiencyMatrix2, "double")
  expect_equal(dim(ExampleEfficiencyMatrix2), c(5,12))

})
