test_that("simulateHubbell", {
  ExampleHubbell <- simulateHubbell(n.species = 8, M = 10, I = 1000, d = 50,
                                        m = 0.02, tend = 100)

  expect_type(ExampleHubbell, "double")
  expect_equal(dim(ExampleHubbell), c(10,100))

  expect_error(Error1 <- simulateHubbell(n.species = 3.4, M = 7.9, I = 6.5, tskip = 4.3, tend = 100))

  #check norm = TRUE
  ExampleHubbell_2 <- simulateHubbell(n.species = 8, M = 10, I = 1000, d = 50,
                                    m = 0.02, tend = 100, norm = TRUE)
})
