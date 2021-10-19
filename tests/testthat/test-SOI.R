test_that("simulateSOI", {
  SOI <- simulateSOI(n.species = 10, I = 1000, A = powerlawA(n.species = 10,
          alpha = 1.2), k=5, com = NULL, tend = 150, norm = TRUE)

  SOI_matrix <- assay(SOI)
  expect_type(SOI_matrix, "double")
  expect_equal(dim(SOI_matrix), c(10,150))
  expect_s4_class(SOI, "SummarizedExperiment")

  #check if(length(com) == n.species+1)
  SOI2 <- simulateSOI(n.species = 10, I = 1000, A = powerlawA(n.species = 10,
          alpha = 1.2), k=5, com = c(1:11) , tend = 150, norm = TRUE)

  #check if(length(com) == n.species)
  SOI3 <- simulateSOI(n.species = 10, I = 1000, A = powerlawA(n.species = 10,
          alpha = 1.2), k=5, com = c(1:10) , tend = 150, norm = TRUE)

  #check else
  expect_error(miaSim:::simulateSOI(n.species = 10, I = 1000,
      A = powerlawA(n.species = 10, alpha = 1.2), k=5, com = c(1:9),
      tend = 150, norm = TRUE), force = TRUE)

})
