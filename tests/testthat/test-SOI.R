test_that("simulateSOI", {
  SOI <- simulateSOI(n_species = 10, I = 1000, A = powerlawA(n_species = 10,
          alpha = 1.2), k=5, com = NULL, tend = 150, norm = TRUE)

  expect_type(SOI, "double")
  expect_equal(dim(SOI), c(10,150))

  #check if(length(com) == n_species+1)
  SOI2 <- simulateSOI(n_species = 10, I = 1000, A = powerlawA(n_species = 10,
          alpha = 1.2), k=5, com = c(1:11) , tend = 150, norm = TRUE)

  #check if(length(com) == n_species)
  SOI3 <- simulateSOI(n_species = 10, I = 1000, A = powerlawA(n_species = 10,
          alpha = 1.2), k=5, com = c(1:10) , tend = 150, norm = TRUE)

  expect_error(Error1 <- simulateSOI(n_species = 0.5, I = 0.6, k = 9.7, tend = 80))

  #check else
  expect_error(miaSim:::simulateSOI(n_species = 10, I = 1000,
      A = powerlawA(n_species = 10, alpha = 1.2), k=5, com = c(1:9),
      tend = 150, norm = TRUE), force = TRUE)

})
