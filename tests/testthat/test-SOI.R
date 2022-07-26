test_that("simulateSOI", {
  SOI <- simulateSOI(n_species = 10, carrying_capacity = 1000, A = powerlawA(n_species = 10,
          alpha = 1.2), k_events=5, x0 = NULL, t_end = 150, norm = TRUE)

  expect_type(SOI$matrix, "double")
  expect_equal(dim(SOI$matrix), c(150, 11))

  #check if(length(x0) == n_species+1)
  SOI2 <- simulateSOI(n_species = 10, carrying_capacity = 1000, A = powerlawA(n_species = 10,
          alpha = 1.2), k_events=5, x0 = c(1:11) , t_end = 150, norm = TRUE)

  #check if(length(x0) == n_species)
  SOI3 <- simulateSOI(n_species = 10, carrying_capacity = 1000, A = powerlawA(n_species = 10,
          alpha = 1.2), k_events=5, x0 = c(1:10) , t_end = 150, norm = TRUE)

  expect_error(Error1 <- simulateSOI(n_species = 0.5, carrying_capacity = 0.6, k_events = 9.7, t_end = 80))

  #check else
  expect_error(miaSim:::simulateSOI(n_species = 10, carrying_capacity = 1000,
      A = powerlawA(n_species = 10, alpha = 1.2), k_events=5, x0 = c(1:9),
      t_end = 150, norm = TRUE), force = TRUE)

})
