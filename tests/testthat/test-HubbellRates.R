test_that("simulateHubbellRates", {
  HubbellWTime <- simulateHubbellRates(community_initial = c(0,5,10),
                                      migration_p = 0.01, metacommunity_p = NULL, k_events = 1,
                                      growth_rates = NULL, norm = FALSE, t.end=1000)
  Abundance_matrix <- assay(HubbellWTime)
  expect_type(Abundance_matrix, "double")
  expect_equal(dim(Abundance_matrix), c(3,1000))
  expect_s4_class(HubbellWTime, "SummarizedExperiment")

  #check errors in inputs
  expect_error(Error1 <- simulateHubbellRates(community_initial = -1))
  expect_error(Error2 <- simulateHubbellRates(community_initial = c(0,5,10),
                        migration_p = '0.01', metacommunity_p = NULL, k_events = '1',
                        growth_rates = NULL, norm = FALSE, t.end=1000))
  expect_error(Error3 <- simulateHubbellRates(community_initial = c(0,5,10),
                        migration_p = 0.01, metacommunity_p = NULL, k_events = 1,
                         growth_rates = NULL, norm = 'FALSE', t.end=1000))

  #check norm = TRUE
  HubbellWTime <- simulateHubbellRates(community_initial = c(0,5,10),
                                       migration_p = 0.01, metacommunity_p = NULL, k_events = 1,
                                       growth_rates = NULL, norm = TRUE, t.end=1000)
})
