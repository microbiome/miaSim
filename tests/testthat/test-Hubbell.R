test_that("simulateHubbell", {
    ExampleHubbell <- simulateHubbell(
        n_species = 8, M = 10,
        carrying_capacity = 1000, k_events = 50, migration_p = 0.02,
        t_end = 100
    )

    expect_type(ExampleHubbell$matrix, "double")
    expect_equal(dim(ExampleHubbell$matrix), c(100, 11))

    expect_error(Error1 <- simulateHubbell(
        n_species = 3.4, M = 7.9,
        carrying_capacity = 6.5, tskip = 4.3, t_end = 100
    ))

    # check norm = TRUE
    ExampleHubbell_2 <- simulateHubbell(
        n_species = 8, M = 10,
        carrying_capacity = 1000, k_events = 50, migration_p = 0.02,
        t_end = 100, norm = TRUE
    )
})
