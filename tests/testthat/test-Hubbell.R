test_that("simulateHubbell", {
    tse <- simulateHubbell(
        n_species = 8, M = 10,
        carrying_capacity = 1000, k_events = 50, migration_p = 0.02,
        t_end = 100
    )

    expect_s4_class(tse, "TreeSummarizedExperiment")
    #expect_type(tse@assays@data@listData[["counts"]], "double")
    expect_type(assay(tse, "counts"), "double")    
    expect_equal(dim(tse), c(8,100))

    expect_error(Error1 <- simulateHubbell(
        n_species = 3.4, M = 7.9,
        carrying_capacity = 6.5, tskip = 4.3, t_end = 100
    ))

    # check norm = TRUE
    tse2 <- simulateHubbell(
        n_species = 8, M = 10,
        carrying_capacity = 1000, k_events = 50, migration_p = 0.02,
        t_end = 100, norm = TRUE
    )
})
