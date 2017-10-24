
test_that("CompoundDb constructor and low level functions", {
    expect_error(CompoundDb())
    cmp <- new("CompoundDb")
    expect_true(is.null(PeakABro:::.dbconn(cmp)))
    expect_true(is.null(dbconn(cmp)))

    ## Create a simple database from the internal files.
    metadata <- data.frame(name = c("source", "url", "source_version",
                                   "source_date", "organism"),
                           value = c("HMDB_tmp", "http://www.hmdb.ca", "4",
                                     "2017", "Hsapiens"),
                           stringsAsFactors = FALSE)
    fl <- system.file("extdata/hmdb/hmdb_sub.xml", package = "PeakABro")
    cmps <- generate_hmdb_tbl(fl)
    db_f <- createCompoundDb(cmps, metadata = metadata, path = tempdir())
    cmp <- CompoundDb(db_f)

    expect_true(!is.null(PeakABro:::.dbconn(cmp)))
    expect_true(PeakABro:::.validCompoundDb(dbconn(cmp)))
    res <- PeakABro:::.metadata(cmp)
    expect_equal(metadata, res[1:nrow(metadata), ])
    res <- PeakABro:::.metadata(dbconn(cmp))
    expect_equal(metadata, res[1:nrow(metadata), ])
    res <- PeakABro:::.metadata_value(cmp, "organism")
    expect_equal(res, "Hsapiens")
    res <- PeakABro:::.metadata_value(dbconn(cmp), "source")
    expect_equal(res, "HMDB_tmp")
})
