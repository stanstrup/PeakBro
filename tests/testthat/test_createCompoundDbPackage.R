test_that(".valid_metadata works", {
    ## Check errors.
    expect_error(PeakABro:::.valid_metadata("something"))
    expect_true(is.character(
        PeakABro:::.valid_metadata("something", error = FALSE)))
    metadata <- data.frame(name = c("source", "url", "source_version",
                                   "source_date", "organism"),
                           value = c("HMDB", "http://www.hmdb.ca", "4", "2017",
                                     "Hsapiens"),
                           stringsAsFactors = FALSE)
    expect_error(PeakABro:::.valid_metadata(metadata[, 1, drop = FALSE]))
    expect_error(PeakABro:::.valid_metadata(metadata[1:4, ]))
    metadata_fail <- metadata
    metadata_fail[1, 2] <- NA
    expect_error(PeakABro:::.valid_metadata(metadata_fail))

    ## Valid one.
    expect_true(PeakABro:::.valid_metadata(metadata))
})

test_that(".db_file_from_metadata works", {
    metadata <- data.frame(name = c("source", "url", "source_version",
                                    "source_date", "organism"),
                           value = c("HMDB", "http://www.hmdb.ca", "v4", "2017",
                                     "Hsapiens"),
                           stringsAsFactors = FALSE)
    db_file <- PeakABro:::.db_file_from_metadata(metadata)
    expect_equal(db_file, "CompoundDb.Hsapiens.HMDB.v4")
})

test_that(".valid_compound works", {
    cmps <- data.frame(id = c("01", "02"), name = c("a", "b"),
                       inchi = c("i1", "i2"), formula = c("some", "thing"),
                       mass = c(1, 3))
    expect_true(PeakABro:::.valid_compound(cmps))
    ## Errors
    expect_error(PeakABro:::.valid_compound("b"))
    expect_true(is.character(PeakABro:::.valid_compound("b", error = FALSE)))
    expect_error(PeakABro:::.valid_compound(data.frame()))
    expect_error(PeakABro:::.valid_compound(cmps[, 1:3]))
    cmps$mass <- c("1", "2")
    expect_error(PeakABro:::.valid_compound(cmps))
})

test_that("createCompoundDb works", {
    fl <- system.file("extdata/hmdb/hmdb_sub.xml", package = "PeakABro")
    cmps <- generate_hmdb_tbl(fl)

    metad <- data.frame(name = c("source", "url", "source_version",
                                 "source_date", "organism"),
                        value = c("HMDB", "http://www.hmdb.ca",
                                  "v4", "2017-08-27", "Hsapiens"),
                        stringsAsFactors = FALSE)
    db_f <- createCompoundDb(cmps, metadata = metad, path = tempdir())

    library(RSQLite)
    con <- dbConnect(dbDriver("SQLite"), db_f)
    db_metad <- dbGetQuery(con, "select * from metadata")
    db_cmps <- dbGetQuery(con, "select * from compound")    
})
