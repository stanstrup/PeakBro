
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

    ## .tables
    tbls <- PeakABro:::.tables(cmp)
    expect_equal(length(tbls), 1)
    expect_equal(names(tbls), "compound")
    tbls <- PeakABro:::.tables(cmp, metadata = TRUE)
    expect_equal(length(tbls), 2)
    expect_equal(names(tbls), c("compound", "metadata"))
    tbls <- PeakABro:::.tables(cmp, name = "not_there")
    expect_equal(length(tbls), 1)
    
    ## .get_property
    prps <- PeakABro:::.get_property(cmp, "tables")
    expect_equal(prps, PeakABro:::.tables(cmp, metadata = TRUE))
    prps <- PeakABro:::.get_property(cmp, "not_there")
    expect_equal(prps, NULL)
})

test_that("compounds works", {
    cmps <- compounds(cmp_db)
    expect_true(is(cmps, "data.frame"))
    cmps_tbl <- compounds(cmp_db, columns = c("id", "name"),
                          return.type = "tibble")
    expect_true(is(cmps_tbl, "tbl"))
    expect_equal(colnames(cmps_tbl), c("id", "name"))

    expect_error(compounds(cmp_db, filter = "something"))
})

test_that("src_compound works", {
    src_cmp <- src_compdb(cmp_db)
    expect_true(is(src_cmp, "src_dbi"))
})
