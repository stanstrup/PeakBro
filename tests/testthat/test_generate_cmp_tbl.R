## Unit tests for generate_cmp_tbl.R

test_that("simple_parse_hmdb_xml works", {
    ## On a single XML file downloaded from HMDB
    x <- system.file("extdata/hmdb/HMDB0000001.xml", package = "PeakABro")
    res <- PeakABro:::simple_parse_hmdb_xml(x)
    expect_equal(nrow(res), 1)
    expect_equal(res$id, c("HMDB0000001"))
    expect_equal(colnames(res), c("id", "name", "inchi", "formula", "mass"))
    
    x <- system.file("extdata/peaklist_pos.tsv", package = "PeakABro")
    expect_error(PeakABro:::simple_parse_hmdb_xml(x))
    
    ## On a subset of the main XML file
    x <- system.file("extdata/hmdb/hmdb_sub.xml", package = "PeakABro")
    res <- PeakABro:::simple_parse_hmdb_xml(x)
    expect_equal(nrow(res), 3)
    expect_equal(res$id, c("HMDB0000001", "HMDB0000002", "HMDB0000005"))
    expect_equal(colnames(res), c("id", "name", "inchi", "formula", "mass"))
})

test_that("generate_hmdb_tbl works", {
    ## Input files
    one <- system.file("extdata/hmdb/HMDB0000001.xml", package = "PeakABro")
    more <- system.file("extdata/hmdb/hmdb_sub.xml", package = "PeakABro")
    ## Exceptions
    suppressWarnings(
        expect_error(generate_hmdb_tbl("I do not exist"))
    )
    suppressWarnings(
        expect_error(generate_hmdb_tbl(c("I do not exist", one)))
    )
    ## Read one file
    res <- generate_hmdb_tbl(one)
    expect_true(is.data.frame(res))
    expect_true(any(class(res) == "tbl"))
    expect_true(nrow(res) == 1)
    expect_equal(colnames(res), c("id" ,"name", "inchi",
                                  "formula", "mass"))
    res <- generate_hmdb_tbl(more)
    expect_true(is.data.frame(res))
    expect_true(any(class(res) == "tbl"))
    expect_true(nrow(res) == 3)
    expect_equal(colnames(res), c("id" ,"name", "inchi",
                                  "formula", "mass"))
    ## Read multiple files
    res <- generate_hmdb_tbl(c(one, more))
    expect_true(is.data.frame(res))
    expect_true(any(class(res) == "tbl"))
    expect_true(nrow(res) == 4)
    expect_equal(colnames(res), c("id" ,"name", "inchi",
                                  "formula", "mass"))

    ## SDF format
    sdf <- system.file("extdata/hmdb/hmdb_sub.sdf", package = "PeakABro")
    res_xml <- generate_hmdb_tbl(more)
    res_sdf <- generate_hmdb_tbl(sdf)
    have_id <- res_sdf$id %in% res_xml$id
    expect_equal(res_xml, res_sdf[have_id, ])
})

test_that("parse_hmdb_sdf works", {
    sdf <- system.file("extdata/hmdb/hmdb_sub.sdf", package = "PeakABro")
    res <- PeakABro:::parse_hmdb_sdf(sdf)
    expect_true(is.numeric(res$mass))
})
