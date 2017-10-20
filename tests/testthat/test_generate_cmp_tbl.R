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
