## show, dbconn.

test_that("show,CompoundDb works", {
    expect_output(show(cmp_db))
})

test_that("dbconn,CompoundDb works", {
    expect_true(!is.null(dbconn(cmp_db)))
    expect_true(is(dbconn(cmp_db), "DBIConnection"))
})
