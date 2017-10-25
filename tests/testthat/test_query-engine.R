test_that(".reduce_tables works", {
    tabs <- list(gene = c("gene_name", "gene_id", "redundant_field"),
                 tx = c("tx_id", "gene_id", "tx_name", "redundant_field",
                        "tx_start"),
                 exon = c("exon_id", "redundant_field", "tx_id"))
    
    res <- .reduce_tables(tabs, columns = c("tx_id", "gene_id"))
    expect_equal(length(res), 1)
    expect_equal(res, list(tx = c("tx_id", "gene_id")))

    res <- .reduce_tables(tabs, columns = c("gene_id", "gene_name"))
    expect_equal(length(res), 1)
    expect_equal(res, list(gene = c("gene_id", "gene_name")))

    res <- .reduce_tables(tabs, columns = c("gene_name", "exon_id",
                                            "redundant_field",
                                            "tx_id"))
    expect_equal(length(res), 2)
    expect_equal(res, list(exon = c("exon_id", "redundant_field", "tx_id"),
                           gene = "gene_name"))
})

test_that(".prefix_columns works", {
    res <- .prefix_columns(list(a = c("b", "c")))
    expect_equal(res, list(a = c("a.b", "a.c")))
    res <- .prefix_columns(list(a = c("b", "c"),
                                b = c("d", "e")))
    expect_equal(res, list(a = c("a.b", "a.c"),
                           b = c("b.d", "b.e")))
})

test_that(".from works", {
    expect_equal(.from("a"), "from a")
    expect_error(.from(c("a", "b")))
})

test_that(".where works", {
    expect_equal(.where(), NULL)
    expect_error(.where("something"))
})

test_that(".select works", {
    expect_equal(.select(c("a", "b")), "select a,b")
    expect_error(.select())
})

test_that(".build_query_CompoundDb works", {
    res <- PeakABro:::.build_query_CompoundDb(
                          cmp_db, columns = c("id", "inchi"))
    expect_equal(res, "select compound.id,compound.inchi from compound")
    expect_error(PeakABro:::.build_query_CompoundDb(
                                cmp_db, columns = c("od", "inchi")))
})
