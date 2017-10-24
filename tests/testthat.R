library("testthat")
library("PeakABro")

## Create a test CompoundDb that can be used across tests.
metadata <- data.frame(name = c("source", "url", "source_version",
                                "source_date", "organism"),
                       value = c("HMDB", "http://www.hmdb.ca", "4", "2017",
                                 "Hsapiens"),
                       stringsAsFactors = FALSE)
fl <- system.file("extdata/hmdb/hmdb_sub.xml", package = "PeakABro")
cmps <- generate_hmdb_tbl(fl)
db_f <- createCompoundDb(cmps, metadata = metadata, path = tempdir())
cmp_db <- CompoundDb(db_f)

test_check("PeakABro")
