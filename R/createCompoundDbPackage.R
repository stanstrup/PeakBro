#' @title Create a CompoundDb database
#'
#' @description
#' 
#' `createCompoundDb` creates a `SQLite`-based `CompoundDb` object/database
#' from a compound resource provided as a `data.frame` or `tbl`.
#' An additional `data.frame` providing metadata information is mandatory.
#' Required columns for the `data.frame` providing the compound information are:
#'
#' 
#' The metadata `data.frame` is supposed to have two columns named `"name"` and
#' `"value"` providing the following minimal information as key-value pairs:
#' + `"source"`: the source from which the data was retrieved (e.g. `"HMDB"`).
#' + `"url"`: the url from which the original data was retrieved.
#' + `"source_version"`: the version from the original data source
#'   (e.g. `"v4"`).
#' + `"source_date"`: the date when the original data source was generated.
#' + `"organism"`: the organism. Should be in the form `"Hsapiens"` or
#'   `"Mmusculus"`.
#'
#' @details
#'
#' Metadata information is also used to create the file name for the database
#' file. The name starts with `"CompoundDb"`, followed by the organism, the
#' data source and its version. A compound database file for HMDB version 4
#' with human metabolites will thus be named: `"CompoundDb.Hsapiens.HMDB.v4"`.
#' 
#' @param x `data.frame` or `tbl` with the compound annotations. See
#'     description for details.
#'
#' @param metadata `data.frame` with metadata information. See description for
#'     details.
#'
#' @param path `character(1)` with the path to the directory where the database
#'     file should be written. Defaults to the current directory.
#'
#' @return A `character` with the database name (invisibly).
#' 
#' @importFrom DBI dbDriver dbWriteTable dbExecute dbDisconnect
#' @importFrom RSQLite dbConnect
#'
#' @export
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Read compounds for a HMDB subset
#' fl <- system.file("extdata/hmdb/hmdb_sub.xml", package = "PeakABro")
#' cmps <- generate_hmdb_tbl(fl)
#'
#' ## Create a metadata data.frame for the compounds.
#' metad <- data.frame(name = c("source", "url", "source_version",
#'     "source_date", "organism"), value = c("HMDB", "http://www.hmdb.ca",
#'     "v4", "2017-08-27", "Hsapiens"))
#'
#' ## Create a SQLite database in the temporary folder
#' db_f <- createCompoundDb(cmps, metadata = metad, path = tempdir())
#'
#' ## connect to the database and query it's tables using RSQlite
#' library(RSQLite)
#' con <- dbConnect(dbDriver("SQLite"), db_f)
#'
#' dbGetQuery(con, "select * from metadata")
#' dbGetQuery(con, "select * from compound")
createCompoundDb <- function(x, metadata, path = ".") {
    .valid_metadata(metadata)
    .valid_compound(x)
    db_file <- paste0(path, "/", .db_file_from_metadata(metadata), ".sqlite")
    con <- dbConnect(dbDriver("SQLite"), dbname = db_file)
    ## Add additional metadata info
    metadata <- rbind(metadata, c("db_creation_date", date()))
    metadata <- rbind(metadata, c("supporting_package", "PeakABro"))
    metadata <- rbind(metadata, c("supporting_object", "CompoundDb"))
    dbWriteTable(con, name = "metadata", metadata, row.names = FALSE)
    dbWriteTable(con, name = "compound", x, row.names = FALSE)
    ## Creating indices
    dbExecute(con, paste0("create index compound_id_idx on compound (id)"))
    dbExecute(con, paste0("create index compound_name_idx on compound (name)"))
    dbDisconnect(con)
    invisible(db_file)
}

.required_metadata_keys <- c("source", "url", "source_version", "source_date",
                             "organism")
.required_compound_columns <- c("id", "name", "inchi", "formula", "mass")

#' @description Create the database file name from the metadata `data.frame`.
#'     The function checks also for the presence of all required metadata fields
#'     ensuring that these are also not `NA` or `NULL`.
#'
#' @noRd
.db_file_from_metadata <- function(x) {
    paste0("CompoundDb.", x$value[x$name == "organism"], ".",
           x$value[x$name == "source"], ".",
           x$value[x$name == "source_version"])
}

#' @description Check the metadata data.frame for required columns.
#'
#' @noRd
.valid_metadata <- function(metadata, error = TRUE) {
    txt <- character()
    if (!is.data.frame(metadata))
        txt <- c(txt, "'metadata' is expected to be a data.frame")
    if (all(c("name", "value") %in% colnames(metadata))) {
        keys <- metadata$name
        vals <- metadata$value
        got_it <- .required_metadata_keys %in% keys
        if (!all(got_it))
            txt <- c(txt, paste0("required fields ",
                                 paste0(.required_metadata_keys[!got_it],
                                        collapse = ", "),
                                 " not found in metadata data.frame"))
        vals <- vals[keys %in% .required_metadata_keys]
        if (length(vals))
            if (any(is.na(vals)) | any(length(vals) == 0))
                txt <- c(txt, paste0("values for metadata data.frame fields ",
                                     "should not be empty or NA"))
    } else {
        txt <- c(txt, paste0("metadata data.frame needs to have columns ",
                             "named 'name' and 'value'"))
    }
    if (length(txt))
        if (error)
            stop(paste(txt, collapse = "\n"))
        else txt
    else TRUE
}

#' @description Check that the compounds table contains all required data.
#'
#' @noRd
.valid_compound <- function(x, error = TRUE) {
    txt <- character()
    if (!is.data.frame(x))
        txt <- c(txt, "'x' is supposed to be a data.frame")
    got_it <- .required_compound_columns %in% colnames(x)
    if (!all(got_it)) {
        txt <- c(txt, paste0("Miss required columns: ",
                             paste0(.required_compound_columns[!got_it],
                                    collapse = ", ")))
    } else {
        if (!is.numeric(x$mass))
            txt <- c(txt, "Column 'mass' should be numeric")
    }
    if (length(txt))
        if (error)
            stop(paste(txt, collapse = "\n"))
        else txt
    else TRUE
}
