#' @include createCompoundDbPackage.R

#' @name CompoundDb
#'
#' @title Simple compound (metabolite) databases
#'
#' @aliases CompoundDb-class show,CompoundDb-method dbconn,CompoundDb-method
#'     show dbconn
#' 
#' @description
#'
#' `CompoundDb` objects provide access to general (metabolite) compound
#' annotations along with *metadata* information such as the annotation's
#' source, date and release version. The data is stored internally in a
#' database (usually an SQLite database).
#'
#' @details
#'
#' `CompoundDb` objects should be created using the constructor function
#' `CompoundDb` providing the name of the (SQLite) database file providing
#' the compound annotation data.
#'
#' @usage
#' dbconn(x)
#' show(object)
#' 
#' @param object For all methods: a `CompoundDb` object.
#'
#' @param x For `CompoundDb`: `character(1)` with the file name of the SQLite
#'     compound database. Alternatively it is possible to provide the connection
#'     to the database with parameter `x`.
#'
#'     For all other methods: a `CompoundDb` object.
#' 
#' @author Johannes Rainer
#'
#' @md
#'
#' @seealso `createCompoundDb()` for the function to create a SQLite compound
#'     database.
#'
#' @examples
#'
#' ## Create a small CompoundDb from a provided HMDB subset
#' cmps <- generate_hmdb_tbl(system.file("extdata/hmdb/hmdb_sub.xml",
#'     package = "PeakABro"))
#' metad <- data.frame(name = c("source", "url", "source_version",
#'     "source_date", "organism"),
#'     value = c("sub_HMDB", "http://www.hmdb.ca", "4", "2017", "Hsapiens"),
#'     stringsAsFactors = FALSE)
#' ## Create the SQLite database:
#' db_file <- createCompoundDb(cmps, metadata = metad, path = tempdir())
#'
#' ## Create a CompoundDb object
#' cmp_db <- CompoundDb(db_file)
#' cmp_db
NULL

#' @exportClass CompoundDb
.CompoundDb <- setClass("CompoundDb",
                        slots = c(dbcon = "DBIConnection",
                                  .properties = "list"),
                        prototype = list(.properties = list(),
                                         dbcon = NULL))

#' @importFrom methods validObject
setValidity("CompoundDb", function(object) {
    if (!is.null(object@dbcon))
        .validCompoundDb(object@dbcon)
    else TRUE
})

#' @importFrom DBI dbListTables dbGetQuery
.validCompoundDb <- function(x) {
    txt <- character()
    tables <- dbListTables(x)
    required_tables <- c("compound", "metadata")
    got <- required_tables %in% tables
    if (!all(got))
        stop("Required tables ", paste0(required_tables[!got]), "not found",
             " in the database")
    ## Check table columns.
    comps <- dbGetQuery(x, "select * from compound limit 3")
    res <- .valid_compound(comps, error = FALSE)
    if (is.character(res))
        txt <- c(txt, res)
    metad <- .metadata(x)
    res <- .valid_metadata(metad, error = FALSE)
    if (is.character(res))
        txt <- c(txt, res)
    if (length(txt)) txt else TRUE
}

#' @description `CompoundDb` *constructs* a `CompoundDb` object by connecting
#'     to the provided database file.
#'
#' @export
CompoundDb <- function(x) {
    if (missing(x))
        stop("Argument 'x' is required")
    if (is.character(x)) {
        ## Assume it's the file name of the SQLite database, open it read only
        x <- dbConnect(dbDriver("SQLite"), dbname = x,
                         flags = RSQLite::SQLITE_RO)
    }
    if (is(x, "DBIConnection")) {
        res <- .validCompoundDb(x)
        if (is.character(res))
            stop(res)
        cdb <- .CompoundDb(dbcon = x)
        return(cdb)
    }
    stop("Can not create a 'CompoundDb' from 'x' of type '", class(x), "'.")
}

.metadata <- function(x) {
    if (!is(x, "DBIConnection"))
        x <- .dbconn(x)
    dbGetQuery(x, "select * from metadata")
}

.metadata_value <- function(x, key) {
    metad <- .metadata(x)
    metad[metad$name == key, "value"]
}

.dbconn <- function(x) {
    x@dbcon
}
