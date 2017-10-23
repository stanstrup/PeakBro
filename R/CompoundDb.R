#' @include createCompoundDbPackage.R

#' @name CompoundDb
#'
#' @title Simple compound (metabolite) databases
#'
#' @description `CompoundDb` objects provide access to general (metabolite)
#'     compound annotations along with *metadata* information such as the
#'     annotation's source, date and release version. The data is stored
#'     internally in a database (usually an SQLite database).
#'
#' @author Johannes Rainer
#'
#' @md
NULL

#' @exportClass CompoundDb
.CompoundDb <- setClass("CompoundDb",
                        slots = c(dbcon = "DBIConnection",
                                  .properties = "list"),
                        prototype = list(.properties = list()))

setValidity("CompoundDb", function(object) {
    txt <- character()
    ## Validity check should be as generic as possible, allowing also for other
    ## implementations (i.e. without SQLite or forcing a certain database
    ## layout. Ideally use the accessor functions here.
    ## 1) metadata should NOT be empty: use .valid_metadata function.
    ## 2) compound data should have required columns.
    if (length(txt)) txt else TRUE
})
