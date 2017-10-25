# Helper functions --------------------------------------------------------

#' Robust helper functions to process structures with rcdk
#' 
#' Some helper functions modified from rcdk to return NA on failure.
#' This makes it easier to process many compounds where some might behave unexpected.
#' 
#' @importFrom rcdk do.aromaticity do.typing do.isotopes get.exact.mass
#' @importFrom purrr possibly
#' @importFrom rinchi parse.inchi
#' @return Returns functions equavalent of \code{\link{do.aromaticity}}, \code{\link{do.typing}}, \code{\link{do.isotopes}}, \code{\link{get.exact.mass}} but made robust with \code{\link{possibly}}.
#' 
#' @keywords internal
#' 
#' @describeIn do.aromaticity_p Robust version of \code{\link{do.aromaticity}}
do.aromaticity_p <- possibly(do.aromaticity, NA)

#' @describeIn do.aromaticity_p Robust version of \code{\link{do.typing}} 
do.typing_p <-  possibly(do.typing, NA)

#' @describeIn do.aromaticity_p Robust version of \code{\link{do.isotopes}}
do.isotopes_p <-  possibly(do.isotopes, NA)

#' @describeIn do.aromaticity_p Robust version of \code{\link{get.exact.mass}}
get.exact.mass_p <- possibly(get.exact.mass, NA)

#' @describeIn do.aromaticity_p Robust version of \code{\link{parse.inchi}}
parse.inchi_p <- possibly(parse.inchi, NA)



#' Helper functions to extract info from LipidBlast json.
#'
#' Some helper functions to extract info from LipidBlast json.
#'
#' @param json LipidBlast json
#'
#' @importFrom dplyr %>%
#' @importFrom rlist list.select list.flatten
#' @return Returns character vector with the relevant info.
#'
#' @keywords internal
#'
#' @describeIn LB_getname Extracts the compound name from LipidBlast json
LB_getname <- function(json){
    # make check happy
     . <- name <- NULL
     
    json %>% 
        {.$compound[[1]]$names} %>% 
        list.select(name) %>% 
        list.flatten %>% 
        paste(collapse="; ")
    }

#' @describeIn LB_getname Extracts the inchi from LipidBlast json
LB_get_inchi <- function(json){
    # make check happy
    . <- NULL
     
    json %>% 
        {.$compound[[1]]$inchi}
}

#' @describeIn LB_getname Extracts the id name from LipidBlast json
LB_get_id <- function(json){
    # make check happy
    . <- NULL
     
    json %>% 
    {.$id}
}



#' Inchi to formula
#'
#' @param inchi Character vector with inchis
#'
#' @return Character vector
# '
#' @export
#'
#' @importFrom stringr str_split_fixed
#'
inchi2formula <- function(inchi){
    str_split_fixed(inchi,"/",3)[,2]
}



#' Inchi to mass
#'
#' @param inchi Character vector with inchis
#'
#' @return Numeric vector
#'
#' @export
#'
#'

inchi2mass <- function(inchi){

    mol <- parse.inchi_p(inchi)[[1]]

    do.aromaticity_p(mol)
    do.typing_p(mol)
    do.isotopes_p(mol)
    mass <- get.exact.mass_p(mol)

    return(mass)
}



# formula2mass  <- . %>% {map_chr(., ~getMass(getMolecule(..1)))} #  rdisop is crashy on windows 10



# Functions to extract a table of data from different database for --------


#' Generate table with LipidMaps compounds
#'
#' @param sdf_path Path to the Lipidmaps SDF file available at: http://www.lipidmaps.org/resources/downloads/
#'
#' @return A \code{\link[tibble]{tibble}} containing the columns: id, name, inchi, formula, and mass.
#'
#' @export
#'
#' @importFrom ChemmineR read.SDFset datablock datablock2ma
#' @importFrom tibble as_data_frame
#' @importFrom tidyr replace_na unite
#' @importFrom dplyr mutate select %>%
#'
#' @family compound table creation functions
#' 
#' @examples
#' \dontrun{
#' lipidmaps_tbl <- generate_lipidmaps_tbl("path/to/LMSDFDownload6Dec16FinalAll.sdf")
#' 
#' # Alternatively get pre-made table from the package
#' lipidmaps_tbl2 <- readRDS(system.file("extdata", "lipidmaps_tbl.rds", package="PeakABro"))
#' }
#'
#'

generate_lipidmaps_tbl <- function(sdf_path){

    # make check happy
    COMMON_NAME <-
    SYNONYMS <-
    SYSTEMATIC_NAME <-
    INCHI <-
    FORMULA <-
    EXACT_MASS <-
    LM_ID <-
    name <-
    mass <-
        NULL

    SDF <- read.SDFset(sdf_path)

    SDF_table <- datablock2ma(datablock(SDF))

    lipidmaps_tbl <- SDF_table %>%
                     as_data_frame %>%
                     replace_na(list(COMMON_NAME = "", SYNONYMS = "", SYSTEMATIC_NAME = "")) %>%
                     unite(name,COMMON_NAME, SYNONYMS, SYSTEMATIC_NAME, sep="; ") %>%
                     mutate(name = sub("^; .{1}", "", name)) %>%
                     mutate(name = sub("; $", "", name)) %>%
                     mutate(name = gsub("; ; ", "; ", name)) %>%
                     select(id = LM_ID, name, inchi = INCHI, formula = FORMULA, mass = EXACT_MASS) %>%
                     mutate(mass=as.numeric(mass))

    return(lipidmaps_tbl)
}



#' Generate table with LipidBlast compounds
#'
#' @param json_path Path to the LipidBlast json file available at: http://mona.fiehnlab.ucdavis.edu/downloads
#'
#' @return A \code{\link[tibble]{tibble}} containing the columns: id, name, inchi, formula, and mass.
# '
#' @export
#'
#' @importFrom jsonlite read_json
#' @importFrom tibble data_frame
#' @importFrom purrr map_chr map_dbl
#' @importFrom dplyr %>% mutate
#' @importFrom magrittr %<>%
#'
#' @family compound table creation functions
#' 
#' @examples
#' \dontrun{
#' lipidblast_tbl <- generate_lipidblast_tbl("path/to/MoNA-export-LipidBlast.json")
#' 
#' # Alternatively get pre-made table from the package
#' lipidblast_tbl2 <- readRDS(system.file("extdata", "lipidblast_tbl.rds", package="PeakABro"))
#' }
#'
#'

generate_lipidblast_tbl <- function(json_path){

    # make check happy
    . <-
        inchi <-
        NULL

    lipidblast <- read_json(json_path)

    lipidblast_tbl <- lipidblast %>% {data_frame(
                                                  id = map_chr(.,LB_get_id),
                                                  name = map_chr(.,LB_getname),
                                                  inchi = map_chr(.,LB_get_inchi)
                                                )
                                     }


    lipidblast_tbl %<>% mutate(formula = inchi2formula(inchi)) %>%
                        # mutate(mass = formula2mass(formula)) %>%
                        mutate(mass = map_dbl(inchi, inchi2mass)) # slow but rdisop is crashy on windows 10

}


#' @title Generate a table with HMDB compounds
#'
#' @description `generate_hmdb_tbl` processes one or more HMDB xml files
#'     (downloaded from http://http://www.hmdb.ca/) or an HMDB file in SDF
#'     format and extracts general compound information.
#'
#' @note At present only a subset of the available data provided by HMDB are
#'     extracted.
#' 
#' @param file `character` with the name(s) of xml files downloaded from HMDB,
#'     or the file name of a HMD file in SDF format.
#'     For xml files: can be a single xml file (containing all compounds from
#'     an HMDB release) or the names of several xml files, each supposed to
#'     provide data for one compound.
#'
#' @return A [tibble::tibble] with general compound information (one row per
#' compound):
#' + `id`: the HMDB ID of the compound.
#' + `name`: the compound's name.
#' + `inchi`: the inchi of the compound.
#' + `formula`: the chemical formula.
#' + `mass`: the compound's mass. The value from HMDB's
#'   `monisotopic_molecular_weight` is returned.
#'
#' @importFrom dplyr as.tbl
#' @importFrom tools file_ext
#' @importFrom purrr map_dfr
#'
#' @family compound table creation functions
#' 
#' @export
#' 
#' @md
#'
#' @export
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Extract compound information from a single xml file from HMDB.
#' fl <- system.file("extdata/hmdb/HMDB0000001.xml", package = "PeakABro")
#' hmdb <- generate_hmdb_tbl(fl)
#' colnames(hmdb)
#' hmdb
#' 
#' ## Extract compound information from a (subset of a) HMDB release xml file.
#' fl <- system.file("extdata/hmdb/hmdb_sub.xml", package = "PeakABro")
#' hmdb <- generate_hmdb_tbl(fl)
#' colnames(hmdb)
#' hmdb
#'
#' ## Extract compound information from a HMDB file in SDF format
#' fl <- system.file("extdata/hmdb/hmdb_sub.sdf", package = "PeakABro")
#' hmdb <- generate_hmdb_tbl(fl)
#' hmdb
#' 
#' ## Extract the masses
#' hmdb$mass
generate_hmdb_tbl <- function(file) {
    if (missing(file))
        stop("'file' is required")
    file <- normalizePath(file)
    if (any(!file.exists(file)))
        stop("Can not find file: ", paste(file[!file.exists(file)],
                                          collapse = ", "))
    ## Check if we've got an xml file or a sdf file...
    if (any(tolower(file_ext(file)) == "xml")) {
        res <- map_dfr(file, simple_parse_hmdb_xml)
    } else {
        ## Assume we've go A SINGLE file in SDF format.
        res <- map_dfr(file, parse_hmdb_sdf)
    }
    as.tbl(res)
}

#' @description Parse one HMDB xml file and return its results as a `data.frame`
#'
#' @details Report the `"monisotopic_molecular_weight"` in column `"mass"`.
#'    The function contains two different parsing approaches: a fast and a slow
#'    one. The fast one works if the file is in the correct format and should
#'    hence always work. The second approach, that iterates through the
#'    individual metabolite definitions within the xml file is much slower, but
#'    is expected to be fail save.
#'
#' @param x `character(1)` with the file name.
#'
#' @return `data.frame` with columns:
#' + `id`: the HMDB ID of the compound.
#' + `name`: the compound's name.
#' + `inchi`: the inchi of the compound.
#' + `formula`: the chemical formula.
#' + `mass`: the compound's mass. The value from HMDB's
#'   `monisotopoc_molecular_weight` is used.
#' 
#' @md
#'
#' @importFrom xml2 read_xml xml_find_all xml_find_first xml_text xml_double
#'     xml_ns
#'
#' @author Johannes Rainer
#' 
#' @noRd
#'
#' @examples
#'
#' x <- system.file("extdata/hmdb/HMDB0000001.xml", package = "PeakABro")
#' simple_parse_hmdb_xml(x)
#' x <- system.file("extdata/hmdb/hmdb_sub.xml", package = "PeakABro")
#' simple_parse_hmdb_xml(x)
simple_parse_hmdb_xml <- function(x) {
    x <- read_xml(x)
    ns <- ""
    if (length(xml_ns(x)))
        ns <- "d1:"
    x_met <- paste0("//", ns, "metabolite")
    metabolites <- xml_find_first(x, x_met)
    if (length(metabolites) == 0)
        stop("Can not find a single metabolite in the XML file")
    ## Now, test if we can get all in one go, otherwise we have to switch to
    ## a (much) slower version
    res <- list(
        id = xml_text(xml_find_all(x, paste0(x_met, "/", ns,
                                             "accession"))),
        name = xml_text(xml_find_all(x, paste0(x_met, "/", ns, "name"))),
        inchi = xml_text(xml_find_all(x, paste0(x_met, "/", ns, "inchi"))),
        formula = xml_text(xml_find_all(x, paste0(x_met, "/", ns,
                                                  "chemical_formula"))),
        mass = xml_double(xml_find_all(x,
                                       paste0(x_met, "/", ns,
                                              "monisotopic_molecular_weight")))
    )
    if (length(unique(lengths(res))) == 1) {
        res <- as.data.frame(res, stringsAsFactors = FALSE)
    } else {
        ## Fail back to a much slower version.
        message("Something went wrong during the processing of the file. Have",
                " to switch to a slower parser.")
        ## Define the function to extract the data
        hmdb_extract_metabolite <- function(met, ns) {
            data.frame(
                id = xml_text(
                    xml_find_first(met, paste0("./", ns, "accession"))),
                name = xml_text(
                    xml_find_first(met, paste0("./", ns, "name"))),
                inchi = xml_text(
                    xml_find_first(met, paste0("./", ns, "inchi"))),
                formula = xml_text(
                    xml_find_first(met, paste0("./", ns, "chemical_formula"))),
                mass = xml_double(
                    xml_find_first(met,
                                   paste0("./", ns,
                                          "monisotopic_molecular_weight"))),
                stringsAsFactors = FALSE
            )
        }
        metabolites <- xml_find_all(x, x_met)
        res <- do.call(rbind,
                       lapply(metabolites, FUN = hmdb_extract_metabolite,
                              ns = ns))
    }
    res
}

#' @description Parse the HMDB data in SDF format.
#'
#' @importFrom ChemmineR read.SDFset datablock datablock2ma
#'
#' @md
#' 
#' @noRd
parse_hmdb_sdf <- function(x) {
    dblock <- datablock(read.SDFset(x))
    full_mat <- datablock2ma(dblock)
    data.frame(id = full_mat[, "DATABASE_ID"],
               name = full_mat[, "GENERIC_NAME"],
               inchi = full_mat[, "INCHI_IDENTIFIER"],
               formula = full_mat[, "FORMULA"],
               mass = as.numeric(full_mat[, "EXACT_MASS"]),
               stringsAsFactors = FALSE)
}

#' Prepare compound table for interactive display
#' 
#' This function modifies some of the columns in a compound table to be displayed better in an interactive table.
#' To the inchi column code to crop very long strings will be added.
#' To the name column "; " separater will be changed to a line break.
#' This is mainly used internally.
#'
#' @param cmp_tbl \code{\link[tibble]{tibble}} of compounds.
#'
#' @return A \code{\link[tibble]{tibble}} containing the same columns as the input table.
#'
#' @export
#'
#' @importFrom dplyr %>% mutate filter select
#' @importFrom magrittr %<>%
#'
#'

cmp_tbl_pretty <- function(cmp_tbl){
        
    # make check happy
        inchi <-
        name <-
        NULL
    
    cmp_tbl %<>% 
                    mutate(inchi = paste0('<div style= "-o-text-overflow: ellipsis; text-overflow: ellipsis;  overflow:hidden;  white-space:nowrap;   width: 20em;">',inchi,'</div>')) %>% 
                    mutate(name = gsub("; ","<br>",name))
}
