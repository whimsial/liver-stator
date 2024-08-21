#' Helper function to load a given Rdata object.
load.rdata <- function(filename) {
    load(filename)
    get(ls()[ls() != "filename"])
}
