#' Gencode prefaces mouse releases with an M, and human releases with nothing.
#'
#' @param release   the release name (e.g. M24 or 33) 
#'
#' @return          the release species (Mus musculus or Homo sapiens)
#' 
#' @export
get_gencode_species <- function(release) { 

  release <- sub("^VM", "M", ignore=TRUE, release) 
  ifelse(toupper(substr(release, 1, 1)) == "M", "Mus musculus", "Homo sapiens")

}
