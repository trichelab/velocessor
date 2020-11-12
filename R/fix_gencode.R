#' helper function to make Gencode version decoding slightly less hellish.
#' 
#' @param release   the release string
#' 
#' @return          a standardized release string ("33" or "M24" or the like)
#'
#' @export
fix_gencode <- function(release) {
  
  release <- toupper(release)
  release <- sub("^VM", "M", ignore=TRUE, release)
  release <- sub("gencode_", "", ignore=TRUE, release)
  release <- sub("gencode", "", ignore=TRUE, release)
  return(release) 

} 
