#' convenience function: split a [Summarized/SingleCell]Experiment by columns
#'
#' @param sce   the SE-like object
#' @param splt  the factor to split by, or a colData column in sce
#' 
#' @return      a list of SE or SCE objects, one per level of splt
#'
#' @import      SingleCellExperiment
#' 
#' @export
split_SE <- function(sce, splt) {
  if (!is(splt, "factor")) {
    if (splt %in% names(colData(sce))) {     
      splt <- factor(colData(sce)[, splt])
    } else {
      splt <- factor(splt)
    }
  } 
  levs <- levels(splt)
  names(levs) <- levs
  lapply(levs, function(lev) sce[, splt == lev])
}
