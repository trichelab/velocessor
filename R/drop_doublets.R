#' @rdname find_doublets
#' @export
drop_doublets <- function(txis, score=NULL, clust=TRUE, ...) {

  dbls <- find_doublets(txis, score=score, clust=clust, ...) 
  if ("scDblFinder.class" %in% names(colData(dbls))) {
    dbls <- dbls[, dbls$scDblFinder.class == "singlet"]
  }
  return(dbls)

} 
