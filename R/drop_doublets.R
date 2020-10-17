#' simple automated doublet marking for pipelining
#' 
#' @param txis  a SingleCellExperiment
#' @param ...   other arguments, currently unused
#' 
#' @return      if successful, a SingleCellExperiment with doublets dropped
#' 
#' @import scDblFinder
#'
#' @export
drop_doublets <- function(txis, ...) { 

  if (!"cluster" %in% names(colData(txis))) {
    txis$cluster <- cluster_velo_txis(txis)
  } 
  
  if (!"sample" %in% names(colData(txis))) {
    txis$sample <- get_sample_from_barcode(txis)
  } 

  # don't assume success here
  dbls <- try(scDblFinder(txis, 
                          clusters=txis$cluster, 
                          samples=txis$sample))

  # and if we fail, return unscathed
  if (inherits(dbls, "try-error")) {
    message("Automated doublet marking could not be performed. Skipping.")
    return(txis)
  } else { 
    # but if we do succeed, drop the doublets
    return(dbls[, dbls$scDblFinder.class == "singlet"])
  }

} 
