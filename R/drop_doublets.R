#' simple automated doublet marking for pipelining
#' 
#' Note: if colLabels(txis) exists, it will be used to supply cluster labels.
#' 
#' @param txis  a SingleCellExperiment
#' @param ...   more params for scDblFinder, e.g. `trajectoryMode` or `BPPARAM`
#' 
#' @return      if successful, a SingleCellExperiment with doublets dropped
#' 
#' @import scDblFinder
#'
#' @export
drop_doublets <- function(txis, ...) { 

  # automated during process_velo_txis
  if (!"sample" %in% names(colData(txis))) {
    txis$sample <- get_sample_from_barcode(txis)
  } 

  # transitioning from colData(.)$cluster to colLabels(.) 
  dbls <- try(scDblFinder(txis, 
                          clusters=colLabels(txis), 
                          samples=txis$sample, 
                          ...))

  # if we fail, return `txis` unscathed
  if (inherits(dbls, "try-error")) {
    message("Automated doublet marking could not be performed. Skipping.")
    return(txis)
  } else { 
    # but if we succeed, return `txis` without doublets
    return(dbls[, dbls$scDblFinder.class == "singlet"])
  }

} 
