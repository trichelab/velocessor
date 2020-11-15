#' simple automated doublet marking for pipelining
#' 
#' Note: if colLabels(txis) exists, it will be used to supply cluster labels.
#' Note: we skip the xgboost option in favor of "weighted" for HPC usage.
#'
#' @param txis  a SingleCellExperiment
#' @param score scoring method (default is "weighted", use "xgb" if you can)
#' @param ...   more params for scDblFinder, e.g. `trajectoryMode` or `BPPARAM`
#' 
#' @return      if successful, a SingleCellExperiment with doublets dropped
#' 
#' @import scDblFinder
#'
#' @export
drop_doublets <- function(txis, score="weighted", ...) { 

  # automated during process_velo_txis
  if (!"sample" %in% names(colData(txis))) {
    txis$sample <- get_sample_from_barcode(txis)
  } 

  # transitioning from colData(.)$cluster to colLabels(.) 
  message("Trying cluster-wise doublet removal with sample labels...")
  dbls <- try(silent=TRUE,
              scDblFinder(txis, 
                          clusters=colLabels(txis), 
                          samples=txis$sample, 
                          score=score, 
                          ...))

  # if we fail, return `txis` unscathed
  if (inherits(dbls, "try-error")) {
    message("Failed. Trying doublet removal *without* sample labels...")
    dbls <- try(silent=TRUE, 
                scDblFinder(txis, clusters=colLabels(txis), score=score, ...))
    if (inherits(dbls, "try-error")) {
      message("Automated doublet marking could not be performed. Skipping.")
      return(txis)
    } 
    message("Done.")
  } else { 
    message("Done.")
  } 
  
  # but if we succeed, return `txis` without doublets
  return(dbls[, dbls$scDblFinder.class == "singlet"])

} 
