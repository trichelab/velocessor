#' simple automated doublet marking for pipelining
#' 
#' Note: if colLabels(txis) exists, it will be used to supply cluster labels.
#' Note: we skip the xgboost option in favor of "weighted" for HPC usage.
#'
#' @param txis  a SingleCellExperiment
#' @param score scoring method (default is "weighted", use "xgb" if you can)
#' @param clust use clusters? (default is TRUE; set to false with many samples)
#' @param drop  drop doublets (instead of just marking them)? (TRUE; drop them) 
#' @param ...   more params for scDblFinder, e.g. `trajectoryMode` or `BPPARAM`
#' 
#' @return      if successful, SingleCellExperiment with doublets marked/dropped
#' 
#' @import scDblFinder
#'
#' @export
drop_doublets <- function(txis, score="weighted", clust=TRUE, drop=TRUE, ...) {

  # automated during process_velo_txis
  if (!"sample" %in% names(colData(txis))) {
    txis$sample <- get_sample_from_barcode(txis)
  } 

  # if colLabels(txis) is NULL, then fix that.
  if (is.null(colLabels(txis))) { 
    txis <- cluster_velo_txis(txis) 
  }

  # transitioning from colData(.)$cluster to colLabels(.) 
  if (clust) { 
    
    message("Trying cluster-wise doublet removal with sample labels...")
    dbls <- try(silent=TRUE,
                scDblFinder(txis, 
                            clusters=colLabels(txis), 
                            samples=txis$sample, 
                            score=score, 
                            ...))

  } else { 

    message("Trying sample-wise doublet removal without cluster labels...")
    dbls <- try(silent=TRUE,
                scDblFinder(txis, 
                            samples=txis$sample, 
                            score=score, 
                            ...))

  }
  
  # if we fail, return `txis` unscathed
  if (inherits(dbls, "try-error")) {
    message("Failed. Trying cluster-wise removal *without* sample labels...")
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
  if (drop == TRUE) { 
    return(dbls[, dbls$scDblFinder.class == "singlet"])
  } else { 
    return(dbls[, colnames(txis)])
  }

} 
