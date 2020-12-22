#' simple automated doublet marking for pipelining, refactored a little bit
#'
#' Note: if colLabels(txis) exists, it will be used to supply cluster labels.
#' Note: if xgboost is installed, we use "xgb", else "weighted" for the score.
#'
#' @name        find_doublets
#' @rdname      find_doublets
#' @aliases     drop_doublets
#' 
#' @param txis  a SingleCellExperiment, ideally with sample and cluster labels
#' @param score scoring method ("xgb" if possible, "weighted" otherwise)
#' @param clust use clusters? (default is TRUE; set to false if many samples)
#' @param ...   more params for scDblFinder, e.g. `trajectoryMode` or `BPPARAM`
#' 
#' @return      a SingleCellExperiment with doublets marked (or dropped)
#' 
#' @details     drop_doublets drops doublets. find_doublets just marks them. 
#' 
#' @import      scDblFinder
#'
#' @export
find_doublets <- function(txis, score=NULL, clust=TRUE, ...) {

  # use xgboost if possible
  if (is.null(score) & requireNamespace("xgboost", quietly=FALSE)) {
    score <- "xgb"
  } else if (is.null(score)) { 
    score <- "weighted"
  } 
  message("Using `", score, "` for doublet scoring.")

  # automated during process_velo_txis
  if (!"sample" %in% names(colData(txis))) {
    txis$sample <- get_sample_from_barcode(txis)
  } 

  # if colLabels(txis) is NULL, then fix that.
  if (is.null(colLabels(txis))) txis <- cluster_velo_txis(txis) 

  # transitioning from colData(.)$cluster to colLabels(.) 
  if (clust) { 
    
    message("Trying cluster-wise doublet removal *with* sample labels...")
    dbls <- try(silent=TRUE,
                scDblFinder(txis, 
                            clusters=colLabels(txis), 
                            samples=txis$sample, 
                            score=score, 
                            ...))

  } else { 

    message("Trying sample-wise doublet removal *without* cluster labels...")
    dbls <- try(silent=TRUE,
                scDblFinder(txis, 
                            samples=txis$sample, 
                            score=score, 
                            ...))

  }
  
  # if we fail once, try cluster-wise doublet finding across samples
  if (inherits(dbls, "try-error")) {
    message("Failed. Trying cluster-wise removal *without* sample labels...")
    dbls <- try(silent=TRUE, 
                scDblFinder(txis, clusters=colLabels(txis), score=score, ...))
  } else { 
    message("Done.")
  }

  # if we fail twice, return `txis` unscathed
  if (inherits(dbls, "try-error")) {
    message("Automated doublet marking could not be performed. Skipping.")
    return(txis)
  } else { 
    message("Done.")
  } 
  
  return(dbls[, colnames(txis)])

} 
