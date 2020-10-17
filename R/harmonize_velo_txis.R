#' test harness for velocity correction, currently just runs on spliced 
#' 
#' Harmony is not a Bioconductor package, so the user needs to install it from 
#' \url{https://github.com/immunogenomics/harmony} prior to using this function.
#' The UMAP projection will be updated to reflect harmonized coordinates unless
#' the user sets UMAP to FALSE when calling the function. 
#' 
#' @param txis      a SingleCellExperiment
#' @param colname   name of colData column with the batch variable ("sample")
#' @param UMAP      update UMAP with runUMAP after harmonizing? (TRUE) 
#' @param how       how to respect spliced/unspliced relation, currently unused
#' @param ...       other arguments to pass to RunHarmony
#' 
#' @return          a SingleCellExperiment with "HARMONY" in reducedDimNames
#' 
#' @import SingleCellExperiment
#' 
#' @export
harmonize_velo_txis <- function(txis, colname="sample", UMAP=TRUE, how=0, ...) {
  
  if (!requireNamespace("harmony")) {
    message("You do not have harmony installed.")
    message("Visit https://github.com/immunogenomics/harmony for instructions.")
    stop("Cannot proceed without (non-Bioconductor) 'harmony' package.")
  } 

  if (!colname %in% names(colData(txis))) {
    if (colname == "sample") txis$sample <- get_sample_from_barcode(txis)
    else stop("The column ", colname, " was not found in names(colData(txis)).")
  }

  message("Harmonizing...")
  txis <- harmony::RunHarmony(txis, colname, ...)

  # short circuit if no UMAP
  if (UMAP) {
    message("Updating UMAP...")
    if ("UMAP" %in% reducedDimNames(txis)) {
      ncomp <- ncol(reducedDim(txis, "UMAP"))
    } else {
      ncomp <- 3 # the minimum that I find useful 
    }
    txis <- runUMAP(txis, ncomponents=ncomp, name="HARMONY")
  }

  message("Done.")
  return(txis)

}
