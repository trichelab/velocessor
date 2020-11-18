#' harness for velocity correction (currently ignorant of spliced/unspliced)
#' 
#' Harmony is not a Bioconductor package, so the user needs to install it from 
#' \url{https://github.com/immunogenomics/harmony} prior to using this function.
#' The UMAP projection will be updated to reflect harmonized coordinates unless
#' the user sets UMAP to FALSE when calling the function. 
#' 
#' @param txis      a SingleCellExperiment
#' @param colname   name of colData column with the batch variable ("sample")
#' @param how       how to respect spliced/unspliced relation, currently unused
#' @param iter      max.iter.harmony, defaults to 20 (harmony default is 10) 
#' @param ...       arguments to pass to RunHarmony, e.g. max.iter.harmony=20
#' 
#' @return          a SingleCellExperiment with "HARMONY" in reducedDimNames
#' 
#' @import SingleCellExperiment
#' 
#' @export
harmonize_velo_txis <- function(txis, colname="sample", how="default", iter=20, ...) {
  
  if (!requireNamespace("harmony")) {
    message("You do not have harmony installed.")
    message("Visit https://github.com/immunogenomics/harmony for instructions.")
    stop("Cannot proceed without (non-Bioconductor) 'harmony' package.")
  } 

  if (!colname %in% names(colData(txis))) {
    if (colname == "sample") txis$sample <- get_sample_from_barcode(txis)
    else stop("The column ", colname, " was not found in names(colData(txis)).")
  }

  if (!"PCA" %in% reducedDimNames(txis)) {
    message("Computing PCA (needed by Harmony)...") 
    txis <- runPCA(txis)
  }

  message("Harmonizing...")
  txis <- harmony::RunHarmony(txis, colname, max.iter.harmony=20, ...)
  message("Done. Consider harmony::HarmonyMatrix for complex experiments.")

  # short circuit if no UMAP
  message("Updating UMAP using HARMONY reduced dimensional representation...")
  if ("UMAP" %in% reducedDimNames(txis)) {
    ncomp <- ncol(reducedDim(txis, "UMAP"))
  } else {
    ncomp <- 3 # the minimum that I find useful 
  }
  txis <- runUMAP(txis, ncomponents=ncomp, dimred="HARMONY")
  message("Done. Consider inspecting eig x eig plots for how many PCs to use.")
  return(txis)

}
