#' rescale a reduced-dimensional representation for plotting purposes
#' 
#' Project the quantiles of a given representation onto a Gaussian and use that.
#' 
#' @param txis        a SingleCellExperiment
#' @param dimred      name of the existing reduced dimensional representation
#' 
#' @return            a SingleCellExperiment with a paste0("scaled_", dimred)
#' 
#' @export 
rescale_dimred <- function(txis, dimred="UMAP") {

  if (dimred == "UMAP") {
    rdn <- reducedDimNames(txis)
    if (!"PCA" %in% rdn) txis <- runPCA(txis)
    if (!"UMAP" %in% rdn) {
      if ("HARMONY" %in% rdn) {
        txis <- runUMAP(txis, dimred="HARMONY")
      } else { 
        txis <- runUMAP(txis) # if harmony was run:
      }
    }
  } else if (!dimred %in% reducedDimNames(txis)) {
    stop(dimred, " is not among reducedDimNames(txis). Cannot proceed.")
  }

  scaled_dr <- paste0("scaled_", dimred)
  reducedDim(txis, scaled_dr) <- apply(reducedDim(txis, dimred), 2, .normit)
  return(txis)

}


# helper fn
.normit <- function(x) qnorm((rank(x) + 1) / (length(x) + 2))
