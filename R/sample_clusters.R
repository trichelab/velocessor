#' downsample by clusters for NMF/NTF exploration
#' 
#' @param txis  a SingleCellExperiment with $cluster set
#' @param n     how many cells per cluster (default is 100)
#'
#' @return      a downsampled SingleCellExperiment
#'
#' @import SingleCellExperiment
#'
#' @export
sample_clusters <- function(txis, n=100) {

  if (!"cluster" %in% names(colData(txis))) stop("You need to cluster first.")

  clusts <- unique(txis$cluster)
  names(clusts) <- clusts

  cells <- do.call(c, lapply(clusts, .block_sample, txis=txis, n=n))
  txis[, cells]

}


# helper fn
.block_sample <- function(clust, txis, n=100) {
 
  message("Sampling cluster ", clust, "...") 
  in_clust <- rownames(subset(colData(txis), txis$cluster == clust))
  sample(in_clust, n)

}
