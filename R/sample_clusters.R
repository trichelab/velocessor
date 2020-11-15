#' simple downsampling for NMF/NTF exploratory work and big datasets
#' 
#' Hint: if it doesn't work, you might have too many clusters. Try different 
#' options to cluster_velo_txis until you get a tranctable representation. Feel
#' free to plot the various representations and/or look at silhouette plots, 
#' modularity, whatever until you feel comfortable that the clusters don't suck.
#'
#' And remember: Clustering Is Difficult Only When It Does Not Matter (ref 1). 
#' 
#' @param txis    a SingleCellExperiment with cluster assignments in colLabels
#' @param n       how many cells per cluster to sample at most? (500)
#' @param batch   batch term (will use $sample if none) 
#' @param block   block term (will omit if none) 
#' @param ...     additional arguments to scuttle::downsampleBatches if using
#' 
#' @return        a downsampled SingleCellExperiment
#' 
#' @references \
#' Daniely A, Linial N, Saks M. (2012) \
#' Clustering is difficult only when it does not matter. \
#' https://arxiv.org/abs/1205.4891
#'
#' @seealso       scuttle::downsampleBatches
#'
#' @import        SingleCellExperiment
#' @import        scuttle
#'
#' @export
sample_clusters <- function(txis, n=100) {

  if (!"cluster" %in% names(colData(txis))) stop("You need to cluster first.")

  clusts <- unique(txis$cluster)
  names(clusts) <- clusts

  cells <- do.call(c, lapply(clusts, .down_sample, txis=txis, n=n))
  txis[, cells]

}


# helper fn
.down_sample <- function(clust, txis, n=100) {
 
  message("Sampling cluster ", clust, "...") 
  in_clust <- rownames(subset(colData(txis), txis$cluster == clust))
  sample(in_clust, n)

}
